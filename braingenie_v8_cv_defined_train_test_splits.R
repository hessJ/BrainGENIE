
# -- require packages
require(data.table)
require(ggplot2)
require(plyr)
require(dplyr)
require(tidyr)
source("~/Google Drive/mac_storage/script_lib/braingenie.functions.R")

######################################################################
######################################################################
#################### custom braingenie functions #####################
######################################################################
######################################################################

# function to run repeated cross-validation
expr_split_for_train_test = function(path_to_expr_data = NULL, brain_structure = NULL,  training_ids = NULL, test_set_ids = NULL){
  
  if(is.null(path_to_expr_data)){stop("Please supply full path to paired blood-brain transcriptome data (instructions for download @ https://github.com/hessJ/BrainGENIE/tree/master/expression")}
  if(is.null(brain_structure)){stop("Please supply name of brain region for model training, must comply with specific GTEx nomenclature (see documentation)!")}
  if(is.null(training_ids)){stop("Please supply GTEx subject IDs for training sample")}
  if(is.null(test_set_ids)){stop("Please supply GTEx subject IDs that comprise the test sample")}
  
  # == loading expression data into memory
  mod_path=paste(path_to_expr_data, brain_structure ,sep="")
  data_files = list.files(mod_path, full.names=T)
  data_files = data_files[!grepl("Correlations", data_files)]
  if(length(data_files) < 1){stop("Warning! Paired blood-brain gene expression data files not found. Note that files need to be in Rdata format.")}
  
  cat("\n")
  cat("\rLoading expression data...")
  blood_counts = data.frame(t(readRDS(data_files[grepl("WholeBlood", data_files)])))
  brain_counts = data.frame(t(readRDS(data_files[!grepl("WholeBlood", data_files)])))
  
  # retain same samples
  common.sids = intersect(colnames(brain_counts), colnames(blood_counts))
  
  blood_counts = blood_counts[,colnames(blood_counts) %in% common.sids]
  brain_counts = brain_counts[,colnames(blood_counts) %in% common.sids]
  
  test.set.sids = test_set_ids
  test.set = gsub("[-]", ".", unique(phenotypes$SUBJID.x[phenotypes$dbGaP_Subject_ID %in% test.set.sids]))
  test.set = gsub("GTEX", "GTEx", test.set)
  
  test.brain = brain_counts[,colnames(brain_counts) %in% test.set]
  test.blood = blood_counts[,colnames(blood_counts) %in% test.set]
  
  # remove test set from cross-validation samples
  brain_counts = brain_counts[,!colnames(brain_counts) %in% test.set]
  blood_counts = blood_counts[,!colnames(blood_counts) %in% test.set]
  
  # create new objects that will store the expression data for the training and test samples
  training_blood <<- blood_counts
  training_brain <<- brain_counts
  testing_blood <<- test.blood
  testing_brain <<- test.brain
}


create.folds = function(nfolds = 10, ids = NULL){
  if(is.null(ids)){stop("Expecting IDs column for creating folds")}
  sample(dplyr::ntile(ids, nfolds)) # random folds
}

run_repeated_cv = function(nfolds = 10, mode = 'lm', brain_structure = NULL, n_repeats = 0, use_n_comps = 20){
  
  list_objects = ls()
  
  save_performance_repeats = list()
  for(r in 1:n_repeats){
    if(n_repeats > 1){
      cat("\n")
      cat("\rCV repeat: ", r)
    }
    # make folds 
    folds = create.folds(nfolds = nfolds, ids=colnames(training_brain))
    
    allmods = list()
    for(fold in 1:max(folds)){
      
      
      # define training/testing splits
      x_train = t(training_blood[,which(folds != fold)])
      y_train = t(training_brain[,which(folds != fold)])
      
      x_test = t(training_blood[,!colnames(training_blood) %in% rownames(x_train)])
      y_test = t(training_brain[,!colnames(training_brain) %in% rownames(y_train)])
      
      # - Run PCA on blood gene expression data
      pca.blood.train = prcomp(x_train)
      n.blood.comps = pca.blood.train$sdev^2 / sum(pca.blood.train$sdev^2)
      # use_n_comps = min(which(cumsum(n.blood.comps) >= pca.prop.min))
      
      
      # apply PCA model to test set
      pca.in.test = predict(pca.blood.train, x_test)
      
      # format PCA predictor matrices based on chosen number of components 
      pca.TRAIN = data.frame(pca.blood.train$x[,1:use_n_comps])
      pca.TEST = data.frame(pca.in.test[,1:use_n_comps])
      
      # if enet
      if(mode =='enet'){
        gene.mods = list()
        nGenes = ncol(y_train)
        # nGenes = 100
        for(g in 1:nGenes){
          cat("\rFitting cv.glmnet:", g)
          fit = cv.glmnet(x = as.matrix(pca.TRAIN), y = y_train[,g], type.measure = 'mse', nfolds = 10, keep = T)
          eval = predict(fit, newx = as.matrix(pca.TRAIN), s = 'lambda.min')[,1]
          eval_cor = cor(eval, y_train[,g])
          
          test.eval = predict(fit, newx = as.matrix(pca.TEST), s = 'lambda.min')[,1]
          test.eval.cor = cor(test.eval, y_test[,colnames(y_test) %in% colnames(y_train)[[g]]])
          
          gene.mods[[g]] = data.frame(gene = colnames(y_train)[[g]], train_cor = eval_cor, train_rsq = eval_cor^2,
                                      test_cor = test.eval.cor, test_rsq = test.eval.cor^2,
                                      N = nrow(pca.TEST))
        }
        allmods[[fold]] = ldply(gene.mods)
      }
      
      if(mode == 'lm'){
        # apply PCA model to predict gene expression in brain
        fitted.model = lm(as.matrix(y_train)  ~ ., data = pca.TRAIN)
        pred.train = predict(fitted.model)
        train.cors = lapply(1:ncol(pred.train), function(x) cor(pred.train[,x], y_train[,x]))
        names(train.cors) = colnames(pred.train)
        train.cors = ldply(train.cors, .id='gene')
        colnames(train.cors)[2] = 'train_cor'
        train.cors$train_rsq = train.cors$train_cor^2
        
        # apply PCA model to test set, evaluate accuracy
        pred.test = predict(fitted.model, newdata = pca.TEST)
        test.cors = lapply(1:ncol(pred.test), function(x) cor(pred.test[,x], y_test[,x]))
        names(test.cors) = colnames(pred.test)
        test.cors = ldply(test.cors, .id='gene')
        colnames(test.cors)[2] = 'test_cor'
        test.cors$test_rsq = test.cors$test_cor^2
        test.cors$N = nrow(pca.TEST)
        allmods[[fold]] = data.frame(cbind(train.cors, test.cors[,-1]))
      }
      
      
    }
    
    stats = ldply(allmods)
    
    # deal with missing values (NAs)
    stats$test_cor[is.na(stats$test_cor)] = 0
    stats$train_cor[is.na(stats$train_cor)] = 0
    stats$test_rsq[is.na(stats$test_rsq)] = 0
    stats$train_rsq[is.na(stats$train_rsq)] = 0
    
    stats$fisherZ = atanh(stats$test_cor)*sqrt(stats$N-3)
    stats = data.table(stats)
    
    
    perf_train = stats[,list(Cor = mean(test_cor, na.rm=TRUE),
                             train.cor = mean(train_cor, na.rm=T),
                             SE = sd(test_cor,  na.rm=TRUE)/sqrt(length(test_cor)),
                             SD = sd(test_cor, na.rm=T),
                             Rsq=mean(test_cor)^2,
                             FisherZ = sum(fisherZ),
                             sample.size = mean(N),
                             Zscore = sum(fisherZ) / sqrt(length(unique(folds)))),by=c("gene")]
    
    perf_train$pval = pnorm(-abs(perf_train$Zscore))
    perf_train$fdr = p.adjust(perf_train$pval, "fdr")
    
    colnames(perf_train)[-c(1:2)] = paste("cv_", colnames(perf_train)[-c(1:2)], sep="")
    perf_train$n_comps = use_n_comps
    
    save_performance_repeats[[r]] = perf_train
  }
  
  save_performance_cv <<- save_performance_repeats
  
}

# compute mean performance over repeated CV folds 
calc_mean_cv_performance = function(nfolds = 10){
  
  names(save_performance_cv) = paste("repeat_", 1:length(save_performance_cv), sep="")
  all_cv = ldply(save_performance_cv, .id='repeats')
  all_cv = all_cv[order(all_cv$gene), ]
  
  
  if(length(save_performance_cv) > 1){
    n_repeats = length(save_performance_cv) # number of repeats that we need to summarize performance over
    
    # summarize performance with ddply
    sum_z = ddply(all_cv, .(gene), summarize, 
                  zscore_est = sum(cv_Zscore) / sqrt(n_repeats), 
                  out_sample_cv_cor = mean(Cor), 
                  in_sample_cv_cor = mean(cv_train.cor))
    sum_z$zscore_pval <- 2*pnorm(abs(sum_z$zscore_est), lower.tail = FALSE)
    sum_z$fdr = p.adjust(sum_z$zscore_pval, 'fdr')
    sum_z$rsq = sum_z$out_sample_cv_cor^2  
  } else {
    sum_z = data.frame(gene = all_cv$gene,
                       zcore_est = all_cv$cv_Zscore, 
                       zscore_pval = all_cv$cv_pval,
                       out_sample_cv_cor = all_cv$Cor,
                       in_sample_cv_cor = all_cv$cv_train.cor,
                       fdr = all_cv$cv_fdr,
                       rsq = all_cv$cv_Rsq)
  }
  
  sum_z_out <<- sum_z
  
}

# apply most generalizable model to test set
apply_models_to_test_set = function(use_n_comps = 20, sig_gene_list = NULL){
  if(is.null(sig_gene_list)){stop("Please provide vector of Ensembl gene IDs for final evaluation in test set")}
  
  cat("\rRetraining models in full training sample (without cross-validation)")
  
  # define training/testing splits
  x_train = t(training_blood)
  y_train = t(training_brain)
  y_train = y_train[,colnames(y_train) %in% sig_gene_list]
  
  x_test = t(testing_blood)
  y_test = t(testing_brain)
  y_test = y_test[,colnames(y_test) %in% sig_gene_list]
  
  # - Run PCA on blood gene expression data
  pca.blood.train = prcomp(x_train)
  n.blood.comps = pca.blood.train$sdev^2 / sum(pca.blood.train$sdev^2)
  
  # apply PCA model to test set
  pca.in.test = predict(pca.blood.train, x_test)
  
  # format PCA predictor matrices based on chosen number of components 
  pca.TRAIN = data.frame(pca.blood.train$x[,1:use_n_comps])
  pca.TEST = data.frame(pca.in.test[,1:use_n_comps])
  
  # apply PCA model to predict gene expression in brain
  fitted.model = lm(as.matrix(y_train)  ~ ., data = pca.TRAIN)
  pred.train = predict(fitted.model)
  train.cors = lapply(1:ncol(pred.train), function(x) cor.test(pred.train[,x], y_train[,x]))
  names(train.cors) = colnames(pred.train)
  est = lapply(train.cors, function(x) x$estimate)
  pval = lapply(train.cors ,function(x) x$p.value)
  train.cors = data.frame(gene = colnames(pred.train), 
                          train_cor = unlist(est),
                          train_pval = unlist(pval))
  train.cors$train_rsq = train.cors$train_cor^2
  
  # apply PCA model to test set, evaluate accuracy
  pred.test = predict(fitted.model, newdata = pca.TEST)
  test.cors = lapply(1:ncol(pred.test), function(x) cor.test(pred.test[,x], y_test[,x]))
  est = lapply(test.cors, function(x) x$estimate)
  pval = lapply(test.cors ,function(x) x$p.value)
  test.cors = data.frame(test_cor = unlist(est), 
                         test_pval = unlist(pval))
  test.cors$test_rsq = test.cors$test_cor^2
  test.cors$N = nrow(pca.TEST)
  
  final_performance_stats <<- data.frame(cbind(train.cors, test.cors))
  predicted_expr_in_test <<- pred.test
  
}
