require(tensorflow)
require(keras)

# re-train PCA models using specified gene lists present in target sample:
retrain_autoencoder_gtex = function(gene_list = NULL, dropout_rate = 0.1, n_batch=25, n_epochs = 100, n_units = 1e3, output = "", tissue = NULL, ncomps = 20, n_folds = 5){
  
  if(base::exists("brain_expr") == FALSE | base::exists("blood_expr") == FALSE){stop("Warning! Paired blood-brain data have not been imported properly. Please run `load_expr_data` first.")}
  if(is.null(n_folds) == T){stop("Specify the number of cross-validation folds")}
  if(is.null(tissue) == T){stop("Specify the name of brain region for model training")}
  if(output == ""){warning("No output directory specified!")}
  
  use_n_comps = ncomps
  
  brain_counts = brain_expr
  blood_counts = blood_expr
  blood_counts = blood_counts[rownames(blood_counts) %in% gene_list, ]
  blood_counts = blood_counts[!rowSums(is.na(blood_counts)), ]
  
 
  
  # set up fold ids for cross-validation
  folds = create.folds(nfolds = n_folds, ids=colnames(brain_counts))
  
  allmods = list();
  for(fold in 1:max(folds)){
    
    cat("\rFold: ", fold);cat("\n")
    
    # define training/testing splits
    x_train = t(blood_counts[,which(folds != fold)])
    y_train = t(brain_counts[,which(folds != fold)])
    
    x_test = t(blood_counts[,!colnames(blood_counts) %in% rownames(x_train)])
    y_test = t(brain_counts[,!colnames(brain_counts) %in% rownames(y_train)])
    model <- keras_model_sequential()
    
    # RELU with dropout layer input
    # Linear activation for output (no activation function specified)
    if(file.exists("model.hdf5") == T){file.remove("model.hdf5")}
    
    model %>%
      layer_dense(units = n_units, input_shape = ncol(x_train)) %>%
      layer_activation_leaky_relu() %>%
      layer_dropout(rate = dropout_rate) %>%
      layer_dense(units = use_n_comps, name = "bottleneck") %>%
      layer_dense(units = n_units) %>%
      layer_dense(units = ncol(x_train))
    
    model %>% compile(
      loss = "mean_squared_error",
      # loss="mean_absolute_error",
      optimizer = "adam"
    )
    
    
    checkpoint <- callback_model_checkpoint(
      filepath = "model.hdf5", 
      save_best_only = TRUE, 
      save_freq  = 'epoch',
      verbose = 0
    )
    
    early_stopping <- callback_early_stopping(patience = 5)
    

    model %>% fit(
      x = x_train,
      y = x_train,
      validation_data = list(x_test, x_test),
      epochs = n_epochs,
      batch_size = n_batch,
      callbacks = list(checkpoint, early_stopping), verbose=1
    )
    
   
    
    pred_train <- predict(model, x_train)
    mse_train <- apply((x_train - pred_train)^2, 1, sum)
    
    pred_test <- predict(model, x_test)
    mse_test <- apply((x_test - pred_test)^2, 1, sum)
    
    # extract hidden layer (bottleneck)
    layer_name <- 'bottleneck'
    intermediate_layer_model <- keras_model(inputs = model$input,
                                            outputs = get_layer(model, layer_name)$output)
    intermediate_output <- predict(intermediate_layer_model, x_train)
    intermediate_output = data.frame(intermediate_output)
    colnames(intermediate_output) = paste("EnC", 1:ncol(intermediate_output), sep="")
    
    # train model to predict brain gene expression using encoder 
    fit = lm(as.matrix(y_train) ~ ., data = intermediate_output)
    
    # prediction accuracy
    pred.train = predict(fit)
    train.cor = lapply(1:ncol(pred.train), function(x) cor(pred.train[,x], y_train[,x]))
    train.cor = unlist(train.cor)
    mean(train.cor)
    
    # predict hidden layer in new data
    bottleneck.test = predict(intermediate_layer_model, x_test)
    bottleneck.test = data.frame(bottleneck.test)
    colnames(bottleneck.test) = paste("EnC", 1:ncol(bottleneck.test), sep="")
    
    pred.test = predict(fit, bottleneck.test)
    test.cor = lapply(1:ncol(pred.test), function(x) cor(pred.test[,x], y_test[,x]))
    test.cor = unlist(test.cor)
    mean(test.cor)
    
    stats = data.frame(gene = colnames(y_train), train_cor = train.cor, test_cor = test.cor)
    stats$N = nrow(x_test)
    allmods[[fold]] = stats
    rm(model)
  }
  
  stats = ldply(allmods)
  
  # deal with missing values (NAs)
  stats$test_cor[is.na(stats$test_cor)] = 0
  stats$train_cor[is.na(stats$train_cor)] = 0
  stats$test_rsq = stats$test_cor^2
  stats$train_rsq = stats$train_cor^2
  stats$test_rsq[is.na(stats$test_rsq)] = 0
  stats$train_rsq[is.na(stats$train_rsq)] = 0
  
  # transform correlation coefficients 
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
  
  perf_train$pval = 2*pnorm(abs(perf_train$Zscore), lower.tail = FALSE) # calculate two-tailed p-value based on combined z-score
  perf_train$fdr = p.adjust(perf_train$pval, "fdr")
  
  all.models.merge = perf_train
  
  all.models.merge$tissue = tissue
  
  # save output as .Rdata file
  saveRDS(all.models.merge, file=paste(output, "/", tissue, "_BrainGENIE_retrain-",ncomps,".Rdata", sep=""))
  
  perf <<- all.models.merge
  
}


# 3. principal component analysis in reference data
recover_latent_space = function(use_n_comps = 20, dropout_rate = 0.1, n_batch=25, n_epochs = 100, n_units = 1e3, newdata = NULL){
  
  if(is.null(newdata) == T){stop("Please supply a data frame with target dataset for recovering latent space")}
  
  
  trained_pca_blood = blood_expr[rownames(blood_expr) %in% colnames(newdata), ]

  return_these_genes = rownames(trained_pca_blood)
  newdata = newdata[,colnames(newdata) %in% return_these_genes]
  
  x_train = data.frame(t(trained_pca_blood))
  newdata = newdata[,match(colnames(x_train), colnames(newdata))]
  
  x_train = as.matrix(x_train)
  newdata = as.matrix(newdata)
  
  if(file.exists("model.hdf5") == T){file.remove("model.hdf5")}
  
  deploy <- keras_model_sequential()
  
  deploy %>%
    layer_dense(units = n_units, input_shape = ncol(x_train)) %>%
    layer_activation_leaky_relu() %>%
    layer_dropout(rate = dropout_rate) %>%
    layer_dense(units = use_n_comps, name = "bottleneck") %>%
    layer_dense(units = n_units) %>%
    layer_dense(units = ncol(x_train))
  
  deploy %>% compile(
    # loss="mean_absolute_error",
    loss="mean_squared_error",
    optimizer = "adam"
  )
  
  summary(deploy)
  
  checkpoint <- callback_model_checkpoint(
    filepath = "model.hdf5", 
    save_best_only = TRUE, 
    save_freq  = 'epoch',
    verbose = 1
  )
  
  early_stopping <- callback_early_stopping(patience = 5)
  
  deploy %>% fit(
    x = x_train, 
    y = x_train,
    validation_data = list(newdata, newdata),
    epochs = n_epochs, 
    batch_size = n_batch,
    callbacks = list(checkpoint, early_stopping), verbose=1
  )
  
  layer_name <- 'bottleneck'
  intermediate_layer_model <- keras_model(inputs = deploy$input,
                                          outputs = get_layer(deploy, layer_name)$output)
  intermediate_output <- predict(intermediate_layer_model, x_train)
  intermediate_output = data.frame(intermediate_output)
  colnames(intermediate_output) = paste("EnC", 1:ncol(intermediate_output), sep="")
  
  # extract latent space model to new data
  latent_space = predict(intermediate_layer_model, newdata)
  latent_space = data.frame(latent_space)
  colnames(latent_space) = colnames(intermediate_output)
  
  return(list(gtex_latent = intermediate_output, newdata_latent = latent_space, genes = return_these_genes))
}




# 5. Fit LR-PCA prediction model to GTEx counts, then apply the resulting weights to the PCA scores derived from new samples
impute_from_latent_space = function(latent_space = NULL){
  
  sig_genes = perf[perf$Cor >= 0.1 & perf$fdr < .05, ]

  # TODO: Linear regression: brain gene ~ blood-based latent variables
  Y = data.frame(t(brain_expr)) # select normalized GTEx brain counts
  # filter Y by genes that are well predicted from 5-fold cross-validation
  
  X = data.frame(latent_space$gtex_latent) # use PCA model derived from GTEx paired blood-brain data
  fit = lm(as.matrix(Y) ~ ., data = X) # fit a LR model per gene
 
  predict_vals = predict(fit, data.frame(latent_space$newdata_latent))
  predict_vals = data.frame(predict_vals)
  predict_vals = predict_vals[,colnames(predict_vals) %in% sig_genes$gene]
  
  return(predict_vals)
}


