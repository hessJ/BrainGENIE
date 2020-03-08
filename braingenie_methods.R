mast="#########################################\n###   BrainGENIE (LR-PCA imputation)  ###\n#########################################"

version = "\n\nversion 0.0.1"

contact = "\n\nJonathan L. Hess, PhD
SUNY Upstate Medical University
hessjo@upstate.edu"

mast_head = c(mast, version, contact)
message(mast_head)

# == load library
require(rtracklayer)
require(plyr)
require(readxl)
require(data.table)
require(rtracklayer)


# 1. Load brain data from GTEx samples
load_expr_data = function(path_to_data = NULL){
  # warning message
  if(is.null(path_to_data) == T){stop("Please provide a directory path")}
  # load data
  data_files = list.files(path_to_data, pattern=".Rdata", full.names=T)
  blood_expr = data.frame(t(readRDS(data_files[grepl("WholeBlood", data_files)])))
  brain_expr = data.frame(t(readRDS(data_files[!grepl("WholeBlood", data_files)])))
  # retain same samples
  common.sids = intersect(colnames(brain_expr), colnames(blood_expr))
  blood_expr <<- blood_expr[,colnames(blood_expr) %in% common.sids]
  brain_expr <<- brain_expr[,colnames(brain_expr) %in% common.sids]
}


# 2. load cross-validation performance 
load_cv_performance = function(file_path = NULL){
  if(is.null(file_path)){stop("Please provide full path to .Rdata file with cross-validation accuracies (if relying on pre-trained models")}
  return(data.frame(readRDS(file_path)))
}

# 3. principal component analysis in reference data
fit_pca = function(gene_list = NULL, autorun = TRUE){
  
  trained_pca_blood = blood_expr
  
  if(autorun == FALSE){
  if(is.null(gene_list)){stop("Please provide a list of genes that are present in new samples in order to run PCA correctly")}
  matched_genes = intersect(rownames(trained_pca_blood), gene_list)
  if(length(matched_genes) < 1){"No genes common to reference and new samples! Please check that IDs are in ENSEMBL gene ID format."}
  trained_pca_blood = trained_pca_blood[rownames(trained_pca_blood) %in% matched_genes, ]
  }
  
  return_these_genes = rownames(trained_pca_blood)
  
  pca_blood = prcomp(t(trained_pca_blood))
  
  return(list(pca = pca_blood, genes = return_these_genes))
}

# 4. predict PCA in new data
predict_pca = function(dat = NULL, pca_model = NULL, mean_imputation = TRUE){
  
  if(is.null(dat)){stop("Please supply data frame for new samples")}
  if(is.null(pca_model)){stop("Please supply fitted PCA model")}
  if(class(pca_model) != "list"){stop("Expecting list object for fitted PCA model")}
  
  new_samples = dat[,colnames(dat) %in% pca_model$genes]
  
  if(mean_imputation == TRUE){
    
    means = colMeans(new_samples,na.rm=TRUE)
    for(n in 1:ncol(new_samples)){
      new_samples[is.na(new_samples[,n]),n] = means[[n]]
    }
    
  }
  
  if(mean_imputation == FALSE){
    
    na_detected = colSums(is.na(new_samples))
    na_detected = na_detected[na_detected > 0]
    if(length(na_detected) > 0){stop("Warning! Missing values (NAs) detected in new samples. Switch mean_imputation to TRUE to resolve this issue.")}
    
  }

  preds = predict(pca_model$pca, newdata = new_samples)
  return(preds)
}


# 5. Fit LR-PCA prediction model to GTEx counts, then apply the resulting weights to the PCA scores derived from new samples
fit_lr_weights_in_gtex = function(pca_model  = NULL, cv_performance = NULL, n_comps = 20){
  if(is.null(pca_model)){stop("Argument required for blood_PCA")}
  # TODO: Linear regression: brain gene ~ blood PCA
  Y = data.frame(t(brain_expr)) # select normalized GTEx brain counts
  # filter Y by genes that are well predicted from 5-fold cross-validation
  matched_genes = intersect(colnames(Y), cv_performance$gene)
  if(length(matched_genes) < 1){stop("Unexpected issue! No matching gene IDs between CV performance file and normalized GTEx counts.")}
  common_genes = intersect(matched_genes, blood.pca$genes)
  if(length(common_genes) < 1){stop("Unexpected issue! No genes available to run LR.")}
  message("\rTraining LR models for: ", length(common_genes), " genes")
  Y = Y[,colnames(Y) %in% cv_performance$gene, ]
  X = data.frame(pca_model$pca$x[,1:n_comps]) # use PCA model derived from GTEx blood counts
  fit = lm(as.matrix(Y) ~ ., data = X) # fit a LR model per gene
  return(fit) # return model
}

# 6. apply brain transcriptome prediction model to blood PCA
impute_gxp = function(pca_in_new_sample = NULL, trained_model = NULL, scale = TRUE){
  
  if(is.null(trained_model)){stop("Please specify object containing pre-trained models")}
  
  if(is.null(pca_in_new_sample)){stop("Please provide blood-based transcriptome PCs for new samples")}
  
  predict_vals = predict(trained_model, data.frame(pca_in_new_sample))
  
  if(scale == TRUE){
    
  predict_vals = scale(predict_vals)
  }
  
  return(data.frame(predict_vals))
}


# additional functions:

# > convert hgnc symbols to ensembl gene ids (gencode v26)
convert_hgnc_to_ensg = function(hgnc = NULL, gtf_path = NULL){
  
  # check if parameters are filled in
  if(is.null(gtf_path)){stop("Please specify path to the gtf file in the BrainGENIE provided repository")}
  if(is.null(hgnc)){stop("Please specify a vector of HGNC gene symbols")}

  # load GTF file
  gtf = data.frame(fread(gtf_path))
  gtf = gtf[,colnames(gtf) %in% c("gene_id", "gene_name")]
  
  symbols = data.frame(gene_name = hgnc)
  convert = merge(symbols, gtf, by='gene_name')
  return(convert)
}
