mast="#####################\n###   BrainGENIE  ###\n#####################"

version = "\n\nversion 0.0.1"

contact = "\n\nJonathan L. Hess, PhD
SUNY Upstate Medical University
hessjo@upstate.edu"

mast_head = c(mast, version, contact)
message(mast_head)

# == load library
require(rtracklayer)
require(plyr)
require(data.table)

# == Functions

# 1. loading trained models
load_models = function(path){
  path_to_models = list.files(path, full.name=T, pattern=".Rdata")
  if(length(path_to_models) < 1) stop("Cannot find files containing trained models! Please make sure path is correctly specified.")
  message("Detected ", length(path_to_models), " files with fitted models!")
  list_models = lapply(path_to_models, function(x) readRDS(x))
  names(list_models) = basename(path_to_models)
  return(list_models)
}

# 2. loading cross-validation performance for genes (R2>0.01, p<0.05 pre-filtered)
load_cv_performance = function(path){
  path_to_models = list.files(path, full.name=T, pattern=".cv.performance.txt")
  cv = data.frame(fread(path_to_models, header = T))
  cv = split(cv, cv$BrainStructure)
  cv.perf <<- cv
}

# 3. impute gene expression using coefficients from GTEx trained models
predict_brain_gxp = function(mod, target, index = NULL, cor = 0.1, pval = 0.05, missing_prop = 0.01){
  
  if(is.null(index) == T) stop("Warning! Value for index is expected.")
  
  # all genes in model
  gene_id  = names(mod)
  
  # genes that have significant cross-validation accuracy
  sig_genes_cv = cv.perf[[index]]
  sig_genes_cv = sig_genes_cv[sig_genes_cv$Cor > cor & sig_genes_cv$Pval < pval, ] # adjustable filter 
  gene_id = gene_id[gene_id %in% sig_genes_cv$gene_id]
  
  message("\rNumber of genes that BrainGENIE will attempt to impute: ", length(gene_id)); cat("\n")
  
  imp.gxp = list() # imputed expression value for gene
  missingness = list()
  for(x in 1:length(gene_id)){
    
    cat("\rImputing expression values:", round(x/length(gene_id) * 100, 1), " %")
    
    # parameters
    intercept = mod[[x]]$intercept
    beta = mod[[x]]$coefs
    names(beta) = gsub("predX", "", names(beta)) # simple patch. TODO: fix the names of LR models later
    
    # missing beta?
    common_beta = intersect(names(beta), colnames(target))
    missing_beta = 1 - length(common_beta)/length(beta)
    
    missingness[[x]] = missing_beta
    names(missingness)[[x]] = gene_id[[x]]
    if(missing_beta > missing_prop) next
    
    # filter target
    sub_target = target[,colnames(target) %in% common_beta]
    sub_target = sub_target[,match(common_beta, colnames(sub_target))]
    weights = beta[names(beta) %in% common_beta]
    weights = weights[match(common_beta, names(weights))]
    
    # sweep function
    beta_hat = sweep(x = sub_target, MARGIN = 2, STATS = weights, FUN = `*`) 
    beta_hat = rowSums(beta_hat,na.rm=T)
    est = beta_hat + intercept
    
    # save imputed gene expression level
    imp.gxp[[x]] = data.frame(est)
    colnames(imp.gxp[[x]]) = gene_id[[x]]
  }
  missing = unlist(lapply(imp.gxp, is.null))
  if(length(which(missing == T)) > 1){
    imp.gxp = imp.gxp[-which(missing == T)]
  }
  imp.gxp.df = do.call(cbind, imp.gxp) # bind column wise
  
  # number of problems genes
  output = length(imp.gxp)
  expected = length(gene_id)
  problems = expected - output
  cat("\n"); message("\rProblems genes: ", problems)
  
  return(list(edat = imp.gxp.df, problem_gene = unlist(missingness)))
  
}
