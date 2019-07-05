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

# 1. Load count data from GTEx samples
load_counts = function(path_to_counts = NULL){
  # warning message
  if(is.null(path_to_counts) == T){stop("Please provide the full directory path to the reference RNA-sequencing counts")}
  temp =  fread(path_to_counts, sep="\t", header = TRUE)
  temp = data.frame(temp)
  rownames(temp) = temp$Name
  temp = temp[,!colnames(temp) %in% c("Name","Description")]
  rnaseq.counts <<- temp
}

# 2. Load sample factors and covariates
load_sample_factors = function(path_to_attributes = NULL, path_to_phenotypes = NULL) {
  # warning message
  if(is.null(path_to_attributes)){stop("Provide path to sample attributes")}
  # warning message
  if(is.null(path_to_phenotypes)){stop("Provide path to sample phenotypes")}
  # == load sample factors files
  
  phenotypes = fread(path_to_attributes)
  phenotypes = data.frame(phenotypes)
  phenotypes$SAMPID = gsub("[-]", ".", phenotypes$SAMPID)
  phenotypes = phenotypes[phenotypes$SAMPID %in% colnames(rnaseq.counts), ]
  
  psid = strsplit(phenotypes$SAMPID, "[.]")
  psid = lapply(psid, function(x) x[[2]])
  phenotypes = data.frame(SID = unlist(psid), phenotypes)
  
  
  dbgap_attributes = data.frame(readxl::read_excel(path_to_phenotypes))
  subid = strsplit(dbgap_attributes$SUBJID, "[-]")
  subid = lapply(subid, function(x) x[[2]])
  subid = unlist(subid)
  dbgap_attributes$SUBJID = subid
  
  colnames(dbgap_attributes)[1] = 'SID'
  
  phens = merge(dbgap_attributes, phenotypes, by='SID')
  return(phens)
}


load_covariates = function(path_to_blood_covariates = NULL, path_to_brain_covariates = NULL){
  # warning message
  if(is.null(path_to_blood_covariates)){stop("Provide path to the covariate file")}
  # warning message
  if(is.null(path_to_brain_covariates)){stop("Provide path to the covariate file")}
  
  # -- load GTEx covariates
  brain_covs = fread(path_to_brain_covariates, header = TRUE)
  brain_covs = data.frame(brain_covs)
  rownames(brain_covs) = brain_covs$ID
  brain_covs = brain_covs[,-1]
  colnames(brain_covs) = gsub("GTEX.", "GTEx.", colnames(brain_covs))

  blood_covs = fread(path_to_brain_covariates, header = TRUE)
  blood_covs = data.frame(blood_covs)
  rownames(blood_covs) = blood_covs$ID
  blood_covs = blood_covs[,-1]
  colnames(blood_covs) = gsub("GTEX.", "GTEx.", colnames(blood_covs))
 
  covariates_from_gtex <<- list(blood_covariates = blood_covs, brain_covariates = brain_covs)
}

# 3. Pair donor ids
paired_donor_ids = function(phens = NULL){
  # warning message
  if(is.null(phens) == T){stop("Provide the object name for GTEX phenotypes")}
  
  # --- blood samples
  blood_samples = phenotypes[phenotypes$SMTSD %in% "Whole Blood", ]
  length(unique(blood_samples$SID))
  
  # -- brain samples (delete unwanted elements from string to match labels in predixcan predictdb data set)
  brain_samples = phenotypes[phenotypes$SMTS %in% "Brain", ]
  brain_samples = brain_samples[!grepl("spinal cord", ignore.case = TRUE, brain_samples$SMTSD), ]
  brain_samples$SMTSD = gsub("Brain - ", "", brain_samples$SMTSD)
  brain_samples$SMTSD = gsub(" ", "_", brain_samples$SMTSD)
  brain_samples$SMTSD = gsub(")", "", brain_samples$SMTSD)
  brain_samples$SMTSD = gsub("\\(", "", brain_samples$SMTSD)
  unique(brain_samples$SMTSD)
  
  # -- number of overlapping samples
  split_names = split(brain_samples[,c("SID", "SMTSD")], brain_samples$SMTSD)
  
  # -- common donor ides
  common_ids <<- lapply(split_names, function(x) intersect(x$SID, blood_samples$SID))
}


load_gtf_file = function(path = NULL, filter_type = "gene"){
  if(is.null(path)){stop("ERROR: Expecting file path for loading the GTF file")}
  if(is.null(filter_type)){stop("Please provide a biotype to filter GTF")}
  loaded = data.frame(rtracklayer::import(path))
  loaded = loaded[loaded$type %in% filter_type,]
  return(loaded)
}


# 4. Normalize transcriptomes
normalize_read_counts  = function(brain_region_index = NULL, min_counts = 5, min_samples = 10, phenotypes = NULL, gtf_table = gtf_table, counts = rnaseq.counts) {
   
  brain_region_index = 1; min_counts = 5; min_samples = 10; counts = rnaseq.counts
  
  # warning message
  if(is.null(brain_region_index) == T){stop("Provide an index for a brain structure to proceed with data normalization properly.")}
  
  b = brain_region_index # which index (brain structure) should be used?
  message("Loading counts from reference data set: ", names(common_ids)[[b]])
  brain_structure = names(common_ids)[[b]]
  
  # filter phenotypes by whole blood and brain structure IDs
  blood_samples = phenotypes[phenotypes$SMTSD %in% "Whole Blood", ]
  brain_samples = phenotypes[phenotypes$SMTS %in% "Brain", ]
  
  # filter blood and brain samples by common ids
  blood_sample = blood_samples[blood_samples$SID %in% common_ids[[b]], ]
  brain_sample = brain_samples[brain_samples$SID %in% common_ids[[b]], ]
  
  # some string munging
  blood_sample$SAMPID = gsub("[-]", ".", blood_sample$SAMPID)
  brain_sample$SAMPID = gsub("[-]", ".", brain_sample$SAMPID)
  brain_sample$SMTSD = gsub("Brain - ", "", brain_sample$SMTSD)
  
  brain_sample = brain_sample[brain_sample$SMTSD %in% brain_structure, ]
  
  # = blood counts
  blood_counts = counts[,colnames(counts) %in% blood_sample$SAMPID]
  blood_counts = blood_counts[,match(blood_sample$SAMPID, colnames(blood_counts))]
  
  # = brain counts
  brain_counts = counts[,colnames(counts) %in% brain_sample$SAMPID]
  brain_counts = brain_counts[,match(brain_sample$SAMPID, colnames(brain_counts))]
  
  # -- remove low counts from blood
  blood_counts = blood_counts[rowSums(blood_counts > min_counts) >= min_samples, ] # min 5 counts in 10 samples
  
  # -- remove low counts from brain
  brain_counts = brain_counts[rowSums(brain_counts > min_counts) >= min_samples, ] # min 5 counts in 10 samples
  
  # -- common genes
  common_genes = intersect(rownames(brain_counts), rownames(blood_counts))
  
  # -- number of genes in common between tissues
  cat("\rNumber of genes in common between tissues: ", length(common_genes))
  
  brain_counts = brain_counts[rownames(brain_counts) %in% common_genes, ]
  blood_counts = blood_counts[rownames(blood_counts) %in% common_genes, ]
  
  # convert COLID to SID
  colnames(brain_counts) = paste("GTEx.", unlist(lapply(strsplit(colnames(brain_counts), "[.]"), function(x) x[[2]])), sep="")
  colnames(blood_counts) = paste("GTEx.", unlist(lapply(strsplit(colnames(blood_counts), "[.]"), function(x) x[[2]])), sep="")
  
  # colnames(blood_counts) = paste("GTEx.", colnames(blood_counts), sep="")
  # colnames(brain_counts) = paste("GTEx.",colnames(brain_counts),sep="")
  
  brain_counts = brain_counts[,match(colnames(blood_counts), colnames(brain_counts))]
  
  # RPKM normalization
  geneSize = gtf_table$width[gtf_table$gene_id %in% rownames(brain_counts)]
  brain_rpkm = edgeR::rpkm(brain_counts, gene.length = geneSize)
  blood_rpkm = edgeR::rpkm(blood_counts, gene.length = geneSize)
  
  # quantile normalization
  brain_rpkm = limma::normalizeQuantiles(brain_rpkm)
  blood_rpkm = limma::normalizeQuantiles(blood_rpkm)
  
  # INT Normalization
  brain_rpkm = data.frame(apply(brain_rpkm, 1, function(x) RNOmni::rankNorm(x)))
  blood_rpkm = data.frame(apply(blood_rpkm, 1, function(x) RNOmni::rankNorm(x)))
  
  # transpose
  blood_counts = data.frame(t(blood_rpkm))
  brain_counts = data.frame(t(brain_rpkm))
  
  # -- load GTEx covariates
  brain_covs <- covariates_from_gtex$brain_covariates
  brain_covs = brain_covs[,colnames(brain_covs) %in% colnames(brain_counts)]
  
  phens = phenotypes
  phens$SID = paste("GTEx.", phens$SID, sep="")
  
  brain_covs = data.frame(t(brain_covs))
  brain_covs = data.frame(SID = rownames(brain_covs), brain_covs)
  brain_covs = merge(brain_covs, phens, by='SID')
  brain_covs = brain_covs[!duplicated(brain_covs$SID), ]
  
  blood_covs = covariates_from_gtex$blood_covariates
  blood_covs = blood_covs[,colnames(blood_covs) %in% colnames(blood_counts)]
  blood_covs = data.frame(t(blood_covs))
  blood_covs = data.frame(SID = rownames(blood_covs), blood_covs)
  blood_covs = merge(blood_covs, phens, by='SID')
  blood_covs = blood_covs[!duplicated(blood_covs$SID), ]
  
  # Common subject IDs between the tissues
  keep_sid = intersect(brain_covs$SID, blood_covs$SID)
  brain_counts = brain_counts[,colnames(brain_counts) %in% keep_sid]
  blood_counts = blood_counts[,colnames(blood_counts) %in% keep_sid]
  brain_counts = brain_counts[,match(keep_sid, colnames(brain_counts))]
  blood_counts = blood_counts[,match(keep_sid, colnames(blood_counts))]
  
  # -- residuals using GTEx supplied covariates (genotype PCs and PEER factors)
  if(nrow(brain_covs) == ncol(brain_counts)){
  brain_counts = resid(lm(t(as.matrix(brain_counts)) ~ C1 + C2 + C3 + AGE + SEX, brain_covs))}
  
  if(nrow(blood_covs) == ncol(blood_counts)){
  blood_counts = resid(lm(t(as.matrix(blood_counts)) ~ C1 + C2 + C3 + AGE + SEX, blood_covs))}
  
  # -- remove effect of age and sex
  brain_counts = data.frame(t(brain_counts))
  blood_counts = data.frame(t(blood_counts))
  
  normalized_counts_gtex <<- list(brain = brain_counts, blood = blood_counts)
}


# 5. principal component analysis in reference data
fit_pca = function(varExplained = 0.8, gene_list = NULL){
  if(is.null(gene_list)){stop("Please provide a list of genes that are present in new samples in order to run PCA correctly")}
  trained_pca_blood = normalized_counts_gtex$blood
  matched_genes = intersect(rownames(trained_pca_blood), gene_list)
  if(length(matched_genes) < 1){"No genes common to reference and new samples! Please check that IDs are in ENSEMBL gene ID format."}
  trained_pca_blood = trained_pca_blood[rownames(trained_pca_blood) %in% matched_genes, ]
  return_these_genes = rownames(trained_pca_blood)
  pca_blood = prcomp(t(trained_pca_blood))
  eigenvalue_blood = (pca_blood$sdev^2)/sum(pca_blood$sdev^2)
  ncomps = cumsum(eigenvalue_blood)
  ncomps_select = length(ncomps[ncomps <= varExplained])
  return(list(PCA = pca_blood, n_comps = ncomps_select, genes = return_these_genes))
}

# 6. predict PCA in new data
predict_pca = function(dat = NULL, pca_model = NULL){
  if(is.null(pca)){stop("Please supply data frame for new samples")}
  if(is.null(pca_model)){stop("Please supply fitted PCA model")}
  if(class(pca_model) != "list"){stop("Expecting list object for fitted PCA model")}
  
  new_samples = dat[,colnames(dat) %in% pca_model$genes]
  new_samples[is.na(new_samples)] = 0
  preds = predict(pca_model$PCA, newdata = new_samples)
  return(dat[,1:pca_model$n_comps])
}

# 7. Fit LR-PCA prediction model to GTEx counts, then apply the resulting weights to the PCA scores derived from new samples
fit_gtex_weights = function(counts = rnaseq.counts, pca_model  = NULL){
  if(is.null(pca_model)){stop("Argument required for blood_PCA")}
  # TODO: Linear regression: brain gene ~ blood PCA
  Y = data.frame(t(normalized_counts_gtex$brain)) # select normalized GTEx brain counts
  X = data.frame(pca_model$PCA$x[,1:pca_model$n_comps]) # use PCA model derived from GTEx blood counts
  fit = lm(as.matrix(Y) ~ ., data = X) # fit a LR model per gene
  return(fit) # return model
}

# 8. apply brain transcriptome prediction model to blood PCA
predict_brain_gxp = function(predicted_pca_in_new_samples = NULL, gtex_model = NULL, scale = TRUE){
  if(is.null(gtex_model)){stop("Please provide weights to gtex_model")}
  if(is.null(blood_pca)){stop("Expecting argument for blood_pca")}
  
  predict_vals = predict(gtex_model, data.frame(predict_pca_new))
  
  if(scale == TRUE){
  predict_vals = scale(predict_vals)}
  
}

