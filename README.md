# Welcome to BrainGENIE!

``` 
Developed using
______ 

R version 3.5.2 (2018-12-20)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS Mojave 10.14.5
______
```

## Clone this repository with:
`git clone https://github.com/hessJ/BrainGENIE`

## Install dependencies for BrainGENIE with the command in console:
```
Rscript install.R
```

## To run braingenie_pca.R, please download the following files from GTEx version 7 release:
```
wget https://storage.googleapis.com/gtex_analysis_v7/rna_seq_data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct.gz
wget https://storage.googleapis.com/gtex_analysis_v7/annotations/GTEx_v7_Annotations_SampleAttributesDS.txt
wget https://storage.googleapis.com/gtex_analysis_v7/annotations/GTEx_v7_Annotations_SubjectPhenotypesDS.txt
wget https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/GTEx_Analysis_v7_eQTL_covariates.tar.gz
```

### Number of trained gene-level prediction models in GTEx verison 7 (elastic net regression models):
`(As of July 5, 2019, prefiltered for Pearson's r ≥ 0.1, p-value < 0.05)`
 |                        Tissue    |  Genes  |  Cor |
 |--------------------------------  | ------- | ---- |
 |                        Amygdala  | 14,666  | 0.47 |
 |  Anterior_cingulate_cortex_BA24  | 13,496  | 0.43 |
 |           Caudate_basal_ganglia  | 11,538  | 0.34 |
 |           Cerebellar_Hemisphere  | 11,350  | 0.37 |
 |                      Cerebellum  | 13,777  | 0.35 |
 |                          Cortex  | 12,386  | 0.36 |
 |              Frontal_Cortex_BA9  | 12,135  | 0.36 |
 |                     Hippocampus  | 11,082  | 0.36 |
 |                    Hypothalamus  | 11,739  | 0.38 |
 | Nucleus_accumbens_basal_ganglia  | 10,842  | 0.34 |
 |           Putamen_basal_ganglia  | 11,010  | 0.36 |
 |                Substantia_nigra  | 13,149  | 0.44 |

### Note: 
`Gene IDs must be in Ensembl gene ID format # ENSG[XXXXXXXXXX].[X]`

### Simple example:
```
bg_models = load_models(path="~/Documents/braingenie/trained_models/gtex_v7/")
load_cv_performance(path="~/Documents/braingenie/trained_models/gtex_v7/")
# eDat (data frame of blood transcriptome profiles; rows = subjects, columns = genes)

bg.output=list() # captures output from BrainGENIE into a list object [1 = imputed transcriptome, 2 = problematic genes]
for(x in 1:length(bg_models)){
bg.output[[x]] = predict_brain_gxp(mod = iter[[x]], target = eDat, index = x, missing_prop = 0.5)
}

```

#### Normalization of transcriptome-wide gene expression data from GTEx v7:
1. Kept gened with > 5 read counts in a minimum of 10 samples
2. RPKM normalization using edgeR 
3. Quantile normalization across samples using limma
4. Inverse normal transformation
5. Used residuals after regressing out age, sex, and top three genotype-based principal components

