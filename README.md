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
`Rscript install.R`

## Download normalized residual gene expression data (paired blood-brain) for GTEx version 8 release:
https://zenodo.org/record/6350240/files/normalized_expression_dat_gtexv8.tar.gz?download=1

### Number of significant gene-level prediction made by BrainGENIE:
`(As of March 6, 2020, prefiltered for cross-validation R^2 ≥ 0.01, cross-valdidation FDRp-value < 0.05, and test-set R^2 ≥ 0.01)`

 |                        Tissue    | # genes|  R^2 |
 | -------------------------------  | ------ | ---- |
 |                        Amygdala  | 3,132  | 0.13 |
 |  Anterior_cingulate_cortex_BA24  | 1,934  | 0.11 |
 |           Caudate_basal_ganglia  | 4,062  | 0.09 |
 |           Cerebellar_Hemisphere  | 4,532  | 0.10 |
 |                      Cerebellum  | 5,691  | 0.10 |
 |                          Cortex  | 4,220  | 0.10 |
 |              Frontal_Cortex_BA9  | 5,554  | 0.10 |
 |                     Hippocampus  | 2,947  | 0.10 |
 |                    Hypothalamus  | 1,565  | 0.12 |
 | Nucleus_accumbens_basal_ganglia  | 3,611  | 0.08 |
 |           Putamen_basal_ganglia  | 6,835  | 0.11 |
 |                Substantia_nigra  | 1,355  | 0.15 |


```

#### Normalization of transcriptome-wide gene expression data from GTEx v8:
1. Kept genes with ≥5 read counts ≥10 samples &  RPKM ≥ 0.1 in ≥10 samples
2. Quantile normalization of RPKM values
3. Inverse normal transformation
4. Regressed covariates out of expression data, used residuals for model training cross-validation
 - Brain: Age, sex, top three genotype-based principal components, PCR method, platform, ischemic time, RIN, and death record by Hardy  scale.
 - Blood: Age, sex, top three genotype-based principal components, top three CIBERSORT principal components, PCR method, platform, ischemic time, RIN, and death record by Hardy  scale.

 
