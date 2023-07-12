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

## Extract the data from the compressed tar archive with the following command in terminal:
`tar -xzvf normalized_expression_dat_gtexv8.tar.gz`

### Vignette for getting started:
https://github.com/hessJ/BrainGENIE/wiki/How-to-run-BrainGENIE

```

#### Normalization of transcriptome-wide gene expression data from GTEx v8:
1. Kept genes with ≥5 read counts ≥10 samples &  RPKM ≥ 0.1 in ≥10 samples
2. Quantile normalization of RPKM values
3. Inverse normal transformation
4. Regressed covariates out of expression data, used residuals for model training cross-validation
 - Brain: Age, sex, top three genotype-based principal components, PCR method, platform, ischemic time, RIN, and death record by Hardy  scale.
 - Blood: Age, sex, top three genotype-based principal components, top three CIBERSORT principal components, PCR method, platform, ischemic time, RIN, and death record by Hardy  scale.

 
