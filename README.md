# Welcome to BrainGENIE!

``` 
Developed using
______ 

R version 3.5.2 (2018-12-20)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS Mojave 10.14.5
______
```

## Please install the following R packages with the commands:
```
install.packages("rtracklayer");
install.packages("plyr");
install.packages("data.table")
```

### Number of trained gene-level prediction models in GTEx verison 7 (elastic net regression models):
`(As of June 24, 2019)`
| Brain region | # of genes |
| -----------  | ---------- |
| Amygdala     | 6,202      |
| ACC          | 9,765      |
| Caudate      | 9,036      |
| Cere. Hem.   | 10,736     |
| Cerebellum   | 9,888      |
| Cortex       | 7,937      |
| FCx (BA9)    | 5,318      |
| Hippocampus  | 5,855      |
| NAcc         | 8,205      |
| Putamen      | 6,326      |
| Subst. Nigra | 5,840      |


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


