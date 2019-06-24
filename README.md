# BrainGENIE

``` 
Developed using
______ 

R version 3.5.2 (2018-12-20)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS Mojave 10.14.5
______
```

# Please install the following R packages with the commands:
```
install.packages("rtracklayer");
install.packages("plyr");
install.packages("data.table")
```

### Note: 
`Gene IDs must be in Ensembl format`

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


