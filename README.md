## ParsimonyC3: prioritize LR pairs that navigate one cell type to others. 

<p align="center">
  <img width="300"  src="logo.png">
</p>
    
## Installation instructions
You can install the package via `devtools::install_github()` function in R
```R
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
devtools::install_github('sheenaseven/ParsimonyC3', dependencies = TRUE)
```
## Usage instructions
```R
library(ParsimonyC3)
cci <- DataTransformCPD('cellphonedb/ccm/out',
                       'cellphonedb/ep/out')
cci <- AddWeights(cci,features=c('means_e','FC_ec'),w = c(1,1))
lr_combination <- FindComb(cci,n=2)
result <- ParsimonyC3(cci,lr_combination)
LollipopPlot(result,top = 10)
NetworkPlot(cci,lr_select='PECAM1_CD38|integrin_a4b7_complex_FN1')
```

## Contact
Please contact us:  
Yingxi Yang: yxyang@ust.hk
