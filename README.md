# iDEA

Integrative Differential expression and gene set Enrichment Analysis Using Summary Statistics for single cell RNAseq studies 

System requirements
--------------------
This package is supported for Windows, macOS and Linux. The package has been tested on the following systems:
- Windows: Home (1903)
- MAC: OSX (10.14.1)
- Linux: Ubuntu (16.04.6)

Dependencies
-------------
- Rcpp (>= 1.0.0) 
- doParallel (>=1.0.15) 
- doSNOW (>=1.0.15) 
- foreach (>=1.4.4) 
- parallel (>=3.6.1) 


Installation
------------
``` r
# install devtools if necessary
install.packages('devtools')

# install the iDEA package
devtools::install_github('xzhoulab/iDEA')

# load package
library(iDEA)
```

How to cite `iDEA`
-------------------
Shiquan Sun, Ying Ma, Xuequn Shang, Evan T. Keller, Mengjie Chen and Xiang Zhou#. *Integrative differential expression and gene set enrichment analysis using summary statistics for single cell RNAseq studies*, 2019. 

How to use `iDEA`
-------------------
Details in [here](https://xzhoulab.github.io/iDEA/)
