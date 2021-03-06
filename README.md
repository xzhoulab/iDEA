# iDEA

Integrative Differential expression and gene set Enrichment Analysis using summary statistics for single cell RNAseq studies 

System requirements
--------------------
This package is supported for Windows 10, MAC and Linux. The package has been tested on the following systems:
- Windows 10: Home (1903)
- MAC: OSX (10.14.1)
- Linux: Ubuntu (16.04.6)

Dependencies (r packages)
-------------
- Rcpp (>= 1.0.0) 
- RcppArmadillo (>=0.9.400.2)
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
Ying Ma, Shiquan Sun, Evan T. Keller, Mengjie Chen and Xiang Zhou#. *Integrative differential expression and gene set enrichment analysis using summary statistics for single cell RNAseq studies*, Nature Communications 2020.(https://www.nature.com/articles/s41467-020-15298-6)

How to use `iDEA`
-------------------
Details in [Tutorial](https://xzhoulab.github.io/iDEA/)
