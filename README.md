# CSsingle

`CSsingle` (Cross-Source SINGLE cell deconvolution) is an R package for the accurate decomposition of bulk transcriptomic data into a set of pre-defined cell types using the scRNA-seq or flow-sorting reference. The key strengths of `CSsingle` include (i) it introduces cell size coefficients in deconvolution, which properly correcting for bias arising from both “technical” library size and “biological” cell size, and (ii) it effectively handles technical and biological variations between individual bulk mixtures and the signature matrix. 

<p align="center">
<img src="./CSsingle_framework.jpg" width="700">
</p>

## Installation

```r
# in R
# install devtools if necessary
if (!"devtools" %in% rownames(installed.packages())) {
  install.packages('devtools')
}

# install the CSsingle package
if (!"CSsingle" %in% rownames(installed.packages())) {
 devtools::install_github('wenjshen/CSsingle')
}
# load
library(CSsingle)
```

## Tutorial

Vignette: [HTML Vignette](http://htmlpreview.github.io/?https://github.com/wenjshen/CSsingle/master/vignettes/CSsingle_vignette.html)
