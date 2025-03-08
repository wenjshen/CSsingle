# CSsingle

CSsingle is a unified tool designed for the accurate and robust decomposition of bulk and spatial transcriptomic data into a set of predefined cell types using the scRNA-seq or flow-sorting reference.

<p align="center">
<img src="./CSsingle_framework.jpg" width="700">
</p>

## Installation

```r
#in R
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

## Help

Vignette: [HTML Vignette](https://github.com/wenjshen/CSsingle/vignettes/CSsingle_vignette.html)
