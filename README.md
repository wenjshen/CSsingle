# CSsingle

CSsingle is a unified tool designed for the accurate and robust decomposition of bulk and spatial transcriptomic data into a set of predefined cell types using the scRNA-seq or flow-sorting reference.

CSsingle addresses core challenges in cellular deconvolution through three integrated capabilities: (i) adaptive weighting that leverages marker-gene sectional linearity to achieve robust cross-source harmonization; (ii) integrated cell size correction that accounts for cell-type-specific RNA content differences using either ERCC spike-ins or a novel computational estimator; and (iii) spatially resolved deconvolution enabling fine-grained dissection of tissue microenvironments at single-spot resolution.

<p align="center">
<img src="./CSsingle_framework.jpg" width="700">
</p>

## Installation

```r
# install devtools if necessary
if (!"devtools" %in% rownames(installed.packages())) {
  install.packages('devtools')
}

# install CSsingle package
if (!"CSsingle" %in% rownames(installed.packages())) {
 devtools::install_github('wenjshen/CSsingle')
}
# load package
library(CSsingle)
```

## Help

Vignette: [HTML Vignette](https://wenjshen.github.io/CSsingle)

## How to cite CSsingle
Please cite the following publication:
Wenjun Shen, Yunfei Hu, Yuanfang Lei, Hau-San Wong, Cheng Liu, Si Wu, Xin Maizie Zhou, CSsingle: a unified tool for robust decomposition of bulk and spatial transcriptomic data across diverse single-cell references, Nucleic Acids Research, Volume 54, Issue 8, 8 May 2026, gkag410, https://doi.org/10.1093/nar/gkag410
