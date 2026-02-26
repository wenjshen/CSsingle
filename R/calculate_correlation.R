calculate_correlation <- function(est.mixture, mixture, genes, method = "spearman") {
  corr <- sapply(1:ncol(mixture), function(j) {
    if (any(is.na(est.mixture[, j]))) {
      return(-1)
    } else {
      return(cor(est.mixture[genes, j], mixture[genes, j],
                 method = method))
    }
  })
  mean(corr, na.rm = TRUE)
}
