#' This function is to identify cell types that are enriched in Spatial transcriptomic (ST) spots
#'
#' @param sigMat numeric matrix, signature matrix
#' @param mixture numeric matrix, spatial trancriptomic gene expression matrix
#' @param markers list of gene markers for each cell type
#' @param enrich.thres threshold for enrichment, default is 0
#'
#' @return a list of cell type enriched in ST spots
#'
#' @author Wenjun Shen
#' @export
#'
enrichmentMarker <- function(sigMat, mixture, markers, enrich.thres = 0){
  rownames(sigMat) <- toupper(rownames(sigMat))
  rownames(mixture) <- toupper(rownames(mixture))
  com.genes <- intersect(rownames(sigMat), rownames(mixture))

  if (length(com.genes) == 0) {
   stop("No common genes between signature matrix and mixture")
 }

  X <- as.matrix(sigMat[com.genes,])
  Y <- as.matrix(mixture[com.genes, , drop = F])
  mean.Y <- rowMeans(Y)
  geneFoldchange <- log2(Y+1) - log2(mean.Y+1)
  prevalent_cells = identify_prevalent_cells(Y, markers)
  if (length(prevalent_cells) > 0){
    Q1.Y <- apply(Y, 1, quantile, probs = 0.25, na.rm = TRUE)
    geneFoldchange_q1 <- log2(Y+1) - log2(Q1.Y+1)
    for (ct in prevalent_cells) {
      ct_markers <- intersect(rownames(Y), markers[[ct]])
      if (length(ct_markers) > 0) {
        geneFoldchange[ct_markers, ] <- geneFoldchange_q1[ct_markers, ]
      }
    }
  }
  mean.fc <- apply(geneFoldchange, 2, mean)
  sd.fc <- apply(geneFoldchange, 2, stats::sd)

  zscore1 <- matrix(NA, nrow = ncol(Y), ncol = ncol(X))
  for (i in 1 : ncol(Y)){
    u <- mean.fc[i]
    sigma <- sd.fc[i]
    y <- Y[, i]
    expressed <- names(y)[y != 0]
    for (j in 1 : ncol(X)){
      ct <- colnames(X)[j]
      id <- names(y) %in% toupper(markers[[ct]])
      markers.j <- names(y)[id]
      m <- length(intersect(markers.j, expressed))
      if (m > 1){
        Sm <- mean(geneFoldchange[markers.j, i])
        zscore1[i,j] <- (Sm - u) * sqrt(m) / sigma
      }
    }
  }
  rownames(zscore1) <- colnames(Y)
  colnames(zscore1) <- colnames(X)
  enrich.ct1 <- apply(zscore1, 1, function(x) colnames(X)[which(x > enrich.thres)], simplify = FALSE)
  enrich.ct2 <- apply(zscore1, 1, function(x) colnames(X)[order(-x)[1]], simplify = FALSE)
  zscore2 <- t(apply(zscore1, 1, function(x) (x-mean(x, na.rm = T))/stats::sd(x, na.rm = T)))

  enrich.ct3 <- apply(zscore2, 1, function(x) colnames(X)[which(x > enrich.thres)], simplify = FALSE)
  enrich.ct <- list()
  for (i in names(enrich.ct1)){
    enrich.ct[[i]] <- union(union(enrich.ct1[[i]], enrich.ct2[[i]]), enrich.ct3[[i]])
  }
  return(enrich.ct)
}
