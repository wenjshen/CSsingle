#' This function is to identify cell types that are enriched in Spatial transcriptomic (ST) spots
#'
#' @param sigMat numeric matrix, signature matrix
#' @param mixture numeric matrix, spatial trancriptomic gene expression matrix
#' @param markers list of gene markers for each cell type
#' @param enrich.thres default is 0
#' @param filepath directory where the results will be saved, default is NULL
#'
#' @return a list of cell type enriched in ST spots
#'
#' @author Wenjun Shen
#' @export
#'
enrichmentMarkers <- function(sigMat, mixture, markers, enrich.thres = 0, filepath = NULL){
  rownames(sigMat) <- toupper(rownames(sigMat))
  rownames(mixture) <- toupper(rownames(mixture))
  com.genes <- intersect(rownames(sigMat), rownames(mixture))
  X <- as.matrix(sigMat[com.genes,])
  Y <- as.matrix(mixture[com.genes, , drop=F])
  mean.Y <- rowMeans(Y)
  geneFoldchange <- log2(Y+1) - log2(mean.Y+1)
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
        zscore1[i,j] <- (Sm - u) * m^(1 / 2) / sigma
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
  if (!is.null(filepath)){
    #save results
    save(zscore1, zscore2, enrich.ct, file = file.path(filepath, "enrichment-Results.rda"))
  }

  return(enrich.ct)
}
