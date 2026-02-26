#' This function is to estimate the cell-type proportions from bulk and spatial transcriptomic (ST) data
#'
#' @param sigMat numeric matrix, signature matrix
#' @param mixture numeric matrix, bulk or spatial trancriptomic gene expression matrix (raw counts for RNA-seq or ST, or RMA normalized intensity values in linear scale for Microarray)
#' @param markers list of gene markers for each cell type
#' @param cellSize users can specify a numeric vector of cell sizes, default is NULL, which assumes that the absolute amount of total mRNA is similar across different cell types in deconvolution. Set to NULL when ERCC spike-ins are unavailable from scRNA-seq reference or for spatial transcriptomics (ST) analysis.
#' @param enrich.ct a list of cell type enriched in ST spots (optional)
#' @param dampened whether to introduce a upper bound constant that limits the maximum value any weight can take on, default is FALSE
#' @param beta numeric vector of [0, 0.5, 1]
#' @param filepath directory where the results will be saved, default is NULL
#'
#' @return a numeric matrix of cell type proportions estimated by CSsingle (add up to 1)
#' @export
#' @import Seurat
#' @importFrom nnls nnls
#'
#' @author Wenjun Shen
#'
CSsingle <- function(sigMat, mixture, markers, cellSize = NULL, enrich.ct = NULL, dampened = FALSE, beta = c(0, 0.5, 1), filepath = NULL){
  # column normalization
  mixture <- apply(mixture, 2, function(x) 1e4*x/sum(x))
  rownames(sigMat) <- toupper(rownames(sigMat))
  rownames(mixture) <- toupper(rownames(mixture))
  com.genes <- intersect(rownames(sigMat),rownames(mixture))
  X <- as.matrix(sigMat[com.genes,])
  Y <- as.matrix(mixture[com.genes, , drop=F])

  if (!is.null(cellSize) & length(unique(cellSize)) > 1){
    scale.factor <- seq(1e4, 1e7, length=1000)
    diff_tmp <- c()
    for (sf in scale.factor){
      X_tmp <- t(t(X) * cellSize[colnames(X)]/sf)
      t_tmp <- rowMeans(X_tmp)
      y_tmp <- rowMeans(Y)
      diff_tmp <- c(diff_tmp, abs(sum(t_tmp) - sum(y_tmp)))
    }
    opt.sf <- scale.factor[which.min(diff_tmp)]
    cellSize <- cellSize/opt.sf
  }

  correlation <- c()
  Est.propS <- list()
  for (b in 1:length(beta)){
    Est.prop <- c()
    for (j in colnames(Y)){
      y <- Y[, j]
      enrich <- enrich.ct[[j]]
      idx.zero <- y == 0
      if (sum(!idx.zero) > 5){
        y <- y[!idx.zero]
        x <- X[!idx.zero,]
        Est <- CSsingle_basis(x, y, markers, cellSize, enrich, dampened, beta[b])
      }else{
        Est <- rep(NA, ncol(X))
      }
      Est.prop <- rbind(Est.prop, Est)
    }
    Est.propS[[b]] <- Est.prop

    if (!is.null(cellSize) & length(unique(cellSize)) > 1){
      Est.Y <- (t(t(X) * cellSize[colnames(X)])) %*% t(Est.prop)
    }else{
      Est.Y <- X %*% t(Est.prop)
    }
    corr <- c()
    for (k in 1:ncol(Est.Y)){
      if (any(is.na(Est.Y[, k]))){
        corr[k] <- -1
      }else{
        corr[k] <- cor(Est.Y[, k], Y[, k], method = "spearman")
      }
    }
    correlation[b] <- round(mean(corr, na.rm = T), 3)
  }
  Est.propS <- Est.propS[order(correlation, decreasing = T)]
  Est.prop <- Est.propS[[1]]
  rownames(Est.prop) <- colnames(Y)

  if (!is.null(filepath)){
    # save results
    saveRDS(Est.prop, file.path(filepath, "CSsingle-Results.rds"))
    write.table(Est.prop, file = file.path(filepath, "CSsingle-Results.txt"), sep="\t", row.names=T, col.names=T, quote=F)
  }
  return(Est.prop)
}
