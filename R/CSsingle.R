#' This function is to estimate the cell-type proportions from bulk gene expression data
#'
#' @param sigMat numeric matrix, signature matrix
#' @param mixture bulk gene expression matrix
#' @param markers list of gene markers for each cell type
#' @param cellSize Users can specify a numeric vector of cell sizes, default is NULL, which means that we assume the absolute amount of total mRNA is similar across different cell types in deconvolution
#' @param dampened whether to introduce a upper bound constant that limits the maximum value any weight can take on, default is FALSE
#' @param filepath directory where the results will be saved, default is NULL
#'
#' @return a numeric matrix of cell type proportions estimated by CSsingle (add up to 1)
#' @export
#' @import Seurat Biobase
#' @importFrom nnls nnls
#'
#' @author Wenjun Shen
#'
CSsingle <- function(sigMat, mixture, markers, cellSize = NULL, dampened = FALSE, filepath = NULL){
  # column normalization
  mixture <- apply(mixture, 2, function(x) 1e4*x/sum(x))

  rownames(sigMat) <- toupper(rownames(sigMat))
  rownames(mixture) <- toupper(rownames(mixture))
  com.genes <- intersect(rownames(sigMat),rownames(mixture))
  X <- as.matrix(sigMat[com.genes,])
  Y <- as.matrix(mixture[com.genes, , drop=F])
  Est.prop <- c()
  for (j in 1 : ncol(Y)){
    y <- Y[, j]
    idx.zero <- y == 0
    y <- y[!idx.zero]
    x <- X[!idx.zero,]
    Est <- CSsingle_basis(x, y, markers, cellSize, dampened)
    Est.prop <- rbind(Est.prop, Est)
  }
  rownames(Est.prop) <- colnames(Y)

  if (!is.null(filepath)){
    #save results
    saveRDS(Est.prop, file.path(filepath, "CSsingle-Results.rds"))
    write.table(Est.prop, file = file.path(filepath, "CSsingle-Results.txt"), sep="\t", row.names=T, col.names=T, quote=F)
  }
  return(Est.prop)
}
