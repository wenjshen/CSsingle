#' This function is to remove sparse features and cells from scRNA-seq data
#'
#' @param sc.eset SeuratObject for scRNA-seq data
#' @param cell.thres include features detected in at least this many cells
#' @param feature.thres include cells where at least this many features are detected
#'
#' @return SeuratObject for scRNA-seq data after quality control
#'
#' @import Seurat
#' @author Wenjun Shen
#'
qualityControl <- function(sc.eset, cell.thres=20, feature.thres=0){
  sc.eset <- sc.eset[which(!duplicated(toupper(rownames(sc.eset)))),]
  exprs <- as.matrix(sc.eset[["RNA"]]$counts)
  cell_sparsity <- apply(exprs != 0, MARGIN = 2, sum)
  keep_cells <- cell_sparsity > cell.thres
  feature_sparsity <- apply(exprs != 0, MARGIN = 1, sum)
  keep_features <- feature_sparsity > feature.thres
  sc.eset <- sc.eset[keep_features, keep_cells]
  return(sc.eset)
}
