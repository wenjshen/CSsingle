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
  exprs <- as.matrix(GetAssayData(sc.eset, assay = "RNA", layer = "counts"))
  cell_sparsity <- apply(exprs != 0, MARGIN = 2, sum)
  keep_cells <- colnames(exprs)[cell_sparsity > cell.thres]
  feature_sparsity <- apply(exprs != 0, MARGIN = 1, sum)
  keep_features <- rownames(exprs)[feature_sparsity > feature.thres]
  DefaultAssay(sc.eset) <- "RNA"
  sc.eset[["SCT"]] <- NULL
  sc.eset <- subset(sc.eset, features = keep_features, cells = keep_cells)
  return(sc.eset)
}
