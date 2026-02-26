#' This function is to exclude cell types that are not interested from scRNA-seq data, remove sparse features and cells and trim bulk and scRNA-seq data to have the same features
#'
#' @param bulk.eset ExpressionSet object for bulk or spatial trancriptomic gene expression data
#' @param sc.eset SeuratObject for scRNA-seq data
#' @param cellType character, name of metadata column indicating cell type
#' @param select.ct character vector of cell types included, default is NULL, which means we include all cell types in the scRNA-seq data
#' @param cell.thres include features detected in at least this many cells
#' @param feature.thres include cells where at least this many features are detected
#'
#' @return a list with two attributes:
#' \item{bulk.eset}{bulk.eset}
#' \item{sc.eset}{SeuratObject for scRNA-seq data}
#'
#' @import Seurat
#' @author Wenjun Shen
#' @export
#'
preprocessData <- function(bulk.eset, sc.eset, cellType, select.ct = NULL, cell.thres = 20, feature.thres = 0){
  # exclude cell types that are not interested from scRNA-seq data
  if (!is.null(select.ct)){
  idx = sc.eset@meta.data[[cellType]] %in% select.ct
  sc.eset <- sc.eset[, idx]
  }
  # remove sparse features and cells
  sc.eset = qualityControl(sc.eset, cell.thres=cell.thres, feature.thres=feature.thres)
  sc.eset@meta.data[[cellType]] <- gsub("[.]", "_", make.names(sc.eset@meta.data[[cellType]]))
  # trim bulk and scRNA-seq data to have the same features
  com.features = intersect(toupper(rownames(bulk.eset)), toupper(rownames(sc.eset)))
  bulk.eset = bulk.eset[match(com.features, toupper(rownames(bulk.eset))), ]
  sc.eset = sc.eset[match(com.features, toupper(rownames(sc.eset))), ]
  return(list(bulk.eset = bulk.eset, sc.eset = sc.eset))
}
