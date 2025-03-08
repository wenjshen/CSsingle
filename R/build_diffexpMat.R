#' This function is to build a differentially expressed matrix
#'
#' @param sc.eset SeuratObject for scRNA-seq data
#' @param cellType character, name of metadata column indicating cell type
#' @param only.pos only return positive gene markers (TRUE by default)
#' @param test.use which statistical test to use
#' @param logfc.threshold limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells
#' @param min.pct only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations
#' @param return.thresh only return gene markers that have a p-value < return.thresh
#' @param filepath directory where the results will be saved, default is NULL
#'
#' @return a list with three attributes:
#' \item{ctDEGs}{data frame, matrix containing a ranked list of differentially expressed genes, and associated statistics}
#' \item{diffexpMat}{numeric matrix, a differential expressed matrix}
#' \item{designMat}{numeric matrix, a gene by cell-type matrix (genome-wide)}
#'
#' @import Seurat
#' @author Wenjun Shen
#' @export
#'
build_diffexpMat <- function(sc.eset, cellType, only.pos = TRUE, test.use = "wilcox", logfc.threshold = 0.25, min.pct = 0.05, return.thresh = 0.01, filepath = NULL, verbose = TRUE, ...){
  # Find marker genes for each cell type
  meta <- sc.eset@meta.data
  cty <- as.character(meta[[cellType]])
  idx <- is.na(cty)
  if (sum(idx)>0){
    sc.eset <- sc.eset[, which(!idx)]
  }
  cty <- as.character(sc.eset@meta.data[[cellType]])
  Idents(object = sc.eset) <- factor(cty)
  # build design matrix
  designMat <- sapply(unique(cty), function(ct) {
    y = as.matrix(sc.eset[["RNA"]]$counts[, cty %in% ct,
                                          drop = FALSE])
    1e4 * rowSums(y)/sum(y)
  })
  # build differentially expressed matrix
  if(verbose){message("Constructing differentially expressed Matrix...")}
  sc.eset <- NormalizeData(sc.eset)
  ctDEGs <- FindAllMarkers(sc.eset, only.pos = only.pos, min.pct = min.pct, logfc.threshold = logfc.threshold, test.use = test.use, return.thresh = return.thresh, ...)
  diffexpMat <- designMat[unique(ctDEGs$gene),]
  res <- list(ctDEGs=ctDEGs, diffexpMat=diffexpMat, designMat=designMat)
  if (!is.null(filepath)) saveRDS(res, file.path(filepath, "diffexpMat.rds"))
  return(res)
}
