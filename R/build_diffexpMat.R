#' This function is to build a differentially expressed matrix
#'
#' @param sc.eset SeuratObject for scRNA-seq data
#' @param cellType character, name of metadata column indicating cell type
#' @param normalization.mode character, normalization strategy to use. Options:
#'   - "pseudobulk": Sum-then-normalize (default). Creates pseudobulk profiles by summing counts across cells, then normalizing.
#'   - "cell-average": Normalize-then-average. Normalizes each cell individually, then averages across cells.
#' @param only.pos only return positive gene markers (TRUE by default)
#' @param test.use which statistical test to use
#' @param logfc.threshold limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells
#' @param min.pct only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations
#' @param return.thresh only return gene markers that have a p-value < return.thresh
#' @param filepath directory where the results will be saved, default is NULL
#' @param verbose logical, whether to print progress messages
#' @param ... additional arguments passed to FindAllMarkers
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
build_diffexpMat <- function(sc.eset, cellType, normalization.mode = c("pseudobulk", "cell-average"),
                             only.pos = TRUE, test.use = "wilcox", logfc.threshold = 0.25,
                             min.pct = 0.05, return.thresh = 0.01, filepath = NULL,
                             verbose = TRUE, ...){

  # Validate normalization mode
  normalization.mode <- match.arg(normalization.mode)

  # Process cell type information
  meta <- sc.eset@meta.data
  cty <- as.character(meta[[cellType]])
  idx <- is.na(cty)

  if (sum(idx) > 0){
    sc.eset <- sc.eset[, which(!idx)]
  }
  cty <- as.character(sc.eset@meta.data[[cellType]])
  Idents(object = sc.eset) <- factor(cty)

  if(normalization.mode == "pseudobulk"){
    # Mode 1: Pseudobulk Normalization (Sum-then-normalize)
    designMat <- sapply(unique(cty), function(ct) {
      # Subset cells of current cell type
      y <- as.matrix(GetAssayData(sc.eset, assay = "RNA", layer = "counts")[
        , cty %in% ct, drop = FALSE])

      # Sum counts across all cells, then normalize to 10,000
      1e4 * rowSums(y) / sum(y)
    })

  } else if(normalization.mode == "cell-average"){
    # Mode 2: Cell-Average Normalization (Normalize-then-average)
    designMat <- sapply(unique(cty), function(ct) {
      # Subset cells of current cell type
      y <- as.matrix(GetAssayData(sc.eset, assay = "RNA", layer = "counts")[
        , cty %in% ct, drop = FALSE])

      # Normalize each cell to 10,000, then average across cells
      # Handle division by zero for cells with no counts
      col_sums <- colSums(y)
      valid_cols <- col_sums > 0
      if(sum(valid_cols) < ncol(y)) {
        # Some cells have zero counts - only use valid ones
        y <- y[, valid_cols, drop = FALSE]
        col_sums <- col_sums[valid_cols]
      }

      # Normalize each cell and take row means
      rowMeans(t(1e4 * t(y) / col_sums))
    })
  }

  # Set row and column names
  rownames(designMat) <- rownames(sc.eset)
  colnames(designMat) <- unique(cty)

  # Build differentially expressed matrix
  if(verbose){message("Constructing differentially expressed matrix...")}

  # Normalize data for DEG analysis
  sc.eset <- NormalizeData(sc.eset)

  # Find marker genes
  ctDEGs <- FindAllMarkers(sc.eset, only.pos = only.pos, min.pct = min.pct,
                           logfc.threshold = logfc.threshold, test.use = test.use,
                           return.thresh = return.thresh, ...)

  # Create differential expression matrix using only DEGs
  diffexpMat <- designMat[unique(ctDEGs$gene), , drop = FALSE]
  # Prepare result list
  res <- list(
    ctDEGs = ctDEGs,
    diffexpMat = diffexpMat,
    designMat = designMat
  )

  # Save results if filepath provided
  if (!is.null(filepath)) {
    saveRDS(res, file.path(filepath, "diffexpMat.rds"))
    if(verbose){
      message(sprintf("Results saved to %s", file.path(filepath, "diffexpMat.rds")))
    }
  }

  return(res)
}
