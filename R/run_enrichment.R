#' This function runs the complete enrichment analysis pipeline. It identifies cell types that are enriched in spatial transcriptomic (ST) spots.
#'
#' @param diffexpMat numeric matrix, a differentially expressed matrix
#' @param mixture numeric matrix, spatial transcriptomic gene expression matrix
#' @param ctDEGs data frame, matrix containing a ranked list of differentially expressed genes
#' @param num.top number of top markers to use for enrichment analysis, default is 50
#' @param sort.by rank the markers in ascending order by their p-values (p.value, default)
#'        or in descending order by their fold changes (fold.change)
#' @param rm.dupli whether to exclude markers shared between two or more cell types, default is FALSE
#' @param enrich.thres threshold for relative enrichment z-score, default is 0
#'
#' @return a list of enriched cell types for each spot (combined from both methods)
#'
#' @export
#'
#'
run_enrichment <- function(diffexpMat, mixture, ctDEGs, num.top = 50,
                          sort.by = "p.value", rm.dupli = FALSE,
                          enrich.thres = 0
                          ) {

  # Validate inputs
  if (is.null(diffexpMat) || is.null(mixture) || is.null(ctDEGs)) {
    stop("diffexpMat, mixture, and ctDEGs must all be provided")
  }

  # Build signature matrix with top markers
  res <- build_sigMat(diffexpMat, ctDEGs, num.top, sort.by, rm.dupli)

  # Run enrichment analysis
  enrich.ct <- enrichmentMarker(res$sigMat, mixture, res$markers, enrich.thres)

  return(enrich.ct)
}
