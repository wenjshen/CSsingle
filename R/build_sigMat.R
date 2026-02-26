#' This function is to build a signature matrix by selecting the top N markers for each cell type
#'
#' @param diffexpMat numeric matrix, a differentially expressed matrix
#' @param ctDEGs data frame, matrix containing a ranked list of differentially expressed genes, and associated statistics
#' @param num.top the number of top markers
#' @param sort.by rank the markers in ascending order by their p-values (p.value, default) or in descending order by their fold changes (fold.change)
#' @param rm.dupli whether to exclude markers shared between two or more cell types, default is TRUE
#' @param filepath directory where the results will be saved, default is NULL
#'
#' @return a list with two attributes:
#' \item{topSigMat}{numeric matrix, signature matrix}
#' \item{topMarkers}{a list of marker genes for each cell type}
#'
#' @author Wenjun Shen
#' @export
#'
build_sigMat <- function(diffexpMat, ctDEGs, num.top, sort.by = "p.value", rm.dupli = TRUE, filepath = NULL){
  rownames(diffexpMat) <- toupper(rownames(diffexpMat))
  # select the top marker genes for each cell type
  cty <- colnames(diffexpMat)
  ctDEGs.cty <- list()
  for (i in cty){
    idx <- ctDEGs[, "cluster"] == i
    cds <- ctDEGs[idx, ]
    if (sort.by == "fold.change"){
      gene.i <- toupper(cds[order(cds[, "avg_log2FC"], decreasing = TRUE), "gene"])
    }else if(sort.by == "p.value"){
      gene.i <- toupper(cds[, "gene"])
    }
    num.top <- as.numeric(num.top)
    ctDEGs.cty[[i]] <- gene.i[1:min(num.top,length(gene.i))]
  }
  if (rm.dupli){
    rm_ctDEGs <- list()
    for (i in 1:length(cty)){
      markers.i <- ctDEGs.cty[[i]]
      idx <- Rfast::rowMaxs(diffexpMat[markers.i, ]) == i
      rm_ctDEGs[[i]] <- markers.i[idx == F]
    }
    rm_ctDEGs <- unique(unlist(rm_ctDEGs))
    topMarkers <- list()
    for (i in cty){
      topMarkers[[i]] <- setdiff(ctDEGs.cty[[i]], rm_ctDEGs)
      if (length(topMarkers[[i]]) < 3){
        topMarkers[[i]] <- ctDEGs.cty[[i]]
      }
    }
  }else topMarkers <- ctDEGs.cty

  num.genes <- sapply(topMarkers, length)
  if (min(num.genes) < 3){
    id <- num.genes < 3
    rm.cty <- cty[id]
    keep.cty <- cty[!id]
    print(paste0("The following cell types have markers less than 3:"))
    print(rm.cty)
  }else keep.cty <- cty

  topMarkers <- topMarkers[keep.cty]
  topSigMat <- diffexpMat[unique(unlist(topMarkers)), keep.cty]
  res <- list("sigMat" = topSigMat, "markers" = topMarkers)
  if (!is.null(filepath)) saveRDS(res, file.path(filepath, "sigMat.rds"))
  return(res)
}
