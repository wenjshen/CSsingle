#' This function is to build a signature matrix by minimizing the condition number
#'
#' @param diffexpMat numeric matrix, a differentially expressed matrix
#' @param ctDEGs data frame, matrix containing a ranked list of differentially expressed genes, and associated statistics
#' @param increment numeric vector, create multiple candidate matrices by specifying the number of top markers from each cell type
#' @param sort.by rank the markers in ascending order by their p-values (p.value, default) or in descending order by their fold changes (fold.change)
#' @param rm.dupli whether to exclude markers shared between two or more cell types, default is TRUE
#' @param filepath directory where the results will be saved, default is NULL
#'
#' @return a list with four attributes:
#' \item{sigMat}{numeric matrix, signature matrix}
#' \item{markers}{a list of marker genes for each cell type}
#' \item{min.cond.num}{the lowest condiction number}
#' \item{num.top}{the number of top markers with the lowest condiction number}
#'
#' @import Seurat
#' @author Wenjun Shen
#' @export
#'
build_sigMat_condNum <- function(diffexpMat, ctDEGs, increment = seq(50, 300, 50), sort.by = "p.value", rm.dupli = TRUE, filepath = NULL){
  cond.num <- c()
  init.cond.num <- Inf
  for (i in 1:length(increment)){
    res <- build_sigMat(diffexpMat, ctDEGs, increment[i], sort.by, rm.dupli)
    smat <- res$sigMat
    mkers <- res$markers
    cond.num[i] <- kappa(smat, exact=T)
    if (cond.num[i] < init.cond.num){
      opt.num <- increment[i]
      init.cond.num <- cond.num[i]
      markers <- mkers
      sigMat <- smat
    }
  }
  res <- list("sigMat" = sigMat, "markers" = markers, "min.cond.num" = init.cond.num, "num.top" = opt.num)
  if (!is.null(filepath)){
    saveRDS(res, file.path(filepath, "sigMat_condNum.rds"))
  }
  return(res)
}
