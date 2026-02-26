#' This function is to identify prevalent cell types from spatial transcriptomic (ST) data
#'
#' @param Y numeric matrix, spatial trancriptomic gene expression matrix
#' @param markers list of gene markers for each cell type
#'
#' @return character vector of prevalent cell types
#' @export
#'
#' @author Wenjun Shen
#'
identify_prevalent_cells <- function(Y, markers) {
  mean.Y <- rowMeans(Y)
  geneFoldchange <- log2(Y+1) - log2(mean.Y+1)
  celltypes <- names(markers)
  absolute = list()
  relative = list()
  median_abs = c()
  median_ref = c()
  for (i in seq_along(celltypes)) {
    ct <- celltypes[i]
    ct_markers <- intersect(rownames(Y), markers[[ct]])
    abs <- Y[ct_markers, , drop = FALSE]
    rel <- geneFoldchange[ct_markers, , drop = FALSE]
    absolute[[ct]] <- colMeans(abs, na.rm = TRUE)
    relative[[ct]] <- colMeans(rel, na.rm = TRUE)
    median_abs[ct] <- median(absolute[[ct]], na.rm = TRUE)
    median_ref[ct] <- median(relative[[ct]], na.rm = TRUE)
  }

  overall_q3_abs <- quantile(unlist(absolute), 0.75, na.rm = TRUE)
  prevalent_cells <- celltypes[median_abs > overall_q3_abs & abs(2^median_ref - 1) < 0.1]
  if (length(prevalent_cells) > 0){
    cat("Prevalent cells:", paste(prevalent_cells, collapse = ", "), "\n")
  }else{
    cat("Prevalent cells: None", "\n")
  }
  return(prevalent_cells)
}
