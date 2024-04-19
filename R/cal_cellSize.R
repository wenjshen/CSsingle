#' This function is to estimate the cell sizes by using ERCC spike-in controls
#'
#' @param sc.eset SeuratObject for scRNA-seq data
#' @param cellType character, name of metadata column indicating cell type
#' @param probs numeric vector of probabilities with values in [0, 1], default is 0.75
#'
#' @return cellSize a vector of calculated cell sizes
#' @export
#' @import Seurat
#'
#' @author Wenjun Shen
#'
cal_cellSize <- function(sc.eset, cellType, probs=0.75){
  is.spikes <- grepl("^ERCC-", toupper(rownames(sc.eset)))
  counts <- as.matrix(sc.eset@assays$RNA@counts)
  spike.counts <- counts[is.spikes, ]
  normFactor <- apply(spike.counts,MARGIN=2,
                     FUN=quantile, probs=probs, na.rm=TRUE)
  norm.counts <- t(t(counts)/(normFactor+1))
  gene.norm.counts <- norm.counts[!is.spikes, ]
  gene.totals <- colSums(gene.norm.counts)
  scf <- aggregate(gene.totals, by=list(sc.eset@meta.data[[cellType]]), FUN=mean)
  cellSize <- scf$x
  names(cellSize) <- gsub("[.]", "_", make.names(scf$Group.1))
  return(cellSize)
}
