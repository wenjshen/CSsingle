#' This function is to combine CSsingle with each signature matrix to get estimates of the cell type proportions. The signature matrices were ordered in descending order based on the Spearman correlation between their inferred and real bulk gene expression data.
#'
#' @param diffexpMat numeric matrix, a differentially expressed matrix
#' @param mixture numeric matrix, bulk or spatial trancriptomic gene expression matrix (raw counts for RNA-seq or ST, RMA normalized intensity values in linear scale for Microarray)
#' @param ctDEGs data frame, matrix containing a ranked list of differentially expressed genes, and associated statistics
#' @param cellSize Users can specify a numeric vector of cell sizes, default is NULL, which means that we assume the absolute amount of total mRNA is similar across different cell types in deconvolution
#' @param enrichment whether to identify cell types enriched in ST spots, default is FALSE
#' @param enrich.thres default is 0
#' @param dampened whether to introduce a upper bound constant that limits the maximum value any weight can take on, default is FALSE
#' @param beta numeric vector of [0, 0.5, 1]
#' @param increment numeric vector, create multiple candidate matrices by specifying the number of top markers from each cell type
#' @param sort.by rank the markers in ascending order by their p-values (p.value, default) or in descending order by their fold changes (fold.change)
#' @param rm.dupli whether to exclude markers shared between two or more cell types, default is TRUE
#' @param filepath directory where the results will be saved, default is NULL
#'
#' @return a list with two attributes:
#' \item{Est.prop}{a numeric matrix of cell type proportions estimated by CSsingle through maximizing the Spearman correlation between the inferred and real bulk gene expression data (add up to 1)}
#' \item{Est.propS}{a ranked list of estimated cell type proportions}
#' \item{corr}{a numeric vector of Spearman's correlation coefficients}
#'
#' @export
#' @import Seurat
#' @importFrom nnls nnls
#'
#' @author Wenjun Shen
#'
CSsingle_optimize <- function(diffexpMat, mixture, ctDEGs, cellSize = NULL, enrichment = FALSE, enrich.thres = 0, dampened = FALSE, beta = c(0, 0.5, 1), increment = seq(50, 200, 50), sort.by = "p.value", rm.dupli = TRUE, filepath = NULL){
  Est.propS <- list()
  topSigMat <- list()
  topMarkers <- list()
  enrich.ct <- NULL
  if (enrichment){
    res <- build_sigMat(diffexpMat, ctDEGs, 50, sort.by, rm.dupli)
    enrich.ct <- enrichmentMarkers(res$sigMat, mixture, res$markers, enrich.thres = enrich.thres)
  }
  for (i in increment){
    res <- build_sigMat(diffexpMat, ctDEGs, i, sort.by, rm.dupli)
    topSigMat[[as.character(i)]] <- res$sigMat
    topMarkers[[as.character(i)]] <- res$markers
    Est.propS[[as.character(i)]] <- CSsingle(topSigMat[[as.character(i)]], mixture, topMarkers[[as.character(i)]], cellSize, enrich.ct = enrich.ct, dampened = dampened, beta = beta)
  }
  ###
  mixture.norm <- apply(mixture, 2, function(x) 1e4*x/sum(x))
  sigmat <- topSigMat[[as.character(increment[length(increment)])]]
  rownames(sigmat) <- toupper(rownames(sigmat))
  rownames(mixture.norm) <- toupper(rownames(mixture.norm))
  com.genes <- intersect(rownames(sigmat), rownames(mixture.norm))
  sigmat <- sigmat[com.genes, ]
  mixture.norm <- mixture.norm[com.genes, , drop = F]

  if (!is.null(cellSize)){
    sigmat <- t(t(sigmat) * cellSize[colnames(sigmat)])
  }

  correlation <- c()
  for (i in 1 : length(Est.propS)){
    est.mixture <- sigmat %*% t(Est.propS[[i]])
    corr <- c()
    for (j in 1 : ncol(mixture.norm)){
      if (any(is.na(est.mixture[, j]))){
        corr[j] <- -1
      }else{
        corr[j] <- cor(est.mixture[, j], mixture.norm[, j], method = "spearman")
      }
    }
    correlation[i] <- round(mean(corr, na.rm = T), 3)
  }
  Est.propS <- Est.propS[order(correlation, decreasing = T)]
  correlation <- correlation[order(correlation, decreasing = T)]
  Est.prop <- Est.propS[[1]]

  if (!is.null(filepath)){
    # save results
    res <- list(Est.prop = Est.prop, Est.propS = Est.propS, correlation = corr)
    saveRDS(res, file.path(filepath, "CSsingle_optimize-Results.rds"))
    write.table(Est.prop, file = file.path(filepath, "CSsingle_optimize-Results.txt"), sep="\t", row.names=T, col.names=T, quote=F)
  }
  return(Est.prop)
}
