#' This function is to combine CSsingle with each signature matrix to get estimates of the cell type proportions. The signature matrices were ordered in descending order based on the Spearman's correlation between their inferred and real bulk gene expression data.
#'
#' @param diffexpMat numeric matrix, a differentially expressed matrix
#' @param mixture bulk gene expression matrix
#' @param ctDEGs data frame, matrix containing a ranked list of differentially expressed genes, and associated statistics
#' @param cellSize Users can specify a numeric vector of cell sizes, default is NULL, which means that we assume the absolute amount of total mRNA is similar across different cell types in deconvolution
#' @param dampened whether to introduce a upper bound constant that limits the maximum value any weight can take on, default is FALSE
#' @param increment numeric vector, create multiple candidate matrices by specifying the number of top markers from each cell type
#' @param sort.by rank the markers in ascending order by their p-values (p.value, default) or in descending order by their fold changes (fold.change)
#' @param rm.dupli whether to exclude markers shared between two or more cell types, default is TRUE
#' @param filepath directory where the results will be saved, default is NULL
#'
#' @return a list with two attributes:
#' \item{Est.prop}{a numeric matrix of cell type proportions estimated by CSsingle through maximizing the Spearman's correlation between the inferred and real bulk gene expression data (add up to 1)}
#' \item{Est.propS}{a ranked list of estimated cell type proportions}
#' \item{corr}{a numeric vector of Spearman's correlation coefficients}
#'
#' @export
#' @importFrom nnls nnls
#'
#' @author Wenjun Shen
#'
CSsingle_optimize <- function(diffexpMat, mixture, ctDEGs, cellSize = NULL, dampened = FALSE, increment = seq(50, 300, 50), sort.by = "p.value", rm.dupli = TRUE, filepath = NULL){
  Est.propS <- list()
  topSigMat <- list()
  topMarkers <- list()
  for (i in increment){
    res <- build_sigMat(diffexpMat, ctDEGs, i, sort.by, rm.dupli)
    topSigMat[[as.character(i)]] <- res$sigMat
    topMarkers[[as.character(i)]] <- res$markers
    Est.propS[[as.character(i)]] <- CSsingle(topSigMat[[as.character(i)]], mixture, topMarkers[[as.character(i)]], cellSize, dampened)
  }
  ###
  # column normalization
  mixture = apply(mixture, 2, function(x) 1e4*x/sum(x))
  sigmat <- topSigMat[[as.character(increment[length(increment)])]]
  rownames(sigmat) <- toupper(rownames(sigmat))
  rownames(mixture) <- toupper(rownames(mixture))
  com.genes <- intersect(rownames(sigmat), rownames(mixture))
  sigmat <- sigmat[com.genes, ]
  if (!is.null(cellSize)){
    sigmat = t(t(sigmat) * cellSize[colnames(sigmat)])
  }
  mixture <- mixture[com.genes, , drop = F]
  log.mixture = log(mixture + 1e-4)
  corr = c()
  for (i in 1 : length(Est.propS)){
    est.mixture <- sigmat %*% t(Est.propS[[i]])
    log.est.mixture = log(est.mixture + 1e-4)
    correlation = c()
    for (j in 1 : ncol(mixture)){
      if (any(is.na(log.est.mixture[, j]))){
        correlation[j] <- 0
      }else{
        correlation[j] <- cor(est.mixture[, j], mixture[, j], method = "spearman", use = "complete.obs")
      }
    }
    corr[i] = round(mean(correlation),3)
  }

  Est.propS <- Est.propS[order(corr, decreasing = T)]
  corr <- corr[order(corr, decreasing = T)]

  Est.prop <- Est.propS[[1]]

  if (!is.null(filepath)){
    #save results
    res <- list(Est.prop = Est.prop, Est.propS = Est.propS, corr = corr)
    saveRDS(res, file.path(filepath, "CSsingle_optimize-Results.rds"))
    write.table(Est.prop, file = file.path(filepath, "CSsingle_optimize-Results.txt"), sep="\t", row.names=T, col.names=T, quote=F)
  }
  return(Est.prop)
}
