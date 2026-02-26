#' This function is to combine CSsingle with each signature matrix to get estimates of the cell type proportions. The signature matrices were ordered in descending order based on the Spearman correlation between their inferred and real bulk gene expression data.
#'
#' @param designMat numeric matrix, a gene by cell-type matrix (genome-wide) from scRNA-seq
#' @param diffexpMat numeric matrix, a differentially expressed matrix
#' @param mixture numeric matrix, bulk or spatial trancriptomic gene expression matrix (raw counts for RNA-seq or ST, RMA normalized intensity values in linear scale for Microarray)
#' @param ctDEGs data frame, matrix containing a ranked list of differentially expressed genes, and associated statistics
#' @param cellSize users can specify a numeric vector of cell sizes, default is NULL, which assumes that the absolute amount of total mRNA is similar across different cell types in deconvolution. Set to NULL when ERCC spike-ins are unavailable from scRNA-seq reference or for spatial transcriptomics (ST) analysis.
#' @param enrich.ct list of pre-computed enriched cell types for each spot (optional)
#' @param dampened whether to introduce a upper bound constant that limits the maximum value any weight can take on, default is FALSE
#' @param beta numeric vector of [0, 0.5, 1]
#' @param increment numeric vector, create multiple candidate matrices by specifying the number of top markers from each cell type
#' @param sort.by rank the markers in ascending order by their p-values (p.value, default) or in descending order by their fold changes (fold.change)
#' @param rm.dupli whether to exclude markers shared between two or more cell types, default is FALSE
#' @param tune.correction when cellSize = NULL and ERCC spike-ins are unavailable, setting tune.correction = TRUE triggers computational estimation of cell size coefficients. Default is FALSE
#' @param tune.threshold threshold for tuning correction, default is 0.25
#' @param n_cores number of cores for parallel processing, default is 1
#' @param filepath directory where results will be saved, default is NULL
#'
#' @return a list with four attributes:
#' \item{Est.prop}{a numeric matrix of cell type proportions estimated by CSsingle through maximizing the Spearman correlation between the inferred and real bulk gene expression data (add up to 1)}
#' \item{opt_num_top}{optimal number of top markers selected}
#' \item{opt_sigmat}{optimal signature matrix}
#' \item{opt_markers}{optimal marker genes}
#'
#' @export
#' @import Seurat
#' @importFrom nnls nnls
#' @importFrom parallel detectCores mclapply
#' @importFrom stats cor
#'
#' @author Wenjun Shen
#'
CSsingle_optimize <- function(designMat, diffexpMat, mixture, ctDEGs, cellSize = NULL, enrich.ct = NULL, dampened = FALSE, beta = c(0, 0.5, 1), increment = seq(50, 200, 50), sort.by = "p.value", rm.dupli = FALSE, tune.correction = FALSE, tune.threshold = 0.25, n_cores = 1, filepath = NULL){
  # Prepare data
  rownames(designMat) <- toupper(rownames(designMat))
  rownames(mixture) <- toupper(rownames(mixture))
  rownames(diffexpMat) <- toupper(rownames(diffexpMat))
  mixture.norm <- apply(mixture, 2, function(x) 1e4 * x / sum(x))
  com.genes <- intersect(rownames(designMat), rownames(mixture))
  designMat <- designMat[com.genes, ]
  mixture.norm <- mixture.norm[com.genes, ]
  non_markers <- setdiff(rownames(designMat), rownames(diffexpMat))
  if (length(non_markers) < 15) non_markers = rownames(designMat)
  if (!is.null(cellSize) & length(unique(cellSize)) > 1){
    designMat <- t(t(designMat) * cellSize[colnames(designMat)])
  }

  # Parallel execution for increment loop
  process_increment <- function(i) {
    res <- build_sigMat(diffexpMat, ctDEGs, i, sort.by, rm.dupli)
    est_prop <- CSsingle(res$sigMat, mixture, res$markers, cellSize,
                         enrich.ct, dampened, beta)
    return(list(sigMat = res$sigMat, markers = res$markers, est_prop = est_prop))
  }

  if (is.null(n_cores)){
    n_cores <- 1
  }
  if (.Platform$OS.type != "windows" && n_cores > 1) {
    message("Running in parallel with ", n_cores, " cores on Unix system")
    results <- parallel::mclapply(increment, process_increment, mc.cores = n_cores)
  } else {
    message("Running sequentially (Windows or n_cores = 1)")
    results <- lapply(increment, process_increment)
  }

  Est.propS <- list()
  topSigMat <- list()
  topMarkers <- list()
  for (i in seq_along(increment)) {
    topSigMat[[as.character(increment[i])]] <- results[[i]]$sigMat
    topMarkers[[as.character(increment[i])]] <- results[[i]]$markers
    Est.propS[[as.character(increment[i])]] <- results[[i]]$est_prop
  }

  # Evaluate initial correlations
  correlation <- sapply(seq_along(Est.propS), function(i) {
    est.mixture <- designMat[, colnames(Est.propS[[i]])] %*% t(Est.propS[[i]])
    calculate_correlation(est.mixture, mixture, non_markers)
  })
  max_cor_init <- max(correlation)

  # Get best initial results
  best_index <- which.max(correlation)
  opt_num_top <- names(Est.propS)[best_index]
  Est.prop <- Est.propS[[best_index]]
  opt_sigmat <- topSigMat[[best_index]]
  opt_markers <- topMarkers[[best_index]]
  Prop_init = Est.prop

  if (!is.null(filepath)){
    # save results
    res <- list(Est.prop = Est.prop, opt_num_top = opt_num_top, opt_sigmat = opt_sigmat, opt_markers = opt_markers)
    saveRDS(res, file.path(filepath, "CSsingle_optimize-Results.rds"))
    write.table(Est.prop, file = file.path(filepath, "CSsingle_optimize-Results.txt"), sep="\t", row.names=T, col.names=T, quote=F)
  }

  if (is.null(cellSize) & tune.correction){
    # Identify potential misestimated cell types by analyzing residual correlations
    Prop_init <- Est.prop
    sigmat <- opt_sigmat
    markers <- opt_markers
    mixture_reconstructed <- designMat[, colnames(Prop_init)] %*% t(Prop_init)
    num.top <- seq(200, 1000, 50)
    cor_mean_top <- matrix(NA, nrow = length(num.top), ncol = ncol(sigmat))
    colnames(cor_mean_top) <- colnames(sigmat)
    rownames(cor_mean_top) <- num.top
    for (n_idx in seq_along(num.top)) {
      n <- num.top[n_idx]
      sigmat.topn <- build_sigMat(diffexpMat, ctDEGs, n, sort.by,
                                  rm.dupli = FALSE)$sigMat
      # Ensure we only use common genes
      common_genes <- intersect(rownames(sigmat.topn), rownames(mixture_reconstructed))
      residual <- mixture.norm[common_genes, ] - mixture_reconstructed[common_genes, ]
      sigmat.topn <- sigmat.topn[common_genes, ]

      # Calculate correlations for each cell type
      cor_mat <- matrix(NA, nrow = ncol(residual), ncol = ncol(sigmat.topn))
      colnames(cor_mat) <- colnames(sigmat.topn)
      for (j in 1:ncol(residual)) {
        for (k in 1:ncol(sigmat.topn)) {
          if (all(!is.na(residual[, j])) && all(!is.na(sigmat.topn[, k]))) {
            cor_mat[j, k] <- abs(cor(sigmat.topn[, k], residual[, j], method = "spearman"))
          }
        }
      }
      # Store mean correlation for each cell type across all samples
      cor_mean_top[n_idx, ] <- colMeans(cor_mat, na.rm = TRUE)
    }

    # Find the best number of top genes
    rmax <- apply(cor_mean_top, 1, max, na.rm = TRUE)
    best_row_index <- which.max(rmax)
    cor_mean <- cor_mean_top[best_row_index, ]

    # Check if correction is needed
    if (max(cor_mean, na.rm = TRUE) <= tune.threshold) {
      cat("No cell size correction needed. All correlations below threshold of",
          tune.threshold, "\n")
      return(list(Est.prop = Prop_init, opt_num_top = opt_num_top, opt_sigmat = opt_sigmat, opt_markers = opt_markers))
    }

    num_celltypes <- length(cor_mean)
    high_cor_mask <- cor_mean > tune.threshold
    high_cor_sorted <- sort(cor_mean[high_cor_mask], decreasing = TRUE)
    max_to_correct <- min(length(high_cor_sorted), num_celltypes - 1)
    high_cor_celltypes <- names(high_cor_sorted)[1:max_to_correct]
    high_prop_celltypes <- colnames(Prop_init)[colMeans(Prop_init) > 0.05]
    candidate_celltypes <- intersect(high_cor_celltypes, high_prop_celltypes)

    # Identify potential overestimated cell types
    lambda <- seq(0.5, 3, 0.5)
    evaluate_candidate_celltypes <- function(ct) {
      scaling_factors <- rep(1, num_celltypes)
      names(scaling_factors) <- names(cor_mean)
      correlations <- sapply(lambda, function(l) {
        norm_cor <- (cor_mean[ct] - tune.threshold) /
          (max(cor_mean, na.rm = TRUE) - tune.threshold)
        scaling_factors[ct] <- 1 + l * (1 - exp(-2 * norm_cor))

        czFactor <- scaling_factors / min(scaling_factors) * 1e5
        Est.propS2 <- CSsingle(sigmat, mixture, markers, czFactor,
                               enrich.ct, dampened, beta)

        curr_sigmat <- t(t(designMat) * czFactor[colnames(designMat)])
        est.mixture <- curr_sigmat[, colnames(Est.propS2)] %*% t(Est.propS2)
        calculate_correlation(est.mixture, mixture, non_markers)
      })

      if (max(correlations, na.rm = TRUE) > max_cor_init) {
        best_lambda_index <- which.max(correlations)
        return(list(
          celltype = ct,
          best_lambda_index = best_lambda_index
        ))
      } else {
        return(NULL)
      }
    }

    if (.Platform$OS.type != "windows" && n_cores > 1) {
      candidate_results <- mclapply(candidate_celltypes, evaluate_candidate_celltypes, mc.cores = n_cores)
    } else {
      candidate_results <- lapply(candidate_celltypes, evaluate_candidate_celltypes)
    }

    celltypes_to_correct <- unlist(sapply(candidate_results, function(x) x$celltype))
    best_lambda_indices <- unlist(sapply(candidate_results, function(x) x$best_lambda_index))

    if (length(celltypes_to_correct) == 0) {
      cat("No overestimated cell types identified - initial estimates maintained\n")
      return(list(Est.prop = Prop_init, opt_num_top = opt_num_top, opt_sigmat = opt_sigmat, opt_markers = opt_markers))
    }

    # Perform correction if potential overestimated cell types identified
    tune_identified_celltypes <- function(i) {
      l <- lambda[i]
      scaling_factors <- rep(1, num_celltypes)
      names(scaling_factors) <- names(cor_mean)

      for (ct in celltypes_to_correct) {
        norm_cor <- (cor_mean[ct] - tune.threshold) /
          (max(cor_mean, na.rm = TRUE) - tune.threshold)
        scaling_factors[ct] <- 1 + l * (1 - exp(-2 * norm_cor))
      }
      czFactor_val <- scaling_factors / min(scaling_factors) * 1e5
      Est.propS2_val <- CSsingle(sigmat, mixture, markers, czFactor_val,
                                 enrich.ct, dampened, beta)

      curr_sigmat <- t(t(designMat) * czFactor_val[colnames(designMat)])
      est.mixture <- curr_sigmat[, colnames(Est.propS2_val)] %*% t(Est.propS2_val)
      correlation_val <- calculate_correlation(est.mixture, mixture, non_markers)

      return(list(czFactor = czFactor_val, Est.propS2 = Est.propS2_val, correlation = correlation_val))
    }

    lambda <- c(0, lambda)
    if (.Platform$OS.type != "windows" && n_cores > 1) {
      lambda_results <- mclapply(seq_along(lambda), tune_identified_celltypes, mc.cores = n_cores)
    } else {
      lambda_results <- lapply(seq_along(lambda), tune_identified_celltypes)
    }

    # Extract results
    czFactor <- lapply(lambda_results, function(x) x$czFactor)
    Est.propS2 <- lapply(lambda_results, function(x) x$Est.propS2)
    correlation_results <- sapply(lambda_results, function(x) x$correlation)

    # Get best correction results
    best_corr_index <- which.max(correlation_results)
    if (best_corr_index == 1) {
      cat("Cell size correction did not improve correlation - initial estimates maintained\n")
      return(list(Est.prop = Prop_init, opt_num_top = opt_num_top, opt_sigmat = opt_sigmat, opt_markers = opt_markers))
    }

    estimate_czf <- czFactor[[best_corr_index]]
    Est.prop <- Est.propS2[[best_corr_index]]

    # Save results
    if (!is.null(filepath)) {
      res <- list(Est.prop = Est.prop, estimate_czf = estimate_czf, opt_num_top = opt_num_top, opt_sigmat = opt_sigmat, opt_markers = opt_markers)
      saveRDS(res, file.path(filepath, "CSsingle_optimize-Results.rds"))
      write.table(Est.prop, file = file.path(filepath, "CSsingle_optimize-Results.txt"), sep="\t", row.names=T, col.names=T, quote=F)
    }
  }
  return(list(Est.prop = Est.prop, opt_num_top = opt_num_top, opt_sigmat = opt_sigmat, opt_markers = opt_markers))
}
