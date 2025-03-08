#' @param X numeric matrix
#' @param y numeric vector
#' @param X_D numeric matrix
#' @param y_D numeric vector
#' @param est.prop numeric vector of estimated cell-type proportions
#' @param markers list of gene markers for each cell type
#' @param dampened whether to introduce a upper bound constant that limits the maximum value any weight can take on, default is FALSE
#' @param beta numeric vector of [0, 0.5, 1]
#'
#' @return a numeric vector of cell type proportions (add up to 1)
#'
#' @importFrom nnls nnls
#'
#' @author Wenjun Shen
#' @export
#'
CSsingle_basisJ <- function(X, y, X_D, y_D, est.prop, markers, dampened = FALSE, beta){
  T_set <- X * rep(est.prop, each=nrow(X))
  t_set <- rowSums(T_set)
  diff <- (abs(y - t_set))^(beta)*t_set^(1-beta)
  gene.weight <- 1/(diff^2 + 1e-4)
  if (dampened){
    limit <- quantile(gene.weight, seq(0.01, 1, 0.05))
    cor.pval <- c()
    correlation <- c()
    for (h in 1:length(limit)){
      gweight <- gene.weight
      gweight[gweight > limit[h]] <- limit[h]
      X.weight <- sweep(X, 1, sqrt(gweight), '*')
      y.weight <- y * sqrt(gweight)
      w <- nnls(X.weight, y.weight)
      wnorm <- w$x/sum(w$x)
      correlation[h] <- cor(X_D %*% matrix(wnorm), y_D, method="spearman")
    }
    idx.opt <- order(correlation, decreasing = T)
    gene.weight[gene.weight > limit[idx.opt[1]]] <- limit[idx.opt[1]]
  }
  ## weighted NNLS
  X.weight <- sweep(X, 1, sqrt(gene.weight), '*')
  y.weight <- y * sqrt(gene.weight)
  w <- nnls(X.weight, y.weight)
  wnorm <- w$x/sum(w$x)
  return(wnorm)
}
