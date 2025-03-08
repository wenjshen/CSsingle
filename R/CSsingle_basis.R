#' This function is to estimate the cell-type proportions from individual bulk mixtures or spatial spots \code{y}
#'
#' @param X numeric matrix, signature matrix
#' @param y numeric vector of bulk or spatial trancriptomic expression
#' @param markers list of gene markers for each cell type
#' @param cellSize Users can specify a numeric vector of cell sizes, default is NULL, which means that we assume the absolute amount of total mRNA is similar across different cell types in deconvolution
#' @param enrich character vector of enriched cell types
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
CSsingle_basis <- function(X, y, markers, cellSize = NULL, enrich = NULL, dampened = FALSE, beta){
  log.X <- log(X + 1e-4)
  log.y <- log(y + 1e-4)
  present.ct.idx <- c()
  intercepts <- c()
  for (j in 1 : ncol(X)){
    ct <- colnames(X)[j]
    if (is.null(enrich) | ct %in% enrich){
      id <- names(y) %in% toupper(markers[[ct]])
      if (sum(id) > 5){
        present.ct.idx <- c(present.ct.idx, j)
        intercept <- mean(log.y[id]-log.X[id, j, drop=F])
        intercepts <- c(intercepts, intercept)
      }
    }
  }

  Est.prop <- matrix(0, 1, ncol(X))
  colnames(Est.prop) <- colnames(X)
  if(length(present.ct.idx) == 0){
    if (!is.null(enrich) & length(enrich) == 1){
      Est.prop[ , enrich] <- 1
    }else{
      Est.prop <- matrix(NA, 1, ncol(Est.prop))
    }
  }else if(length(present.ct.idx) == 1){
    Est.prop[ , colnames(X)[present.ct.idx]] <- 1
  }else{
    present.ct <- colnames(X)[present.ct.idx]
    eb <- exp(intercepts)
    delta <- eb/sum(eb)
    T_set <- X[, present.ct] * rep(delta, each=nrow(X))
    t_set <- rowSums(T_set)

    marker.genes <- toupper(unique(unlist(markers[present.ct])))
    marker.genes <- marker.genes[marker.genes %in% rownames(X)]
    X_D <- X[, present.ct]
    y_D <- y
    X <- X[marker.genes, present.ct]
    y <- y[marker.genes]
    t_set <- t_set[marker.genes]

    diff <- (abs(y - t_set))^(beta)*t_set^(1-beta)
    gene.weight <- 1/(diff^2 + 1e-4)
    if (!is.null(cellSize)){
      cellSize <- cellSize[present.ct]
      X <- t(t(X) * cellSize)
    }

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
        correlation[h] <- cor(X_D %*% matrix(wnorm), y_D, method = "spearman")
      }
      idx.opt <- order(correlation, decreasing = T)
      gene.weight[gene.weight > limit[idx.opt[1]]] <- limit[idx.opt[1]]
    }

    ## weighted NNLS
    X.weight <- sweep(X, 1, sqrt(gene.weight), '*')
    y.weight <- y * sqrt(gene.weight)
    w <- nnls(X.weight, y.weight)
    wnorm <- w$x/sum(w$x)
    num.iter <- 0
    change <- 1
    changes <- c()
    while(change>0.01 & num.iter<1000){
      wnormNew <- CSsingle_basisJ(X, y, X_D, y_D, wnorm, markers, dampened, beta)
      wnormNew <- colMeans(rbind(wnormNew, wnorm))
      change <- norm(as.matrix(wnormNew - wnorm))
      wnorm <- wnormNew
      num.iter <- num.iter + 1
      changes <- c(changes, change)
    }
    Est.prop[, colnames(X)] <- wnorm
  }
  return(Est.prop)
}
