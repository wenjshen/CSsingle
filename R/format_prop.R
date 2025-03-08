format_prop <- function(data.list){
  uni.colnames = Reduce(union, lapply(data.list, colnames))
  num.col <- length(uni.colnames)
  int.rownames <- Reduce(intersect, lapply(data.list, rownames))
  num.row <- length(int.rownames)
  new.data <- matrix(0.0, num.row, num.col)
  rownames(new.data) <- int.rownames
  colnames(new.data) <- uni.colnames
  new.list <- list()
  for (i in 1:length(data.list)){
    new.data.i <- new.data
    com.row <- intersect(rownames(data.list[[i]]), int.rownames)
    new.data.i[com.row, colnames(data.list[[i]])] <- data.list[[i]][com.row, , drop=F]
    new.list[[i]] <- new.data.i
  }
  return(new.list)
}
