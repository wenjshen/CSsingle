% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CSsingle.R
\name{CSsingle}
\alias{CSsingle}
\title{This function is to estimate the cell-type proportions from bulk gene expression data}
\usage{
CSsingle(
  sigMat,
  mixture,
  markers,
  cellSize = NULL,
  dampened = FALSE,
  filepath = NULL
)
}
\arguments{
\item{sigMat}{numeric matrix, signature matrix}

\item{mixture}{bulk gene expression matrix}

\item{markers}{list of gene markers for each cell type}

\item{cellSize}{Users can specify a numeric vector of cell sizes, default is NULL, which means that we assume the absolute amount of total mRNA is similar across different cell types in deconvolution}

\item{dampened}{whether to introduce a upper bound constant that limits the maximum value any weight can take on, default is FALSE}

\item{filepath}{directory where the results will be saved, default is NULL}
}
\value{
a numeric matrix of cell type proportions estimated by CSsingle (add up to 1)
}
\description{
This function is to estimate the cell-type proportions from bulk gene expression data
}
\author{
Wenjun Shen
}
