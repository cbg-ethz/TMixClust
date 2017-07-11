#' @title Simulated time-series gene expression data
#' @description This data set contains a toy example of time-series gene
#' expression data.
#' @format A \code{data frame} with 91 rows and 6 columns.
#' The columns correspond to different time points, while the rownames of the
#' data frame correspond to gene names.
#' @author Monica Golumbeanu, \email{monica.golumbeanu@bsse.ethz.ch}
#' @references Golumbeanu M, Desfarges S, Hernandez C, Quadroni M, Rato S,
#' Mohammadi P, Telenti A, Beerenwinkel N, Ciuffi A. (2017) Dynamics of
#' Proteo-Transcriptomic Response to HIV-1 Infection.
#' @examples
#' # Load the toy time series data provided with the TMixClust package
#' data("toy_data_df")
#'
#' # Print the first lines of the toy data frame
#' head(toy_data_df)
#'
#' @return toy data
#'
"toy_data_df"

#' @title TMixClust object containing the optimal clustering solution for the
#'  toy data with 3 clusters.
#' @description This object contains the result of clustering and
#' stability analysis
#' corresponding to the clustering solution with the highest likelihood among 10
#' different runs of clustering on the toy data with K=3 clusters.
#' @format A \code{TMixClust} object.
#' @author Monica Golumbeanu, \email{monica.golumbeanu@bsse.ethz.ch}
#' @references Golumbeanu M, Desfarges S, Hernandez C, Quadroni M, Rato S,
#' Mohammadi P, Telenti A, Beerenwinkel N, Ciuffi A. (2017) Dynamics of
#' Proteo-Transcriptomic Response to HIV-1 Infection.
#'
#' @examples
#' # Load the optimal clustering solution for the toy data
#' # provided with the TMixClust package
#' data("best_clust_toy_obj")
#'
#' # Print the first lines of the toy clustering object
#' head(best_clust_toy_obj)
#'
#' @return optimal clustering solution for the toy data
"best_clust_toy_obj"

#' @title TMixClust object containing the optimal clustering solution for the
#' yeast data.
#' @description This object contains the result of clustering and
#' stability analysis
#' corresponding to the clustering solution with the highest likelihood among 10
#' different runs of clustering on the yeast data with K=4 clusters.
#' @format A \code{TMixClust} object.
#' @author Monica Golumbeanu, \email{monica.golumbeanu@bsse.ethz.ch}
#' @references Golumbeanu M, Desfarges S, Hernandez C, Quadroni M, Rato S,
#' Mohammadi P, Telenti A, Beerenwinkel N, Ciuffi A. (2017) Dynamics of
#' Proteo-Transcriptomic Response to HIV-1 Infection.
#'
#' @examples
#' # Load the optimal clustering solution for the yeast data
#' # provided with the TMixClust package
#' data("best_clust_yeast_obj")
#'
#' # Print the first lines of the yeast clustering object
#' head(best_clust_yeast_obj)
#'
#' @return optimal clustering solution for the yeast data
"best_clust_yeast_obj"
