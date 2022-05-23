#' Get Laplacian Matrix for a given distance matrix
#'
#' @param dist_mat Distance matrix. Can be obtained from GetDistanceMatrix.
#' @param cutoff dbl. The maximal distance two residues are consider adjacent.
#' @param normalize logic. whether normalized Laplacian is returned
#' @return The Laplacian matrix
#' @description Compute Laplacian Matrix from given distance matrix. Distance matrix
#' can be obtained from GetDistanceMatrix function.
#' @export

GetLaplacianMatrix <- function(dist_mat, cutoff = 4, normalize = FALSE) {
  A_mat <- (dist_mat <= cutoff) + 0
  diag(A_mat) <- 0
  D_mat <- diag(x = rowSums(A_mat))
  L_mat <- D_mat - A_mat
  if (normalize == TRUE) {
    D_mat_2 <- diag(x = rowSums(A_mat)^-0.5)
    L_mat <- D_mat_2 %*% L_mat %*% D_mat_2
  }
  return(L_mat)
}

