#' Compute smoothness of residue ranks over a given structure.
#'
#' @param dist_mat Distance matrix. Can be obtained from GetDistanceMatrix.
#' @param vari a numerical vector that stores the ranks/scores of the residues. The
#' length of vari should be the same as the dim of dist_mat. The order should also be
#' matched.
#' @param cutoff dbl. The maximal distance two residues are consider adjacent.
#' @param variable_normalization Normallization type of the variable. Default is L2.
#' @param seed dbl. Seeding for randomization.
#' @return A list that stores the smoothness of the vari vector, 10000 random smoothness,
#' the z score and p value.
#' @description Compute Laplacian smoothness of vector vari over the protein structure.
#' Smoothness background is computed for 10000 random vari vector which are generated
#' by shuffling the vari vector. Z score and p values (normal, one sided) are calculated
#' for the vari smoothness comparing to background.
#' @export

GetLaplacianSmoothness <- function(dist_mat, vari, cutoff = 4, variable_normalization = c("L2", "L1", "none"),
                                   seed = 100) {
  check_POS <- nrow(dist_mat) == length(vari)
  if (check_POS == FALSE) {
    stop("Length don't match in dist_mat and vari")
  }

  mat_L <- GetLaplacianMatrix(dist_mat, cutoff = cutoff, normalize = FALSE)

  variable_normalization <- match.arg(variable_normalization)
  if (variable_normalization == "L2") {
    scale_fun <- function(x) {x / sqrt(sum(x^2))}
  } else if (variable_normalization == "L1") {
    scale_fun <- function(x) {x / sum(abs(x))}
  } else {
    scale_fun <- function(x) {x}
  }
  scaled_vari <- scale_fun(vari)

  GetSmoothness <- function(x, mat_L) {
    output <- (t(x) %*% mat_L %*% x)[1,1]
    return(output)
  }
  smoothness <- GetSmoothness(scaled_vari, mat_L)

  size <- length(scaled_vari)
  set.seed(seed)
  random_vari <- replicate(10000, sample(1:size, size = size))
  random_vari <- matrix(sort(scaled_vari)[random_vari], nrow = size)

  smoothness_random <- apply(random_vari, MARGIN = 2, GetSmoothness, mat_L = mat_L)
  z <- (smoothness - mean(smoothness_random))/sd(smoothness_random)
  pval <- pnorm(z)
  output <- list(smoothness=smoothness, smoothness_random=smoothness_random,
                 z=z, pval=pval)
  return(output)
}
