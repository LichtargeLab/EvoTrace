#' Apply 1.5 IQR filter for a given numeric vector
#'
#' @param x A given numeric vector to apply the filter
#' @param type If "bool" is used, them a boolean vector with the same length as x is returned.
#' If "ends" is used, then the lower and upper boundaries are returned
#' @param fold Fold of IQR when computing the upper and lower boundaries.
#' @return A boolean vector or the lower and upper boundaries
#' @description Apply IQR to vector x.
#' @export

IQR_filter <- function(x, type = c("bool", "ends"), fold = 1.5) {
  type <- match.arg(type)
  IQR <- quantile(x, 0.75) - quantile(x, 0.25)
  ul <- quantile(x, 0.75) + fold*IQR
  ll <- quantile(x, 0.25) - fold*IQR
  if(type != "bool") {
    output <- c(ll, ul)
    names(output) <- NULL
  } else {
    output <- (x <= ul) & (x >= ll)
  }
  return(output)
}
