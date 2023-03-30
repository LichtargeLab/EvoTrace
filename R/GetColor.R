#' Assign colors for a given numeric vector
#'
#' @param x A numeric vector
#' @param lower_bound The lower bound of the scale
#' @param upper_bound The upper bound of the scale
#' @param color_type The color ranges to choose from.
#' @return A vector of hex colors
#' @description There are 100 colors stored in each color_type.
#' The color value assigned for each element in x is based on its relative position (percent ranking)
#' regarding the lower and upper bounds. The assignment is upper bound included. A value equal the lower
#' bound will not be assigned correctly, and an error will be given.
#' @export


GetColor <- function(x, lower_bound, upper_bound,
                     color_type = c("ET", "red_white", "red_white_blue", "white_red",
                                    "white_blue", "alphafold")) {
  if (min(x, na.rm = TRUE) <= lower_bound) {
    stop("Lower bound should be smaller than the minimum value of x.")
  }
  if (max(x, na.rm = TRUE) > upper_bound) {
    stop("Upper bound should be larger than or equal to the maximum value of x.")
  }
  scaled_x <- (x - lower_bound)/(upper_bound - lower_bound) * 100
  colorRange <- SelectColor(color_type = color_type)
  output <- colorRange[ceiling(scaled_x)]
  return(output)
}


