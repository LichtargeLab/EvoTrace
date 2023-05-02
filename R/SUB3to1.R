#' Convert AA substitutions from 3 letter to 1 letter
#'
#' @param SUB substitution using 3 letter representation
#' @return substitution using 1 letter representation
#' @description Convert amino acids from 3 letter to 1 letter in a substitution
#' @export
SUB3to1 <- function(SUB) {
  ref <- str_sub(SUB, 1, 3) %>%
    AA3to1()
  pos <- str_sub(SUB, 4, -4)
  alt <- str_sub(SUB, -3, -1) %>%
    AA3to1()
  output <- paste0(ref, pos, alt)
  return(output)
}
