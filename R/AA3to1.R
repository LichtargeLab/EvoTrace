#' Convert amino acids from 3 letter to 1 letter
#'
#' @param AA a string of 3 letter AA
#' @return A string of 1 letter AA
#' @description Convert amino acids from 3 letter to 1 letter
#' @export
AA3to1 <- function(AA) {
  AA3 <- c("HIS", "PRO", "GLU", "THR", "LEU", "VAL", "LYS", "ASP", "ALA",
           "GLN", "GLY", "ARG", "TYR", "ILE", "ASN", "SER", "PHE", "MET",
           "CYS", "TRP", "MSE", "TER")
  AA1 <- c("H", "P", "E", "T", "L", "V", "K", "D", "A",
           "Q", "G", "R", "Y", "I", "N", "S", "F", "M",
           "C", "W", "M", "X")
  names(AA1) <- AA3
  output <- AA1[toupper(AA)]
  names(output) <- NULL
  return(output)
}
