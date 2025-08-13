#' Convert AA substitutions from 3 letter to 1 letter
#'
#' @param SUB substitution using 3 letter representation
#' @return substitution using 1 letter representation
#' @description Convert amino acids from 3 letter to 1 letter in a substitution
#' @export
SUB3to1 <- function(SUB) {
  ref <- str_sub(SUB, 1, 3)
  pos_loc <- str_locate(SUB, "[:digit:]+")
  pos <- str_sub(SUB, pos_loc[,1], pos_loc[,2])
  alt <- str_sub(SUB, pos_loc[,2] + 1, pos_loc[,2] + 3)
  # %3D are synonymous_variant labels in VEP
  alt[alt == "%3D"] <- ref[alt == "%3D"]
  ref1 <- AA3to1(ref)
  alt1 <- AA3to1(alt)
  # If 3 letter AA to 1 letter AA is not present (del),
  # copy the letters intead.
  alt1[is.na(alt1)] <- alt[is.na(alt1)]
  trailing <- str_sub(SUB, pos_loc[,2] + 4, -1)
  output <- paste0(ref1, pos, alt1, trailing)
  return(output)
}
