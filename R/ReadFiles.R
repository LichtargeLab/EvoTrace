#' Read legacy ET format
#'
#' @param ET_file Path to ET file
#' @return a dataframe(tibble) with all the columns from ET output
#' @description Read legacy ET file and output a tibble
#' @export
ReadET <- function(ET_file) {
  ET.lines <- read_lines(ET_file)
  if(length(ET.lines) == 0) {
    return(NA)
  } else {
    ET.df <- ET.lines %>%
      .[!str_detect(., "%")] %>%
      I(.) %>%
      read_fwf(fwf_cols(align.num = 10, POS = 10, AA = 10, rank = 10,
                        vari.n = 10, vari = 22, rho = 10, coverage = 10),
               col_types = "ndcddcdd")
    return(ET.df)
  }
}

#' Read fasta file
#'
#' @param fasta_file Path to fasta file
#' @param output_type vector or df
#' @return If vector is returned, then all sequences are stored in a vector, ids are
#' stored as the names. If a dataframe is returned, then ids are stored in id column and
#' sequences are stored in seq column.
#' @description Read fasta file
#' @export
ReadFasta <- function(fasta_file, output_type = c("vector", "df")) {
  all.seq <- Biostrings::readAAStringSet(fasta_file)
  output_type <- match.arg(output_type)
  if (output_type == "vector") {
    output <- as.character(all.seq)
    names(output) <- names(all.seq)
  } else {
    output <- data.frame(id = names(all.seq),
                         seq = as.character(all.seq),
                         stringsAsFactors = FALSE)
  }
  return(output)
}

