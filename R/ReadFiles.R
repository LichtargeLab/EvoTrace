#' Read legacy ET format
#'
#' @param ET_file Path to ET file
#' @return a dataframe(tibble) with all the columns from ET output
#' @description Read legacy ET file and output a tibble
#' @export
ReadET <- function(ET_file) {
  ET.df <- read_lines(ET_file) %>%
    str_replace_all(" +", " ") %>%
    .[!str_detect(., "%")] %>%
    read_delim(col_names = c("align.num", "POS", "AA", "rank","vari.n",
                             "vari", "rho", "coverage"),
               col_types = "ndcddcdd",
               delim = " ")
  return(ET.df)
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
  all.seq <- Biostrings::readAAStringSet(glue("{new.dir}/{gisaid.fasta}"))
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

