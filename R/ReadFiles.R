#' Read legacy ET format
#'
#' @param ET_file Path to ET file
#' @param ET_format EA format or UET format.
#' @return a dataframe(tibble) with all the columns from ET output
#' @description Read legacy ET file and output a tibble. Support output
#' from EA calculation and UET server.
#' @export
ReadET <- function(ET_file, ET_format = c("EA", "UET")) {
  ET_format <- match.arg(ET_format)
  ET_lines <- readLines(ET_file)
  if(length(ET_lines) == 0) {
    return(NA)
  } else {
    empty_lines <- ET_lines == ""
    header_lines <- grepl("%", ET_lines)
    skip_lines <- empty_lines | header_lines
    entry_line_count <- sum(1 - skip_lines)
    skip_lines <- which(skip_lines == FALSE)[1]
    if (ET_format == "EA") {
      ET_df <- read.fwf(file = ET_file,
                        skip = skip_lines - 1,
                        col.names = c("align.num", "POS", "AA",
                                      "rank", "vari.n", "vari",
                                      "rho", "coverage"),
                        widths = c(10,10,10,10,10,22,10,10),
                        strip.white = TRUE, nrow = entry_line_count)
    } else {
      ET_df <- read.fwf(file = ET_file,
                        skip = skip_lines - 1,
                        col.names = c("align.num", "POS", "AA",
                                      "coverage", "vari.n", "vari",
                                      "rho"),
                        widths = c(10,10,10,10,10,22,10),
                        strip.white = TRUE, nrow = entry_line_count)
      ET_df$rank <- rank(ET_df$rho, ties.method = "max")
      ET_df <- ET_df[,c("align.num", "POS", "AA",
                        "rank", "vari.n", "vari",
                        "rho", "coverage")]
    }
    return(ET_df)
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

