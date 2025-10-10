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
#' @description Read files with fasta like format. The function allows most characters (no space)
#' in the sequence and preserves upper/lower cases.
#' @export
ReadFasta <- function(fasta_file, output_type = c("vector", "df")) {
  output_type <- match.arg(output_type)

  fasta_lines <- readLines(fasta_file, warn = FALSE)
  if (length(fasta_lines) == 0) stop("File is empty.")

  fasta_lines <- trimws(fasta_lines)
  fasta_lines <- fasta_lines[nchar(fasta_lines) > 0]

  header_positions <- which(substr(fasta_lines, 1, 1) == ">")
  if (length(header_positions) == 0) stop("No FASTA headers ('>') found.")
  header_positions <- c(header_positions, length(fasta_lines) + 1)

  id_vec  <- character(length(header_positions) - 1)
  seq_vec <- character(length(header_positions) - 1)

  for (i in seq_len(length(header_positions) - 1)) {
    start_idx <- header_positions[i]
    end_idx   <- header_positions[i + 1] - 1

    seq_id <- sub("^>", "", fasta_lines[start_idx])  # remove '>' before ID

    seq_lines <- fasta_lines[(start_idx + 1):end_idx]
    seq_combined <- gsub("\\s+", "", paste(seq_lines, collapse = ""))

    id_vec[i]  <- seq_id
    seq_vec[i] <- seq_combined
  }

  if (output_type == "vector") {
    output <- seq_vec
    names(output) <- id_vec
  } else {
    output <- data.frame(id = id_vec,
                         seq = seq_vec,
                         stringsAsFactors = FALSE)
  }
  return(output)
}

#' Read Clustal MSA file
#'
#' @param clustal_path Path to MSA file
#' @return A tibble format of the MSA. Each sequence is in one column, each cell is a residue.
#' The position of each residue in each sequence is also returned. NA is returned in the position
#' columns for gaps.
#' @description Convert Clustal MSA file into a tibble
#' @export
ReadClustalMSA <- function(clustal_path) {
  # Read each line as string then remove empty lines, skip header
  clustal_raw <- read_lines(clustal_path, skip = 1)
  clustal_raw <- clustal_raw[clustal_raw != ""]

  # The trailing positions are separted from the sequences with tab in the example.
  # The delim could also be space. Remove those
  clustal_raw <- str_remove(clustal_raw, "(\t| )[:digit:]+")

  # The remaining part is a fix width structure
  # Add id to the symbol lines
  str_sub(clustal_raw[str_starts(clustal_raw, "      ")], 1, 3) <- "symbol"

  # Identify the position of the last space
  last_space <- str_locate(clustal_raw, " +")[1, "end"]

  # Read fixed width lines into dataframe
  clustal_df <- read_fwf(I(clustal_raw),
                         fwf_cols(id = c(1, last_space), seq = c(last_space + 1, nchar(clustal_raw[1]))),
                         col_types = "cc",
                         trim_ws = FALSE) %>%
    mutate(id = str_remove_all(id, " ")) %>%
    # Collapse sequences for each id
    group_by(id) %>%
    summarize(seq = paste0(seq, collapse = ""),
              .groups = "drop")

  # Get the squence names
  seq_names <- clustal_df$id[1:(nrow(clustal_df) - 1)]

  # Convert to vector
  clustal_vec <- clustal_df$seq
  names(clustal_vec) <- clustal_df$id

  # Split vector by charachters, and convert to dataframe
  output <- str_split(clustal_vec, pattern = "", simplify = TRUE) %>%
    t() %>%
    as.data.frame()
  names(output) <- names(clustal_vec)

  # This function takes a broken down sequence and return the position for
  # each residue. Gaps return NAs.
  GetPosition <- function(vec) {
    # mark non-gap positions
    vec_bool <- vec != "-"
    # residue position +1 if for non-gap
    vec_out <- cumsum(vec_bool)
    vec_out[vec_bool == FALSE] <- NA
    return(vec_out)
  }

  output <- output %>%
    # Get the residue positions for each sequence
    mutate(across(all_of(seq_names), GetPosition, .names = "{.col}_POS")) %>%
    # Rearrange the column orders
    relocate(ends_with("POS")) %>%
    relocate(starts_with(seq_names)) %>%
    # Replace space as NA in the symbol column.
    mutate(symbol = ifelse(symbol == " ",
                           NA,
                           symbol))
  return(output)
}
