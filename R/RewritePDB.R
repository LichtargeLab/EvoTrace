#' Rewrite PDB file according to linear sequence
#'
#' @param pdb_file path to PDB file
#' @param chain A character. Chain that needs to be rewritten
#' @param linear_seq A character. The linear sequence.
#' @param fix_positions TRUE or FALSE. If TRUE, the positions in the PDB file are repleaced
#' with linear sequence positions.
#' @param keep_other_chains TRUE or FALSE. Whether other chains are kept.
#' @param output_file path to output PDB file
#' @return a messege reporting the output file is written
#' @description The PDB files are often messy. The sequences positions might not match the linear
#' sequence positions. The pdb file might contain insertions, like His tag. This function reads the
#' squence in PDB file and aligns it with the linear sequence. Any positions in the pdb file that
#' don't align to the linear sequence will be removed. The position number can be rewritten to match
#' linear positions. Water molecules are also removed. Only coordinates fields are kept in the output
#' pdb file.
#' @export

RewritePDB <- function(pdb_file, chain, linear_seq, fix_positions = TRUE,
                       keep_other_chains = TRUE, output_file) {
  pdb <- read_lines(pdb_file)
  pdb_row_type <- str_sub(pdb, 1, 6)

  seq_map <- CompareSeqs(pdb_file = pdb_file, chain = chain, seq = linear_seq, pos.only = TRUE)

  atom_data <- pdb[str_detect(pdb_row_type, "ATOM  |HETATM|TER   ")] %>%
    read_fwf(., col_positions = fwf_widths(widths = c(6, 5, 1, 4, 1, 3, 1, 1, 4, 1, 3, 8, 8, 8, 6, 6, 6, 4, 2, 2)),
             trim_ws = TRUE)

  other_chains <- atom_data %>%
    filter(X8 != chain)

  mapped_chain <- atom_data %>%
    filter(X8 == chain) %>%
    filter(X9 %in% seq_map$AA.POS.pdb)

  if(fix_positions == TRUE) {
    mapped_chain <- mapped_chain %>%
      left_join(seq_map, by = c("X9" = "AA.POS.pdb")) %>%
      mutate(X9 = AA.POS.seq) %>%
      filter(!is.na(X9)) %>%
      select(-AA.POS.seq)
  }

  if(keep_other_chains == TRUE) {
    output <- bind_rows(mapped_chain, other_chains) %>%
      arrange(X2) %>%
      filter(X6 != "HOH") %>%
      as.data.frame()
  } else {
    output <- mapped_chain %>%
      as.data.frame()
  }


  gdata::write.fwf(output, file = output_file, append = FALSE, width = c(6, 5, 1, 4, 1, 3, 1, 1, 4, 1, 3, 8, 8, 8, 6, 6, 6, 4, 2, 2),
                   colnames = FALSE, justify="left", sep = "")
  write_lines("END", output_file, append = TRUE)
  return(paste0("write ", output_file))
}
