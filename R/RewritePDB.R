#' Rewrite PDB file according to linear sequence
#'
#' @param pdb_file path to PDB file
#' @param chains A character vector. Chains that need to be rewritten.
#' @param linear_seq A character. The linear sequence.
#' @param fix_positions TRUE or FALSE. If TRUE, the positions in the PDB file are repleaced
#' with linear sequence positions.
#' @param keep_other_chains TRUE or FALSE. Whether other chains are kept.
#' @param output_file path to output PDB file
#' @return a message reporting the output file is written
#' @description The PDB files are often messy. The sequences positions might not match the linear
#' sequence positions. The pdb file might contain insertions, like His tag. This function reads the
#' sequence in PDB file and aligns it with the linear sequence. Any positions in the pdb file that
#' don't align to the linear sequence will be removed. The position number can be rewritten to match
#' linear positions. Water molecules are also removed. Only coordinates fields are kept in the output
#' pdb file.
#' @export

RewritePDB <- function(pdb_file, chains, linear_seq, fix_positions = TRUE,
                       keep_other_chains = TRUE, output_file) {
  pdb <- read_lines(pdb_file)
  pdb_row_type <- str_sub(pdb, 1, 6)

  atom_data <- pdb[str_detect(pdb_row_type, "ATOM  |HETATM|TER   ")] %>%
    I() %>%
    read_fwf(., col_positions = fwf_widths(widths = c(6, 5, 1, 4, 1, 3, 1, 1, 4, 1, 3, 8, 8, 8, 6, 6, 6, 4, 2, 2)),
             trim_ws = FALSE,
             show_col_types = FALSE) %>%
    # Convert AA positions into numeric
    mutate(X9 = as.numeric(X9))

  other_chains <- atom_data %>%
    filter(!X8 %in% chains)

  rewrite_chains <- atom_data %>%
    filter(X8 %in% chains)

  seq_map_list <- lapply(chains, CompareSeqs, pdb_file = pdb_file,
                         seq = linear_seq)
  names(seq_map_list) <- chains

  FixChain <- function(rewrite_chains, chain_to_fix, fix_positions = fix_positions,
                       seq_map) {
    ignore_chains <- filter(rewrite_chains, X8 != chain_to_fix)
    working_chain <- filter(rewrite_chains, X8 == chain_to_fix) %>%
      filter(X9 %in% seq_map$AA.POS.pdb)
    if(fix_positions == TRUE) {
      working_chain <- working_chain %>%
        left_join(seq_map, by = c("X9" = "AA.POS.pdb")) %>%
        mutate(X9 = AA.POS.seq) %>%
        filter(!is.na(X9)) %>%
        select(-AA.POS.seq)
    }
    output <- bind_rows(working_chain, ignore_chains) %>%
      arrange(X2)
    return(output)
  }

  for (i in 1:length(chains)) {
    rewrite_chains <- FixChain(rewrite_chains, chain_to_fix = chains[i], fix_positions = fix_positions,
                               seq_map = seq_map_list[[chains[i]]])
  }

  if(keep_other_chains == TRUE) {
    output <- bind_rows(rewrite_chains, other_chains) %>%
      arrange(X2) %>%
      filter(X6 != "HOH")
  } else {
    output <- rewrite_chains
  }
  # Pad white space to the left of residue number to make it to width of 4
  output <- output %>%
    mutate(X9 = str_pad(X9, width = 4, side = "left", pad = " ")) %>%
    as.data.frame()

  gdata::write.fwf(output, file = output_file, append = FALSE, width = c(6, 5, 1, 4, 1, 3, 1, 1, 4, 1, 3, 8, 8, 8, 6, 6, 6, 4, 2, 2),
                   colnames = FALSE, justify="left", sep = "")
  write_lines("END", output_file, append = TRUE)
  return(paste0("write ", output_file))
}
