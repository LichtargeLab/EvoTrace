#' Get per-residue confidence metric (pLDDT) from AlphaFold structures
#'
#' @param pdb_file path to AlphaFold PDB file
#' @param chain Character vector. Chains of interest. Usually AlphaFold structures are stored
#' in chain A.
#' @return A tibble that contains the chains, amino acid positions and pLDDT scores.
#' @description Amino acid positions and pLDDT scores are extracted from AlphaFold pdb
#' file. pLDDT values are stored in b factor column. So this function can also be used
#' to extract b factor from experimental pdb files. Use cation when applying to experimental
#' data due to isoforms have the same AA position.
#' @export
GetpLDDT <- function(pdb_file, chain = "A") {
  pdb <- read_lines(pdb_file)
  pdb_row_type <- stringr::str_sub(pdb, 1, 6)
  output <- pdb[stringr::str_detect(pdb_row_type, "ATOM  |HETATM")] %>%
    I() %>%
    readr::read_fwf(., col_positions = readr::fwf_widths(widths = c(6, 5, 1, 4, 1, 3, 1, 1, 4, 1, 3, 8, 8, 8, 6, 6, 6, 4, 2, 2)),
                    trim_ws = TRUE,
                    show_col_types = FALSE) %>%
    dplyr::filter(X8 %in% chain) %>%
    dplyr::filter(X6 %in% c("HIS", "PRO", "GLU", "THR", "LEU", "VAL", "LYS", "ASP", "ALA",
                            "GLN", "GLY", "ARG", "TYR", "ILE", "ASN", "SER", "PHE", "MET",
                            "CYS", "TRP", "MSE")) %>%
    dplyr::filter(X4 == "CA") %>%
    dplyr::select(chain = X8, POS = X9, pLDDT = X16) %>%
    dplyr::arrange(chain, POS)
  return(output)
}
