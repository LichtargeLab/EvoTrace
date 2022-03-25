#' Extract coordinates from PDB file
#'
#' @param pdb_file path to PDB file
#' @param chain Character string. Chains that needs to be extracted
#' @param CA_only Logic. Whether only C alpha are extracted
#' @return A tibble that contains chain, POS (AA position), AA, ATOM, ATOM_id, x, y, z
#' @description Extract coordinates from PDB file. Output will be a tibble
#' that contains protein chain, AA position, AA type (one letter), atom type,
#' atom id and coordinates (x, y, z)
#' @export
GetCoordinates <- function(pdb_file, chain, CA_only = TRUE) {
  pdb <- read_lines(pdb_file)
  pdb_row_type <- str_sub(pdb, 1, 6)
  output <- pdb[str_detect(pdb_row_type, "ATOM  |HETATM")] %>%
    I(.) %>%
    read_fwf(., col_positions = fwf_widths(widths = c(6, 5, 1, 4, 1, 3, 1, 1, 4, 1, 3, 8, 8, 8, 6, 6, 6, 4, 2, 2)),
             col_types = "cdccccccdccdddddcccc",
             trim_ws = TRUE) %>%
    filter(X8 %in% chain) %>%
    filter(X6 %in% c("HIS", "PRO", "GLU", "THR", "LEU", "VAL", "LYS", "ASP", "ALA",
                     "GLN", "GLY", "ARG", "TYR", "ILE", "ASN", "SER", "PHE", "MET",
                     "CYS", "TRP", "MSE"))
  if(CA_only == TRUE) {
    output <- filter(output, X4 == "CA")
  }
  output <- output %>%
    select(chain = X8, POS = X9, AA = X6, ATOM = X4, ATOM_id = X2, x = X12, y = X13, z = X14, X15) %>%
    group_by(chain, POS, AA, ATOM) %>%
    # If there are 2 conformation for a residue, the dominant one is picked
    arrange(X15) %>%
    filter(row_number()==1) %>%
    select(-X15) %>%
    ungroup() %>%
    mutate(AA = AA3to1(AA))
  return(output)
}


#' Extract coordinates from PDB file without AA filter
#'
#' @param pdb_file path to PDB file
#' @param chain Character string. Chains that needs to be extracted
#' @return A tibble that contains chain, POS (AA position), AA, ATOM, ATOM_id, x, y, z
#' @description Extract coordinates from PDB file. Output will be a tibble
#' that contains protein chain, AA position, AA type (one letter), atom type,
#' atom id and coordinates (x, y, z)
GetCoordinates_no_filter <- function(pdb_file, chain) {
  pdb <- read_lines(pdb_file)
  pdb_row_type <- str_sub(pdb, 1, 6)
  output <- pdb[str_detect(pdb_row_type, "ATOM  |HETATM")] %>%
    read_fwf(., col_positions = fwf_widths(widths = c(6, 5, 1, 4, 1, 3, 1, 1, 4, 1, 3, 8, 8, 8, 6, 6, 6, 4, 2, 2)),
             trim_ws = TRUE) %>%
    filter(X8 %in% chain)
  output <- output %>%
    select(chain = X8, POS = X9, MOLECUE = X6, ATOM = X4, ATOM_id = X2, x = X12, y = X13, z = X14)
  return(output)
}

