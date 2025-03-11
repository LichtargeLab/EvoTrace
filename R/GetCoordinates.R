#' Extract coordinates from PDB file
#'
#' @param pdb_file path to PDB file
#' @param chain Character string. Chains that needs to be extracted
#' @param CA_only Logic. Whether only C alpha are extracted
#' @param remove_insertions Logic. Some residues are labeled as insertions in the pdb
#' file and they have the same POS index as the residue that they insert after. If True,
#' these insertion residues are removed.
#' @param pdb_format Use "pdb" if the input is traditional pdb. Use "pdbx" for PDBx/mmCIF.
#' @return A tibble that contains chain, POS (AA position), AA, ATOM, ATOM_id, x, y, z
#' @description Extract coordinates from PDB file. Output will be a tibble
#' that contains protein chain, AA position, AA type (one letter), atom type,
#' atom id and coordinates (x, y, z)
#' @export
GetCoordinates <- function(pdb_file, chain, CA_only = TRUE, remove_insertions = TRUE,
                           pdb_format = c("pdb", "pdbx")) {
  pdb_format <- match.arg(pdb_format)
  if (pdb_format == "pdb") {
    output <- GetCoordinates_pdb(pdb_file = pdb_file,
                                 chain = chain,
                                 CA_only = CA_only,
                                 remove_insertions = remove_insertions)
  } else {
    output <- GetCoordinates_cif(pdb_file = pdb_file,
                                 chain = chain,
                                 CA_only = CA_only,
                                 remove_insertions = remove_insertions)
  }
  return(output)
}

#' Extract coordinates from PDB file
#'
#' @param pdb_file path to PDB file
#' @param chain Character string. Chains that needs to be extracted
#' @param CA_only Logic. Whether only C alpha are extracted
#' @param remove_insertions Logic. Some residues are labeled as insertions in the pdb
#' file and they have the same POS index as the residue that they insert after. If True,
#' these insertion residues are removed.
#' @return A tibble that contains chain, POS (AA position), AA, ATOM, ATOM_id, x, y, z
#' @description Extract coordinates from PDB file. Output will be a tibble
#' that contains protein chain, AA position, AA type (one letter), atom type,
#' atom id and coordinates (x, y, z)
GetCoordinates_pdb <- function(pdb_file, chain, CA_only = TRUE, remove_insertions = TRUE) {
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
                     "CYS", "TRP", "MSE")) %>%
    filter(X1 == "ATOM"|X6 == "MSE")
  if(CA_only == TRUE) {
    output <- filter(output, X4 == "CA")
  }
  if(remove_insertions == TRUE) {
    output <- output %>%
      filter(is.na(X10)) %>%
      select(chain = X8, POS = X9, AA = X6, ATOM = X4, ATOM_id = X2, x = X12, y = X13, z = X14, X15) %>%
      group_by(chain, POS, AA, ATOM) %>%
      # If there are 2 conformation for a residue, the dominant one is picked
      arrange(X15) %>%
      filter(row_number()==1) %>%
      select(-X15) %>%
      ungroup() %>%
      mutate(AA = AA3to1(AA)) %>%
      arrange(ATOM_id)
  } else {
    output <- output %>%
      select(chain = X8, POS = X9, insertion_code = X10, AA = X6, ATOM = X4, ATOM_id = X2, x = X12, y = X13, z = X14, X15) %>%
      group_by(chain, POS, insertion_code, AA, ATOM) %>%
      # If there are 2 conformation for a residue, the dominant one is picked
      arrange(X15) %>%
      filter(row_number()==1) %>%
      select(-X15) %>%
      ungroup() %>%
      mutate(AA = AA3to1(AA)) %>%
      arrange(ATOM_id)
  }

  return(output)
}

#' Extract coordinates from PDBx/mmcif file
#'
#' @param pdb_file path to PDBx/mmcif file
#' @param chain Character string. Chains that needs to be extracted
#' @param CA_only Logic. Whether only C alpha are extracted
#' @param remove_insertions Logic. Some residues are labeled as insertions in the pdb
#' file and they have the same POS index as the residue that they insert after. If True,
#' these insertion residues are removed.
#' @return A tibble that contains chain, POS (AA position), AA, ATOM, ATOM_id, x, y, z
#' @description Extract coordinates from PDB file. Output will be a tibble
#' that contains protein chain, AA position, AA type (one letter), atom type,
#' atom id and coordinates (x, y, z)
GetCoordinates_cif <- function(pdb_file, chain, CA_only = TRUE, remove_insertions = TRUE) {
  pdb <- read_lines(pdb_file)
  info_lines <- which(str_sub(pdb, 1, 1) %in% c("_", "#"))
  info_lines <- c(info_lines, length(pdb))
  header_start <- which(str_starts(pdb, "_atom_site.group_PDB"))
  header_end <- header_start
  while((header_end + 1) %in% info_lines) {
    header_end <- header_end + 1
  }

  header <- pdb[header_start:header_end] %>%
    str_sub(., 12, -1) %>%
    # trim white spaces for the column names
    str_trim("both")

  coord_start <- header_end + 1
  coord_end <- min(info_lines[info_lines > coord_start]) - 1


  output <- pdb[coord_start:coord_end] %>%
    # replace tandem spaces to one, so the lines
    # can be easily read from read_delim
    str_replace_all(., " +", " ") %>%
    str_trim(., "right") %>%
    I(.) %>%
    read_delim(., delim = " ", col_names = header,
               col_types = cols(.default = col_character())) %>%
    filter(auth_asym_id %in% chain)

  # When subseting pdb in pymol, auth_comp_id, auth_seq_id,
  # auth_atom_id columns get removed
  if (!"auth_comp_id" %in% names(output)) {
    output[["auth_comp_id"]] <- output[["label_comp_id"]]
  }
  if (!"auth_seq_id" %in% names(output)) {
    output[["auth_seq_id"]] <- output[["label_seq_id"]]
  }
  if (!"auth_atom_id" %in% names(output)) {
    output[["auth_atom_id"]] <- output[["label_atom_id"]]
  }

  output <- output %>%
    filter(auth_comp_id %in% c("HIS", "PRO", "GLU", "THR", "LEU", "VAL", "LYS", "ASP", "ALA",
                               "GLN", "GLY", "ARG", "TYR", "ILE", "ASN", "SER", "PHE", "MET",
                               "CYS", "TRP", "MSE")) %>%
    filter(group_PDB == "ATOM"|auth_comp_id == "MSE")

  if(CA_only == TRUE) {
    output <- filter(output, auth_atom_id == "CA")
  }

  # check if insertion code column is used in the pdbx file
  if((remove_insertions == FALSE) & ("pdbx_PDB_ins_code" %in% header)) {
    output <- output %>%
      select(chain = auth_asym_id,
             POS = auth_seq_id,
             insertion_code = pdbx_PDB_ins_code,
             AA = auth_comp_id,
             ATOM = auth_atom_id,
             ATOM_id = id, # unique num id for each atom
             x = Cartn_x,
             y = Cartn_y,
             z = Cartn_z,
             element = type_symbol) %>%
      group_by(chain, POS, insertion_code, AA, ATOM) %>%
      # If there are 2 conformation for a residue, the dominant one is picked
      arrange(element) %>%
      slice(1) %>%
      select(-element) %>%
      ungroup() %>%
      mutate(AA = AA3to1(AA)) %>%
      mutate(across(c(POS, x, y, z, ATOM_id), as.numeric)) %>%
      arrange(ATOM_id)
  } else {
    output <- output %>%
      filter(pdbx_PDB_ins_code == "?") %>%
      select(chain = auth_asym_id,
             POS = auth_seq_id,
             AA = auth_comp_id,
             ATOM = auth_atom_id,
             ATOM_id = id, # unique num id for each atom
             x = Cartn_x,
             y = Cartn_y,
             z = Cartn_z,
             element = type_symbol) %>%
      group_by(chain, POS, AA, ATOM) %>%
      # If there are 2 conformation for a residue, the dominant one is picked
      arrange(element) %>%
      slice(1) %>%
      select(-element) %>%
      ungroup() %>%
      mutate(AA = AA3to1(AA)) %>%
      mutate(across(c(POS, x, y, z, ATOM_id), as.numeric)) %>%
      arrange(ATOM_id)
  }
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
    I(.) %>%
    read_fwf(., col_positions = fwf_widths(widths = c(6, 5, 1, 4, 1, 3, 1, 1, 4, 1, 3, 8, 8, 8, 6, 6, 6, 4, 2, 2)),
             trim_ws = TRUE) %>%
    filter(X8 %in% chain)
  output <- output %>%
    select(chain = X8, POS = X9, MOLECUE = X6, ATOM = X4, ATOM_id = X2, x = X12, y = X13, z = X14)
  return(output)
}

