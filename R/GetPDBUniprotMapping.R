#' Map PDB residue positions to UniProt positions
#'
#' @param pdb_id A character scalar with the PDB identifier to query.
#' @return A tibble with residue-level PDBe, PDB, and UniProt mappings. The
#' returned columns are `POS_PDBe`, `AA_PDBe`, `CHAIN_PDB`, `POS_PDB`,
#' `AA_PDB`, `Uniprot_ID`, `POS_Uniprot`, and `AA_Uniprot`.
#' @description Downloads the SIFTS XML mapping file for a PDB structure from
#' the EBI FTP server and extracts residue-level mappings among PDBe residue
#' records, PDB chain residue numbers, and UniProt residue numbers. PDBe and
#' PDB amino acid residue names are converted from three-letter to one-letter
#' codes using `AA3to1`.
#' @export
#' @examples
#' \dontrun{
#' GetPDBUniprotMapping("1a2b")
#' }
GetPDBUniprotMapping <- function(pdb_id) {
  if (!is.character(pdb_id) || length(pdb_id) != 1 || is.na(pdb_id)) {
    stop("pdb_id must be a non-empty character scalar")
  }

  pdb_id <- tolower(trimws(pdb_id))
  if (!nzchar(pdb_id)) {
    stop("pdb_id must be a non-empty character scalar")
  }

  url <- glue::glue("https://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/{pdb_id}.xml.gz")
  temp_file <- tempfile(fileext = ".xml.gz")
  on.exit(unlink(temp_file), add = TRUE)

  download_status <- tryCatch(
    utils::download.file(url, destfile = temp_file, mode = "wb", quiet = TRUE),
    warning = function(w) {
      stop(
        glue::glue("Failed to download SIFTS XML for PDB ID '{pdb_id}': {conditionMessage(w)}"),
        call. = FALSE
      )
    },
    error = function(e) {
      stop(
        glue::glue("Failed to download SIFTS XML for PDB ID '{pdb_id}': {conditionMessage(e)}"),
        call. = FALSE
      )
    }
  )

  if (!isTRUE(download_status == 0) ||
      !file.exists(temp_file) ||
      file.info(temp_file)$size == 0) {
    stop(
      glue::glue("Failed to download SIFTS XML for PDB ID '{pdb_id}' from {url}"),
      call. = FALSE
    )
  }

  xml_data <- tryCatch(
    xml2::read_xml(temp_file),
    error = function(e) {
      stop(
        glue::glue("Failed to parse SIFTS XML for PDB ID '{pdb_id}': {conditionMessage(e)}"),
        call. = FALSE
      )
    }
  )

  xml_data <- xml2::xml_ns_strip(xml_data)

  residues <- xml2::xml_find_all(xml_data, ".//residue[@dbSource='PDBe']")

  PDB_nodes <- xml2::xml_find_first(
    residues,
    "./crossRefDb[@dbSource='PDB']"
  )

  Uniprot_nodes <- xml2::xml_find_first(
    residues,
    "./crossRefDb[@dbSource='UniProt']"
  )

  output <- tibble::tibble(
    POS_PDBe = xml2::xml_attr(residues, "dbResNum"),
    AA_PDBe = xml2::xml_attr(residues, "dbResName"),

    CHAIN_PDB = xml2::xml_attr(PDB_nodes, "dbChainId"),
    POS_PDB = xml2::xml_attr(PDB_nodes, "dbResNum"),
    AA_PDB = xml2::xml_attr(PDB_nodes, "dbResName"),

    Uniprot_ID = xml2::xml_attr(Uniprot_nodes, "dbAccessionId"),
    POS_Uniprot = xml2::xml_attr(Uniprot_nodes, "dbResNum"),
    AA_Uniprot = xml2::xml_attr(Uniprot_nodes, "dbResName")
  )

  output <- dplyr::mutate(
    output,
    AA_PDBe = AA3to1(AA_PDBe),
    POS_PDBe = suppressWarnings(as.integer(POS_PDBe)),
    POS_PDB = suppressWarnings(as.integer(POS_PDB)),
    POS_Uniprot = suppressWarnings(as.integer(POS_Uniprot)),
    AA_PDB = AA3to1(AA_PDB)
  )

  return(output)
}
