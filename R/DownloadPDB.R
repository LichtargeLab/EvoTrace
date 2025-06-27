#' Download PDB file from rcsb
#'
#' @param pdb The pdb id that needs to be downloaded
#' @param out_dir The directory that the pdb will be saved
#' @param pdb_format Use "pdb" if the input is traditional pdb. Use "pdbx" for PDBx/mmCIF.
#' @return a message reporting the pdb is downloaded
#' @description This function intake a pdb id and download the pdb file to
#' the path location as {pdb.id}.pdb
#' @export


DownloadPDB <- function(pdb, out_dir = ".", pdb_format = c("pdb", "pdbx")) {
  pdb_format <- match.arg(pdb_format)
  if (pdb_format == "pdb") {
    pdb_url <- paste0("https://files.rcsb.org/download/", pdb, ".pdb")
    save_location <- file.path(out_dir, paste0(pdb, ".pdb"))
  } else {
    pdb_url <- paste0("https://files.rcsb.org/download/", pdb, ".cif")
    save_location <- file.path(out_dir, paste0(pdb, ".cif"))
  }
  download.file(pdb_url, save_location)
  return(paste0(pdb, " saved at ", save_location))
}
