#' Rewrite PDB file according to linear sequence
#'
#' @param pdb The pdb id that needs to be downloaded
#' @param out_dir The directory that the pdb will be saved
#' @return a message reporting the pdb is downloaded
#' @description This function intake a pdb id and download the pdb file to
#' the path location as {pdb.id}.pdb
#' @export


DownloadPDB <- function(pdb, out_dir = ".") {
  pdb_url <- paste0("https://files.rcsb.org/download/", pdb, ".pdb")
  save_location <- file.path(out_dir, paste0(pdb, ".pdb"))
  download.file(pdb_url, save_location)
  return(paste0(pdb, " saved at ", save_location))
}
