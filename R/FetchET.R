#' Fetch prestored ET file
#'
#' @param prot_id Protein id. Protein identification. Use full ENSP for human proteins, eg. ENSP00000388638.1.
#' Use locus tag for E. coli MG1655 proteins, eg. b4112.
#' @return A dataframe(tibble) with prestored ET. The ET file is reduced to only include POS, AA and coverage.
#' @description Fetch prestored ET file. Only supports human and E. coli MG1655 proteins
#' @export
FetchET <- function(prot_id) {
  # Check if there is ET for this protein
  if (!prot_id %in% id_map$prot_id) {
    stop("No prestored ET for this protein")
  }
  ET_storage <- id_map$storage[which(id_map$prot_id == prot_id)]
  ET_df <- read_rds(file.path(system.file("extdata", package = "EvoTrace"), ET_storage))
  ET_df <- ET_df$data[ET_df$prot_id == prot_id]%>%
    .[[1]] %>%
    mutate(POS = 1:n()) %>%
    select(POS, AA, coverage)
  return(ET_df)
}
