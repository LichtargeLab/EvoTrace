#' Get residues at protein interface of two chains from a PDB file
#'
#' @param pdb_file path to pdb file
#' @param chain1 the first chain
#' @param chain2 the second chain
#' @param dist_cutoff the distance cutoff that consider two atoms are adjacent
#' @param output_type "paired" or "stacked"
#' @return a dataframe(tibble) contains all residues at the protein interface
#' @description Residues are at the interface when any of their heavy atoms are
#' within the cutoff range of a heavy atom from the other chain. The paired output
#' contains columns chain_i, chain_j, POS_i, POS_j, and dist. The stacked output
#' contains columns chain and POS.
#' @export
GetInterface <- function(pdb_file, chain1, chain2, dist_cutoff = 4,
                         output_type = c("paired", "stacked")) {
  coord <- GetCoordinates_no_filter(pdb_file, chain = c(chain1, chain2)) %>%
    mutate(POS = paste0(chain, "_", POS))
  df1 <- coord %>%
    filter(chain == chain1) %>%
    select(POS_i = POS, ATOM_i = ATOM, x_i = x, y_i = y, z_i = z)
  df2 <- coord %>%
    filter(chain == chain2) %>%
    select(POS_j = POS, ATOM_j = ATOM, x_j = x, y_j = y, z_j = z)
  interface <- expand_grid(df1, df2, .name_repair = "minimal") %>%
    mutate(dist = sqrt((x_i - x_j)^2 + (y_i - y_j)^2 + (z_i - z_j)^2)) %>%
    group_by(POS_i, POS_j) %>%
    summarize(dist = min(dist)) %>%
    ungroup() %>%
    separate(POS_i, into = c("chain_i", "POS_i")) %>%
    separate(POS_j, into = c("chain_j", "POS_j")) %>%
    mutate(POS_i = as.numeric(POS_i),
           POS_j = as.numeric(POS_j)) %>%
    filter(chain_i == chain1,
           chain_j == chain2,
           dist <= dist_cutoff)
  if (output_type == "paired") {
    output <- interface
  } else {
    output <- bind_rows(select(interface, chain = chain_i, POS = POS_i),
                        select(interface, chain = chain_j, POS = POS_j)) %>%
      unique() %>%
      arrange(chain, POS)
  }
  return(output)
}
