#' Compare protein sequence in a pdb file with linear seq
#'
#' @param pdb_file path to PDB file
#' @param chain Character string. Chains of interest
#' @param seq Character string. Linear sequence of the protein
#' @param pos.only If TRUE, then only matched positions between the pdb sequence and
#' the linear sequence are returned.
#' @return A tibble that contains the matching information of the pdb sequence and
#' linear sequence.
#' @description Sequence from pdb_file is extracted and compared (aligned) to linear
#' sequence. The matching positions in the pdb entry and in the linear sequence are
#' returned.
#' @export
CompareSeqs <- function(pdb_file, chain, seq, pos.only = TRUE) {
  pdb.df <- GetCoordinates(pdb_file = pdb_file, chain = chain, CA_only = TRUE) %>%
    select(AA.pdb = AA, AA.POS.pdb = POS)
  pdb.str <- pdb.df$AA.pdb %>%
    paste0(., collapse = "")

  alignment <- Biostrings::pairwiseAlignment(pattern = pdb.str, subject = seq,
                                             substitutionMatrix = "BLOSUM62",
                                             gapOpening = 11, gapExtension = 1)

  align.df <- c(Biostrings::alignedPattern(alignment),
                Biostrings::alignedSubject(alignment)) %>%
    as.character() %>%
    str_split(., pattern = "", simplify = TRUE) %>%
    t() %>%
    as.data.frame() %>%
    dplyr::rename(pdb.AA.align = V1, seq.AA.align = V2) %>%
    mutate(align.POS = 1:n()) %>%
    mutate(pdb.gap = pdb.AA.align == "-",
           seq.gap = seq.AA.align == "-")

  align.df <- align.df %>%
    arrange(align.POS) %>%
    group_by(pdb.gap) %>%
    nest() %>%
    arrange(pdb.gap)
  align.df$data[[1]]$AA.POS.pdb <- pdb.df$AA.POS.pdb
  align.df <- align.df %>%
    unnest(cols = c(data)) %>%
    ungroup()

  align.df <- align.df %>%
    arrange(align.POS) %>%
    group_by(seq.gap) %>%
    nest() %>%
    arrange(seq.gap)
  align.df$data[[1]]$AA.POS.seq <- 1:nrow(align.df$data[[1]])
  align.df <- align.df %>%
    unnest(cols = c(data)) %>%
    ungroup()
  if (pos.only == TRUE) {
    align.df <- select(align.df, AA.POS.pdb, AA.POS.seq) %>%
      filter(!is.na(AA.POS.pdb), !is.na(AA.POS.seq))
  }
  return(align.df)
}
