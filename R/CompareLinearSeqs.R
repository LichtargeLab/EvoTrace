#' Compare two protein linear sequences
#' @param seq1 Character string. Linear sequence of protein 1
#' @param seq2 Character string. Linear sequence of protein 2
#' @param pos.only If TRUE, then only matched positions between the pdb sequence and
#' the linear sequence are returned.
#' @param penalty Gap penalties used for the pairwise alignment.
#' @return A tibble that contains the matching information of the pdb sequence and
#' linear sequence.
#' @description Compared two linear sequences. The matching positions in the two
#' linear sequences are returned.
#' @export
CompareLinearSeqs <- function(seq1, seq2, pos.only = TRUE, penalty = c("blastp", "pyETv", "DNA")) {
  penalty <- match.arg(penalty)
  if(penalty == "blastp") {
    gap.open = 11
    gap.extend = 1
  } else if (penalty == "pyETV"){
    gap.open = 10
    gap.extend = 0.5
  } else {
    gap.open = 1
    gap.extend = 1
  }

  if (penalty == "DNA") {
    mat <- Biostrings::nucleotideSubstitutionMatrix(match = 1, mismatch = -2, baseOnly = TRUE)
    alignment <- Biostrings::pairwiseAlignment(pattern = seq1, subject = seq2,
                                               substitutionMatrix = mat,
                                               gapOpening = gap.open, gapExtension = gap.extend)
  } else {
    alignment <- Biostrings::pairwiseAlignment(pattern = seq1, subject = seq2,
                                               substitutionMatrix = "BLOSUM62",
                                               gapOpening = gap.open, gapExtension = gap.extend)
  }

  align.df <- c(Biostrings::alignedPattern(alignment),
                Biostrings::alignedSubject(alignment)) %>%
    as.character() %>%
    str_split(., pattern = "", simplify = TRUE) %>%
    t() %>%
    as.data.frame() %>%
    dplyr::rename(seq1.AA.align = V1, seq2.AA.align = V2) %>%
    mutate(align.POS = 1:n()) %>%
    mutate(seq1.gap = seq1.AA.align == "-",
           seq2.gap = seq2.AA.align == "-")

  align.df <- align.df %>%
    arrange(align.POS) %>%
    group_by(seq1.gap) %>%
    nest() %>%
    arrange(seq1.gap)
  align.df$data[[1]]$AA.POS.seq1 <- 1:nrow(align.df$data[[1]])
  align.df <- align.df %>%
    unnest(cols = c(data)) %>%
    ungroup() %>%
    arrange(align.POS)

  align.df <- align.df %>%
    arrange(align.POS) %>%
    group_by(seq2.gap) %>%
    nest() %>%
    arrange(seq2.gap)
  align.df$data[[1]]$AA.POS.seq2 <- 1:nrow(align.df$data[[1]])
  align.df <- align.df %>%
    unnest(cols = c(data)) %>%
    ungroup() %>%
    arrange(align.POS) %>%
    select(seq1.gap, seq2.gap, seq1.AA.align, seq2.AA.align, align.POS, AA.POS.seq1, AA.POS.seq2)

  if (pos.only == TRUE) {
    align.df <- select(align.df, AA.POS.seq1, AA.POS.seq2) %>%
      filter(!is.na(AA.POS.seq1), !is.na(AA.POS.seq2))
  }
  if (penalty == "DNA") {
    names(align.df) <- str_replace(names(align.df), "AA", "DNA")
  }
  return(align.df)
}
