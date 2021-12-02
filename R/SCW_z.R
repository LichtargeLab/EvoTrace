#' Compute the background stats for SCW z score calculation
#'
#' @param pdb_file path to pdb file
#' @param chain the protein chain SCW z score is calculated on
#' @param dist_cutoff the distance cutoff that consider two atoms are adjacent
#' @param resi a vector of residues numbers that are included in SCW bg calculation.
#' If NULL, then all residues in that chain is used.
#' @return a list with parameters that are needed for SCW z scores calculation
#' @description Calculates the required stats for w_hat, w^2_hat calculation in
#' SCW z scores. Residues are consider adjacent when any of their heavy atoms are
#' within the cutoff range. The output is feed into ComputeSCWzscore function for
#' z score calculation.
#' @export
GetSCWBackgound <- function(pdb_file, chain, dist_cutoff = 4, resi = NULL) {
  # This computes class counts used in calculation of
  # w_hat, w^2_hat in SCW z score.
  # It will produce values for both unbiased (1) and biased (j-i) approach
  # Args:
  #   pdb.file: path to pdb file
  #   chain: chain of interest
  #   dist_cutoff: the cutoff to determine adjacency of two residues in the structure
  # output:
  #   list contains the distance df and the class count df
  cord <- GetCoordinates(pdb_file, chain, CA_only = FALSE)
  if(!is.null(resi)) {
    cord <- cord %>%
      filter(POS %in% resi)
  }
  dist <- GetDistanceMatrix(cord, output_type = "df") %>%
    mutate(A = ifelse(dist < dist_cutoff, 1, 0))

  dist_filt <- dist %>%
    filter(POS_i < POS_j) %>%
    filter(A == 1)

  AssignType <- function(i,j,k,l) {
    output <- rep("B", length(i))
    output[which((i == k & j == l)|(i == l & j == k))] <- "A"
    output[which((i != k) & (i != l) & (j != k) & (j != l))] <- "C"
    return(output)
  }

  w2_class <- select(dist_filt, i = POS_i, j = POS_j) %>%
    crossing(., select(dist_filt, k = POS_i, l = POS_j)) %>%
    mutate(bias_term = (j-i)*(l-k),
           unbias_term = 1) %>%
    mutate(type = AssignType(i,j,k,l))

  w_ave_df <- dist_filt %>%
    mutate(bias_term = POS_j - POS_i,
           unbias_term = 1)

  class_count <- w2_class %>%
    group_by(type) %>%
    summarize(bias = sum(bias_term),
              unbias = sum(unbias_term)) %>%
    ungroup() %>%
    add_row(type = "SS_hat", bias = sum(w_ave_df$bias_term), unbias = sum(w_ave_df$unbias_term))
  pdb_length = length(unique(dist$POS_i))

  return(list("dist_filt" = dist_filt, "class_count" = class_count, "pdb_length" = pdb_length))
}



#' Compute SCW z scores
#'
#' @param SCW_background A list with parameters required for SCW z score calculation.
#' This is obtained through GetSCWBackground function.
#' @param resi numeric vector. The selection of residue positions. Residue positions
#' should be based on pdb structure not the sequence positions
#' @param output_df Whether a df should be returned instead of a vector.
#' @return a vector or df with biased and unbiased SCW z scores.
#' @description This function calculates the biased and unbiased SCW z scores
#' for a given selection of residues in the structure.
#' @export
ComputeSCWzscore <- function(SCW_background, resi, output_df = FALSE) {
  # This function calculates the biased and unbiased SCW z scores for a given
  # selection of residues in the structure.
  # Args:
  #   SCW_background: list. Output from GetSCWBackgound()
  #   resi: numeric vector. The selection of residue positions. Residue positions should
  #   be based on pdb structure not the sequence positions
  #   output_df: logi. Whether a df should be returned instead of a vector.
  # output:
  #   biased and unbiased SCW z scores
  w_df <- SCW_background[["dist_filt"]] %>%
    filter(POS_i %in% resi,
           POS_j %in% resi) %>%
    mutate(bias_term = POS_j - POS_i,
           unbias_term = 1)

  w <- c(sum(w_df$bias_term), sum(w_df$unbias_term))
  names(w) <- c("bias", "unbias")
  m <- length(resi)
  l <- SCW_background[["pdb_length"]]

  pi1 <- m * (m - 1.0) / (l * (l - 1.0))
  pi2 <- pi1 * (m - 2.0) / (l - 2.0)
  pi3 <- pi2 * (m - 3.0) / (l - 3.0)
  w_ave <- as.numeric(SCW_background[["class_count"]][4,2:3]) * pi1
  w2_ave <- (pi1 * as.numeric(SCW_background[["class_count"]][1,2:3]) +
               (pi2 * as.numeric(SCW_background[["class_count"]][2,2:3])) +
               (pi3 * as.numeric(SCW_background[["class_count"]][3,2:3])))
  sigma <- sqrt(w2_ave - w_ave * w_ave)
  z <- (w - w_ave)/sigma
  if (output_df == TRUE) {
    output <- tibble(bias = z[1], unbias = z[2])
  } else {
    output <- z
  }
  return(output)
}



#' Compute SCW z scores with multiple selection of residues
#'
#' @param SCW_background A list with parameters required for SCW z score calculation.
#' This is obtained through GetSCWBackground function.
#' @param resi_list a list of numeric vector. The selection of residue positions.
#' Residue positions should be based on pdb structure not the sequence positions.
#' @return df with biased and unbiased SCW z scores for each residue selection in
#' the resi_list. Names of the resi_list is used as the condition column. If no names
#' are provided, list numbering is used.
#' @description This function calculates the biased and unbiased SCW z scores
#' for several selections of residues in the structure.
#' @export
ComputeSCWzscore_list <- function(SCW_background, resi_list) {
  if(is.null(names(resi_list))) {
    names(resi_list) <- seq(1:length(resi_list))
  }
  output <- tibble(condition = names(resi_list), resi = resi_list) %>%
    mutate(z = map(resi, ~ComputeSCWzscore(SCW_background, resi = .x, output_df = TRUE))) %>%
    unnest(cols = c(z))
  return(output)
}


#' Compute SCW z scores for ET file
#'
#' @param pdb_file Path to pdb file, needs to be specified if SCW_bg is not provided
#' @param chain Chain in the pdb file to use, needs to be specified if SCW_bg is not provided
#' @param SCW_background A list with parameters required for SCW z score calculation.
#' This is obtained through GetSCWBackground function. Needs to be specified if pdb_file/chain
#' are not provided
#' @param ET_file Path to ranks file
#' @param coverage A numeric (0-1) value to specify the top fraction of ET residues are used
#' in SCW z score calculation
#' @param adjust_gap Logic. When False, this function behaves like PyET viewer, all residues
#' in the pdb structure are used to calculate SCW background. When TRUE, only aligned positions
#' in the pdb structure are used to calculate this background. Should set to TRUE when there is
#' a large insertion in the pdb comparing to the ET sequence. Adjust_gap only works when using
#' pdb_file/chain as input.
#' @return a df with biased and unbiased SCW z scores.
#' @description This function calculates the biased and unbiased SCW z scores
#' for top residues in the ET ranks file. Can use pdb_file/chain or pre calculated SCW_background
#' as input.
#' @export
ET_SCWzscore <- function(pdb_file = NULL, chain = NULL, SCW_background = NULL, ET_file,
                     coverage = 0.3, adjust_gap = FALSE) {
  ET <- ReadET(ET_file)
  ET_seq <- paste0(ET$AA, collapse = "")
  if (!is.null(SCW_background)) {
    scw_bg <- SCW_background
  } else {
    if (is.null(pdb_file) | is.null(chain)) {
      stop("Provide either pdb_file/chain or SCW_background")
    }
    resi_map <- CompareSeqs(pdb_file, chain, ET_seq, pos.only = TRUE)
    if (adjust_gap == TRUE) {
      scw_bg <- GetSCWBackgound(pdb_file = pdb_file, chain = chain, resi = resi_map$AA.POS.pdb)
    } else {
      scw_bg <- GetSCWBackgound(pdb_file = pdb_file, chain = chain)
    }
  }
  # When pdb file is not provided, use residue numbers in the ET output to map to strucutre,
  # if pdb file is provided, pdb seq is aligned with ET seq to map the ET top residues to
  # structure.
  if (!is.null(pdb_file)) {
    resi_select <- ET %>%
      left_join(rename(resi_map, POS = AA.POS.seq, pdb.POS = AA.POS.pdb), by = "POS")
  } else {
    resi_select <- ET %>%
      mutate(pdb.POS = POS)
  }
  resi_select <- resi_select %>%
    mutate(cov = rank(rho, ties.method = "max")/n()) %>%
    filter(cov <= coverage) %>%
    filter(!is.na(pdb.POS)) %>%
    .$pdb.POS
  output <- ComputeSCWzscore(SCW_background = scw_bg, resi = resi_select, output_df = TRUE)
  return(output)
}
