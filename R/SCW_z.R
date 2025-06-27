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
  # It will produce values for both unbiased (1), biased (j-i) and
  # adj_dist (cutoff - dist) approach
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
    filter(A == 1) %>%
    mutate(bias_term = POS_j - POS_i,
           unbias_term = 1,
           adj_dist_term = dist_cutoff - dist)

  AssignType <- function(i,j,k,l) {
    output <- rep("C", length(i))
    check <- (i == k) + (i == l) + (j == k) + (j == l)
    output[check == 2] <- "A"
    output[check == 1] <- "B"
    return(output)
  }

  w2_class <- select(dist_filt, i = POS_i, j = POS_j, dist1 = dist) %>%
    crossing(., select(dist_filt, k = POS_i, l = POS_j, dist2 = dist)) %>%
    mutate(bias_term = (j-i)*(l-k),
           unbias_term = 1,
           adj_dist_term = (dist_cutoff - dist1)*(dist_cutoff - dist2)) %>%
    mutate(type = AssignType(i,j,k,l))

  w_ave_df <- dist_filt %>%
    mutate(bias_term = POS_j - POS_i,
           unbias_term = 1,
           adj_dist_term = dist_cutoff - dist)

  class_count <- w2_class %>%
    group_by(type) %>%
    summarize(bias = sum(bias_term),
              unbias = sum(unbias_term),
              adj_dist = sum(adj_dist_term)) %>%
    ungroup() %>%
    add_row(type = "SS_hat",
            bias = sum(w_ave_df$bias_term),
            unbias = sum(w_ave_df$unbias_term),
            adj_dist = sum(w_ave_df$adj_dist_term))
  pdb_length = length(unique(dist$POS_i)) + 1

  return(list("dist_filt" = select(dist_filt, POS_i, POS_j, dist),
              "class_count" = class_count,
              "pdb_length" = pdb_length,
              "dist_cutoff" = dist_cutoff))
}




#' Compute SCW z scores
#'
#' @param SCW_background A list with parameters required for SCW z score calculation.
#' This is obtained through GetSCWBackground function.
#' @param resi numeric vector. The selection of residue positions. Residue positions
#' should be based on pdb structure not the sequence positions
#' @param output_df Whether a df should be returned instead of a vector.
#' @return a vector or df with biased, unbiased, and adj_dist SCW z scores.
#' @description This function calculates the biased, unbiased, and adj_dist SCW z scores
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
  #   biased, unbiased, and adj_dist SCW z scores
  w_df <- SCW_background[["dist_filt"]] %>%
    mutate(bias_term = POS_j - POS_i,
           unbias_term = 1,
           adj_dist_term = SCW_background[["dist_cutoff"]] - dist) %>%
    filter(POS_i %in% resi,
           POS_j %in% resi)

  w <- c(sum(w_df$bias_term), sum(w_df$unbias_term), sum(w_df$adj_dist_term))
  names(w) <- c("bias", "unbias", "adj_dist")
  m <- length(resi)
  l <- SCW_background[["pdb_length"]]

  pi1 <- m * (m - 1.0) / (l * (l - 1.0))
  pi2 <- pi1 * (m - 2.0) / (l - 2.0)
  pi3 <- pi2 * (m - 3.0) / (l - 3.0)
  w_ave <- as.numeric(SCW_background[["class_count"]][4,2:4]) * pi1
  w2_ave <- (pi1 * as.numeric(SCW_background[["class_count"]][1,2:4]) +
               (pi2 * as.numeric(SCW_background[["class_count"]][2,2:4])) +
               (pi3 * as.numeric(SCW_background[["class_count"]][3,2:4])))
  sigma <- sqrt(w2_ave - w_ave * w_ave)
  z <- (w - w_ave)/sigma
  if (output_df == TRUE) {
    output <- tibble(bias = z[1], unbias = z[2], adj_dist = z[3])
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
  # When pdb file is not provided, use residue numbers in the ET output to map to strucutre,
  # if pdb file is provided, pdb seq is aligned with ET seq to map the ET top residues to
  # structure.
  if (!is.null(pdb_file)) {
    resi_map <- CompareSeqs(pdb_file, chain, ET_seq, pos.only = TRUE)
    resi_select <- ET %>%
      left_join(dplyr::rename(resi_map, POS = AA.POS.seq, pdb.POS = AA.POS.pdb), by = "POS")
  } else {
    resi_select <- ET %>%
      mutate(pdb.POS = POS)
  }
  if (!is.null(SCW_background)) {
    scw_bg <- SCW_background
  } else {
    if (is.null(pdb_file) | is.null(chain)) {
      stop("Provide either pdb_file/chain or SCW_background")
    }
    if (adjust_gap == TRUE) {
      scw_bg <- GetSCWBackgound(pdb_file = pdb_file, chain = chain, resi = resi_map$AA.POS.pdb)
    } else {
      scw_bg <- GetSCWBackgound(pdb_file = pdb_file, chain = chain)
    }
  }

  resi_select <- resi_select %>%
    mutate(cov = rank(rho, ties.method = "max")/n()) %>%
    select(-coverage) %>%
    filter(cov <= coverage) %>%
    filter(!is.na(pdb.POS)) %>%
    .$pdb.POS
  output <- ComputeSCWzscore(SCW_background = scw_bg, resi = resi_select, output_df = TRUE)
  return(output)
}


#' Compute the background stats for SCW z score calculation
#'
#' @param pdb_file path to pdb file
#' @param chain the protein chain SCW z score is calculated on
#' @param dist_cutoff the distance cutoff that consider two atoms are adjacent
#' @param resi a vector of residues numbers that are included in SCW bg calculation.
#' If NULL, then all residues in that chain is used.
#' @param temp_dir temp dir to store intermediate files. Default is to create a temp dir
#' in the current working directory.
#' @param cores number of cores to parallel computing.
#' @return a list with parameters that are needed for SCW z scores calculation
#' @description Calculates the required stats for w_hat, w^2_hat calculation in
#' SCW z scores. Residues are consider adjacent when any of their heavy atoms are
#' within the cutoff range. Intermediate files are produced to resolve memory issue.
#' The output is feed into ComputeSCWzscore function for
#' z score calculation.
#' @export
GetSCWBackgound_reduce_memory <- function(pdb_file, chain, dist_cutoff = 4, resi = NULL,
                                          temp_dir = "temp",
                                          cores = 1) {
  # This computes class counts used in calculation of
  # w_hat, w^2_hat in SCW z score.
  # It will produce values for both unbiased (1), biased (j-i) and
  # adj_dist (cutoff - dist) approach
  # Args:
  #   pdb.file: path to pdb file
  #   chain: chain of interest
  #   dist_cutoff: the cutoff to determine adjacency of two residues in the structure
  #   resi: a vector of residue numbers to compute SCW background with. If NULL, all residues in the
  #     chain are used.
  #   temp_dir: temp dir to store intermediate files
  #   cores: number of cores for parallel computing
  # output:
  #   list contains the distance df and the class count df
  if (!dir.exists(temp_dir)) {
    dir.create(temp_dir)
  }

  if (!dir.exists(paste0(temp_dir, "/class1/"))) {
    dir.create(paste0(temp_dir, "/class1/"))
  }
  if (!dir.exists(paste0(temp_dir, "/class2/"))) {
    dir.create(paste0(temp_dir, "/class2/"))
  }


  cord <- GetCoordinates(pdb_file, chain, CA_only = FALSE)
  if(!is.null(resi)) {
    cord <- cord %>%
      filter(POS %in% resi)
  }
  dist <- GetDistanceMatrix(cord, output_type = "df") %>%
    mutate(A = ifelse(dist < dist_cutoff, 1, 0))

  dist_filt <- dist %>%
    filter(A == 1) %>%
    mutate(bias_term = POS_j - POS_i,
           unbias_term = 1,
           adj_dist_term = dist_cutoff - dist)

  AssignType <- function(i,j,k,l) {
    output <- rep("C", length(i))
    check <- (i == k) + (i == l) + (j == k) + (j == l)
    output[check == 2] <- "A"
    output[check == 1] <- "B"
    return(output)
  }

  Compute_w2_class_1 <- function(dist_filt, POS, temp_file = NULL) {
    output <- select(dist_filt, i = POS_i, j = POS_j,
                     dist1 = dist,
                     bias_term1 = bias_term,
                     unbias_term1 = unbias_term,
                     adj_dist_term1 = adj_dist_term) %>%
      filter(i == POS) %>%
      crossing(., select(dist_filt, k = POS_i, l = POS_j,
                         dist2 = dist,
                         bias_term2 = bias_term,
                         unbias_term2 = unbias_term,
                         adj_dist_term2 = adj_dist_term) %>%
                 filter(k > POS)) %>%
      mutate(bias_term = (j-i)*(l-k),
             unbias_term = 1,
             adj_dist_term = (dist_cutoff - dist1)*(dist_cutoff - dist2)) %>%
      mutate(type = AssignType(i,j,k,l)) %>%
      group_by(type) %>%
      summarize(bias = sum(bias_term),
                unbias = sum(unbias_term),
                adj_dist = sum(adj_dist_term))
    if (!is.null(temp_file)) {
      saveRDS(output, file = temp_file, compress = FALSE)
    }
    return(output)
  }


  Compute_w2_class_2 <- function(dist_filt, POS, temp_file = NULL) {
    filt_df <- filter(dist_filt, POS_i == POS)
    output <- select(filt_df, i = POS_i, j = POS_j,
                     dist1 = dist,
                     bias_term1 = bias_term,
                     unbias_term1 = unbias_term,
                     adj_dist_term1 = adj_dist_term) %>%
      crossing(., select(filt_df, k = POS_i, l = POS_j,
                         dist2 = dist,
                         bias_term2 = bias_term,
                         unbias_term2 = unbias_term,
                         adj_dist_term2 = adj_dist_term)) %>%
      mutate(bias_term = (j-i)*(l-k),
             unbias_term = 1,
             adj_dist_term = (dist_cutoff - dist1)*(dist_cutoff - dist2)) %>%
      mutate(type = AssignType(i,j,k,l)) %>%
      group_by(type) %>%
      summarize(bias = sum(bias_term),
                unbias = sum(unbias_term),
                adj_dist = sum(adj_dist_term))
    if (!is.null(temp_file)) {
      saveRDS(output, file = temp_file, compress = FALSE)
    }
    return(output)
  }

  POS_filt <- sort(unique(c(dist_filt$POS_i, dist_filt$POS_j)))


  class1_temp_df <- tibble(POS = POS_filt,
                           class1_temp = paste0(temp_dir, "/class1/", POS, ".rds")) %>%
    filter(!file.exists(class1_temp))

  class2_temp_df <- tibble(POS = POS_filt,
                           class2_temp = paste0(temp_dir, "/class2/", POS, ".rds")) %>%
    filter(!file.exists(class2_temp))



  if (nrow(class2_temp_df) > 0) {
    parallel::mcmapply(Compute_w2_class_2, POS = class2_temp_df$POS,
             temp_file = class2_temp_df$class2_temp,
             MoreArgs = list(dist_filt = dist_filt), mc.cores = cores,
             SIMPLIFY = FALSE, USE.NAMES = FALSE)
  }

  if (nrow(class1_temp_df) > 0) {
    parallel::mcmapply(Compute_w2_class_1, POS = class1_temp_df$POS,
             temp_file = class1_temp_df$class1_temp,
             MoreArgs = list(dist_filt = dist_filt), mc.cores = cores,
             SIMPLIFY = FALSE, USE.NAMES = FALSE)
  }


  w2_class_1 <- tibble(file = list.files(paste0(temp_dir, "/class1/"), full.names = TRUE)) %>%
    mutate(data = map(file, readRDS, .progress = TRUE)) %>%
    select(-file) %>%
    unnest(cols = c(data)) %>%
    group_by(type) %>%
    summarize(bias1 = sum(bias) * 2,,
              unbias1 = sum(unbias) * 2,
              adj_dist1 = sum(adj_dist) * 2)

  w2_class_2 <- tibble(file = list.files(paste0(temp_dir, "/class2/"), full.names = TRUE)) %>%
    mutate(data = map(file, readRDS, .progress = TRUE)) %>%
    select(-file) %>%
    unnest(cols = c(data)) %>%
    group_by(type) %>%
    summarize(bias2 = sum(bias),
              unbias2 = sum(unbias),
              adj_dist2 = sum(adj_dist))

  w_ave_df <- dist_filt %>%
    mutate(bias_term = POS_j - POS_i,
           unbias_term = 1,
           adj_dist_term = dist_cutoff - dist)

  class_count <- full_join(w2_class_2, w2_class_1, by = "type") %>%
    replace_na(list(unbias1 = 0, unbias2 = 0,
                    bias1 = 0, bias2 = 0,
                    adj_dist1 = 0, adj_dist2 = 0)) %>%
    mutate(bias = bias1 + bias2,
           unbias = unbias1 + unbias2,
           adj_dist = adj_dist1 + adj_dist2) %>%
    select(type, bias, unbias, adj_dist) %>%
    add_row(type = "SS_hat",
            bias = sum(w_ave_df$bias_term),
            unbias = sum(w_ave_df$unbias_term),
            adj_dist = sum(w_ave_df$adj_dist_term))

  pdb_length = length(unique(dist$POS_i)) + 1

  return(list("dist_filt" = select(dist_filt, POS_i, POS_j, dist),
              "class_count" = class_count,
              "pdb_length" = pdb_length,
              "dist_cutoff" = dist_cutoff))
}


