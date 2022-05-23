#' Get residue distance matrix/dataframe
#'
#' @param df atom coordinates dataframe from GetCoordinates function.
#' It should contain these columns POS, ATOM, x, y, and z
#' @param output_type output as dataframe/tibble or matrix
#' @return a tibble with POS_i, POS_j and distance or a distance matrix
#' @description Calculates distance matrix for residues from atom distances.
#' Distance between two residues is defined as the minimal distance
#' between any two atoms within those residues.
#' To calculate distance with only C-alpha, only include C-alpha atoms
#' in the input df.
#' @export
GetDistanceMatrix <- function(df, output_type = c("df", "matrix")) {
  # Calculates distance matrix for residues from atom distances
  # Distance between two residues is defined as the minimal distance
  # between any two atoms within those residues.
  # To calculate distance with only C-alpha, only include C-alpha atoms
  # in the input df
  # Args:
  #   df: tibble from GetCoordinates function. It should contain these columns
  #   POS, ATOM, x, y, and z
  #   output_type: "df" or "matrix. The df form display all pairwise distance
  #   between two residues
  # output: tibble or matrix
  output_type <- match.arg(output_type)
  df_i <- df %>%
    select(POS_i = POS, ATOM_i = ATOM, x_i = x, y_i = y, z_i = z)
  df_j <- df %>%
    select(POS_j = POS, ATOM_j = ATOM, x_j = x, y_j = y, z_j = z)
  workdf <- expand_grid(df_i, df_j, .name_repair = "minimal") %>%
    mutate(dist = sqrt((x_i - x_j)^2 + (y_i - y_j)^2 + (z_i - z_j)^2)) %>%
    # mutate(dist = round(dist, 3)) %>%
    group_by(POS_i, POS_j) %>%
    summarize(dist = min(dist)) %>%
    ungroup() %>%
    arrange(POS_i, POS_j)
  if (output_type == "matrix") {
    resn <- length(unique(df$POS))
    output <- workdf %>%
      select(POS_i, POS_j, dist) %>%
      pivot_wider(names_from = POS_j, values_from = dist) %>%
      column_to_rownames("POS_i") %>%
      as.matrix()
  } else {
    output <- workdf
  }
  return(output)
}
