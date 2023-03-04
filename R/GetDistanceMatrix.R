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
  df <- mutate(df, id = 1:n())
  workdf <- tibble(id_i = head(df$id, -1)) %>%
    mutate(id_j = map(id_i, ~seq(.+1, nrow(df)))) %>%
    unnest(cols = id_j)
  workdf$POS_i <- df$POS[workdf$id_i]
  workdf$POS_j <- df$POS[workdf$id_j]
  workdf <- filter(workdf, POS_i != POS_j)
  workdf$x_i <- df$x[workdf$id_i]
  workdf$y_i <- df$y[workdf$id_i]
  workdf$z_i <- df$z[workdf$id_i]
  workdf$x_j <- df$x[workdf$id_j]
  workdf$y_j <- df$y[workdf$id_j]
  workdf$z_j <- df$z[workdf$id_j]
  workdf <- workdf %>%
    mutate(dist = sqrt((x_i - x_j)^2 + (y_i - y_j)^2 + (z_i - z_j)^2)) %>%
    # mutate(dist = round(dist, 3)) %>%
    group_by(POS_i, POS_j) %>%
    summarize(dist = min(dist), .groups = "drop") %>%
    ungroup() %>%
    arrange(POS_i, POS_j)
  if (output_type == "matrix") {
    resi <- unique(c(workdf$POS_i, workdf$POS_j))
    output <- bind_rows(select(workdf, POS_i, POS_j, dist),
                        select(workdf, POS_i = POS_j, POS_j = POS_i, dist),
                        tibble(POS_i = resi, POS_j = resi, dist = 0)) %>%
      arrange(POS_i, POS_j) %>%
      pivot_wider(names_from = POS_j, values_from = dist) %>%
      column_to_rownames("POS_i") %>%
      as.matrix()
  } else {
    output <- workdf
  }
  return(output)
}
