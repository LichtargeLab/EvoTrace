#' Assign groups with minimal overlaps
#'
#' @param start start positions
#' @param end end positions
#' @param min_dist minimal distance to merge two groups
#' @return a vector of the group ids
#' @description Group genomic features into non overlapping groups. Each entry/row is one feature.
#' Start and end describe the start end positions of the features. Features with in each group have to
#' be at least min_dist away from each other.
#' @export

AssignNonOverlapGroup <- function(start, end, min_dist = 0) {
  # Returns the single linkage (min distance) between two sets.
  SetDist <- function(x, y) {
    output <- sapply(x, function(x) min(abs(x - y)))
    output <- min(output)
    return(output)
  }
  span <- map2(start, end, seq)
  group_list <- list(span[[1]])
  group <- rep(1, length(start))
  if (length(start) > 1) {
    for (i in 2:length(start)) {
      set_dist <- map_dbl(group_list, ~SetDist(., span[[i]]))
      if (max(set_dist) <= min_dist) {
        group_list <- append(group_list, list(span[[i]]))
        group[i] <- length(group_list)
      } else {
        group[i] <- which.max(set_dist)
        group_list[[group[i]]] <- union(group_list[[group[i]]], span[[i]])
      }
    }
  }
  return(group)
}
