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
  IntervalDist <- function(interval1, interval2){
    a <- interval1[1]
    b <- interval1[2]
    c <- interval2[1]
    d <- interval2[2]
    if (b < c) {
      return(c - b)
    } else if (d < a) {
      return(a - d)
    } else {
      return(0)  # The ranges overlap
    }
  }

  span <- map2(start, end, c)
  # Get pairwise distances for all segments
  dist_df <- crossing(id1 = 1:length(span),
                      id2 = 1:length(span)) %>%
    filter(id1 < id2) %>%
    mutate(range1 = map(id1, ~span[[.]]),
           range2 = map(id2, ~span[[.]])) %>%
    mutate(dist = map2_dbl(range1, range2, IntervalDist))
  # Assign segment 1 to group 1
  group_list <- list(dist_df[dist_df$id1 == 1,])
  group <- rep(1, length(start))
  if (length(start) > 1) {
    # Starting from segment 2, check the distance
    # between each group and target segment
    for (i in 2:length(start)) {
      set_dist <- map_dbl(group_list, ~min(filter(., id2 == i)$dist))
      # Assign to new group if all the distances are smaller than min_dist
      # Otherwise, add that segment to the most distant group
      if (max(set_dist) <= min_dist) {
        group_list <- append(group_list, list(dist_df[dist_df$id1 == i,]))
        group[i] <- length(group_list)
      } else {
        group[i] <- which.max(set_dist)
        group_list[[group[i]]] <- bind_rows(group_list[[group[i]]], dist_df[dist_df$id1 == i,])
      }
    }
  }
  return(group)
}

# AssignNonOverlapGroup <- function(start, end, min_dist = 0) {
#   # Returns the single linkage (min distance) between two sets.
#   SetDist <- function(x, y) {
#     output <- sapply(x, function(x) min(abs(x - y)))
#     output <- min(output)
#     return(output)
#   }
#   span <- map2(start, end, seq)
#   group_list <- list(span[[1]])
#   group <- rep(1, length(start))
#   if (length(start) > 1) {
#     for (i in 2:length(start)) {
#       set_dist <- map_dbl(group_list, ~SetDist(., span[[i]]))
#       if (max(set_dist) <= min_dist) {
#         group_list <- append(group_list, list(span[[i]]))
#         group[i] <- length(group_list)
#       } else {
#         group[i] <- which.max(set_dist)
#         group_list[[group[i]]] <- union(group_list[[group[i]]], span[[i]])
#       }
#     }
#   }
#   return(group)
# }
