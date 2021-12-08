#' @export
PymolColorChainByResidue <- function(chain, position, color, object = NULL) {
  if(is.null(object)) {
    object_cmd <- ""
  } else {
    object_cmd <- paste0(object, " and ")
  }
  workdf <- tibble(chain, position, color) %>%
    mutate(pymol_cmd = glue::glue("color {color}, {object_cmd}chain {chain} and resi {position}"))
  return(workdf$pymol_cmd)
}

#' @export
PymolSelectResidues <- function(chain, position, selection_name) {
  position_cmd <- paste0("chain ", chain, " & resi ", position, collapse = " OR ")
  output_cmd <- glue::glue("select {selection_name}, ({position_cmd})")
  return(output_cmd)
}

#' @export
PymolLoadFile <- function(input_file) {
  return(glue::glue("load {input_file}"))
}

#' @export
PymolColorChainSingleColor <- function(chain = NULL, color = "white") {
  if(is.null(chain)) {
    output <- glue::glue("color white")
  } else {
    output <- glue::glue("color {color}, chain {chain}")
  }
  return(output)
}

#' @export
PymolSaveScene <- function(scene_name, color = 1, view = 0, active = 0, rep = 0, frame = 0, animate = 0) {
  output <- glue::glue("scene {scene_name}, store, color = {color}, view = {view}, active = 0, rep = {rep}, frame = {frame}, animate = {animate}")
  return(output)
}

#' @export
PymolExecuteAndSave <- function(pymol_cmd, output_file = NULL, pymol_cmd_output = NULL) {
  if(is.null(output_file)) {
    pymol_cmd <- c(pymol_cmd)
  } else {
    pymol_cmd <- c(pymol_cmd,
                   glue::glue("save {output_file}"))
  }
  if(is.null(pymol_cmd_output)) {
    temp_pml <- tempfile("pymol", fileext = ".pml")
    write_lines(pymol_cmd, file = temp_pml)
    system(glue::glue("pymol -c {temp_pml}"))
    file.remove(temp_pml)
    rm(temp_pml)
    return("PSE file written")
  } else {
    if(str_sub(pymol_cmd_output, -4, -1) == ".pml") {
      write_lines(pymol_cmd, file = pymol_cmd_output)
      system(glue::glue("pymol -c {pymol_cmd_output}"))
      return("PSE file written")
    } else {
      return("pymol_cmd_output should end with .pml")
    }
  }
}

#' @export
PymolColorPairs <- function(chain1, position1,
                            chain2, position2,
                            pair_coverage,
                            color_type = c("ET", "red_white", "red_white_blue"),
                            group_name = "pairs",
                            pair_coverage_cutoff = NULL) {
  colorRange <- SelectColor(color_type = color_type, prefix = "0x")
  if (is.null(pair_coverage_cutoff)) {
    pair_coverage_cutoff <- max(pair_coverage)
  }
  workdf <- tibble(chain1, position1,
                   chain2, position2) %>%
    mutate(pair_label = glue::glue("{group_name}_{chain1}_{position1}_{chain2}_{position2}"),
           coverage = pair_coverage) %>%
    filter(coverage <= pair_coverage_cutoff) %>%
    mutate(coverage = 100*coverage,
           color = colorRange[ceiling(coverage)]) %>%
    mutate(pymol_cmd1 = glue::glue("distance {pair_label}, chain {chain1} and resi {position1} and name CA, chain {chain2} and resi {position2} and name CA"),
           pymol_cmd2 = glue::glue("set dash_color, {color}, {pair_label}"),
           pymol_cmd3 = glue::glue("group {group_name}, {pair_label}"))
  output_cmd <- c(
    workdf$pymol_cmd1,
    workdf$pymol_cmd2,
    workdf$pymol_cmd3,
    "hide labels",
    "set dash_gap, 0.0",
    "set dash_radius, 0.1"
  )
  return(output_cmd)
}
