#' Pymol color residues in a given chain
#'
#' @param chain The chain in the pymol session that needs to be colored
#' @param poistion A numeric vector specifying the residues that needs to be
#' colored
#' @param color A string vector with hex colors specifying the color. Should be
#' equal length of position. Color format "0xFFFFFF".
#' @param object A given object in the pymol session that needs to be colored. If is
#' NULL, then all object in the pymol session will be colored.
#' @return A string vector with pymol commands
#' @description Produce pymol commands to color residues in a given chain. This function should
#' be used in combination with PymolLoadFile and PymolExecuteAndSave. Combine PymolLoadFile and
#' all other pymol cmds with c(), then pipe into PymolExecuteAndSave.
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

#' Pymol select residues
#'
#' @param chain The chain in the pymol session that needs to be selected
#' @param poistion A numeric vector specifying the residues that needs to be
#' selected
#' @return A string vector with pymol commands
#' @description Produce pymol commands to select residues in a given chain as a selection.
#' This function should be used in combination with PymolLoadFile and PymolExecuteAndSave.
#' Combine PymolLoadFile and all other pymol cmds with c(), then pipe into PymolExecuteAndSave.
#' @export
PymolSelectResidues <- function(chain, position, selection_name) {
  position_cmd <- paste0("chain ", chain, " & resi ", position, collapse = " OR ")
  output_cmd <- glue::glue("select {selection_name}, ({position_cmd})")
  return(output_cmd)
}

#' Load pdb/pse file
#'
#' @param input_file Path to the pdb/pse file of interest
#' @return A string vector with pymol commands
#' @description Produce pymol commands open a pdb/pse file. This function should be used in
#' combination with PymolExecuteAndSave and other pymol functions. Combine PymolLoadFile and
#' all other pymol cmds with c(), then pipe into PymolExecuteAndSave.
#' @export
PymolLoadFile <- function(input_file) {
  return(glue::glue("load {input_file}"))
}

#' Pymol color residues in a given chain with one color
#'
#' @param chain The chain in the pymol session that needs to be colored
#' @param color A string vector with hex colors specifying the color. Should be
#' equal length of position. Color format "0xFFFFFF" or color name. Default is white.
#' @return A string vector with pymol commands
#' @description Produce pymol commands to color a given chain with on color. This function should
#' be used in combination with PymolLoadFile and PymolExecuteAndSave. Combine PymolLoadFile and
#' all other pymol cmds with c(), then pipe into PymolExecuteAndSave.
#' @export
PymolColorChainSingleColor <- function(chain = NULL, color = "white") {
  if(is.null(chain)) {
    output <- glue::glue("color white")
  } else {
    output <- glue::glue("color {color}, chain {chain}")
  }
  return(output)
}

#' Pymol save scene
#'
#' @param scene_name Name of the scene in pymol
#' @param color 0 or 1. Whether color in the current pymol session is saved in the scene.
#' @param view 0 or 1. Whether view in the current pymol session is saved in the scene.
#' @param active 0 or 1. Whether active objects in the current pymol session are saved in
#' the scene.
#' @param rep 0 or 1. Whether representations in the current pymol session are saved in
#' the scene.
#' @param frame 0 or 1. Whether frame in the current pymol session is saved in the scene.
#' Not sure what this does in pymol, probably relate to movie.
#' @param animate 0 or 1. Whether animation in the current pymol session is saved in the scene.
#' Not sure what this does in pymol, probably relate to movie.
#' @return A string vector with pymol commands
#' @description Produce pymol commands to save current state of pymol objects as a scene. This
#' is very useful to store different coloring of an object in pymol.This function should be
#' used in combination with PymolLoadFile and PymolExecuteAndSave. Combine PymolLoadFile and
#' all other pymol cmds with c(), then pipe into PymolExecuteAndSave.
#' @export

PymolSaveScene <- function(scene_name, color = 1, view = 0, active = 0, rep = 0, frame = 0, animate = 0) {
  output <- glue::glue("scene {scene_name}, store, color = {color}, view = {view}, active = 0, rep = {rep}, frame = {frame}, animate = {animate}")
  return(output)
}

#' Execute pymol commands
#'
#' @param pymol_cmd a string vector of all pymol command that need to be run
#' @param output_file file name of the pymol session to be saved. If is null, pymol
#' commands will be run, but there is no additional save step.
#' @param pymol_cmd_output file name of the pymol script, should end with .pml. Default is
#' null, the pymol script will be saved as a temporary file, and will be removed after the
#' run.
#' @param pymol_path path to pymol. Default is null, default pymol is called.
#' @return Nother is returned.
#' @description pymol commands in pymol_cmd are written as a pymol script, then pymol is evoked
#' to run that script. Pymol session and pymol script can be saved.
#' @export
PymolExecuteAndSave <- function(pymol_cmd, output_file = NULL, pymol_cmd_output = NULL, pymol_path = NULL) {
  if(is.null(pymol_path)) {
    pymol_path <- "pymol"
  }
  if(is.null(output_file)) {
    pymol_cmd <- c(pymol_cmd)
  } else {
    pymol_cmd <- c(pymol_cmd,
                   glue::glue("save {output_file}"))
  }
  if(is.null(pymol_cmd_output)) {
    temp_pml <- tempfile("pymol", fileext = ".pml")
    write_lines(pymol_cmd, file = temp_pml)
    system(glue::glue("{pymol_path} -c {temp_pml}"))
    file.remove(temp_pml)
    rm(temp_pml)
    return("PSE file written")
  } else {
    if(str_sub(pymol_cmd_output, -4, -1) == ".pml") {
      write_lines(pymol_cmd, file = pymol_cmd_output)
      system(glue::glue("{pymol_path} -c {pymol_cmd_output}"))
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
