#' Pymol color pdb structure using ET
#'
#' @param pdb_file Path to pdb or pse file. Can be path to local file or url ("https://files.rcsb.org/view/1xb7.pdb").
#' @param chain The chain in the pymol session that needs to be colored. Can provide multiple chains as a vector.
#' @param ET_format Specify the ET file format. "UET" for ET server, "EA" for ET in EA output, and "ENSP"
#' for precalculated ET scores from EA VEP output.
#' @param ET_file Path to ET file. If "ENSP" is used in ET_format, provide ENSP id here, including version.
#' @param color_type The type of color range to use. Available selections: "ET",
#' "red_white", "red_white_blue", and "gray_scale".
#' @param coverage_cutoff Value between (0,1]. Only residue below this cutoff will be colored.
#' @param other_color Color for the remaining chain. Default "white". If set to NULL, color
#' are not adjusted for those chains. To color heterodimer, color the first chain with other_color as "white" and
#' save as pse. Then color the second chain using the pse file with other_color = NULL to retain the first chain color.
#' @param pml_output Path to store the pymol command file.
#' @param output_file Path to write the output pymol file.
#' @param output_format Output pymol script (pml) or pymol session (pse). Pse output requires pymol to
#' be accessible from terminal with "pymol".
#' @return A tibble with the sequence alignment between the linear sequence in ET file and the sequence in
#' the target chain in PDB.
#' @description The function first align the linear sequence from the ET file to the sequence in PDB. PDB
#' structure is then colored according the ET scores. A pymol output file is written to provided path.
#' The function returns an alignment table. The percent matching between sequence in the ET file and PDB
#' is printed in the console.
#' @export
#' @examples
#' ColorPDBByET(pdb_file = "https://files.rcsb.org/view/2pjl.pdb", chain = "A", ET_format = "ENSP",
#'              ET_file = "ENSP00000000442.6", color_type = "ET", coverage_cutoff = 0.3,
#'              output_format = "pml", output_file = "test.pml")
ColorPDBByET <- function(pdb_file, chain, ET_format = c("EA", "UET", "ENSP"),
                         ET_file, color_type = c("ET", "red_white", "red_white_blue", "gray_scale"),
                         coverage_cutoff = 1, remaining_color = "white",
                         output_file, output_format = c("pml", "pse")) {
  ET_format <- match.arg(ET_format)
  color_type <- match.arg(color_type)
  output_format <- match.arg(output_format)

  if (ET_format == "ENSP") {
    ET <- FetchET(ET_file)
  } else {
    ET <- ReadET(ET_file, ET_format = ET_format)
  }

  if (ET_format == "UET") {
    ET <- ET %>%
      mutate(coverage = coverage*100)
  }

  align_df <- tibble(chain = chain) %>%
    mutate(data = map(chain, ~CompareSeqs(pdb_file = pdb_file, chain = .,
                                          seq = paste0(ET$AA, collapse = ""), pos.only = FALSE))) %>%
    unnest(cols = data) %>%
    left_join(select(ET, AA.POS.seq = POS, coverage), by = join_by(AA.POS.seq)) %>%
    select(chain, ends_with("AA.align"), align.POS, starts_with("AA.POS"), coverage)

  align_df_filter <- align_df %>%
    filter(!is.na(AA.POS.seq)) %>%
    filter(!is.na(AA.POS.pdb))

  pymol_cmd <- c(PymolLoadFile(pdb_file))

  if (!is.null(remaining_color)) {
    pymol_cmd <- c(pymol_cmd,
                   PymolColorChainSingleColor(color = remaining_color))
  }

  for (i in chain) {
    align_df_filter_chain <- align_df_filter[align_df_filter$chain == i,]
    pymol_cmd <- c(pymol_cmd,
                   PymolColorChainSingleColor(chain = i),
                   # ET coverage is rescaled to 0-1 to be consistent with UET output
                   PymolColorChainByET(chain = i, position = align_df_filter_chain$AA.POS.pdb, coverage = align_df_filter_chain$coverage/100,
                                       coverage_cutoff = coverage_cutoff, color_type = color_type))
    aligned_count <- sum(align_df_filter_chain$seq.AA.align == align_df_filter_chain$pdb.AA.align)
    print(glue::glue("{nrow(align_df_filter_chain)}/{nrow(ET)} = {round(nrow(align_df_filter_chain)/nrow(ET),3)*100}% of the positions in ET aligned to chain {i}"))
    print(glue::glue("{aligned_count}/{nrow(align_df_filter_chain)} = {round(aligned_count/nrow(align_df_filter_chain),3)*100}% of the aligned positions in chain {i} are exact match"))
  }

  if (output_format == "pml") {
    writeLines(pymol_cmd, con = output_file)
  } else if (output_format == "pse") {
    PymolExecuteAndSave(pymol_cmd = pymol_cmd, output_file = output_file)
  }
  return(align_df)
}
