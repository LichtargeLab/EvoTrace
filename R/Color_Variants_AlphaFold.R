#' Map variants and ET scores to AlphaFold structures
#'
#' @param variants_case A dataframe that stores the mutations in the cases in target gene.
#' The dataframe should contain at least 2 columns: SUB and EA. Single letter animo acid code should
#' be used in SUB, eg. L10P. EA can be replaced with other values, but have to scale to the
#' range of (0, 100].
#' @param variants_ctrl A dataframe that stores the mutations in the controls in target gene.
#' @param pml_output Path to store the pymol command file.
#' @param prot_id Protein identification. Use full ENSP for human proteins, eg. ENSP00000388638.1.
#' Use locus tag for E. coli MG1655 proteins, eg. b4112.
#' @param output_format Output pymol script (pml) or pymol session (pse). Pse output requires pymol to
#' be accessble from terminal with "pymol".
#' @return A tibble that contains the coloring information. The pml file will
#' be written to the location indicated by pml_output.
#' @description Map variants and ET scores to AlphaFold structures. Pymol is required to view the pml
#' file. Network connection is also required to load the pdb file from AlphaFold.
#' @export

Color_Variants_AlphaFold <- function(variants_case, variants_ctrl = NULL,
                                     prot_id, pml_output, output_format = c("pml", "pse")) {
  output_format <- match.arg(output_format)
  if ((sum(c("SUB", "EA") %in% names(variants_case))) < 2) {
    stop("variants_case df should contain these cols: SUB, EA")
  }

  # Check if there is ET for this protein
  if (!prot_id %in% id_map$prot_id) {
    return("No ET for this protein")
  }
  # When ctrl variants are not provided, ctrl panel is first generated
  # using case variants, then gets removed.
  rm_ctrl <- is.null(variants_ctrl)

  AF_url <- id_map$AF_url[which(id_map$prot_id == prot_id)]
  if (is.na(AF_url)) {
    return("No AlphaFold structure for this protein")
  }
  AF_pdb <- read_lines(AF_url)

  ET_df <- FetchET(prot_id) %>%
    dplyr::rename(AA_POS = POS)

  variants_case <- variants_case %>%
    mutate(AA_REF = str_sub(SUB, 1,1),
           AA_POS = as.numeric(str_extract(SUB, "[[:digit:]]+"))) %>%
    group_by(AA_POS, AA_REF) %>% # If multiple mutations are reported in the same position, the max EA is used.
    summarize(EA = max(EA), .groups = "drop") %>%
    mutate(AA_ET = ET_df$AA[AA_POS])

  if (sum(variants_case$AA_REF != variants_case$AA_ET) > 0) {
    unmatch_POS <- variants_case$AA_POS[variants_case$AA_REF != variants_case$AA_ET]
    unmatch_POS <- paste0(unmatch_POS, collapse = ",")
    stop(paste0("These positions have unmatched residues between the case mutations and ET: ", unmatch_POS))
  }

  if (!is.null(variants_ctrl)) {
    if ((sum(c("SUB", "EA") %in% names(variants_ctrl))) < 2) {
      stop("variants_ctrl df should contain these cols: SUB, EA")
    }
    variants_ctrl <- variants_ctrl %>%
      mutate(AA_REF = str_sub(SUB, 1,1),
             AA_POS = as.numeric(str_extract(SUB, "[[:digit:]]+"))) %>%
      group_by(AA_POS, AA_REF) %>% # If multiple mutations are reported in the same position, the max EA is used.
      summarize(EA = max(EA), .groups = "drop") %>%
      mutate(AA_ET = ET_df$AA[AA_POS])
    if (sum(variants_ctrl$AA_REF != variants_ctrl$AA_ET) > 0) {
      unmatch_POS <- variants_ctrl$AA_POS[variants_ctrl$AA_REF != variants_ctrl$AA_ET]
      unmatch_POS <- paste0(unmatch_POS, collapse = ",")
      stop(paste0("These positions have unmatched residues between the control mutations and ET: ", unmatch_POS))
    }
  } else {
    # copy case to control if no control is provided
    variants_ctrl <- variants_case
  }

  pymol_color_df <- ET_df %>%
    left_join(GetpLDDT(I(AF_pdb), chain = "A"), by = c("AA_POS" = "POS")) %>%
    left_join(select(variants_case, AA_POS, EA_case = EA), by = "AA_POS") %>%
    left_join(select(variants_ctrl, AA_POS, EA_ctrl = EA), by = "AA_POS") %>%
    # due to rounding issues, some missense variants have EA = 0
    # adjust them to EA = 0.01 for easier color assignment.
    mutate(EA_case = ifelse(EA_case == 0, 0.01, EA_case)) %>%
    mutate(EA_ctrl = ifelse(EA_ctrl == 0, 0.01, EA_ctrl)) %>%
    dplyr::rename(ET = coverage) %>%
    dplyr::mutate(ET_color = GetColor(ET, lower_bound = 0, upper_bound = 100, color = "ET"),
                  pLDDT_color = GetColor(pLDDT, lower_bound = 0, upper_bound = 100, color = "alphafold"),
                  EA_case_color = GetColor((100.0001-EA_case), lower_bound = 0, upper_bound = 100, color = "ET"),
                  EA_ctrl_color = GetColor((100.0001-EA_ctrl), lower_bound = 0, upper_bound = 100, color = "ET")
    ) %>%
    dplyr::mutate(dplyr::across(dplyr::ends_with("color"),
                                ~stringr::str_replace(., "#", "0x"))) %>%
    replace_na(list(EA_case_color = "white", EA_ctrl_color = "white"))

  pdb_name <- str_split(AF_url, "/|.pdb", simplify = TRUE)
  pdb_name <- pdb_name[length(pdb_name)-1]
  pymol_cmd <- c(PymolLoadFile(AF_url),
                 paste0("set_name ", pdb_name,", ET"),
                 "copy pLDDT, ET",
                 "copy EA_case, ET",
                 "copy EA_ctrl, ET",
                 "color white",
                 PymolColorChainByResidue(chain = "A", position = pymol_color_df$AA_POS,
                                          color = pymol_color_df$ET_color, object = "ET"),
                 PymolColorChainByResidue(chain = "A", position = pymol_color_df$AA_POS,
                                          color = pymol_color_df$pLDDT_color, object = "pLDDT"),
                 PymolColorChainByResidue(chain = "A", position = pymol_color_df$AA_POS,
                                          color = pymol_color_df$EA_case_color, object = "EA_case"),
                 PymolColorChainByResidue(chain = "A", position = pymol_color_df$AA_POS,
                                          color = pymol_color_df$EA_ctrl_color, object = "EA_ctrl"),
                 "set grid_mode, 1",
                 "set seq_view, 1",
                 paste0("select case_mut, chain A and (ET or pLDDT or EA_case) and resi ", paste0(variants_case$AA_POS, collapse = "+")),
                 paste0("select ctrl_mut, chain A and (ET or pLDDT or EA_ctrl) and resi ", paste0(variants_ctrl$AA_POS, collapse = "+")),
                 "show sphere, case_mut",
                 "show sphere, ctrl_mut",
                 "deselect")
  if (rm_ctrl == TRUE) {
    pymol_cmd <- c(pymol_cmd,
                   "delete ctrl_mut",
                   "delete EA_ctrl")
    pymol_color_df <- pymol_color_df %>%
      mutate(EA_ctrl = NA, EA_ctrl_color = NA)
  }
  if (output_format == "pml") {
    writeLines(pymol_cmd, con = pml_output)
  } else if (output_format == "pse") {
    PymolExecuteAndSave(pymol_cmd = pymol_cmd, output_file = pml_output)
  }
  return(pymol_color_df)
}

