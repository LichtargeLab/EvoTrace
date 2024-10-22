#' Plot lollipop plot
#'
#' @param variants A dataframe with the variants in the cases. "SUB", "EA" and "AC" columns
#' are required. SUB should be 1 letter format, e.g. S10L. EA should be numerical between 0-100.
#' AC is allele count.
#' @param prot_id Protein identification. Use ENSP id for human proteins and b number for E. coli
#' proteins.
#' @param prot_color Linear protein color. Default set as ET.
#' @param plot_domain logical. When domain annotations are available, whether plot domain
#' track.
#' @param AC_scale Whether to use linear or log scale for allele count.
#' @param show_EA_bin When TRUE, the circle sizes in the lollipop plot are adjusted based on
#' EA bins. High: (70,100], medium: (30,70] and low (0,30].
#' @param EA_color Coloring style for EA (circle in the lollipop plots). "prismatic": higher EA
#' mutations are warmer. "gray_scale": higher EA mutations are darker. "EA_bin": based on 3 EA bins, high,
#' medium and low. "black": uniform black color.
#' @param pad_ratio Controls the extra space between the ET (center) track and the mutation tracks.
#' Use larger values for larger space.
#' @param domain_min_dist Controls the minimum distances between two domains in the domain track.
#' If the distance between two domains are less than domain_min_dist, they will be plotted in separate lines.
#' If domain annotations overlap, set to a larger value.
#' @param title Title for the plot.
#' @param add_legend Whether to include ET/EA legends at the bottom of the plot.
#' @param return_individual_plots If TRUE, individual tracks from the lollipop plot is return. This is useful
#' when aligning extract track to the plot.
#' @return lollipop plot
#' @description This function graphs a lollipop plot to show mutational profile in a given gene. The center ET track is colored as prismatic style, with the most important
#' ET positions as red. The height of the lollipops reflects allele count. The color and/or size of the
#' circles reflects EA scores for the mutations. This function fetch pre stored ET scores, which only
#' works for human and E. coli reference proteins. See LollipopPlot2 for ploting mutations in both case and
#' control.
#' @export
#' @examples
#' # Prepare variant data
#' set.seed(755)
#' mut_case <- read_tsv(file.path(system.file("extdata", package = "EvoTrace"),
#'                               "basS_muts.tsv"),
#'                      show_col_types = FALSE) %>%
#'   mutate(AC = runif(11, 0, 1000)) %>%
#'   mutate(AC = round(AC)) %>%
#'   mutate(AF = AC/1000)
#'
#' # Make log scale plot using prismatic coloring
#' LollipopPlot(variants = mut_case, prot_id = "b4112",
#'              AC_scale = "log", plot_domain = TRUE,
#'              show_EA_bin = TRUE, EA_color = "prismatic")
#'
#' # Make linear scale plot using EA_bin coloring
#' LollipopPlot(variants = mut_case, prot_id = "b4112",
#'              AC_scale = "linear", plot_domain = FALSE,
#'              show_EA_bin = TRUE, EA_color = "EA_bin",
#'              pad_ratio = 0)


LollipopPlot <- function(variants,
                         prot_id, prot_color = "ET",
                         plot_domain = TRUE,
                         AC_scale = c("log", "linear"),
                         show_EA_bin = TRUE,
                         EA_color = c("prismatic", "gray_scale", "EA_bin", "black"),
                         fix_scale = TRUE,
                         domain_min_dist = 0,
                         pad_ratio = 0.1,
                         title = NULL,
                         add_legend = TRUE,
                         return_individual_plots = FALSE) {
  AC_scale <- match.arg(AC_scale)
  EA_color <- match.arg(EA_color)
  if ((sum(c("SUB", "EA", "AC") %in% names(variants))) < 3) {
    stop("variants df should contain these cols: SUB, EA and AC")
  }

  # Check if there is ET for this protein
  if (!prot_id %in% id_map$prot_id) {
    stop("No ET for this protein")
  }

  if (prot_color == "ET") {
    # Fetch ET
    ET <- FetchET(prot_id) %>%
      dplyr::rename(AA_POS = POS, ET = coverage) %>%
      mutate(color = SelectColor("ET")[ceiling(ET)])
  } else {
    ET <- FetchET(prot_id) %>%
      dplyr::rename(AA_POS = POS, ET = coverage) %>%
      mutate(color = prot_color)
  }

  variants <- variants %>%
    mutate(AA_REF = str_sub(SUB, 1,1),
           AA_POS = as.numeric(str_extract(SUB, "[[:digit:]]+"))) %>%
    mutate(AA_ET = ET$AA[AA_POS])

  if (sum(variants$AA_REF != variants$AA_ET) > 0) {
    unmatch_POS <- variants$AA_POS[variants$AA_REF != variants$AA_ET]
    unmatch_POS <- paste0(unmatch_POS, collapse = ",")
    stop(paste0("These positions have unmatched residues between the variants and ET: ", unmatch_POS))
  }




  # Generate ET plot
  ET_plot <- ggplot(ET) +
    geom_col(aes(x = AA_POS, y = 1, fill = color), width = 1) +
    annotate("segment", x=0.5, xend= max(ET$AA_POS)+0.5, y=0, yend=0, linewidth = 1) +
    annotate("segment", x=0.5, xend= max(ET$AA_POS)+0.5, y=1, yend=1, linewidth = 1) +
    scale_fill_identity() +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(position = "bottom", limits = c(0, max(ET$AA_POS) + 1),
                       breaks = scales::pretty_breaks(n = 10)) +
    xlab("Amino acid position") +
    scale_size(range = c(2, 4)) +
    theme_classic(base_size = 12) +
    theme(line = element_line(linewidth = 1),
          plot.margin = margin(t = 0, r = 0, b = 10, l = 10),
          axis.line.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none"
    )

  if (prot_color != "ET") {
    ET_plot <- ET_plot +
      annotate("segment", x=0.5, xend= 0.5, y=0, yend=1, linewidth = 0.5) +
      annotate("segment", x=max(ET$AA_POS)+0.5, xend= max(ET$AA_POS)+0.5, y=0, yend=1, linewidth = 0.5)
  }
  # Generate variant plots
  mut_case <- PrepareMuts(variants, y_var = AC_scale, EA_color = EA_color)

  max_lim_case <- max(c(mut_case$y_var))
  pad_case <- pad_ratio * max_lim_case


  y_label <- switch(AC_scale,
                    log = ylab(bquote('Log'[10]~'(Allele Count)')),
                    linear = ylab("Allele Count"))

  if (show_EA_bin == TRUE) {
    pop <- geom_point(aes(x=.data[["AA_POS"]], y=.data[["y_var"]],
                          color = .data[["color"]], size = .data[["EA_bin"]]))
  } else {
    pop <- geom_point(aes(x=.data[["AA_POS"]], y=.data[["y_var"]],
                          color = .data[["color"]]), size = 2)
  }

  mut_case_plot <- ggplot(mut_case) +
    geom_segment(aes(x=AA_POS, xend=AA_POS, y=0-pad_case, yend=y_var), linewidth = 0.5) +
    pop +
    xlim(0, max(ET$AA_POS) + 1) +
    scale_y_continuous(expand = expansion(mult = c(0, 0)), limits = c(-pad_case, max_lim_case*1.2),
                       breaks = pretty(x = c(0, max_lim_case*1.2), n = 5)) +
    y_label +
    scale_size(range = c(2, 4)) +
    scale_color_identity() +
    theme_classic(base_size = 12) +
    theme(line = element_line(linewidth = 1),
          plot.margin = margin(t = 10, r = 0, b = 0, l = 10),
          axis.line.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.length.x = unit(0, "pt"),
          panel.spacing.x = element_blank(),
          legend.position = "none"
    )

  title_grob <- grid::textGrob(title, gp=grid::gpar(fontsize=20, fontface="bold"), x = unit(0.1, "npc"),
                               just = "left")
  if (plot_domain == TRUE) {
    domain <- id_map[id_map$prot_id == prot_id,]
    if (is.na(domain$domain) == TRUE) {
      cat("No domain information for this protein.")
      domain_plot <- NA
      domain_plot_height <- NA
    } else {
      domain <- id_map[id_map$prot_id == prot_id,] %>%
        mutate(domain_df = str_split(domain, "; ", simplify = FALSE)) %>%
        mutate(domain_df = map(domain_df, ParseDomain)) %>%
        select(domain_df) %>%
        unnest(cols = c(domain_df)) %>%
        mutate(group = AssignNonOverlapGroup(start, end, min_dist = domain_min_dist))

      domain_plot_height <- (max(domain$group+0.5) - 1) * 3

      domain_plot <- ggplot(domain) +
        geom_rect(aes(xmin = start, xmax = end, ymin = group, ymax = group + 0.5), fill = "gray90", color = "gray30") +
        geom_text(aes(x = (start+end)/2, y = group + 0.25, label = domain), size = 3, color = "red") +
        xlim(0, max(ET$AA_POS) + 1) +
        ylim(1, max(domain$group+0.5)) +
        theme_nothing() +
        theme(plot.margin = margin(t = 10, r = 0, b = 0, l = 10))
    }
  } else {
    domain_plot <- NA
    domain_plot_height <- NA
  }
  c("prismatic", "gray_scale", "EA_bin", "black")
  if (add_legend == TRUE) {
    EA_legend <- switch(EA_color,
                        "prismatic" = lolliplot_legend$legend.EA.prismatic,
                        "gray_scale" = lolliplot_legend$legend.EA.gray,
                        "EA_bin" = lolliplot_legend$legend.EA.bin,
                        "black" = NA)
    EA_legend_height <- ifelse(EA_color == "black",
                               NA,
                               0.8)
    if (prot_color == "ET") {
      ET_legend <- lolliplot_legend$legend.ET
      ET_legend_height <- 0.8
    } else {
      ET_legend <- NA
      ET_legend_height <- NA
    }
  } else {
    ET_legend <- NA
    ET_legend_height <- NA
    EA_legend <- NA
    EA_legend_height <- NA
  }

  plot_list <- list(domain_plot, mut_case_plot, ET_plot, ET_legend, EA_legend)
  plot_height <- c(domain_plot_height,10,1,ET_legend_height, EA_legend_height)

  plot_list <- plot_list[!is.na(plot_list)]
  plot_height <- plot_height[!is.na(plot_height)]

  if (return_individual_plots == FALSE) {
    output <- egg::ggarrange(plots = plot_list,
                             ncol = 1, heights = plot_height, draw = FALSE,
                             top = title_grob)
  } else {
    output <- list(plots = plot_list, heights = plot_height, top = title_grob)
  }
  return(output)
}
