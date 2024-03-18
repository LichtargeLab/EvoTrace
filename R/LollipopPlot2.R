#' Plot lollipop plot for case and control
#'
#' @param variants_case A dataframe with the variants in the cases. "SUB", "EA" and "AC" columns
#' are required. SUB should be 1 letter format, e.g. S10L. EA should be numerical between 0-100.
#' AC is allele count.
#' @param variants_ctrl A dataframe with the variants in the controls. Same requirements as
#' variants_case.
#' @param prot_id Protein identification. Use ENSP id for human proteins and b number for E. coli
#' proteins.
#' @param plot_domain logical. When domain annotations are available, whether plot domain
#' track.
#' @param AC_scale Whether to use linear or log scale for allele count.
#' @param show_EA_bin When TRUE, the circle sizes in the lollipop plot are adjusted based on
#' EA bins. High: (70,100], medium: (30,70] and low (0,30].
#' @param EA_color Coloring style for EA (circle in the lollipop plots). "prismatic": higher EA
#' mutations are warmer. "gray_scale": higher EA mutations are darker. "EA_bin": based on 3 EA bins, high,
#' medium and low. "black": uniform black color.
#' @param fix_scale Whether axis for AC in case and control should be the same scale or not.
#' @param pad_ratio Controls the extra space between the ET (center) track and the mutation tracks.
#' Use larger values for larger space.
#' @param domain_min_dist Controls the minimum distances between two domains in the domain track.
#' If the distance between two domains are less than domain_min_dist, they will be plotted in separate lines.
#' If domain annotations overlap, set to a larger value.
#' @param title Title for the plot.
#' @return lollipop plot
#' @description This function graphs a lollipop plot to compare mutational profile between cases and
#' controls in a given gene. The center ET track is colored as prismatic style, with the most important
#' ET positions as red. The height of the lollipops reflects allele count. The color and/or size of the
#' circles reflects EA scores for the mutations. This function fetch pre stored ET scores, which only
#' works for human and E. coli reference proteins. See LollipopPlot for ploting mutations only from the
#' case.
#' @export
#' @examples
#' # Prepare variant data
#' set.seed(940)
#' mut_case <- read_tsv(file.path(system.file("extdata", package = "EvoTrace"),
#'                               "basS_muts.tsv"),
#'                      show_col_types = FALSE) %>%
#'   mutate(AC = runif(11, 0, 1000)) %>%
#'   mutate(AC = round(AC)) %>%
#'   mutate(AF = AC/1000)
#'
#' mut_ctrl <- mut_case %>%
#'   mutate(AC = runif(11, 0, 1000)) %>%
#'   mutate(AC = round(AC)) %>%
#'   mutate(AF = AC/1000)
#'
#' # Make log scale plot using prismatic coloring
#' LollipopPlot2(variants_case = mut_case, variants_ctrl = mut_ctrl,
#'               prot_id = "b4112", AC_scale = "log",
#'               plot_domain = TRUE, show_EA_bin = TRUE,
#'               EA_color = "prismatic")
#'
#' # Make linear scale plot using EA_bin coloring
#' LollipopPlot2(variants_case = mut_case, variants_ctrl = mut_ctrl,
#'               prot_id = "b4112", AC_scale = "linear",
#'               plot_domain = FALSE, show_EA_bin = TRUE,
#'               EA_color = "EA_bin", pad_ratio = 0)


LollipopPlot2 <- function(variants_case, variants_ctrl,
                          prot_id, plot_domain = TRUE,
                          AC_scale = c("log", "linear"),
                          show_EA_bin = TRUE,
                          EA_color = c("prismatic", "gray_scale", "EA_bin", "black"),
                          fix_scale = TRUE,
                          pad_ratio = 0.05,
                          domain_min_dist = 0,
                          title = NULL) {
  AC_scale <- match.arg(AC_scale)
  EA_color <- match.arg(EA_color)
  if ((sum(c("SUB", "EA", "AC") %in% names(variants_case))) < 3) {
    stop("variants_case df should contain these cols: SUB, EA and AC")
  }

  # Check if there is ET for this protein
  if (!prot_id %in% id_map$prot_id) {
    stop("No ET for this protein")
  }

  # Fetch ET
  ET <- FetchET(prot_id) %>%
    dplyr::rename(AA_POS = POS, ET = coverage) %>%
    mutate(color = SelectColor("ET")[ceiling(ET)])

  variants_case <- variants_case %>%
    mutate(AA_REF = str_sub(SUB, 1,1),
           AA_POS = as.numeric(str_extract(SUB, "[[:digit:]]+"))) %>%
    mutate(AA_ET = ET$AA[AA_POS])

  if (sum(variants_case$AA_REF != variants_case$AA_ET) > 0) {
    unmatch_POS <- variants_case$AA_POS[variants_case$AA_REF != variants_case$AA_ET]
    unmatch_POS <- paste0(unmatch_POS, collapse = ",")
    stop(paste0("These positions have unmatched residues between the case mutations and ET: ", unmatch_POS))
  }

  if ((sum(c("SUB", "EA", "AC") %in% names(variants_ctrl))) < 3) {
    stop("variants_ctrl df should contain these cols: SUB, EA and AC")
  }
  variants_ctrl <- variants_ctrl %>%
    mutate(AA_REF = str_sub(SUB, 1,1),
           AA_POS = as.numeric(str_extract(SUB, "[[:digit:]]+"))) %>%
    mutate(AA_ET = ET$AA[AA_POS])
  if (sum(variants_ctrl$AA_REF != variants_ctrl$AA_ET) > 0) {
    unmatch_POS <- variants_ctrl$AA_POS[variants_ctrl$AA_REF != variants_ctrl$AA_ET]
    unmatch_POS <- paste0(unmatch_POS, collapse = ",")
    stop(paste0("These positions have unmatched residues between the control mutations and ET: ", unmatch_POS))
  }



  # Generate ET plot
  ET_plot <- ggplot(ET) +
    geom_col(aes(x = AA_POS, y = 1, fill = color), width = 1) +
    annotate("segment", x=0.5, xend= max(ET$AA_POS)+0.5, y=0, yend=0, linewidth = 1) +
    annotate("segment", x=0.5, xend= max(ET$AA_POS)+0.5, y=1, yend=1, linewidth = 1) +
    scale_fill_manual(values = GetManualColor(ET$color)) +
    scale_y_continuous(expand = c(0, 0)) +
    xlim(0, max(ET$AA_POS) + 1) +
    theme_nothing() +
    theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 10))


  # Generate variant plots
  mut_case <- PrepareMuts(variants_case, y_var = AC_scale, EA_color = EA_color)

  mut_ctrl <- PrepareMuts(variants_ctrl, y_var = AC_scale, EA_color = EA_color)

  if (fix_scale == TRUE) {
    max_lim_case <- max(c(mut_case$y_var, mut_ctrl$y_var))
    if (max_lim_case == 0) {
      max_lim_case <- 0.01
    }
    max_lim_ctrl <- max_lim_case
    pad_case <- pad_ratio * max_lim_case
    pad_ctrl <- pad_case
  } else {
    max_lim_case <- max(mut_case$y_var)
    if (max_lim_case == 0) {
      max_lim_case <- 0.01
    }
    max_lim_ctrl <- max(mut_ctrl$y_var)
    if (max_lim_ctrl == 0) {
      max_lim_ctrl <- 0.01
    }
    pad_case <- pad_ratio * max_lim_case
    pad_ctrl <- pad_ratio * max_lim_ctrl
  }

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
    # geom_point(aes(x=AA_POS, y=y_var, color = color), size = 2) +
    # geom_point(aes(x=AA_POS, y=y_var, alpha = EA), size = 2, color = "black", fill = "black") +
    xlim(0, max(ET$AA_POS) + 1) +
    scale_y_continuous(expand = expansion(mult = c(0, 0)), limits = c(-pad_case, max_lim_case*1.2),
                       breaks = pretty(x = c(0, max_lim_case*1.2), n = 5)) +
    y_label +
    scale_size(range = c(2, 4)) +
    scale_color_manual(values = GetManualColor(mut_case$color)) +
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
    ) +
    annotate("text", x = max(ET$AA_POS)*0.05, y = max_lim_case*1.1, label = "Case",
             fontface = "bold.italic")

  mut_ctrl_plot <- ggplot(mut_ctrl) +
    geom_segment(aes(x=AA_POS, xend=AA_POS, y=0-pad_ctrl, yend=y_var), linewidth = 0.5) +
    pop +
    # geom_point(aes(x=AA_POS, y=y_var, color = color), size = 2) +
    # geom_point(aes(x=AA_POS, y=y_var, alpha = EA), size = 2, color = "black", fill = "black") +
    scale_y_reverse(expand = expansion(mult = c(0, 0)), limits = c(max_lim_ctrl*1.2, -pad_ctrl),
                    breaks = pretty(x = c(0, max_lim_ctrl*1.2), n = 5)) +
    scale_x_continuous(position = "bottom", limits = c(0, max(ET$AA_POS) + 1),
                       breaks = scales::pretty_breaks(n = 10)) +
    y_label +
    xlab("Amino acid poistion") +
    scale_size(range = c(2, 4)) +
    scale_color_manual(values = GetManualColor(mut_ctrl$color)) +
    theme_classic(base_size = 12) +
    theme(line = element_line(linewidth = 1),
          plot.margin = margin(t = 0, r = 0, b = 10, l = 10),
          # axis.line.x = element_blank(),
          # axis.text.x = element_blank(),
          # axis.ticks.x = element_blank(),
          # axis.title.x = element_blank(),
          legend.position = "none"
    ) +
    annotate("text", x = max(ET$AA_POS)*0.05, y = max_lim_ctrl*1.1, label = "Control",
             fontface = "bold.italic")

  title_grob <- grid::textGrob(title, gp=grid::gpar(fontsize=20, fontface="bold"), x = unit(0.1, "npc"),
                               just = "left")

  if (plot_domain == TRUE) {
    domain <- id_map[id_map$prot_id == prot_id,]
    if (is.na(domain$domain) == TRUE) {
      cat("No domain information for this protein.")
      output <- egg::ggarrange(mut_case_plot, ET_plot, mut_ctrl_plot,
                               ncol = 1, heights = c(10,1,10), draw = FALSE,
                               top = title_grob)
    } else {
      domain <- id_map[id_map$prot_id == prot_id,] %>%
        mutate(domain_df = str_split(domain, "; ", simplify = FALSE)) %>%
        mutate(domain_df = map(domain_df, ParseDomain)) %>%
        select(domain_df) %>%
        unnest(cols = c(domain_df)) %>%
        mutate(group = AssignNonOverlapGroup(start, end, min_dist = domain_min_dist))

      domian_plot_height <- (max(domain$group+0.5) - 1) * 3

      domain_plot <- ggplot(domain) +
        geom_rect(aes(xmin = start, xmax = end, ymin = group, ymax = group + 0.5), fill = "gray90", color = "gray30") +
        geom_text(aes(x = (start+end)/2, y = group + 0.25, label = domain), size = 3, color = "red") +
        xlim(0, max(ET$AA_POS) + 1) +
        ylim(1, max(domain$group+0.5)) +
        theme_nothing() +
        theme(plot.margin = margin(t = 10, r = 0, b = 0, l = 10))

      output <- egg::ggarrange(domain_plot, mut_case_plot, ET_plot, mut_ctrl_plot,
                               ncol = 1, heights = c(domian_plot_height,10,1,10), draw = FALSE,
                               top = title_grob)
    }
  } else {
    output <- egg::ggarrange(mut_case_plot, ET_plot, mut_ctrl_plot,
                             ncol = 1, heights = c(10,1,10), draw = FALSE,
                             top = title_grob)
  }
  return(output)
}



# internal functions
GetManualColor <- function(color_vec) {
  output <- color_vec %>%
    unique() %>%
    sort()
  return(output)
}

PrepareMuts <- function(df, y_var = c("log", "linear"), EA_color) {
  y_var <- match.arg(y_var)
  output <- df %>%
    mutate(AA_POS = str_extract(SUB, "[[:digit:]]+"),
           AA_POS = as.numeric(AA_POS)) %>%
    mutate(logAC = log10(AC)) %>%
    mutate(EA_bin = ifelse(EA >= 70, 3,
                           ifelse(EA >= 30, 2, 1)))
  if (EA_color == "prismatic") {
    output$color <- GetColor((100.001-output$EA), lower_bound = 0, upper_bound = 100, color = "ET")
  } else if (EA_color == "gray_scale") {
    output$color <- GetColor((output$EA+0.001), lower_bound = 0, upper_bound = 100, color = "gray_scale")
  } else if (EA_color == "EA_bin") {
    output$color <- GetColor((output$EA+0.001), lower_bound = 0, upper_bound = 100, color = "EA_bin")
  } else if (EA_color == "black") {
    output$color <- "black"
  }

  if (y_var == "log") {
    output$y_var <- output$logAC
  } else {
    output$y_var <- output$AC
  }
  return(output)
}

ParseDomain <- function(x) {
  domain <- x[str_starts(x, "/note")]
  domain <- str_sub(domain, 7, -1)
  domain <- str_remove_all(domain, '"')
  coord <- x[str_starts(x, "DOMAIN")]
  coord <- str_sub(coord, 8, -1)
  output <- tibble(domain = domain, coord = coord) %>%
    separate(coord, into = c("start", "end"), sep = "\\..") %>%
    mutate(start = as.numeric(str_extract(start, "[[:digit:]]+")),
           end = as.numeric(str_extract(end, "[[:digit:]]+"))) %>%
    arrange(start, end)
  return(output)
}


# Import theme_nothing from cowplot
theme_nothing <- function (font_size = 14, font_family = "", rel_small = 12/14)
{
  theme_void(base_size = font_size, base_family = font_family) %+replace%
    theme(line = element_blank(), rect = element_blank(),
          text = element_text(family = font_family, face = "plain",
                              color = "black", size = font_size, lineheight = 0.9,
                              hjust = 0.5, vjust = 0.5, angle = 0, margin = margin(),
                              debug = FALSE), axis.line = element_blank(),
          axis.line.x = NULL, axis.line.y = NULL, axis.text = element_blank(),
          axis.text.x = NULL, axis.text.x.top = NULL, axis.text.y = NULL,
          axis.text.y.right = NULL, axis.ticks = element_blank(),
          axis.ticks.length = unit(0, "pt"), axis.title = element_blank(),
          axis.title.x = NULL, axis.title.x.top = NULL, axis.title.y = NULL,
          axis.title.y.right = NULL, legend.background = element_blank(),
          legend.spacing = unit(font_size, "pt"), legend.spacing.x = NULL,
          legend.spacing.y = NULL, legend.margin = margin(0,
                                                          0, 0, 0), legend.key = element_blank(), legend.key.size = unit(1.1 *
                                                                                                                           font_size, "pt"), legend.key.height = NULL, legend.key.width = NULL,
          legend.text = element_text(size = rel(rel_small)),
          legend.text.align = NULL, legend.title = element_text(hjust = 0),
          legend.title.align = NULL, legend.position = "none",
          legend.direction = NULL, legend.justification = "center",
          legend.box = NULL, legend.box.margin = margin(0,
                                                        0, 0, 0), legend.box.background = element_blank(),
          legend.box.spacing = unit(font_size, "pt"), panel.background = element_blank(),
          panel.border = element_blank(), panel.grid = element_blank(),
          panel.grid.major = NULL, panel.grid.minor = NULL,
          panel.spacing = unit(font_size/2, "pt"), panel.spacing.x = NULL,
          panel.spacing.y = NULL, panel.ontop = FALSE, strip.background = element_blank(),
          strip.text = element_blank(), strip.text.x = NULL,
          strip.text.y = NULL, strip.placement = "inside",
          strip.placement.x = NULL, strip.placement.y = NULL,
          strip.switch.pad.grid = unit(0, "cm"), strip.switch.pad.wrap = unit(0,
                                                                              "cm"), plot.background = element_blank(), plot.title = element_blank(),
          plot.subtitle = element_blank(), plot.caption = element_blank(),
          plot.tag = element_text(face = "bold", hjust = 0,
                                  vjust = 0.7), plot.tag.position = c(0, 1), plot.margin = margin(0,
                                                                                                  0, 0, 0), complete = TRUE)
}
