#' Plot lollipop plot for case and control
#'
#' @param variants_case A dataframe with the variants in the cases. "SUB", "EA" and "AC" columns
#' are required. SUB should be 1 letter format, e.g. S10L. EA should be numerical between 0-100.
#' AC is allele count.
#' @param variants_ctrl A dataframe with the variants in the controls. Same requirements as
#' variants_case.
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
#' @param fix_scale Whether axis for AC in case and control should be the same scale or not.
#' @param pad_ratio Controls the extra space between the ET (center) track and the mutation tracks.
#' Use larger values for larger space.
#' @param domain_min_dist Controls the minimum distances between two domains in the domain track.
#' If the distance between two domains are less than domain_min_dist, they will be plotted in separate lines.
#' If domain annotations overlap, set to a larger value.
#' @param title Title for the plot.
#' @param y_lab Label for y axis. Default "Allele Count".
#' @param add_legend Whether to include ET/EA legends at the bottom of the plot.
#' @param AA_start The starting residue position for the lollipop plot. Default NA, 1 is used.
#' AA_start - 1 is used as the xlim start.
#' @param AA_end The ending residue position for the lollipop plot. Default NA, max(AA_POS) is used.
#' AA_end + 1 is used as the xlim end.
#' @param return_individual_plots If TRUE, individual tracks from the lollipop plot is return. This is useful
#' when aligning extract track to the plot.
#' @param interactive If TRUE, plotly is used to render an interactive plot. Interactive plot does
#' not support ET/EA legends. AA_start and AA_end are also set to default.
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


LollipopPlot2_plotly <- function(variants_case, variants_ctrl,
                                 prot_id, prot_color = "ET",
                                 plot_domain = TRUE,
                                 AC_scale = c("log", "linear"),
                                 show_EA_bin = TRUE,
                                 EA_color = c("prismatic", "gray_scale", "EA_bin", "black"),
                                 fix_scale = TRUE,
                                 pad_ratio = 0.05,
                                 domain_min_dist = 0,
                                 title = NULL,
                                 y_lab = "Allele Count") {
  AC_scale <- match.arg(AC_scale)
  EA_color <- match.arg(EA_color)
  if ((sum(c("SUB", "EA", "AC") %in% names(variants_case))) < 3) {
    stop("variants_case df should contain these cols: SUB, EA and AC")
  }

  # Check if there is ET for this protein
  if (!prot_id %in% id_map$prot_id) {
    stop("No ET for this protein")
  }

  if (prot_color == "ET") {
    # Fetch ET
    ET <- FetchET(prot_id) %>%
      dplyr::rename(AA_POS = POS, ET = coverage) %>%
      mutate(Residue = paste0(AA, AA_POS)) %>%
      mutate(color = SelectColor("ET")[ceiling(ET)])
  } else {
    ET <- FetchET(prot_id) %>%
      dplyr::rename(AA_POS = POS, ET = coverage) %>%
      mutate(Residue = paste0(AA, AA_POS)) %>%
      mutate(color = prot_color)
  }

  AA_start <- 1
  AA_end <- max(ET$AA_POS)
  # if NA, use 1 for AA_start and max(ET$AA_POS) for AA_end
  if (AA_start < 1) {
    stop("AA_start should be larger than AA_end")
  }
  if (AA_start >= AA_end) {
    stop("AA_start should be smaller than AA_end")
  }
  if (AA_start > max(ET$AA_POS)) {
    stop("AA_start should be smaller than protein length")
  }

  if (AA_end > max(ET$AA_POS) ) {
    stop("AA_end should be smaller than protein length")
  }


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

  variants_case <- variants_case %>%
    left_join(select(variants_ctrl, SUB, AC_ctrl = AC)) %>%
    mutate(AC_ctrl = ifelse(is.na(AC_ctrl), 0, AC_ctrl)) %>%
    mutate(AC_case = AC)

  variants_ctrl <- variants_ctrl %>%
    left_join(select(variants_case, SUB, AC_case = AC)) %>%
    mutate(AC_case = ifelse(is.na(AC_case), 0, AC_case)) %>%
    mutate(AC_ctrl = AC)


  # Generate ET plot
  ET_plot <- plotly::plot_ly(data = ET) %>%

    # Add bars (equivalent to geom_col)
    plotly::add_bars(
      x = ~AA_POS,
      y = rep(1, nrow(ET)),  # Set y=1 for all bars
      marker = list(color = ~color),  # Use color column
      text = NULL,
      hovertext = ~paste("Residue:", Residue, "<br>ET:", ET),
      hoverinfo = "text",  # Display only custom text in tooltip
      width = 1,  # Set bar width
      showlegend = FALSE  # Hide legend
    ) %>%
    # Add top and bottom boundary lines (equivalent to annotate("segment"))
    plotly::add_segments(
      x = AA_start - 0.5, xend = AA_end + 0.5,
      y = 0, yend = 0,
      line = list(width = 1, color = "black"),
      hoverinfo = "none",
      showlegend = FALSE
    ) %>%
    plotly::add_segments(
      x = AA_start - 0.5, xend = AA_end + 0.5,
      y = 1, yend = 1,
      line = list(width = 1, color = "black"),
      hoverinfo = "none",
      showlegend = FALSE
    ) %>%
    # Customize axes
    plotly::layout(
      xaxis = list(
        showgrid = FALSE,
        zeroline = FALSE,
        # tickvals = NULL,  # Remove y-axis ticks
        title = ""
      ),
      yaxis = list(
        range = c(0, 1),  # Keep y-axis range tight
        showgrid = FALSE,
        zeroline = FALSE,
        title = "",
        tickvals = NULL,  # Remove y-axis ticks
        showticklabels = FALSE  # Hide y-axis labels
      ),
      margin = list(t = 0, r = 0, b = 0, l = 10),  # Set margins
      plot_bgcolor = "white",  # Match ggplot theme_nothing()
      paper_bgcolor = "white"
    )


  if (prot_color != "ET") {
    ET_plot <- ET_plot %>%
      plotly::add_segments(
        x = AA_start-0.5, xend = AA_start-0.5,
        y = 0, yend = 1,
        line = list(width = 1, color = "black"),
        hoverinfo = "none",
        showlegend = FALSE
      ) %>%
      plotly::add_segments(
        x = AA_end+0.5, xend = AA_end+0.5,
        y = 0, yend = 1,
        line = list(width = 1, color = "black"),
        hoverinfo = "none",
        showlegend = FALSE
      )
  }

  # Generate variant plots
  mut_case <- PrepareMuts(variants_case, y_var = AC_scale, EA_color = EA_color) %>%
    left_join(select(ET, AA_POS, ET), by = "AA_POS")

  mut_ctrl <- PrepareMuts(variants_ctrl, y_var = AC_scale, EA_color = EA_color) %>%
    left_join(select(ET, AA_POS, ET), by = "AA_POS")

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
                    log = paste0("Log<sub>10</sub> (", y_lab, ")"),
                    linear = y_lab)

  if (show_EA_bin == TRUE) {
    pop <- function(p) {
      output <- plotly::add_trace(p,
                          x = ~AA_POS, y = ~y_var,
                          type = 'scatter',
                          mode = 'markers',
                          text = ~paste("SUB:", SUB, "<br>EA:", EA, "<br>ET:", ET, "<br>AC_case:", AC_case,
                                        "<br>AC_ctrl:", AC_ctrl),
                          hoverinfo = "text",  # Display only custom text in tooltip
                          marker = list(size = ~EA_bin*3,
                                        color = ~color,
                                        sizemode = "diameter",
                                        opacity = 1,
                                        line = list(color = "black",
                                                    width = 1)))
      return(output)
    }
  } else {
    pop <- function(p) {
      output <- plotly::add_trace(p,
                          x = ~AA_POS, y = ~y_var,
                          type = 'scatter',
                          mode = 'markers',
                          text = ~paste("SUB:", SUB, "<br>EA:", EA, "<br>ET:", ET, "<br>AC_case:", AC_case,
                                        "<br>AC_ctrl:", AC_ctrl),
                          hoverinfo = "text",  # Display only custom text in tooltip
                          marker = list(size = 10,
                                        color = ~color,
                                        sizemode = "diameter",
                                        opacity = 1,
                                        line = list(color = "black",
                                                    width = 1)))
      return(output)
    }

  }

  mut_case_plot <- plotly::plot_ly(data = mut_case) %>%
    # Add segment lines (geom_segment)
    plotly::add_segments(x = ~AA_POS, xend = ~AA_POS,
                 y = 0 - pad_case, yend = ~y_var,
                 line = list(width = 1, color = "black"),
                 hoverinfo = "none",
                 showlegend = FALSE) %>%
    pop() %>%
    plotly::layout(
      xaxis = list(
        # showticklabels = FALSE,  # Hide x-axis text
        zeroline = FALSE,
        showline = FALSE,  # Hide x-axis line
        showgrid = FALSE,  # Remove grid lines
        ticks = "",
        title = ""
      ),

      # Adjust y-axis
      yaxis = list(
        range = c(-pad_case, max_lim_case * 1.2),
        title = y_label,  # y-axis label
        tickvals = pretty(c(0, max_lim_case * 1.2), n = 5),
        showgrid = FALSE,
        showline = TRUE,
        zeroline = FALSE,
        ticks = "outside"
      ),

      # Styling
      # margin = list(t = 10, r = 0, b = 0, l = 10),
      showlegend = FALSE,

      # Add annotation (Case label)
      annotations = list(
        list(
          x = max(mut_case$AA_POS) * 0.05,
          y = max(mut_case$y_var) * 1.1,
          text = "<b><i>Case</i></b>",
          showarrow = FALSE,
          font = list(size = 15)
        )
      )
    )


  mut_ctrl_plot <- plotly::plot_ly(data = mut_ctrl) %>%
    # Add segment lines (geom_segment)
    plotly::add_segments(x = ~AA_POS, xend = ~AA_POS,
                 y = 0 - pad_ctrl, yend = ~y_var,
                 line = list(width = 1, color = "black"),
                 hoverinfo = "none",
                 showlegend = FALSE) %>%
    pop() %>%
    plotly::layout(
      xaxis = list(
        tickmode = "array",
        ticks = "outside",
        tickvals = scales::pretty_breaks(n = 10)(c(0, max(ET$AA_POS))),
        showline = TRUE,
        showgrid = FALSE,  # Remove grid lines
        zeroline = FALSE,
        side = "bottom",
        title = "Amino acid position"
      ),

      # Adjust y-axis
      yaxis = list(
        range = c(max_lim_ctrl*1.2, -pad_ctrl),
        title = y_label,
        tickvals = pretty(c(0, max_lim_ctrl * 1.2), n = 5),
        showgrid = FALSE,
        showline = TRUE,
        zeroline = FALSE,
        autorange = FALSE,
        # autorange = "reversed",
        ticks = "outside"
      ),

      # Styling
      # margin = list(t = 0, r = 0, b = 10, l = 10),
      showlegend = FALSE,

      # Add annotation (Case label)
      annotations = list(
        list(
          x = max(mut_case$AA_POS) * 0.05,
          y = max(mut_case$y_var) * 1.1,
          text = "<b><i>Ctrl</i></b>",
          showarrow = FALSE,
          font = list(size = 15)
        )
      )
    )


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
        mutate(group = AssignNonOverlapGroup(start, end, min_dist = domain_min_dist)) %>%
        mutate(x0 = start, x1 = end, y0 = group, y1 = group + 0.5) %>%
        mutate(layout_list = pmap(list(x0, x1, y0, y1), ~list(type = "rect",
                                                              fillcolor = "#e5e5e5",
                                                              line = list(color = "black"),
                                                              opacity = 1,
                                                              x0 = ..1, x1 = ..2, xref = "x",
                                                              y0 = ..3, y1 = ..4, yref = "y",
                                                              hoverinfo = "text")))

      domain_plot_height <- (max(domain$group+0.5) - 1) * 3

      domain_plot <- plotly::plot_ly(data = domain) %>%
        plotly::add_trace(type = "scatter", mode = "none") %>%
        plotly::layout(
          shapes = domain$layout_list
        ) %>%
        plotly::add_annotations(
          x = ~(start + end) / 2,  # Midpoint of rectangle
          y = ~group + 0.25,  # Adjust y-position
          text = ~domain,  # Domain names
          showarrow = FALSE,
          hovertext = ~paste0(domain, ", ", start, "-", end),
          font = list(color = "red", size = 12)
        ) %>%
        plotly::layout(
          xaxis = list(
            showgrid = FALSE,
            zeroline = FALSE,
            title = ""
          ),
          yaxis = list(
            showgrid = FALSE,
            zeroline = FALSE,
            title = "",
            tickvals = NULL,  # Remove y-axis ticks
            showticklabels = FALSE  # Hide y-axis labels
          ),
          margin = list(t = 0, r = 0, b = 0, l = 10),  # Set margins
          plot_bgcolor = "white",  # Match ggplot theme_nothing()
          paper_bgcolor = "white"
        )

    }
  } else {
    domain_plot <- NA
    domain_plot_height <- NA
  }

  plot_list <- list(domain_plot, mut_case_plot, ET_plot, mut_ctrl_plot)
  plot_height <- c(domain_plot_height,10,1,10)

  plot_list <- plot_list[!is.na(plot_list)]
  plot_height <- plot_height[!is.na(plot_height)]


  output <- plotly::subplot(plot_list,
                            shareX = TRUE, shareY = TRUE, margin = 0,
                            nrows = length(plot_height),
                            heights = plot_height/sum(plot_height))


  return(output)
}

