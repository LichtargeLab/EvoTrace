library(tidyverse)

legend.df <- tibble(x = 1:100, y = 1,
                    ET_color = (SelectColor("ET")),
                    EA_prismatic_color = rev(SelectColor("ET")),
                    EA_gray_color = SelectColor("gray_scale"),
                    EA_bin_color = SelectColor("EA_bin"))

# ET legend
legend.ET <- ggplot(legend.df) +
  geom_col(aes(x, y, fill = ET_color), show.legend = FALSE, width = 1) +
  scale_fill_identity() +
  scale_x_continuous(sec.axis = dup_axis(breaks = c(12, 88),
                                 labels = c("More important positions", "Less important positions")),
                     breaks = c(0:10) * 10) +
  geom_polygon(data = tibble(x = c(0.5, 100.5, 100.5, 0.5), y = c(0,0,1,1)),
               aes(x = x, y = y), color = "black", fill = NA, linewidth = 0.5) +
  cowplot::theme_cowplot() +
  theme(axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 10, hjust = 0.05, face = "plain")) +
  ggtitle("ET color:")

legend.ET


# EA prismatic
legend.EA.prismatic <- ggplot(legend.df) +
  geom_col(aes(x, y, fill = EA_prismatic_color), show.legend = FALSE, width = 1) +
  scale_fill_identity() +
  scale_x_continuous(sec.axis = dup_axis(breaks = c(10, 90),
                                         labels = c("Less impactful", "More impactful")),
                     breaks = c(0:10) * 10) +
  geom_polygon(data = tibble(x = c(0.5, 100.5, 100.5, 0.5), y = c(0,0,1,1)),
               aes(x = x, y = y), color = "black", fill = NA, linewidth = 0.5) +
  cowplot::theme_cowplot() +
  theme(axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 10, hjust = 0.05, face = "plain")) +
  ggtitle("EA color:")

legend.EA.prismatic


# EA gray
legend.EA.gray <- ggplot(legend.df) +
  geom_col(aes(x, y, fill = EA_gray_color), show.legend = FALSE, width = 1) +
  scale_fill_identity() +
  scale_x_continuous(sec.axis = dup_axis(breaks = c(10, 90),
                                         labels = c("Less impactful", "More impactful")),
                     breaks = c(0:10) * 10) +
  geom_polygon(data = tibble(x = c(0.5, 100.5, 100.5, 0.5), y = c(0,0,1,1)),
               aes(x = x, y = y), color = "black", fill = NA, linewidth = 0.5) +
  cowplot::theme_cowplot() +
  theme(axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 10, hjust = 0.05, face = "plain")) +
  ggtitle("EA color:")

legend.EA.gray



# EA bin
legend.EA.bin <- ggplot(legend.df) +
  geom_col(aes(x, y, fill = EA_bin_color), show.legend = FALSE, width = 1) +
  scale_fill_identity() +
  scale_x_continuous(sec.axis = dup_axis(breaks = c(15, 50, 85),
                                         labels = c("Low impact",
                                                    "Medium impact",
                                                    "High impact")),
                     breaks = c(0, 30, 70, 100)) +
  geom_polygon(data = tibble(x = c(0.5, 100.5, 100.5, 0.5), y = c(0,0,1,1)),
               aes(x = x, y = y), color = "black", fill = NA, linewidth = 0.5) +
  cowplot::theme_cowplot() +
  theme(axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 10, hjust = 0.05, face = "plain")) +
  ggtitle("EA color:")

legend.EA.bin

lolliplot_legend <- list(legend.ET = legend.ET,
                         legend.EA.prismatic = legend.EA.prismatic,
                         legend.EA.gray = legend.EA.gray,
                         legend.EA.bin = legend.EA.bin)

save(lolliplot_legend, file = "data/lolliplot_legend.rda")
