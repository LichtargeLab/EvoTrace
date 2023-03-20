## ---- echo = FALSE, message = FALSE-------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
options(tibble.print_min = 4, tibble.print_max = 4)

## ----setup, message = FALSE---------------------------------------------------
library(EvoTrace)
library(dplyr)

## -----------------------------------------------------------------------------
scales::show_col(SelectColor("ET"), labels = FALSE)
ET <- ReadET(ET_file = "b0140.cov") 
b0140.color <- ET %>%
  select(POS, coverage) %>%
  mutate(ET.color = SelectColor(color_type = "ET", prefix = "0x")[ceiling(coverage)])
b0140.color

## -----------------------------------------------------------------------------
scales::show_col(SelectColor("alphafold"), labels = FALSE)

## -----------------------------------------------------------------------------
pLDDT <- GetpLDDT("b0140.pdb", chain = "A")
b0140.color <- b0140.color %>%
  left_join(pLDDT) %>%
  mutate(pLDDT.color = SelectColor(color_type = "alphafold", prefix = "0x")[ceiling(pLDDT)])
b0140.color

## -----------------------------------------------------------------------------
custom.color <- tibble(POS = c(139, 141), color = c("red", "cyan"))

## -----------------------------------------------------------------------------
pymol.cmd <- PymolLoadFile("b0140.pdb")

## -----------------------------------------------------------------------------
pymol.cmd <- c(pymol.cmd,
               "set_name b0140, ET", 
               "copy pLDDT, ET",
               "copy custom_color, ET")

## -----------------------------------------------------------------------------
pymol.cmd <- c(pymol.cmd,
               PymolColorChainSingleColor(color = "white"))

## -----------------------------------------------------------------------------
pymol.cmd <- c(pymol.cmd,
               PymolColorChainByResidue(chain = "A", position = b0140.color$POS,
                                        color = b0140.color$ET.color,
                                        object = "ET"),
               PymolColorChainByResidue(chain = "A", position = b0140.color$POS,
                                        color = b0140.color$pLDDT.color,
                                        object = "pLDDT"),
               PymolColorChainByResidue(chain = "A", position = custom.color$POS,
                                        color = custom.color$color,
                                        object = "custom_color"))

## -----------------------------------------------------------------------------
pymol.cmd <- c(pymol.cmd,
               PymolSelectResidues("A", position = c(139, 141),
                                   selection_name = "highlight"),
               "show spheres, custom_color and highlight")

## -----------------------------------------------------------------------------
pymol.cmd <- c(pymol.cmd,
               "set seq_view, 1",
               "set grid_mode, 1",
               "cmd._cmd._draw(cmd._COb)", # This pymol cmd is necessary for png to work with grid view in pymol batch mode.
               "png b0140.png, width=5cm, height=5cm, dpi=300, ray=1")

## -----------------------------------------------------------------------------
PymolExecuteAndSave(pymol_cmd = pymol.cmd, output_file = "b0140_color.pse")

