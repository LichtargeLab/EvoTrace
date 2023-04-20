#' Return a preset color range (hex code) with 100 colors
#'
#' @param color_type The type of color range to use. Available selections: "ET",
#' "red_white", "red_white_blue", "white_red", "white_blue, "alphafold", "gray_sacle", "EA_bin".
#' @param prefix "#" or "0x". Use "0x" when passing color values to pymol.
#' @return a character string with length of 100. Contains hex code of 100 colors.
#' @description Return a color scheme of 100 colors. For ET color, the smaller the
#' index, the warmer the color. For alphafold color, the larger the index, the more
#' confident of the structure prediction.
#' @export
SelectColor <- function(color_type = c("ET", "red_white", "red_white_blue", "white_red",
                                       "white_blue", "alphafold", "gray_scale", "EA_bin"),
                        prefix = "#") {
  if (color_type == "ET") {
    colorRange <- c("ff0000", "ff0c00", "ff1800", "ff2400", "ff3000", "ff3d00", "ff4900", "ff5500", "ff6100", "ff6e00",
                    "ff7a00", "ff8600", "ff9200", "ff9f00", "ffab00", "ffb700", "ffc300", "ffd000", "ffdc00", "ffe800",
                    "fff400", "fcff00", "f0ff00", "e4ff00", "d8ff00", "cbff00", "bfff00", "b3ff00", "a7ff00", "9bff00",
                    "8eff00", "82ff00", "76ff00", "6aff00", "5dff00", "51ff00", "45ff00", "39ff00", "2cff00", "20ff00",
                    "14ff00", "08ff00", "00ff04", "00ff10", "00ff1c", "00ff28", "00ff35", "00ff41", "00ff4d", "00ff59",
                    "00ff66", "00ff72", "00ff7e", "00ff8a", "00ff96", "00ffa3", "00ffaf", "00ffbb", "00ffc7", "00ffd4",
                    "00ffe0", "00ffec", "00fff8", "00f8ff", "00ecff", "00e0ff", "00d4ff", "00c7ff", "00bbff", "00afff",
                    "00a3ff", "0096ff", "008aff", "007eff", "0072ff", "0066ff", "0059ff", "004dff", "0041ff", "0035ff",
                    "0028ff", "001cff", "0010ff", "0004ff", "0800ff", "1400ff", "2000ff", "2c00ff", "3900ff", "4500ff",
                    "5100ff", "5d00ff", "6a00ff", "7600ff", "8200ff", "8e00ff", "9b00ff", "a700ff", "b300ff", "bf00ff") %>%
      paste0("#", .)
  } else if (color_type == "alphafold") {
    colorRange <- c(rep("#f17c42", 50),
                    rep("#fcdb4b", 20),
                    rep("#65cbf3", 20),
                    rep("#0153d6", 10))
  } else if (color_type == "red_white") {
    colorRange <- colorRampPalette(c("red", "white"))(100)
  } else if (color_type == "red_white_blue") {
    colorRange <- colorRampPalette(c("red", "white", "blue"))(100)
  } else if (color_type == "white_red") {
    colorRange <- colorRampPalette(c("white", "red"))(100)
  } else if (color_type == "white_blue") {
    colorRange <- colorRampPalette(c("white", "blue"))(100)
  } else if (color_type == "gray_scale") {
    colorRange <- colorRampPalette(c("gray90", "black"))(100)
  } else if (color_type == "EA_bin") {
    colorRange <- c(rep("#90d1e2", 30),
                    rep("#ffba57", 40),
                    rep("#ff0000", 30))
  }
  if (prefix == "0x") {
    colorRange <- colorRange %>%
      `str_sub<-`(start = 1, end = 1, value = "0x")
  }
  return(colorRange)
}
