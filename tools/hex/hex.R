# tools/hex/hex.R
library(hexSticker)

sticker(
  subplot   = "man/figures/stage_icon.png",
  package   = "stage",
  p_size    = 30,
  p_color   = "white",
  p_y       = 1.5,
  s_x       = 1,
  s_y       = 0.92,
  s_width   = 0.72,
  s_height  = 0.72,
  h_fill    = "#051838",
  h_color   = NA,
  filename  = "man/figures/logo.png",
  dpi       = 300
)
