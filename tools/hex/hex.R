# tools/hex/hex.R
library(hexSticker)
library(fs)

# --- 1. Build the sticker -----------------------------------------------

sticker(
  subplot   = "tools/hex/stage_icon.png",
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
  filename  = "inst/figures/logo.png",
  dpi       = 300
)


# --- 2. Ensure pkgdown gets copies --------------------------------------

dir_create("man/figures")

file_copy("inst/figures/logo.png",  "man/figures/logo.png",  overwrite = TRUE)
file_copy("tools/hex/stage_icon.png", "man/figures/stage_icon.png", overwrite = TRUE)

message("Logo and icon copied into man/figures for pkgdown.")
