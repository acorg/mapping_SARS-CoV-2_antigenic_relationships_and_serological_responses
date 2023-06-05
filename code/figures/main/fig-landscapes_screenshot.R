
# Setup workspace
rm(list = ls())
source("functions/screenshot.R")

# Screenshot landscapes
lndscps_dir <- "figures/main/lndscps"
lndscp_pages <- list.files(lndscps_dir, pattern = "\\.html$")
lndscp_screenshots <- gsub("\\.html$", ".png", lndscp_pages)

for (i in seq_along(lndscp_pages)) {
  
  # Do screenshot
  screenshotWebpage(
    url = file.path(lndscps_dir, lndscp_pages[i]),
    file = file.path(lndscps_dir, lndscp_screenshots[i]),
    vwidth = 1600,
    vheight = 880,
    cliprect = c(200, 0, 1200, 880)
  )
  
}

