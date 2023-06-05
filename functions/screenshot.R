
screenshotWebpage <- function(
    url,
    file,
    vwidth = 600,
    vheight = 600,
    cliprect = NULL
  ) {
  try({
    
    # Set up paths
    chrome <- "/Applications/Google\\ Chrome.app/Contents/MacOS/Google\\ Chrome"
    
    # Run the command
    file.create(file)
    command <- paste0(
      chrome, " --headless --screenshot='", normalizePath(file),"' --use-angle --window-size=", vwidth,",", vheight," '", normalizePath(url), "'"
    )
    system(command, ignore.stdout = TRUE, ignore.stderr = TRUE)
    
    # Crop the image
    if (!is.null(cliprect)) {
    
      image <- magick::image_read(file)
      image <- magick::image_crop(
        image = image,
        geometry = magick::geometry_area(
          width = cliprect[3]*2,
          height = cliprect[4]*2,
          x_off = cliprect[1]*2,
          y_off = cliprect[2]*2
        )
      )
      magick::image_write(image, file)
      
    }
    
    message(sprintf("Successfully screenshotted '%s'", file))
    
  })
}
