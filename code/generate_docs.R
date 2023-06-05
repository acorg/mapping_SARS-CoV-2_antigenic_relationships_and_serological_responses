
# Clear generated data directories
unlink("docs/", recursive = TRUE)
dir.create("docs")

rmd_files <- list.files("code/docs/", pattern = "Rmd$")
lapply(rmd_files, function(file){
  rmarkdown::render(
    input = file.path("code", "docs", file),
    output_file = gsub("Rmd", "html", file),
    output_format = flexdashboard::flex_dashboard(),
    output_dir = "docs"
  )
})
