#' Fusion visualisation
#'
#' Fusion output visualisation.
#'
#' @param arriba_file Output fusions from Arriba.
#' @param arriba_results PDF files.
#' @param results_dir Results directory.
#'
#' @return PNG images.
#' @export
arriba_plots <- function(arriba_file, arriba_results, results_dir) {

  ##### Get path to fusion visualisation  pdf file
  arriba_dir <- unlist(strsplit(arriba_file, split='/', fixed=TRUE))
  arriba_plots.pdf <- list.files(paste(arriba_dir[1:length(arriba_dir)-1], collapse = "/"), pattern="\\.pdf$")
  arriba_dir <- paste(arriba_dir[1:length(arriba_dir)-1], collapse = "/")
  arriba_plots.pdf <- paste(arriba_dir, arriba_plots.pdf, sep = "/")

  ##### Create directory for results
  if ( !file.exists(results_dir) ) {
    dir.create(results_dir, recursive=TRUE)
  }

  ##### Export pdf images to png
  for ( i in 1:nrow(arriba_results) ) {
    arriba_plots.png <- gsub(":", ".", paste0(results_dir, "/", make.names(paste(arriba_results$X.gene1[i], arriba_results$gene2[i], sep = "__")), "_", arriba_results$breakpoint1[i], "-", arriba_results$breakpoint2[i], ".png"))
    fusion <- pdftools::pdf_render_page(arriba_plots.pdf, page = i, dpi = 300, numeric = TRUE, opw = "", upw = "")
    png::writePNG(fusion, arriba_plots.png)
  }

  ##### Clean the space
  rm(arriba_plots.pdf, arriba_plots.png, fusion)

  #### Clear plots to free up some memory
  if(!is.null(grDevices::dev.list())) invisible(grDevices::dev.off())
}
