##### Function to generate a full-resolution pdf image before generating a small image in the chunk (from https://stackoverflow.com/questions/37834053/what-is-a-simple-way-to-thumbnail-some-plots-in-r-markdown-knitr )
#' Generate a full-resolution pdf image before generating a smaller image
#'
#' @param x A pdf image
#' @param options
#'
#' @return A full resolution pdf image
#' @export
#'

allow_thumbnails <- function(x, options) {
  if (!is.null(options$thumb)) {
    filename <- sprintf("%s.full.pdf", strsplit(basename(x), "\\.")[[1]][1])
    absolute_path <- file.path(dirname(x), filename)

    ##### Generate the full resolution pdf
    grDevices::pdf(absolute_path, width = options$thumb$width, height = options$thumb$height)
    eval(parse(text = options$code))
    grDevices::dev.off()

    ##### Add an html link to the low resolution png
    options$fig.link = absolute_path
  }

  knitr:::hook_plot_md_base(x, options)
}
