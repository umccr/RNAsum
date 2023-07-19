#' Generate UMCCR RNAsum Report
#'
#' Generates a UMCCR RNAsum HTML report. It does so with the following steps:
#' 1. copy the rmd into tmp/rnasum.Rmd
#' 2. render the rmd inside tmp/
#' 3. return the path to the output HTML
#'
#' @param out_file Path to output HTML file (needs '.html' suffix).
#' @param quiet Suppress printing during rendering.
#' @param pars List of named parameters that override custom params specified within the Rmd's YAML front-matter.
#'
#' @return Path to rendered HTML report.
#' @export
rnasum_rmd <- function(out_file = NULL, quiet = FALSE, pars) {
  assertthat::assert_that(
    length(out_file) == 1, !is.null(out_file),
    is.logical(quiet),
    is.list(pars),
    tools::file_ext(out_file) == "html", length(pars) > 0
  )
  tmp_dir <- tempdir()
  rmd_dir <- system.file("rmd", package = "RNAsum")
  fs::dir_copy(rmd_dir, tmp_dir, overwrite = TRUE)
  rmd_file <- file.path(tmp_dir, "rnasum.Rmd")
  out_dir <- dirname(out_file)
  fs::dir_create(out_dir)
  rmarkdown::render(
    input = rmd_file,
    params = pars,
    output_dir = out_dir,
    output_file = I(out_file),
    quiet = quiet
  )
}
