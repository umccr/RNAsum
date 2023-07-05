#' Wrapper for saveWidget bug fix
#'
#' A wrapper to saveWidget which compensates for bug where the temporary
#' directory with widget dependencies is not deleted - see
#' https://github.com/ramnathv/htmlwidgets/issues/296.
#'
#' @param widget Widget to save.
#' @param file Full path to save HTML into.
#' @param selfcontained Whether to save the HTML as a single self-contained file
#' (with external resources base64 encoded) or a file with external resources
#' placed in an adjacent directory.
#' @param ... Additional arguments passed to [htmlwidgets::saveWidget].
#'
#' @return Change work dir to compensate for saveWidget bug
#' @export
saveWidgetFix <- function(widget, file, selfcontained = TRUE, ...) {
  file_dir <- base::dirname(file)
  fs::dir_create(file_dir)
  file_bname <- base::basename(file)
  file_bname_noext <- tools::file_path_sans_ext(file_bname)
  dir2del <- file.path(file_dir, glue::glue("{file_bname_noext}_files"))
  htmlwidgets::saveWidget(widget = widget, file = file, selfcontained = selfcontained, ...)
  base::unlink(dir2del, recursive = TRUE)
}
