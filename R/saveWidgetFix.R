#' Wrapper for saveWidget bug fix
#'
#' A wrapper to saveWidget which compensates for arguable BUG in saveWidget which requires
#' `file` to be in current working directory (see https://github.com/ramnathv/htmlwidgets/issues/299).
#'
#' @param widget Input widget.
#' @param file Input file.
#' @param ... Additional arguments passed to [htmlwidgets::saveWidget].
#'
#' @return Change work dir to compensate for saveWidget bug
#' @export
saveWidgetFix <- function ( widget, file, ...) {
  wd<-getwd()
  on.exit(setwd(wd))
  outDir<-dirname(file)
  file<-basename(file)
  setwd(outDir);
  htmlwidgets::saveWidget(widget,file=file,...)
}
