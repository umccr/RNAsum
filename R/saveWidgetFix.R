##### A wrapper to saveWidget which compensates for arguable BUG in saveWidget which requires `file` to be in current working directory (see post https://github.com/ramnathv/htmlwidgets/issues/299 )
saveWidgetFix <- function ( widget, file, ...) {
  wd<-getwd()
  on.exit(setwd(wd))
  outDir<-dirname(file)
  file<-basename(file)
  setwd(outDir);
  htmlwidgets::saveWidget(widget,file=file,...)
}
