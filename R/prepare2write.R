##### Prepare object to write into a file
#' Prepare object to write into a file
#'
#' @param x Input object.
#'
#' @return Object for writing to a file.
#' @export
#'
prepare2write <- function (x) {

  x2write <- BiocGenerics::cbind(BiocGenerics::rownames(x), x)
  colnames(x2write) <- c("",BiocGenerics::colnames(x))

  ##### Clean the space and return output
  rm(x)
  return(x2write)
}
