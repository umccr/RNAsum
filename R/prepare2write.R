#' Prepare object to write into a file
#'
#' Prepare object to write into a file
#'
#' @param x Input object.
#'
#' @return Object for writing to a file.
#' @export
prepare2write <- function(x) {
  x2write <- base::cbind(base::rownames(x), x)
  colnames(x2write) <- c("", base::colnames(x))

  return(x2write)
}
