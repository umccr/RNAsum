##### Check for nearest value in a vector
#' Check for nearest value in a vector
#'
#' @param vector Input vector
#' @param x Value to search
#'
#' @return Nearest value in a vector
#' @export
#'
nearest_position <- function(vector, x) {

  y <- BiocGenerics::which.min(abs(vector - x))

  ##### Clean the space and return output
  rm(vector, x)
  return( y )
}
