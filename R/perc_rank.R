#' Convert a vector of numbers into corresponding vector of their percentiles
#'
#' Convert a vector of numbers into corresponding vector of their percentiles
#'
#' @param x Input vector
#'
#' @return Vector of percentiles
#' @export
perc_rank <- function(x) {
  base::trunc(base::rank(x))*100/base::length(x)
}
