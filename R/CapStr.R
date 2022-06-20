#' Capitalise input string
#'
#' Capitalises input string.
#'
#' @param x Input string.
#'
#' @return Capitalised string.
#' @export
#'
CapStr <- function(x) {
  c <- base::strsplit(x, " ")[[1]]
  paste(base::toupper(base::substring(c, 1,1)), base::substring(c, 2),
        sep="", collapse=" ")
}
