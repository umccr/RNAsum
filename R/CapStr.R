##### Capitalize input string
#' Capitalize input string
#'
#' @param y Input
#'
#' @return Capitalized string.
#' @export
#'
CapStr <- function(y) {
  c <- base::strsplit(y, " ")[[1]]
  paste(toupper(base::substring(c, 1,1)), base::substring(c, 2),
        sep="", collapse=" ")
}
