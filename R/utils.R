#' Create Empty Tibble Given Column Names
#'
#' From https://stackoverflow.com/a/62535671/2169986. Useful for handling
#' edge cases with empty data.
#'
#' @param ctypes Character vector of column types corresponding to `cnames`.
#' @param cnames Character vector of column names to use.
#'
#' @return A tibble with 0 rows and the given column names.
#' @export
empty_tbl <- function(cnames, ctypes = readr::cols(.default = "c")) {
  readr::read_csv("\n", col_names = cnames, col_types = ctypes)
}
