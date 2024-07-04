#' Tidy eval helpers
#'
#' <https://www.tidyverse.org/blog/2019/06/rlang-0-4-0/>
#'
#' @name tidyeval
#' @keywords internal
#' @importFrom rlang .data := %||%
NULL

#' @keywords internal
"_PACKAGE"

#' @noRd
dummy1 <- function() {
  # Solves R CMD check: Namespaces in Imports field not imported from
  EDASeq::plotRLE
  conflicted::conflict_scout
  edgeR::DGEList
  ggforce::geom_sina
  ggplot2::ggplot
  limma::removeBatchEffect
  manhattanly::manhattanly
  optparse::make_option
  preprocessCore::normalize.quantiles
  ragg::agg_png
}
