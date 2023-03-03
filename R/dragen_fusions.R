#' Read DRAGEN Fusions
#'
#' Reads the `fusion_candidates.final` file output by DRAGEN.
#' @param x Path to the `fusion_candidates.final` file.
#' @return A tibble with the contents of the input TSV file, or NULL if input is NULL.
#' @examples
#' x <- system.file("rawdata/test_data/dragen",
#'   "test_sample_WTS.fusion_candidates.final",
#'   package = "RNAsum"
#' )
#' (a <- dragen_fusions_read(x))
#' @testexamples
#' expect_equal(colnames(a)[ncol(a)], "ReadNames")
#' expect_equal(colnames(a)[1], "gene1")
#' expect_null(dragen_fusions_read(NULL))
#' @export
dragen_fusions_read <- function(x = NULL) {
  if (is.null(x)) {
    return(NULL)
  }
  ctypes <- readr::cols(
    .default = "c", Score = "d", NumSplitReads = "d",
    NumSoftClippedReads = "d", NumPairedReads = "d"
  )
  readr::read_tsv(x, col_types = ctypes) |>
    dplyr::rename(FusionGene = "#FusionGene") |>
    tidyr::separate(.data$FusionGene, into = c("gene1", "gene2"), sep = "--")
}
