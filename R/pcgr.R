#' Read PCGR Tiers TSV File
#'
#' Reads the `snvs_indels.tiers.tsv` file output by PCGR.
#' @param x Path to PCGR `snvs_indels.tiers.tsv` file.
#' @return A tibble with the contents of the input TSV file, or NULL if x is NULL.
#'
#' @examples
#' x <- system.file(
#'   "rawdata/test_data/umccrised/test_sample_WGS/small_variants",
#'   "TEST-somatic.pcgr.snvs_indels.tiers.tsv",
#'   package = "RNAsum"
#' )
#' (ptsv <- pcgr_tiers_tsv_read(x))
#' @testexamples
#' expect_null(pcgr_tiers_tsv_read(x = NULL))
#' @export
pcgr_tiers_tsv_read <- function(x = NULL) {
  if (is.null(x)) {
    return(NULL)
  }
  d <- x |>
    readr::read_tsv(
      col_types = readr::cols(
        .default = "c"
      )
    )

  assertthat::assert_that(
    all(c("CONSEQUENCE", "TIER", "AF_TUMOR") %in% colnames(d))
  )

  d |>
    dplyr::mutate(
      CONSEQUENCE = gsub("_variant", "", .data$CONSEQUENCE),
      CONSEQUENCE = gsub("_", " ", .data$CONSEQUENCE),
      TIER = gsub("TIER ", "", .data$TIER),
      AF_TUMOR = round(as.numeric(.data$AF_TUMOR), digits = 2)
    )
}
