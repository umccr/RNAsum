#' Read Sample Data
#'
#' Reads sample data, including Arriba fusions, Arriba plots, Salmon counts, and
#' DRAGEN fusions.
#'
#' @param p RNAsum params list.
#' @param results_dir Directory to output extracted Arriba PNGs to (created
#'        automatically if it does not already exist).
#' @param tx2gene data.frame with tx_name and gene_id columns. Required for salmon.
#'        See [tximport::tximport].
#'
#' @return A list of the input sample data.
#' @examples
#' p <- list(
#'   dragen_rnaseq = system.file("rawdata/test_data/dragen", package = "RNAsum"),
#'   arriba_pdf = system.file("rawdata/test_data/dragen/arriba/fusions.pdf", package = "RNAsum"),
#'   arriba_tsv = system.file("rawdata/test_data/dragen/arriba/fusions.tsv", package = "RNAsum"),
#'   salmon = NULL,
#'   dragen_fusions = system.file(
#'     "rawdata/test_data/dragen/test_sample_WTS.fusion_candidates.final",
#'     package = "RNAsum"
#'   )
#' )
#' res <- read_sample_data(p, tempdir())
#' @testexamples
#' expect_equal(length(res), 4)
#' expect_null(res$salmon)
#' @export
read_sample_data <- function(p, results_dir, tx2gene = NULL) {
  assertthat::assert_that(
    all(c("arriba_pdf", "arriba_tsv", "salmon", "dragen_fusions") %in% names(p))
  )

  arriba_tsv <- arriba_read_tsv(p[["arriba_tsv"]])
  arriba_pdf <- arriba_read_pdf(p[["arriba_pdf"]],
    fusions = arriba_tsv,
    outdir = file.path(results_dir, "arriba")
  )
  salmon <- salmon_counts(p[["salmon"]], tx2gene = tx2gene)
  dragen_fusions <- dragen_fusions_read(p[["dragen_fusions"]])
  list(
    arriba_tsv = arriba_tsv,
    arriba_pdf = arriba_pdf,
    salmon = salmon,
    dragen_fusions = dragen_fusions
  )
}
