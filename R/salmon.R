#' Get Salmon Counts
#'
#' Get salmon counts via tximport.
#' @param x Path to `quant.sf` or `quant.genes.sf` file with abundances. See [tximport::tximport].
#' @param tx2gene data.frame with tx_name and gene_id columns. See [tximport::tximport].
#' @return Tibble with the counts per gene transcript, or NULL if any of the
#'         input params are NULL.
#' @examples
#' x <- system.file("rawdata/test_data/dragen/TEST.quant.sf", package = "RNAsum")
#' tx2gene <- NULL
#' (sc <- salmon_counts(x, tx2gene)) # NULL since no tx2gene specified
#' @testexamples
#' expect_null(sc)
#' @export
salmon_counts <- function(x, tx2gene = NULL) {
  if (is.null(x) || is.null(tx2gene)) {
    return(NULL)
  }
  assertthat::assert_that(
    file.exists(x), length(x) == 1,
    msg = glue::glue("File {x} does not exist!")
  )
  assertthat::assert_that(
    all(colnames(tx2gene) == c("tx_name", "gene_id")),
    msg = "The tx2gene df object has incorrect column names."
  )
  # check if gene level counts are provided or transcript level - we filter gene names with "PAR_Y" to avoid
  # duplicate row.names issues when combining dataset
  if (grepl("genes.sf", basename(x), fixed = TRUE)) {
    counts <- readr::read_tsv(x, col_types = readr::cols(.default = "c", NumReads = "d")) |>
      dplyr::select("Name", "NumReads") |>
      dplyr::rename(rowname = "Name", count = "NumReads") |>
      dplyr::filter(!grepl("PAR_Y", .data$rowname))
  } else {
    txi_salmon <- tximport::tximport(files = x, type = "salmon", tx2gene = tx2gene)
    counts <- txi_salmon[["counts"]] |>
      tibble::as_tibble(rownames = "rowname", .name_repair = make.names) |>
      dplyr::rename(count = .data$X)
  }
  counts <- counts |>
    dplyr::mutate(rowname = sub("\\..*", "", .data$rowname))
  return(counts)
}
