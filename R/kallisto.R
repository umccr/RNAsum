#' Get Kallisto Counts
#'
#' Get kallisto counts via tximport.
#' @param x Path to `abundance.tsv` file with abundances. See [tximport::tximport].
#' @param tx2gene data.frame with tx_name and gene_id columns. See [tximport::tximport].
#' @return Tibble with the counts per gene transcript, or NULL if any of the
#'         input params are NULL.
#' @examples
#' x <- system.file("rawdata/test_data/quant/abundance.tsv", package = "RNAsum")
#' tx2gene <- NULL
#' (kc <- kallisto_counts(x, tx2gene)) # NULL since no tx2gene specified
#' @testexamples
#' expect_null(sc)
#' @export
kallisto_counts <- function(x, tx2gene = NULL) {
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
  # kallisto outputs transcript level counts - we covert these to gene level and
  # filter gene names with "PAR_Y" to avoid duplicate row.names issues when combining dataset
  txi_kallisto <- tximport::tximport(files = x, type = "kallisto", tx2gene = tx2gene)
  counts <- txi_kallisto[["counts"]] |>
    tibble::as_tibble(rownames = "rowname", .name_repair = make.names) |>
    dplyr::rename(count = X) |>
    dplyr::filter(!grepl("PAR_Y", .data$rowname))

  counts <- counts |>
    dplyr::mutate(
      rowname = sub("\\..*", "", .data$rowname),
      count = as.integer(.data$count)
      )
  return(counts)
}
