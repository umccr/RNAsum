#' Get Salmon Counts
#'
#' Get salmon counts via tximport.
#' @param x Path to `quant.sf` file with abundances. See [tximport::tximport].
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
  # check if list.files returns a value or not based on 'gene' pattern in salmon output name
  if (grepl("genes.sf", basename(x), fixed = TRUE)) {
    counts <- read.table(x, sep="\t", as.is=TRUE, header=TRUE) |>
      dplyr::select(Name, NumReads) |>
      dplyr::rename(rowname = Name, count = NumReads) |>
      dplyr::mutate(rowname = gsub("\\..*", "", rowname))
  } else {
    txi_salmon <- tximport::tximport(files = x, type = "salmon", tx2gene = tx2gene)
    counts <- txi_salmon[["counts"]] |>
      tibble::as_tibble(rownames = "rowname", .name_repair = make.names) |>
      dplyr::rename(count = .data$X) |>
      dplyr::mutate(rowname = gsub("\\..*", "", rowname))
  }
  counts
}
