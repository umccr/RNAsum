#' Get Salmon Counts
#'
#' Returns salmon counts.
#' @param x Path to `quant.sf` file with abundances. See [tximport::tximport].
#' @param tx2gene See [tximport::tximport]
#' @return Tibble with the counts per gene transcript.
#' @export
salmon_counts <- function(x, tx2gene) {
  assertthat::assert_that(file.exists(x), length(x) == 1,
    msg = glue::glue("File {x} does not exist!")
  )
  txi.salmon <- tximport::tximport(files = x, type = "salmon", tx2gene = tx2gene)
  counts <- tibble::as_tibble(txi.salmon[["counts"]], rownames = "rowname", .name_repair = make.names) |>
    dplyr::rename(count = .data$X)
  counts
}
