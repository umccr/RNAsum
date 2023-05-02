#' Process Manta Object
#'
#' @param manta_tsv_obj Manta list object read via `gpgr::process_sv`.
#'
#' @return List with:
#' - total Manta (unmelted) variants
#' - tibble with melted variants
#' - genes involved in multi-gene events
#' @export
manta_process <- function(manta_tsv_obj) {
  assertthat::assert_that(all(c("melted", "unmelted") %in% names(manta_tsv_obj)))
  total_variants <- nrow(manta_tsv_obj[["unmelted"]])
  melted <- manta_tsv_obj[["melted"]]
  multigenes <- melted |>
    dplyr::mutate(multigene = grepl(", ", .data$Genes)) |>
    dplyr::filter(.data$multigene)
  if (nrow(multigenes) > 0) {
    multigenes <- multigenes |>
      dplyr::rowwise() |>
      dplyr::mutate(g = list(unlist(strsplit(.data$Genes, split = ", ", fixed = TRUE)[[1]]))) |>
      dplyr::pull("g") |>
      unlist() |>
      unique() |>
      tibble::as_tibble_col(column_name = "Genes")
  } else {
    multigenes <- empty_tbl(cnames = "Genes")
  }

  list(
    total_variants = total_variants,
    melted_variants = melted,
    multigenes = multigenes
  )
}
