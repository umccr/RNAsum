#' Process CGI Known Translocations
#'
#' Processes parsed CGI known translocations tibble, removing duplicated
#' translocations that have a different source.
#' @param kt_cgi Parsed CGI tibble.
#'
#' @return Tibble with dups removed.
#' @export
known_translocations_cgi_process <- function(kt_cgi) {
  # merge duplicated fusions by source (https://github.com/umccr/RNAsum/issues/89)
  dup_id_cols <- c("translocation", "effector_gene", "cancer_acronym")
  assertthat::assert_that(
    inherits(kt_cgi, "data.frame"),
    all(dup_id_cols %in% colnames(kt_cgi))
  )
  kt_cgi_dup <- kt_cgi |>
    dplyr::group_by(.data$translocation) |>
    dplyr::filter(n() > 1) |>
    dplyr::ungroup() |>
    tidyr::pivot_wider(
      id_cols = dplyr::all_of(dup_id_cols),
      names_from = "source", values_from = "source"
    )
  # handle potential no dup findings
  if (nrow(kt_cgi_dup) > 0) {
    kt_cgi_dup <- kt_cgi_dup |>
      tidyr::unite(!dplyr::all_of(dup_id_cols), col = "source", sep = ";")
  } else {
    kt_cgi_dup <- RNAsum::empty_tbl(cnames = colnames(kt_cgi), ctypes = readr::spec(kt_cgi))
  }
  # now filter out dups and re-bind with clean tibble
  kt_cgi |>
    dplyr::filter(!.data$translocation %in% kt_cgi_dup[["translocation"]]) |>
    dplyr::bind_rows(kt_cgi_dup) |>
    dplyr::mutate(
      cancer_acronym = gsub(";", ", ", .data$cancer_acronym),
      source = gsub(";", ", ", .data$source),
      translocation = gsub("__", "_", .data$translocation)
    ) |>
    dplyr::arrange(.data$translocation)
}


#' @export
known_trans <- function(kt_fusiongdb, kt_cgi) {
  ##### Flag known fusions based on info from:
  ##### FusionGDB (https://ccsm.uth.edu/FusionGDB)
  ##### Cancer Biomarkers database (CGI) (https://www.cancergenomeinterpreter.org/biomarkers)
  kt <- dplyr::full_join(
    kt_fusiongdb, kt_cgi,
    by = c("FGname" = "translocation")
  ) |>
    tidyr::separate_wider_delim(
      cols = "FGname", delim = "_", names = c("geneA", "geneB"),
      cols_remove = FALSE, too_few = "align_start", too_many = "merge"
    ) |>
    tidyr::unite(col = "FGname2", "geneA", "geneB", sep = "-", remove = FALSE) |>
    dplyr::select(
      "FGname", "FGname2", "Hgene", "HgeneID", "Tgene", "TgeneID",
      "FGID", "effector_gene", "cancer_acronym", "source", "geneA", "geneB"
    )
  # NOTE: when geneB is missing we'll get NA e.g. COX6C-NA instead of COX6C-COX6C
  # NOTE: FGname2 is basically the old trans.pairs
  kt
}
