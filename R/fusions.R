#' @export
known_trans <- function(kt_fusiongdb, kt_cgi) {
  kt <- dplyr::full_join(
    kt_fusiongdb, kt_cgi,
    by = c("FGname" = "translocation")
  ) |>
    tidyr::separate_wider_delim(
      cols = "FGname", delim = "_", names = c("geneA", "geneB"),
      cols_remove = FALSE, too_few = "align_start", too_many = "merge"
    ) |>
    dplyr::select(
      "FGname", "Hgene", "HgeneID", "Tgene", "TgeneID",
      "FGID", "effector_gene", "cancer_acronym", "source", "geneA", "geneB"
    )
  # NOTE: when geneB is missing we'll get NA e.g. COX6C-NA instead of COX6C-COX6C
  tpairs <- kt |>
    tidyr::unite(col = "tpair", "geneA", "geneB", sep = "-") |>
    dplyr::select("tpair") # keep as single-col for now
  # dplyr::pull("tpair")

  list(
    known_translocations = kt,
    trans.pairs = tpairs
  )
}
