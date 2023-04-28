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
    dplyr::filter(dplyr::n() > 1) |>
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


#' Get Known Translocations
#'
#' Flag known fusions based on info from:
#' - FusionGDB (https://ccsm.uth.edu/FusionGDB)
#' - Cancer Biomarkers database (CGI) (https://www.cancergenomeinterpreter.org/biomarkers)
#'
#' @param kt_fusiongdb FusionGDB parsed tibble.
#' @param kt_cgi CGI parsed tibble.
#' @return Tibble with known translocations.
#' @export
known_trans <- function(kt_fusiongdb, kt_cgi) {
  # need to remove dup trans from cgi
  kt_cgi <- known_translocations_cgi_process(kt_cgi)
  kt <- dplyr::full_join(
    kt_fusiongdb, kt_cgi,
    by = c("FGname" = "translocation")
  ) |>
    tidyr::separate_wider_delim(
      cols = "FGname", delim = "_", names = c("geneA", "geneB"),
      cols_remove = FALSE, too_few = "align_start", too_many = "merge"
    ) |>
    # use a comparator gene column where the genes are sorted
    dplyr::mutate(
      gene_cmp1 = dplyr::if_else(.data$geneA < .data$geneB, .data$geneA, .data$geneB),
      gene_cmp2 = dplyr::if_else(.data$geneA < .data$geneB, .data$geneB, .data$geneA),
      gene_pair = glue::glue("{gene_cmp1}--{gene_cmp2}")
    ) |>
    dplyr::select(
      "FGname", "gene_pair", "Hgene", "HgeneID", "Tgene", "TgeneID",
      "FGID", "effector_gene", "cancer_acronym", "source", "geneA", "geneB"
    )
  # NOTE: when geneB is missing we'll get NA e.g. COX6C-NA instead of COX6C-COX6C
  kt
}

#' Pre-Process Arriba or DRAGEN Fusions
#'
#'
#' @param d Parsed tibble with Arriba or DRAGEN fusions.
#' @param known_translocations Tibble with known translocations.
#' @param genes_cancer Character vector of cancer genes.
#' @return Tibble with post-processed fusion calls.
#' @export
fusions_preprocess <- function(d, known_translocations, genes_cancer) {
  k <- known_translocations
  assertthat::assert_that(
    inherits(k, "data.frame"),
    inherits(d, "data.frame"),
    is.character(genes_cancer),
    all(c("gene1", "gene2") %in% colnames(d)),
    all(c("FGname", "gene_pair", "geneA", "geneB") %in% colnames(k))
  )
  colnames(d) <- sub("1", "A", colnames(d))
  colnames(d) <- sub("2", "B", colnames(d))
  fgid_from_fgname <- k |>
    dplyr::select("gene_pair", "FGID") |>
    tibble::deframe()
  d |>
    # use a comparator gene column where the genes are sorted
    dplyr::mutate(
      gene_cmp1 = dplyr::if_else(.data$geneA < .data$geneB, .data$geneA, .data$geneB),
      gene_cmp2 = dplyr::if_else(.data$geneA < .data$geneB, .data$geneB, .data$geneA),
      gene_pair = glue::glue("{gene_cmp1}--{gene_cmp2}")
    ) |>
    dplyr::mutate(
      reported_fusion = .data$gene_pair %in% k[["gene_pair"]],
      FGID = fgid_from_fgname[.data$gene_pair],
      FGID = dplyr::if_else(is.na(.data$FGID), fgid_from_fgname[.data$gene_pair], .data$FGID),
      reported_fusion_geneA = .data$geneA %in% c(k[["geneA"]], k[["geneB"]]),
      reported_fusion_geneB = .data$geneB %in% c(k[["geneA"]], k[["geneB"]]),
      effector_gene = any(c(.data$geneA, .data$geneB) %in% k[["effector_gene"]]),
      fusions_cancer = any(c(.data$geneA, .data$geneB) %in% genes_cancer)
    ) |>
    dplyr::select(-c("gene_cmp1", "gene_cmp2"))
}


#' Annotate Fusions with Ensembl Coordinates
#'
#' @param fusions Tibble with DRAGEN and Arriba fusions.
#' @param gene_ann Dataframe with gene annotations.
#'
#' @return Tibble with annotated genes.
#' @export
fusions_annot <- function(fusions, gene_ann) {
  sel_cols <- c("ENSEMBL", "SYMBOL", "SEQNAME", "GENESEQSTART", "GENESEQEND")
  assertthat::assert_that(all(sel_cols %in% colnames(gene_ann)))
  a <- gene_ann |>
    dplyr::select(dplyr::all_of(sel_cols)) |>
    tibble::as_tibble()
  fusions |>
    dplyr::left_join(a, by = c("geneA" = "SYMBOL")) |>
    dplyr::left_join(a, by = c("geneB" = "SYMBOL"))
}
