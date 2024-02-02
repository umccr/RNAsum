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

#' Get PCGR genes with hits
#'
#' @param pcgr Path to parsed PCGR tiered TSV tibble if available.
#' @param expr_data.z Tibble with gene symbol and z-score expr values.
#' @param expr_data.perc Tibble with gene symbol and percentile expr values.
#'
#' @return List.
#' @export
pcgr_expr <- function(pcgr = NULL, expr_data.z, expr_data.perc) {
  if (is.null(pcgr)) {
    # just return the input expr named vectors
    return(
      list(
        gene.mut = NULL,
        expr_data.perc = expr_data.perc,
        expr_data.z = expr_data.z
      )
    )
  }
  pcgr <- pcgr |>
    dplyr::select("SYMBOL", "CONSEQUENCE") |>
    dplyr::filter(!.data$SYMBOL == "")
  # named vector with gene and csq
  pcgr_csq_hits <- pcgr |>
    dplyr::group_by(.data$SYMBOL) |>
    dplyr::mutate(
      nr = dplyr::n(),
      res = dplyr::if_else(
        .data$nr == 1,
        .data$CONSEQUENCE,
        "multiple hits"
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::select("SYMBOL", "res") |>
    dplyr::distinct() |>
    tibble::deframe()

  gene.mut1 <- tibble::tibble(
    gene = expr_data.z[["SYMBOL"]]
  ) |>
    dplyr::mutate(
      csq_hits = pcgr_csq_hits[.data$gene],
      alterations = dplyr::if_else(is.na(.data$csq_hits), "None", paste0("Mutation: ", .data$csq_hits))
    ) |>
    dplyr::select("gene", "alterations")

  ## If there is no expression value for a specific gene then assume it's not
  ## expressed at all and assign the lowest value observed in that sample.
  ## NOTE: these will actually get removed during the intersect call in the Rmd...
  # - get pcgr genes that are not in gene.mut
  # - assign them the min expr_data value
  # - add to expr_data tbls and gene.mut
  gene.mut2 <- pcgr |>
    dplyr::filter(!.data$SYMBOL %in% gene.mut1$gene) |>
    dplyr::select(gene = "SYMBOL") |>
    dplyr::distinct() |>
    dplyr::mutate(alterations = pcgr_csq_hits[.data$gene])

  # use second column of expr_data for the values
  e <- gene.mut2 |>
    dplyr::select("gene") |>
    dplyr::mutate(
      expr_min_z = min(expr_data.z[[2]]),
      expr_min_perc = min(expr_data.perc[[2]])
    )

  e_perc <- e |>
    dplyr::select("gene", "expr_min_perc") |>
    rlang::set_names(colnames(expr_data.perc))
  e_z <- e |>
    dplyr::select("gene", "expr_min_z") |>
    rlang::set_names(colnames(expr_data.z))
  gene.mut <- dplyr::bind_rows(gene.mut1, gene.mut2)
  expr_data.perc <- dplyr::bind_rows(expr_data.perc, e_perc)
  expr_data.z <- dplyr::bind_rows(expr_data.z, e_z)

  list(
    gene.mut = gene.mut,
    expr_data.perc = expr_data.perc,
    expr_data.z = expr_data.z
  )
}

#' Copy Number Subset Table
#'
#' @param gene.mut Tibble with gene and alterations, if PCGR was run.
#' @param cn_data Tibble with gene and MeanCopyNumber
#' @param expr_data.z Tibble with SYMBOL and z-score expr values.
#' @param expr_data.perc Tibble with SYMBOL and percentile expr values.
#'
#' @return A tibble with Gene, CN, Diff_Perc, Diff_Z_score, and Alterations if
#' PCGR was run.
#' @export
cn_subset <- function(gene.mut = NULL, cn_data, expr_data.z, expr_data.perc) {
  assertthat::assert_that(
    is.data.frame(cn_data),
    is.data.frame(expr_data.perc),
    is.data.frame(expr_data.z),
    all(c("gene", "MeanCopyNumber") %in% colnames(cn_data)),
    "SYMBOL" %in% colnames(expr_data.perc),
    "SYMBOL" %in% colnames(expr_data.z)
  )
  genes.intersect <- intersect(cn_data[["gene"]], expr_data.perc[["SYMBOL"]])
  gene.mut.sub <- NULL
  if (!is.null(gene.mut)) {
    assertthat::assert_that(
      all(colnames(gene.mut) == c("gene", "alterations"))
    )
    genes.intersect <- intersect(genes.intersect, gene.mut[["gene"]])
    gene.mut.sub <- gene.mut |> dplyr::filter(.data$gene %in% genes.intersect)
  }
  cn_data.sub <- cn_data |> dplyr::filter(.data$gene %in% genes.intersect)
  expr_data.perc.sub <- expr_data.perc |> dplyr::filter(.data$SYMBOL %in% genes.intersect)
  expr_data.z.sub <- expr_data.z |> dplyr::filter(.data$SYMBOL %in% genes.intersect)
  d <- cn_data.sub |>
    dplyr::left_join(expr_data.perc.sub, by = c("gene" = "SYMBOL")) |>
    dplyr::left_join(expr_data.z.sub, by = c("gene" = "SYMBOL"), suffix = c("_Perc", "_Z_score")) |>
    dplyr::rename("CN" = "MeanCopyNumber", "Gene" = "gene", "Diff_Perc" = 3, "Diff_Z_score" = 4)
  if (!is.null(gene.mut)) {
    d <- d |>
      dplyr::left_join(gene.mut.sub, by = c("Gene" = "gene")) |>
      dplyr::rename("Alterations" = "alterations")
  }
  return(d)
}
