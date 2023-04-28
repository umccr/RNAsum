#' Read DRAGEN Fusions
#'
#' Reads the `fusion_candidates.final` file output by DRAGEN.
#' @param x Path to the `fusion_candidates.final` file.
#' @return A tibble with the contents of the input TSV file, or NULL if input is NULL or TSV has no data rows.
#' @examples
#' x <- system.file("rawdata/test_data/dragen",
#'   "test_sample_WTS.fusion_candidates.final",
#'   package = "RNAsum"
#' )
#' (a <- dragen_fusions_read(x))
#' @testexamples
#' expect_equal(colnames(a)[ncol(a)], "ReadNames")
#' expect_equal(colnames(a)[1], "gene1")
#' expect_null(dragen_fusions_read(NULL))
#' @export
dragen_fusions_read <- function(x = NULL) {
  if (is.null(x)) {
    return(NULL)
  }
  ctypes <- readr::cols(
    .default = "c", Score = "d", NumSplitReads = "d",
    NumSoftClippedReads = "d", NumPairedReads = "d"
  )
  d <- readr::read_tsv(x, col_types = ctypes)
  if (nrow(d) == 0) {
    return(NULL)
  }
  # Old Dragen version (e.g. 3.7.5) had only:
  # #FusionGene, Score, Left/RightBreakpoint, ReadNames columns.
  # New version (e.g. 3.9.3) also has:
  # Gene1/2Location, Gene1/2Sense, Gene1/2Id, NumSplitReads, NumSoftClippedReads, NumPairedReads
  # If you try to read an old file with the above readr spec, you'll just get a
  # warning about the missing columns provided in ctypes (which is okay!).
  d |>
    dplyr::rename(FusionGene = "#FusionGene") |>
    tidyr::separate_wider_delim(
      cols = "FusionGene", delim = "--", names = c("gene1", "gene2"),
      cols_remove = TRUE, too_few = "align_start", too_many = "merge"
    )
}

#' Process DRAGEN Fusions
#'
#' @param dragen.fusions Parsed tibble with DRAGEN fusions.
#' @param known_translocations Tibble with known translocations.
#' @param genes_cancer Character vector of cancer genes.
#' @return Tibble with post-processed DRAGEN fusion calls.
#' @export
dragen_fusions_process <- function(dragen.fusions, known_translocations, genes_cancer) {
  d <- fusions_preprocess(d = dragen.fusions, known_translocations = known_translocations, genes_cancer = genes_cancer)
  d |>
    dplyr::mutate(fusion_caller = "dragen") |>
    dplyr::arrange(
      .data$reported_fusion, .data$Score, .data$fusions_cancer,
      .data$reported_fusion_geneA, .data$reported_fusion_geneB
    ) |>
    # just reverse order
    dplyr::arrange(-dplyr::row_number()) |>
    dplyr::select(
      "geneA", "geneB", "Score", "LeftBreakpoint", "RightBreakpoint",
      "GeneALocation", "GeneBLocation", "NumSplitReads",
      "NumSoftClippedReads", "FGID", "reported_fusion", "reported_fusion_geneA",
      "reported_fusion_geneB", "effector_gene", "fusions_cancer", "fusion_caller", "gene_pair"
    ) |>
    dplyr::rename(
      breakpointA = "LeftBreakpoint",
      breakpointB = "RightBreakpoint",
      siteA = "GeneALocation",
      siteB = "GeneBLocation",
      split_reads = "NumSplitReads",
      soft_clipped_reads = "NumSoftClippedReads",
      score = "Score"
    )
}
