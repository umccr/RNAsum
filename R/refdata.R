#' Get Reference Data File Paths
#'
#' Get a list of paths to internal and external reference data, based on a
#' selected dataset name.
#' @param dataset Reference RNAsum dataset of interest.
#'
#' @return List with paths to internal and external reference datasets
#'         (Counts, Targets and Name).
#'
#' @examples
#' x <- get_refdata(dataset = "TEST")
#' @export
get_refdata <- function(dataset) {
  assertthat::assert_that(dataset %in% REFERENCE_DATASETS)
  refdata_dir <- system.file("rawdata", package = "RNAsum")
  d_clean <- base::strsplit(dataset, split = "-", fixed = TRUE)[[1]][1]
  list(
    "ext_ref" = c(
      file.path(refdata_dir, "ref_data", paste0("TCGA_", d_clean, "_Counts.exp.gz")),
      file.path(refdata_dir, "ref_data", paste0("TCGA_", dataset, "_Target.txt")),
      paste0(d_clean, " (TCGA)")
    ),
    "int_ref" = c(
      file.path(refdata_dir, "ref_data", "UMCCR_PDAC_Counts.exp.gz"),
      file.path(refdata_dir, "ref_data", "UMCCR_PDAC_Target.txt"),
      "PAAD (UMCCR)"
    )
  )
}

#' Get Reference Genes
#'
#' Get a list of reference genes.
#' @param p RNAsum params list.
#' @return List of tibbles corresponding to cancer, oncokb, immune and hrd genes
#'
#' @examples
#' p <- list(
#'   genes_cancer = system.file("rawdata/genes/umccr_cancer_genes.2019-03-20.tsv", package = "RNAsum"),
#'   genes_oncokb = system.file("rawdata/OncoKB/CancerGenesList.txt", package = "RNAsum"),
#'   genes_immune_markers = system.file("rawdata/genes/Genes_immune_markers.txt", package = "RNAsum"),
#'   genes_hrd = system.file("rawdata/genes/Genes_HRD.txt", package = "RNAsum")
#' )
#' x <- get_refgenes(p)
#' @testexamples
#' expect_equal(length(x), 4)
#' expect_null(x$foo)
#' @export
get_refgenes <- function(p) {
  .read <- function(x, ...) {
    if (is.null(x)) {
      return(NULL)
    }
    x |>
      readr::read_tsv(col_types = readr::cols(...))
  }
  genes_cancer <- p[["genes_cancer"]] |>
    .read(.default = "l", driver = "d", n = "i", symbol = "c", sources = "c")
  genes_oncokb <- p[["genes_oncokb"]] |>
    .read(.default = "c", "# of occurence within resources" = "i")
  genes_hrd <- .read(p[["genes_hrd"]], SYMBOL = "c")
  genes_immune_markers <- .read(p[["genes_immune_markers"]], .default = "c")
  genes_immunogram <- .read(p[["genes_immunogram"]], .default = "c")

  list(
    genes_cancer = genes_cancer,
    genes_oncokb = genes_oncokb,
    genes_hrd = genes_hrd,
    genes_immune = list(
      genes_immune_markers = genes_immune_markers,
      genes_immunogram = genes_immunogram
    )
  )
}

#' RNAsum Reference Datasets
#'
#' Reference datasets available in RNAsum.
#'
#' @format A list of lists with the project abbreviation code, name,
#'         tissue code and number of samples for each of the 40 datasets.
"REFERENCE_DATASETS"
