#' Create AnnotationHub object to extract ensmbldb
#'
#' Create AnnotationHub object to extract ensmbldb
#' @param x AnnotationHub object
#' @return Vector of gene and transcript ids
#' @examples
#' x <- get_ensembl_db(x = ah)
#' @export
get_ensembl_db <- function(edb) {
  # Find keys
  keys_ah <- AnnotationDbi::keys(edb_105, keytype = "GENEID")
  # Extract gene and transcript ids
  tx_gene_id <- edb |>
    AnnotationDbi::select(keys = keys_ah, columns = c("TXID", "GENEID"), keytype = "GENEID") |>
    dplyr::rename(tx_name = "TXID", gene_id = "GENEID")
  tx_gene_id
}

#' Extract AnnotationHub keys
#'
#' Extract AnnotationHub keys
#' @param record AnnotationHub record number
#' @return Vector of AnnotationHub keys
#' @examples
#' x <- get_ah_keys(record = "AH98047")
#' @export
get_edb <- function(record) {
  ah <- AnnotationHub()
  # example of query command
  #ahDb <- AnnotationHub::query(x, pattern = c("Homo Sapiens", "EnsDb"))
  # Download db - this takes sometime - needs to be written under inst/extdata/
  edb_105 <- ah[[record, force=TRUE]]
}

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
  assertthat::assert_that(dataset %in% names(RNAsum::REFERENCE_DATASETS))
  refdata_dir <- system.file("rawdata", package = "RNAsum")
  d_clean <- base::strsplit(dataset, split = "-", fixed = TRUE)[[1]][1]

  fn <- list(
    "ext_ref" = list(
      counts = file.path(refdata_dir, "ref_data", paste0("TCGA_", d_clean, "_Counts.exp.gz")),
      target = file.path(refdata_dir, "ref_data", paste0("TCGA_", dataset, "_Target.txt")),
      dataset = paste0(d_clean, " (TCGA)")
    ),
    "int_ref" = list(
      counts = file.path(refdata_dir, "ref_data", "UMCCR_PDAC_Counts.exp.gz"),
      target = file.path(refdata_dir, "ref_data", "UMCCR_PDAC_Target.txt"),
      dataset = "PAAD (UMCCR)"
    )
  )
  fn
}

#' Get Reference Genes
#'
#' Get a list of reference genes.
#' @param p RNAsum params list.
#' @return List of tibbles corresponding to cancer, oncokb, immune and hrd genes
#'
#' @examples
#' p <- list(
#'   genes_cancer = system.file("rawdata/genes/umccr_cancer_genes.2019-03-20.tsv",
#'     package = "RNAsum"
#'   ),
#'   genes_oncokb = system.file(
#'     "rawdata/OncoKB/CancerGenesList.txt",
#'     package = "RNAsum"
#'   ),
#'   civic_var_summaries = system.file(
#'     "rawdata/CIViC/01-Oct-2018-VariantSummaries.tsv",
#'     package = "RNAsum"
#'   )
#' )
#' x <- get_refgenes(p)
#' @testexamples
#' expect_equal(nrow(x[["genes_cancer"]]), 1248)
#' @export
get_refgenes <- function(p) {
  .read <- function(x, backup, ...) {
    x <- x %||% backup
    x |>
      readr::read_tsv(col_types = readr::cols(...))
  }
  genes_cancer <- system.file("rawdata/genes/umccr_cancer_genes.2019-03-20.tsv", package = "RNAsum")
  genes_oncokb <- system.file("rawdata/OncoKB/CancerGenesList.txt", package = "RNAsum")
  genes_immune_markers <- system.file("rawdata/genes/Genes_immune_markers.txt", package = "RNAsum")
  genes_immunogram <- system.file("rawdata/genes/Genes_immunogram.txt", package = "RNAsum")
  genes_hrd <- system.file("rawdata/genes/Genes_HRD.txt", package = "RNAsum")
  oncokb_clin_vars <- system.file("rawdata/OncoKB/allActionableVariants.txt", package = "RNAsum")
  oncokb_all_vars <- system.file("rawdata/OncoKB/allAnnotatedVariants.txt", package = "RNAsum")
  civic_var_summaries <- system.file("rawdata/CIViC/01-Oct-2018-VariantSummaries.tsv", package = "RNAsum")
  civic_clin_evid <- system.file("rawdata/CIViC/01-Oct-2018-ClinicalEvidenceSummaries.tsv", package = "RNAsum")
  cancer_biomarkers_trans <- system.file("rawdata/cancer_biomarkers_database/cancer_genes_upon_trans.tsv", package = "RNAsum")
  FusionGDB <- system.file("rawdata/FusionGDB/TCGA_ChiTaRS_combined_fusion_ORF_analyzed_gencode_h19v19_fgID.txt", package = "RNAsum")

  genes_cancer <- p[["genes_cancer"]] |>
    .read(backup = genes_cancer, .default = "l", driver = "d", n = "i", symbol = "c", sources = "c")
  genes_oncokb <- p[["genes_oncokb"]] |>
    .read(backup = genes_oncokb, .default = "c", "# of occurence within resources" = "i") |>
    dplyr::rename(Hugo_Symbol = "Hugo Symbol")
  genes_hrd <- p[["genes_hrd"]] |>
    .read(backup = genes_hrd, SYMBOL = "c")
  genes_immune_markers <- p[["genes_immune_markers"]] |>
    .read(backup = genes_immune_markers, .default = "c")
  genes_immunogram <- p[["genes_immunogram"]] |>
    .read(backup = genes_immunogram, .default = "c")
  oncokb_clin_vars <- p[["oncokb_clin_vars"]] |>
    .read(backup = oncokb_clin_vars, .default = "c")
  oncokb_all_vars <- p[["oncokb_all_vars"]] |>
    .read(backup = oncokb_all_vars, .default = "c")
  civic_var_summaries <- p[["civic_var_summaries"]] |>
    .read(
      backup = civic_var_summaries, .default = "c",
      start = "i", stop = "i", start2 = "i", stop2 = "i"
    )
  civic_clin_evid <- p[["civic_clin_evid"]] |>
    .read(
      backup = civic_clin_evid, .default = "c",
      start = "i", stop = "i", start2 = "i", stop2 = "i"
    )
  cancer_biomarkers_trans <- p[["cancer_biomarkers_trans"]] |>
    .read(backup = cancer_biomarkers_trans, .default = "c")
  FusionGDB <- FusionGDB |>
    readr::read_tsv(
      col_names = c("Hgene", "HgeneID", "Tgene", "TgeneID", "FGname", "FGID"),
      col_types = readr::cols(.default = "c")
    )

  list(
    genes_cancer = genes_cancer,
    genes_oncokb = genes_oncokb,
    genes_hrd = genes_hrd,
    genes_immune = list(
      genes_immune_markers = genes_immune_markers,
      genes_immunogram = genes_immunogram
    ),
    oncokb_clin_vars = oncokb_clin_vars,
    oncokb_all_vars = oncokb_all_vars,
    civic_var_summaries = civic_var_summaries,
    civic_clin_evid = civic_clin_evid,
    cancer_biomarkers_trans = cancer_biomarkers_trans,
    FusionGDB = FusionGDB
  )
}

#' RNAsum Reference Datasets
#'
#' Reference datasets available in RNAsum.
#'
#' @usage data(REFERENCE_DATASETS)
#' @docType data
#' @format A list of lists with the project abbreviation code, name,
#'         tissue code and number of samples for each of the 40 datasets.
"REFERENCE_DATASETS"
