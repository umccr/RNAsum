#' Read Sample Data
#'
#' Reads sample data, including Arriba fusions, Arriba plots, Salmon counts,
#' DRAGEN fusions, DRAGEN mapping metrics, SVs, PURPLE CNVs, and PCGR SNVs.
#'
#' @param p RNAsum params list.
#' @param results_dir Directory to output extracted Arriba PNGs to (created
#'        automatically if it does not already exist).
#' @param tx2gene data.frame with tx_name and gene_id columns. Required for salmon.
#'        See [tximport::tximport].
#'
#' @return A list of the input sample data.
#' @examples
#' p <- list(
#'   dragen_wts_dir = system.file("rawdata/test_data/dragen", package = "RNAsum"),
#'   arriba_dir = system.file("rawdata/test_data/dragen/arriba", package = "RNAsum"),
#'   umccrise = system.file(
#'     "rawdata/test_data/umccrised/test_sample_WGS",
#'     package = "RNAsum"
#'   )
#' )
#' res <- read_sample_data(p, tempdir())
#' @testexamples
#' expect_equal(length(res), 6)
#' expect_null(res$salmon) # because tx2gene in NULL
#' @export
read_sample_data <- function(p, results_dir, tx2gene = NULL) {
  #---- Arriba ----#
  # if any of arriba tsv/pdf missing and arriba_dir is provided,
  # construct the paths to files from that.
  arriba_tsv <- p[["arriba_tsv"]]
  arriba_pdf <- p[["arriba_pdf"]]
  arriba_dir <- p[["arriba_dir"]]
  if (!is.null(arriba_dir) && (list(NULL) %in% list(arriba_tsv, arriba_pdf))) {
    arriba_tsv <- file.path(arriba_dir, "fusions.tsv")
    arriba_pdf <- file.path(arriba_dir, "fusions.pdf")
  }
  arriba_tsv <- arriba_tsv_read(x = arriba_tsv)
  arriba_pdf <- arriba_pdf_read(pdf = arriba_pdf, fusions = arriba_tsv, outdir = file.path(results_dir, "arriba"))

  #---- DragenWTS ----#
  # if any of salmon, fusions, mapping metrics missing and dragen_wts_dir is provided,
  # construct the paths to files from that.
  salmon <- p[["salmon"]]
  dragen_fusions <- p[["dragen_fusions"]]
  dragen_mapping_metrics <- p[["dragen_mapping_metrics"]]
  dragen_wts_dir <- p[["dragen_wts_dir"]]
  if (!is.null(dragen_wts_dir) && (list(NULL) %in% list(salmon, dragen_fusions, dragen_mapping_metrics))) {
    salmon <- list.files(dragen_wts_dir, pattern = "quant\\.genes\\.sf", full.names = TRUE)
    if (length(salmon) != 1) {
      # check for quant.sf file
      salmon <- list.files(dragen_wts_dir, pattern = "quant\\.sf", full.names = TRUE)
    }
    dragen_mapping_metrics <- list.files(dragen_wts_dir, pattern = "mapping_metrics\\.csv", full.names = TRUE)
    dragen_fusions <- list.files(dragen_wts_dir, pattern = "fusion_candidates\\.final", full.names = TRUE)
    if (length(salmon) != 1) {
      salmon <- NULL
    }
    if (length(dragen_mapping_metrics) != 1) {
      dragen_mapping_metrics <- NULL
    }
    if (length(dragen_fusions) != 1) {
      dragen_fusions <- NULL
    }
  }
  dragen_fusions <- dragen_fusions_read(dragen_fusions)
  dragen_mapping_metrics <- dragen_mapping_metrics_read(dragen_mapping_metrics)
  #---- Kallisto ----#
  kallisto <- p[["kallisto"]]

  # check which quant input is provided
  if (is.null(salmon)) {
    counts <- salmon_counts(salmon, tx2gene = tx2gene)
  } else if ((is.null(kallisto))) {
    counts <- kallisto_counts(kallisto, tx2gene = tx2gene)
  }

  #---- WGS ----#
  wgs <- read_wgs_data(p)
  list(
    arriba_tsv = arriba_tsv,
    arriba_pdf = arriba_pdf,
    salmon = salmon,
    dragen_fusions = dragen_fusions,
    dragen_mapping_metrics = dragen_mapping_metrics,
    wgs = wgs
  )
}

#' Read WGS Data
#'
#' Reads WGS data, including PCGR `tiers.tsv`, PURPLE `cnv.gene.tsv`, and
#' `sv.prioritised.tsv`. If the file path has been specified in the RNAsum params and is
#' valid, it is returned. As a fallback, if the umccrise directory param has
#' been specified, then there is an attempt to detect the file pattern in there.
#'
#' @param p RNAsum params list.
#' @return A list of the input sample data.
#' @examples
#' p <- list(
#'   umccrise = system.file("rawdata/test_data/umccrised/test_sample_WGS", package = "RNAsum"),
#'   pcgr_tiers_tsv = system.file(
#'     "rawdata/test_data/umccrised/test_sample_WGS/small_variants",
#'     "TEST-somatic.pcgr.snvs_indels.tiers.tsv",
#'     package = "RNAsum"
#'   ),
#'   sash_tsv = system.file(
#'     "rawdata/test_data/test_sample_WGS/structural/TEST.sv.prioritised.tsv",
#'     package = "RNAsum"
#'   )
#' )
#' (res <- read_wgs_data(p))
#' @testexamples
#' expect_equal(length(res), 3)
#' @export
read_wgs_data <- function(p) {
  um_dir <- p[["umccrise"]]
  .read <- function(p, subdir, pat, nm, func, ...) {
    # - If file nm provided directly in the params:
    #   - If exists, return it.
    #   - If doesn't exist or not provided, check umccrise param
    # - If umccrise has been provided, look in there for the file.
    #   - If found, return it.
    #   - If not found, check if provided directly in the params.
    #     - If found and exists, return it.
    # - Return NULL if all else fails.
    if (!is.null(p[[nm]]) && file.exists(p[[nm]])) {
      return(func(p[[nm]], ...))
    }
    if (!is.null(um_dir)) {
      x <- list.files(
        file.path(um_dir, subdir),
        pattern = pat,
        recursive = TRUE, full.names = TRUE
      )
      if (length(x) == 1) {
        return(func(x, ...))
      }
    }
    return(NULL)
  }

  pcgr_tiers_tsv <- .read(
    p = p,
    subdir = "small_variants", pat = "pcgr\\.snvs_indels\\.tiers\\.tsv$",
    nm = "pcgr_tiers_tsv", func = pcgr_tiers_tsv_read
  )

  purple_gene_tsv <- .read(
    p = p,
    subdir = "purple", pat = "purple\\.cnv\\.gene\\.tsv$",
    nm = "purple_gene_tsv", func = ppl_cnv_som_gene_read
  )

  sv_tsv <- .read(
    p = p,
    subdir = "structural", pat = "prioritised\\.tsv$",
    nm = "sv_tsv", func = sv_prioritize
  )

  list(
    pcgr_tiers_tsv = pcgr_tiers_tsv,
    purple_gene_tsv = purple_gene_tsv,
    sv_tsv = sv_tsv
  )
}

#' Get Mutated Genes from PCGR
#'
#' @param pcgr_tbl Parsed PCGR TSV tibble.
#' @param tiers PCGR tiers to keep. Default: 1, 2, 3, 4.
#' @param splice_vars Include non-coding splice region variants reported in PCGR?
#'
#' @return Character vector of genes.
#' @export
genes_pcgr_summary <- function(pcgr_tbl, tiers = c(1:4), splice_vars = TRUE) {
  assertthat::assert_that(
    all(tiers %in% 1:4),
    all(c("TIER", "SYMBOL", "CONSEQUENCE") %in% names(pcgr_tbl)),
    is.logical(splice_vars)
  )
  res1 <- pcgr_tbl |>
    dplyr::filter(.data$TIER %in% tiers) |>
    dplyr::pull("SYMBOL")
  res2 <- NULL
  if (splice_vars) {
    res2 <- pcgr_tbl |>
      dplyr::mutate(tier_csq = glue::glue("{TIER}.{CONSEQUENCE}")) |>
      dplyr::filter(
        grepl("NONCODING.*splice region", .data$tier_csq, fixed = FALSE)
      ) |>
      dplyr::pull("SYMBOL")
  }
  res <- unique(c(res1, res2)) |> stats::na.omit()
  if (length(res) == 0) {
    return(NULL)
  }
  res
}

#' Get Arriba/DRAGEN Fusions Summary
#'
#' @param tbl Tibble with fusions in columns 'gene1' and 'gene2'.
#'
#' @return Character vector of the 'gene1' and 'gene2' genes.
#' @export
fusions_summary <- function(tbl = NULL) {
  if (is.null(tbl)) {
    return(NULL)
  }
  assertthat::assert_that(
    inherits(tbl, "data.frame"), all(c("gene1", "gene2") %in% names(tbl))
  )
  if (nrow(tbl) == 0) {
    return(NULL)
  }
  res <- tbl |>
    dplyr::select("gene1", "gene2") |>
    unlist(use.names = FALSE) |>
    as.character()
  res
}

#' Get SV Summary
#'
#' @param tbl Tibble with melted SVs from sash, containing 'Genes' column.
#' @return Character vector of Genes.
#' @export
sv_summary <- function(tbl) {
  assertthat::assert_that(
    inherits(tbl, "data.frame"), (c("Genes") %in% names(tbl))
  )
  res <- tbl |>
    dplyr::filter(.data$Genes != "") |>
    dplyr::pull("Genes") |>
    base::strsplit(split = "&", fixed = TRUE) |>
    base::unlist() |>
    stats::na.omit() |>
    base::unique()

  if (length(res) == 0) {
    return(NULL)
  }
  res
}

#' Get PURPLE CNV Summary
#'
#' @param tbl PURPLE gene CNV tibble
#' @param cancer_genes_symbol Character vector of gene symbols to filter on.
#' @param cn_bottom CN threshold value to classify genes within lost regions.
#' @param cn_top CN threshold value to classify genes within gained regions.
#'
#' @return List with genes, top and bottom copy number thresholds.
#' @export
purple_cnv_summary <- function(tbl, cancer_genes_symbol, cn_bottom, cn_top) {
  assertthat::assert_that(
    inherits(tbl, "data.frame"),
    "gene" %in% names(tbl),
    is.character(cancer_genes_symbol)
  )
  dat <- tbl |>
    dplyr::mutate(
      MeanCopyNumber = base::rowMeans(dplyr::select(tbl, c("minCopyNumber", "maxCopyNumber"))),
      MeanCopyNumber = dplyr::if_else(.data$MeanCopyNumber < 0, 0, .data$MeanCopyNumber)
    ) |>
    # keep only cancer genes
    dplyr::filter(
      .data$gene %in% cancer_genes_symbol
    )
  quants <- stats::quantile(
    dat[["MeanCopyNumber"]],
    probs = seq(0, 1, 0.05), na.rm = TRUE
  )
  # Keep only altered genes with CN values below loss threshold
  # (default 5th percentile) and above gain threshold (default 95th percentile)
  if (cn_bottom == 5 && cn_top == 95) {
    cn_bottom <- quants["5%"] |>
      base::unname() |>
      base::round(digits = 2)
    cn_top <- quants["95%"] |>
      base::unname() |>
      base::round(digits = 2)
  }
  # If their diff is 0 then increase/decrease threshold by 1
  if (abs(cn_top - cn_bottom) == 0) {
    cn_top <- cn_top + 1
    cn_bottom <- cn_bottom - 1
  }
  dat <- dat |>
    dplyr::filter(.data$MeanCopyNumber <= cn_bottom | .data$MeanCopyNumber >= cn_top) |>
    dplyr::pull("gene") |>
    unique() |>
    stats::na.omit()
  if (length(dat) == 0) {
    dat <- NULL
  }
  list(
    dat = dat,
    cn_bottom = cn_bottom,
    cn_top = cn_top
  )
}

#' Get Immune Gene Summary
#'
#' @param tbl_imarkers Tibble with immune response marker genes.
#' @param tbl_igram Tibble with immunogram genes.
#' @param igram_param Immunogram inclusion TRUE/FALSE parameter.
#'
#' @return Character vector of genes.
#' @export
immune_summary <- function(tbl_imarkers, tbl_igram = NULL, igram_param = TRUE) {
  assertthat::assert_that(
    inherits(tbl_imarkers, "data.frame"),
    "SYMBOL" %in% names(tbl_imarkers)
  )
  if (igram_param) {
    assertthat::assert_that(
      !is.null(tbl_igram),
      inherits(tbl_igram, "data.frame"),
      "SYMBOL" %in% names(tbl_igram)
    )
  }
  # Immune reponse markers
  res1 <- tbl_imarkers |>
    dplyr::pull("SYMBOL")

  res2 <- NULL
  if (igram_param) {
    res2 <- tbl_igram[["SYMBOL"]]
  }
  res <- unique(c(res1, res2)) |> stats::na.omit()
  res
}

#' Read PURPLE CNV Gene File
#'
#' Reads the `purple.cnv.gene.tsv` file, which summarises copy number
#' alterations of each gene in the HMF panel
#' (see https://github.com/hartwigmedical/hmftools/tree/master/purple#gene-copy-number-file).
#'
#' @param x Path to `purple.cnv.gene.tsv` file.
#'
#' @return The input file as a tibble.
#'
#' @export
ppl_cnv_som_gene_read <- function(x) {
  ct <- list(
    "chromosome" = "c", "start" = "d", "end" = "d", "gene" = "c",
    "minCopyNumber" = "d", "maxCopyNumber" = "d",
    "unused" = "c", "somaticRegions" = "d", "germlineHomDeletionRegions" = "d",
    "germlineHetToHomDeletionRegions" = "d",
    "transcriptId" = "c", "transcriptVersion" = "c", "chromosomeBand" = "c",
    "minRegions" = "d", "minRegionStart" = "d", "minRegionEnd" = "d",
    "minRegionStartSupport" = "c", "minRegionEndSupport" = "c",
    "minRegionMethod" = "c", "minMinorAlleleCopyNumber" = "d"
  )
  # PURPLE as of at least v3.9.2 (2023-09-05) has some different columns
  ct2 <- list(
    "chromosome" = "c", "start" = "d", "end" = "d", "gene" = "c",
    "minCopyNumber" = "d", "maxCopyNumber" = "d",
    "somaticRegions" = "d", "transcriptId" = "c", "isCanonical" = "c",
    "chromosomeBand" = "c",
    "minRegions" = "d", "minRegionStart" = "d", "minRegionEnd" = "d",
    "minRegionStartSupport" = "c", "minRegionEndSupport" = "c",
    "minRegionMethod" = "c", "minMinorAlleleCopyNumber" = "d",
    "depthWindowCount" = "d"
  )

  hdr <- readr::read_tsv(x, n_max = 0, show_col_types = FALSE)
  if (all(names(ct2) %in% colnames(hdr))) {
    ct <- ct2
  }
  d <- readr::read_tsv(x, col_types = ct)
  d[]
}
