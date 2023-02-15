#' Read Sample Data
#'
#' Reads sample data, including Arriba fusions, Arriba plots, Salmon counts,
#' DRAGEN fusions, and clinical info.
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
#'   dragen_rnaseq = system.file("rawdata/test_data/dragen", package = "RNAsum"),
#'   arriba_pdf = system.file("rawdata/test_data/dragen/arriba/fusions.pdf", package = "RNAsum"),
#'   arriba_tsv = system.file("rawdata/test_data/dragen/arriba/fusions.tsv", package = "RNAsum"),
#'   dragen_fusions = system.file(
#'     "rawdata/test_data/dragen/test_sample_WTS.fusion_candidates.final",
#'     package = "RNAsum"
#'   ),
#'   umccrise = system.file(
#'     "rawdata/test_data/umccrised/test_sample_WGS",
#'     package = "RNAsum"
#'   ),
#'   manta_tsv = system.file(
#'     "rawdata/test_data/umccrised/test_sample_WGS/structural/manta.tsv",
#'     package = "RNAsum"
#'   )
#' )
#' res <- read_sample_data(p, tempdir())
#' @testexamples
#' expect_equal(length(res), 6)
#' expect_null(res$salmon)
#' @export
read_sample_data <- function(p, results_dir, tx2gene = NULL) {
  arriba_tsv <- arriba_tsv_read(p[["arriba_tsv"]])
  arriba_pdf <- p[["arriba_pdf"]] |>
    arriba_pdf_read(fusions = arriba_tsv, outdir = file.path(results_dir, "arriba"))
  salmon <- salmon_counts(p[["salmon"]], tx2gene = tx2gene)
  dragen_fusions <- dragen_fusions_read(p[["dragen_fusions"]])
  clininfo <- clinical_info_read(p[["clinical_info"]])
  wgs <- read_wgs_data(p)
  list(
    arriba_tsv = arriba_tsv,
    arriba_pdf = arriba_pdf,
    salmon = salmon,
    dragen_fusions = dragen_fusions,
    clinical_info = clininfo,
    wgs = wgs
  )
}

#' Read WGS Data
#'
#' Reads WGS data, including PCGR `tiers.tsv`, PURPLE `cnv.gene.tsv`, and Manta
#' `manta.tsv`. If the file path has been specified in the RNAsum params and is
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
#'   manta_tsv = system.file(
#'     "rawdata/test_data/umccrised/test_sample_WGS/structural/TEST-prioritize-manta.tsv",
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
    nm = "purple_gene_tsv", func = gpgr::purple_cnv_som_gene_read
  )

  manta_tsv <- .read(
    p = p,
    subdir = "structural", pat = "manta\\.tsv$",
    nm = "manta_tsv", func = gpgr::process_sv
  )

  list(
    pcgr_tiers_tsv = pcgr_tiers_tsv,
    purple_gene_tsv = purple_gene_tsv,
    manta_tsv = manta_tsv
  )
}
