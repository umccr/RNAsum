#' Read Sample Data
#'
#' Reads sample data, including Arriba fusions, Arriba plot, Salmon counts, and
#' DRAGEN fusions.
#' @param p RNAsum params list.
#' @param results_dir Directory to output extracted Arriba PNGs to (created
#' automatically if it does not already exist).
#'
#' @return A list of the input sample data.
#'
#'
#' @examples
#' @export
read_sample_data <- function(p, results_dir) {
  assertthat::assert_that(
    all(c("arriba_pdf", "arriba_tsv", "salmon", "dragen_fusions") %in% names(p))
  )

  # if the param hasn't been provided, check inside the dragen
  # directory for the file pattern and return it, else return NULL.
  .find_file <- function(x, pat, d) {
    if (is.null(x)) {
      x <- list.files(d, pattern = pat, full.names = TRUE, recursive = TRUE)
      if (length(x) == 1) {
        return(x)
      }
    } else {
      # Not NULL
      if (file.exists(x)) {
        return(x)
      }
    }
    return(NULL) # catch-all
  } # .find_file end

  .read_file <- function(x, pat, d, func, ...) {
    f <- .find_file(x, pat, d)
    if (!is.null(f)) {
      return(func(f, ...))
    }
    return(NULL)
  }

  d <- p[["dragen_rnaseq"]]
  # returns a tibble with fusions
  arriba_tsv <- .read_file(p[["arriba_tsv"]], "^fusions\\.tsv$", d, RNAsum::arriba_read_tsv)
  # returns a tibble with paths to individual png plots
  arriba_pdf <- .read_file(p[["arriba_pdf"]], "^fusions\\.pdf$", d, RNAsum::arriba_read_pdf, arriba_tsv, file.path(results_dir, "arriba"))
  # returns a tibble with counts
  salmon <- .read_file(p[["salmon"]], "\\.quant\\.sf$", d, RNAsum::salmon_counts, tx2gene = tx2ensembl)
  # returns a tibble with fusions
  dragen_fusions <- .read_file(p[["dragen_fusions"]], "\\.fusion_candidates\\.final$", d, RNAsum::dragen_fusions_read)
  list(
    arriba_tsv = arriba_tsv,
    arriba_pdf = arriba_pdf,
    salmon = salmon,
    dragen_fusions = dragen_fusions
  )
}
