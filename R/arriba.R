#' Read Arriba TSV File
#'
#' Reads the TSV file output by Arriba.
#' @param x Path to Arriba TSV file.
#' @return A tibble with the contents of the input TSV file, or NULL if x is NULL.
#'
#' @examples
#' x <- system.file("rawdata/test_data/dragen/arriba/fusions.tsv", package = "RNAsum")
#' (a <- arriba_tsv_read(x))
#' @testexamples
#' expect_equal(colnames(a)[ncol(a)], "read_identifiers")
#' expect_null(arriba_tsv_read(x = NULL))
#' @export
arriba_tsv_read <- function(x = NULL) {
  if (is.null(x)) {
    return(NULL)
  }
  x |>
    readr::read_tsv(
      col_types = readr::cols(
        .default = "c", discordant_mates = "d",
        split_reads1 = "d", split_reads2 = "d",
        coverage1 = "d", coverage2 = "d"
      )
    ) |>
    dplyr::rename(gene1 = .data$`#gene1`)
}

#' Read Arriba PDF File
#'
#' Reads the PDF output by Arriba containing one fusion plot per page and
#' converts each page to a separate PNG file.
#'
#' @param pdf Output PDF file from Arriba.
#' @param fusions Tibble with fusions from Arriba.
#' @param outdir Directory to write output PNGs to.
#' @return Tibble with paths to the created PNG images and their clean names.
#'         If the PDF (or any of the required params) is NULL, returns NULL.
#' @examples
#' pdf <- system.file("rawdata/test_data/dragen/arriba/fusions.pdf", package = "RNAsum")
#' tsv <- system.file("rawdata/test_data/dragen/arriba/fusions.tsv", package = "RNAsum")
#' fusions <- arriba_tsv_read(tsv)
#' (pngs <- arriba_pdf_read(pdf, fusions, tempdir()))
#' @testexamples
#' expect_equal(nrow(pngs), 4)
#' expect_null(arriba_pdf_read(pdf = NULL))
#' @export
arriba_pdf_read <- function(pdf = NULL, fusions = NULL, outdir = NULL) {
  if (is.null(pdf) || is.null(fusions) || is.null(outdir)) {
    return(NULL)
  }
  mkdir(outdir)
  # construct png filenames
  clean_fusions <- fusions |>
    dplyr::select(
      g1 = "gene1", g2 = "gene2",
      bp1 = "breakpoint1", bp2 = "breakpoint2"
    ) |>
    dplyr::mutate(
      nm = glue::glue("{.data$g1}__{.data$g2}_{.data$bp1}_{.data$bp2}"),
      nm = make.names(.data$nm),
      png = file.path(outdir, glue::glue("{.data$nm}.png")),
    ) |>
    dplyr::select("nm", "png")
  # Export pdf images to png
  for (i in seq_len(nrow(clean_fusions))) {
    png <- pdftools::pdf_render_page(
      pdf = pdf, page = i, dpi = 300, numeric = TRUE, opw = "", upw = ""
    )
    png::writePNG(png, clean_fusions[["png"]][i])
  }
  # Return paths to PNGs in a tibble
  return(clean_fusions)
}

#' Write Arriba Fusion Summary
#'
#' Writes Arriba fusion summary information (`gene1__gene2_bp1_bp2`).
#'
#' @param x Tibble with Arriba fusions or NULL.
#' @param file File to write the summary info to.
#'
#' @return The input invisibly.
#' @export
arriba_summary_write <- function(x, file) {
  if (is.null(x)) {
    tibble::tibble(empty = character()) |>
      readr::write_tsv(file = file, col_names = FALSE)
    return()
  }
  x |>
    dplyr::select("nm") |>
    readr::write_tsv(file = file, col_names = FALSE)
}
