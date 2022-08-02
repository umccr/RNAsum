#' Read Arriba TSV File
#'
#' Reads the TSV file output by Arriba.
#' @param x Path to Arriba TSV file.
#'
#' @return A tibble with the contents of the input TSV file.
#' @examples
#' x <- system.file("rawdata/test_data/dragen/arriba/fusions.tsv", package = "RNAsum")
#' (a <- arriba_read_tsv(x))
#' @testexamples
#' expect_equal(colnames(a)[ncol(a)], "read_identifiers")
#' @export
arriba_read_tsv <- function(x) {
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
#' @return Single-column tibble with paths to the created PNG images.
#' @examples
#' pdf <- system.file("rawdata/test_data/dragen/arriba/fusions.pdf", package = "RNAsum")
#' tsv <- system.file("rawdata/test_data/dragen/arriba/fusions.tsv", package = "RNAsum")
#' fusions <- arriba_read_tsv(tsv)
#' (pngs <- arriba_read_pdf(pdf, fusions, tempdir()))
#' @testexamples
#' expect_equal(nrow(pngs), 4)
#' @export
arriba_read_pdf <- function(pdf, fusions, outdir) {
  mkdir(outdir)
  # construct png filenames
  clean_fusions <- fusions |>
    dplyr::select(
      g1 = .data$gene1, g2 = .data$gene2,
      bp1 = .data$breakpoint1, bp2 = .data$breakpoint2
    ) |>
    dplyr::mutate(
      nm = glue::glue("{.data$g1}__{.data$g2}_{.data$bp1}_{.data$bp2}.png"),
      nm = file.path(outdir, make.names(.data$nm))
    ) |>
    dplyr::select(.data$nm)
  # Export pdf images to png
  for (i in seq_len(nrow(clean_fusions))) {
    png <- pdftools::pdf_render_page(pdf,
      page = i, dpi = 300, numeric = TRUE,
      opw = "", upw = ""
    )
    png::writePNG(png, clean_fusions[["nm"]][i])
  }
  # Return paths to PNGs in a tibble
  return(clean_fusions)
}
