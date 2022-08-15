#' Read Clinical Information Spreadsheet
#'
#' Reads clinical information spreadsheet.
#' @param x Path to spreadsheet containing clinical information.
#'
#' @return A data.frame representation of the input spreadsheet, or NULL if input
#'         is NULL.
#' @examples
#' x <- system.file("rawdata/test_data/test_clinical_data.xlsx", package = "RNAsum")
#' (d <- clinical_info_read(x))
#' @testexamples
#' expect_equal(colnames(d)[1], "Subject.ID")
#' expect_null(clinical_info_read(NULL))
#' @export
clinical_info_read <- function(x = NULL) {
  if (is.null(x)) {
    return(NULL)
  }

  openxlsx::read.xlsx(
    xlsxFile = x, sheet = 1, colNames = TRUE, rowNames = FALSE,
    detectDates = TRUE, skipEmptyRows = TRUE, skipEmptyCols = TRUE, check.names = TRUE
  )
}
