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
