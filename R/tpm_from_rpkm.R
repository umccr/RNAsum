#' Calculate TPM from RPKM
#'
#' Calculate TPM from RPKM (from http://luisvalesilva.com/datasimple/rna-seq_units.html )
#'
#' @param x RPKM values
#' @return TPM values
#' @export
tpm_from_rpkm <- function(x) {
  rpkm.sum <- base::colSums(x)
  return(base::t(base::t(x) / (1e-06 * rpkm.sum)))
}
