##### Calculate TPM from RPKM (from http://luisvalesilva.com/datasimple/rna-seq_units.html )
tpm_from_rpkm <- function(x) {
  rpkm.sum <- colSums(x)
  return(t(t(x) / (1e-06 * rpkm.sum)))
}
