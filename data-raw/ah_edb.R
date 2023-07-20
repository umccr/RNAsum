# This script is used to download a specific ensembl database version using Annotation Hub and extract tarnscript and gene IDS for a specific version
require(AnnotationHub)
require(ensembldb)
require(here)

ah <- AnnotationHub()
# Use query(ah, "EnsDb") to find specific id for a ensembl database record
# Downloading v105
edb_105 <- ah[["AH98047", force=TRUE]]

# Find keys
keys_ah <- AnnotationDbi::keys(edb_105, keytype = "GENEID")
# Extract gene and transcript ids
tx_gene_id <- edb_105 |>
  AnnotationDbi::select(
    keys = keys_ah,
    columns = c(
      "TXID", "GENEID", "GENENAME", "GENEBIOTYPE", "SEQNAME", "GENESEQSTART", "GENESEQEND", "SEQLENGTH"
    ), keytype = "GENEID")

# Extract gene lengths from ensbldb and assign values to dataframe
gene.length <- ensembldb::lengthOf(edb_105, of="gene")
tx_gene_id$GENELENGTH <- gene.length[tx_gene_id$GENEID]

# Save gene and transcript ids to a RDS file
saveRDS(tx_gene_id, file = here::here("inst/extdata/tx_gene_id_105.rds"))
