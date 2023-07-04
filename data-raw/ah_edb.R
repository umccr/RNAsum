# This script is used to download a specific ensembl database version using Annotation Hub.
require(AnnotationHub)
require(here)

ah <- AnnotationHub()
# Downloading v105
edb_105 <- ah[["AH98047", force=TRUE]]
# Save RDS object to a file
saveRDS(edb_105, file = here::here("inst/extdata/edb_105.rds"))
