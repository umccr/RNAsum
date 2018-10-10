# Adapted from Matthew Bashton 2017
# Load filtered tsv output from pizzly add in gene co-ordinates + compute distance when on same chr

# CRAN
#install.packages("jsonlite")
#install.packages("dplyr")

# Bioconductor
#source("https://bioconductor.org/biocLite.R")
#biocLite()
#biocLite("ensembldb")
#biocLite("EnsDb.Hsapiens.v79")

# Set up
# Change to your own output location
setwd("~/Documents/UMCCR/data/fusions/pizzly-validation/grolar")

# CRAN libs
library(jsonlite)
library(dplyr)

# Assuming you have used GRCh38 gene models
library(EnsDb.Hsapiens.v79)
edb <- EnsDb.Hsapiens.v79
listColumns(edb)
supportedFilters(edb)


# Suffix which apears after sample id in output file name
suffix <- ".tsv"

# Lets get a list of all the files we want to process
tsv_files <- list.files(path = ".", pattern = paste0("*",suffix))

# Make a list of Ids
Ids <- gsub(suffix, "", tsv_files)

# This function will flatten out the JSON giving you a list of gene A and gene B
# sorted by splitcount then paircount, it will also get co-ordinates for each
# gene, and if on same chr caculate distance between the two
GetFusionz_and_namez <- function(sample, suffix) {

  # To test
  #sample <- Ids[1]
  #suffix = "_fusion_GRCh38.json"

  tsv_file <- paste0(sample, suffix)

  cat(paste("### Working on sample:", sample, "input file:", tsv_file, "###", "\n"))
  cat("Reading TSV..\n")
  TSV <- read.table(file =  paste(getwd(), tsv_file, sep = "/"), header = TRUE)

  # Add uniq-key
  output <- cbind(rownames(TSV), TSV)
  tmp <- colnames(output)[-1]
  colnames(output) <- c("ID", tmp)
  colnames(output)
  #head(output)

  # Now extract geneA geneB locations
  # get all geneAs
  cat("Getting co-ordinates for gene A\n")
  geneAs <- TSV$geneA.id
  geneAs_info <- genes(edb,
        columns = c("gene_id","seq_name","gene_seq_start","gene_seq_end","seq_strand"),
        filter = GeneIdFilter(geneAs),
        return.type = "DataFrame")

  # get all geneBs
  cat("Getting co-ordinates for gene B\n")
  geneBs <- TSV$geneB.id
  geneBs_info <- genes(edb,
                       columns = c("gene_id","seq_name","gene_seq_start","gene_seq_end","seq_strand"),
                       filter = GeneIdFilter(geneBs),
                       return.type = "DataFrame")

  # Rename cols
  colnames(geneAs_info) <- paste0("geneA.", colnames(geneAs_info))
  colnames(geneBs_info) <- paste0("geneB.", colnames(geneBs_info))
  #colnames(geneAs_info)
  #colnames(geneBs_info)

  # Merge in geneA_info
  cat("Merging DFs\n")
  tmp1 <- merge(output, geneAs_info, by.x = "geneA.id", by.y = "geneA.gene_id", all.x = TRUE)
  #excluding transcripts columns, so that it is not repeated when merging tmp1 and tmp2
  tmp2 <- merge(output[-8], geneBs_info, by.x = "geneB.id", by.y = "geneB.gene_id", all.x = TRUE)

  # Check we have the same length
  stopifnot(nrow(tmp1) == nrow(tmp2))

  # Merge these DFs
  output <- merge(tmp1, tmp2[,c(2,8:11)], by = "ID")[,c(1,6,7,2,3,9,10,11,12,5,4,13,14,15,16,8)]
  colnames(output)
  output <- as.data.frame(output)

  # Sort by splitcount then paircount
  cat("Sorting by number of events\n")
  idx <- order(output$splitcount, output$paircount, decreasing = TRUE)
  output <- output[idx,]

  # Write out output
  cat(paste0("Writing out table: ", sample, "_fusions_filt_sorted.txt", "\n"))
  write.table(output, file = paste0(sample, "_fusions_filt_sorted.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  cat("\n")
  return(TRUE)
}

# Use above funciton on all output files
lapply(Ids, function(x) GetFusionz_and_namez(x, suffix = ".tsv") )


