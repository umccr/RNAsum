#A script to select fusions from pizzly output which have a decent TPM support in the quantification file.

#read in the pizzly fusion calls
pizzly.fusions <- read.table(file = '~/Documents/UMCCR/data/fusions/pizzly-validation/MH17T001P013-oncofuse-test-flat-filtered.tsv', header = TRUE)

#read in the transcripts quantification file
quant <- read.table(file = '~/Documents/UMCCR/data/fusions/pizzly-validation/abundance.tsv', header = TRUE)

#sort and filter quantification file on tpm values
quant.sorted.filtered <- filter(arrange(quant, desc(quant$tpm)), tpm >= (mean(quant$tpm)))

#pizzly.fusions$transcripts.list is a factor. Coerce the argument to character first to be able to use sapply
pizzly.fusions$transcripts.list <- as.character(pizzly.fusions$transcripts.list)


quant.transcripts <- sapply(pizzly.fusions$transcripts.list, function(x){
  y <- strsplit(x, ";")
})
quant.transcripts.mod <- unname(quant.transcripts)
