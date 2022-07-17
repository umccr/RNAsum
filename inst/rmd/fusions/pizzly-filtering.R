# Script to select fusions from pizzly output which have a decent TPM support in the quantification file.
# This analysis is done to support adding fusion calling results to RNAseq report

# read in the pizzly fusion calls
pizzly.fusions <- read.table(file = "~/Documents/UMCCR/data/fusions/pizzly-validation/MH17T001P013-oncofuse-test-flat-filtered.tsv", header = TRUE)

# read in the transcripts quantification file
quant <- read.table(file = "~/Documents/UMCCR/data/fusions/pizzly-validation/abundance.tsv", header = TRUE)

# sort and filter quantification file on tpm values. Currently filtering on quantiles. Selected 0.997 because that reduces the final fusion
# calls to the value we are interested in (~15)

quant.sorted.filtered <- filter(arrange(quant, desc(quant$tpm)), tpm >= (quantile(quant$tpm, 0.999)))

# pizzly.fusions$transcripts.list is a factor. Coerce the argument to character first to be able to use sapply
# pizzly.fusions$transcripts.list <- as.character(pizzly.fusions$transcripts.list)

# initialize an empty dataframe
result <- data.frame()
result.list <- vector("list", length = nrow(pizzly.fusions))

# let's try using for loop for iterating over pizzly.fusions dataframe and get transcriptID and fusion gene pair information.
# can also filter quant.sorted.filtered$target_id to have only fusion gene target ids (that is two tracscripts instead of one-
# this will increase speed

for (row in 1:nrow(pizzly.fusions)) {
  y <- strsplit(as.character(pizzly.fusions[row, 7]), "\\;")
  y <- unname(y)
  for (i in 1:length(y[[1]])) {
    if (y[[1]][i] %in% quant.sorted.filtered$target_id) {
      # val <- c(y[[1]][i], as.character(pizzly.fusions[row,1]), as.character(pizzly.fusions[row, 3]))
      # creating a new dataframe for the filtered pizzly results
      result <- rbind(result, data.frame(pizzly.fusions[row, ]))
      # result.list[[row]] <- val
    }
  }
}

# write final output to a file
write.csv(result, file = "pizzly-filtered-MH17T001P013.csv", row.names = FALSE)

# remove nulls from result list
# result.list <- result.list[lapply(result.list, length)>0]

# remove duplicated values from result (as multiple transcripts might support fusion between same gene)
deduped.result <- unique(result)

# apply function on every row in a df, but selected columns. Split the transcript list on ';' and
# check for occurence of every split element in the quantification file. If it exists, extract the correponding fusion gene pair
# apply screws out the result list, only giving me the final value evaluated in the condition
apply(pizzly.fusions[, c("geneA.name", "geneB.name", "transcripts.list")], 1, function(x) {
  y <- strsplit(x["transcripts.list"], ";")
  y <- unname(y)
  for (i in 1:length(y[[1]])) {
    if (y[[1]][i] %in% quant.sorted.filtered$target_id) {
      val <- c(y[[1]][i], x["geneA.name"], x["geneB.name"])
      result.list[[i]] <- val
      print(result.list[[i]])
      inter <- inter + 1
      # result <<- rbind(result, data.frame(transcriptID = y[[1]][i], geneA = x['geneA.name'], geneB = x['geneB.name'], stringsAsFactors=FALSE))
    }
  }
})
