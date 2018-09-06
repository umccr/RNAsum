#A script to select fusions from pizzly output which have a decent TPM support in the quantification file.

#read in the pizzly fusion calls
pizzly.fusions <- read.table(file = '~/Documents/UMCCR/data/fusions/pizzly-validation/MH17T001P013-oncofuse-test-flat-filtered.tsv', header = TRUE)

#read in the transcripts quantification file
quant <- read.table(file = '~/Documents/UMCCR/data/fusions/pizzly-validation/abundance.tsv', header = TRUE)

#sort and filter quantification file on tpm values. Currently filtering on transcripts that have higher value than
#the mean transcript value.
quant.sorted.filtered <- filter(arrange(quant, desc(quant$tpm)), tpm >= (mean(quant$tpm)))

#pizzly.fusions$transcripts.list is a factor. Coerce the argument to character first to be able to use sapply
#pizzly.fusions$transcripts.list <- as.character(pizzly.fusions$transcripts.list)

#initialize an empty dataframe
result <- data.frame()

#apply function on every row in a df, but selected columns. Split the transcript list on ';' and check for occurence of every split element
#in the quantification file. If it exists, extract the correponding fusion gene pair
apply(pizzly.fusions[,c('geneA.name','geneB.name','transcripts.list')],1, function(x){
  y <- strsplit(x['transcripts.list'], ";")
  y <- unname(y)
  for (i in 1:length(y[[1]])){
    if (y[[1]][i] %in% quant.sorted.filtered$target_id){
      result <<- rbind(result, data.frame(geneA = x['geneA.name'], geneB = x['geneB.name'], stringsAsFactors=FALSE))
    }
  }
})

#remove duplicated values from result (as multiple transcripts might support fusion between same gene)
deduped.result <- unique(result)



