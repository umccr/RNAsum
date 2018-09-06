#A script to select fusions from pizzly output which have a decent TPM support in the quantification file.

#read in the pizzly fusion calls
pizzly.fusions <- read.table(file = '~/Documents/UMCCR/data/fusions/pizzly-validation/MH17T001P013-oncofuse-test-flat-filtered.tsv', header = TRUE)

#read in the transcripts quantification file
quant <- read.table(file = '~/Documents/UMCCR/data/fusions/pizzly-validation/abundance.tsv', header = TRUE)

#sort and filter quantification file on tpm values
quant.sorted.filtered <- filter(arrange(quant, desc(quant$tpm)), tpm >= (mean(quant$tpm)))

#pizzly.fusions$transcripts.list is a factor. Coerce the argument to character first to be able to use sapply
#pizzly.fusions$transcripts.list <- as.character(pizzly.fusions$transcripts.list)

data <- data.frame()
df <- apply(pizzly.fusions[,c('geneA.name','geneB.name','transcripts.list')],1, function(x){
  y <- strsplit(x['transcripts.list'], ";")
  y <- unname(y)
  for (i in 1:length(y[[1]])){
    #print(length(y[[1]]))
    #print(y[[1]][i])
    if (y[[1]][i] %in% quant.sorted.filtered$target_id){
      #print("yes")
      #print(x['geneA.name'])
      df <- data.frame(geneA = x['geneA.name'], geneB = x['geneB.name'], stringsAsFactors=FALSE)
      data <- rbind(data, df)
    }
  }
  print(data)
})


