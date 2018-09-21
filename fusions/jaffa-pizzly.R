#required packages
library(dplyr)
library(tidyr)

#read in input files
pizzly <- read.table(file = '~/Documents/UMCCR/data/fusions/comparison/MH17T001P013-oncofuse-test-flat-filtered.tsv', header = TRUE)
jaffa <- read.table(file = '~/Documents/UMCCR/data/fusions/comparison/jaffa/jaffa_results.csv', sep = ",", header = TRUE)

#extract fusion gene pairs from both files
pizzly.fusiongenes <- pizzly %>% select(geneA.name, geneB.name)
jaffa.fusiongenes <- jaffa %>% dplyr::select(fusion.genes)
jaffa.fusiongenes <- tidyr::separate(data = jaffa.fusiongenes, col = fusion.genes, into = c("geneA", "geneB"), sep = ":")
