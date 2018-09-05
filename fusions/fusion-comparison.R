library("dplyr")

#read in the input files
oncofuse <- read.table(file = '~/Documents/UMCCR/data/fusions/comparison/MH17T001P013-oncofuse-edited', header = TRUE)
pizzly <- read.table(file = '~/Documents/UMCCR/data/fusions/comparison/MH17T001P013-oncofuse-test-flat-filtered.tsv', header = TRUE)

#sort oncofuse output on driver probability of fusions AND filter the sorted oncofuse file to extract
#rows with driver probability above a certain threshhold (currently using 0.6)
oncofuse.sorted.filtered <- filter(arrange(oncofuse, desc(oncofuse$DRIVER_PROB)), DRIVER_PROB >= 0.6)

#sort pizzly output on paircount for the fusion genes AND filter the sorted pizzly file to extract
#rows with paircount values above a certain threshhold (currently using 5)
pizzly.sorted.filtered <- filter(arrange(pizzly, desc(pizzly$paircount)), paircount >= 5)

#applying semi_join on oncofuse.sorted.filter and pizzly.sorted.filtered
joint.fusion.calls <- semi_join(oncofuse.sorted.filtered, pizzly.sorted.filtered,
                                by = c("X5_FPG_GENE_NAME" = "geneA.name", "X3_FPG_GENE_NAME" = "geneB.name"))

#in case pizzly pizzly gene naming is not according to 5 -> 3 convention. Couldn't find a clear answer to this online
joint.fusion.calls.2 <- semi_join(oncofuse.sorted.filtered, pizzly.sorted.filtered,
                                by = c("X5_FPG_GENE_NAME" = "geneB.name", "X3_FPG_GENE_NAME" = "geneA.name"))

#selectiong all columns from the two dataframes
joint.fusion.calls.inner <- inner_join(oncofuse.sorted.filtered, pizzly.sorted.filtered,
                                by = c("X5_FPG_GENE_NAME" = "geneA.name", "X3_FPG_GENE_NAME" = "geneB.name"))


