library("dplyr")

#Read in the input files
oncofuse <- read.table(file = '~/Documents/UMCCR/data/fusions/comparison/MH17T001P013-oncofuse-edited', header = TRUE)
pizzly <- read.table(file = '~/Documents/UMCCR/data/fusions/comparison/MH17T001P013-oncofuse-test-flat-filtered.tsv', header = TRUE)

#sort oncofuse output on driver probability of fusions AND filter the sorted oncofuse file to extract
#rows with driver probability above a certain threshhold (currently using 0.6)
oncofuse.sorted <- filter(arrange(oncofuse, desc(oncofuse$DRIVER_PROB)), DRIVER_PROB >= 0.6)


