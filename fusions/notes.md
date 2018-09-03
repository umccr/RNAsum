# Introduction

Gene fusions are hybrid genes formed when two previously independent genes become juxtaposed. The fusion can result from structural rearrangements like translocations and deletions, transcription read-through of neighboring genes, or the trans- and cis-splicing of pre-mRNAs [ref](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4889949/). 

![Alt text](./images/fusion.jpeg)

Many gene fusions are associated with oncogenic properties, and often act as driver mutations in a wide array of cancer types.
We are interested in including fusion information into RNAseq report and validate it against genomic evidence. 

## Explore Pizzly VS Oncofuse

bcbio-RNAseq pipeline support a couple of tools for fusion calling, [pizzly](https://github.com/pmelsted/pizzly) and [oncofuse](http://www.unav.es/genetica/oncofuse.html). Ideally, in the report, we would like to use output from one of the tools that we are confident in. 

### Pizzly

Pizzly produces a `tsv` file of the genes with the breakpoints indicated relative to the transcript, e.g. 

```geneA.name geneA.id geneB.name geneB.id paircount splitcount transcripts.list```

For comparison purposes, I sorted pizzly output on the paircount value:

```sort -k5 -n -r MH17T001P013-oncofuse-test-flat-filtered.tsv > MH17T001P013-oncofuse-test-flat-filtered-sorted.tsv```

### Oncofuse

Oncofuse is a post-processing tool that takes fusions called by other programs (such as  STAR) as an input and predict oncogenic potential i.e. the probability of being 'driver' events to fusion sequences.

Oncofuse produces a 36 columns output. For comparison purpose, I have selected the following columns (and copied those to a new life). 

```5_FPG_GENE_NAME 3_FPG_GENE_NAME SPANNING_READS ENCOMPASSING_READS GENOMIC P_VAL_CORR DRIVER_PROB```

* `FPG` stands for fusion pair gene.
* Encompassing mate pairs refer to those in which each read aligns to an independent transcript, thereby encompassing the fusion junction. 
* Spanning mate pairs refer to those in which one sequence read aligns to a gene and its mate spans the fusion junction.

Sorted the oncofuse output on the encompassing reads: 

```sort -k4 -n -r MH17T001P013-oncofuse-edited >MH17T001P013-oncofuse-edited-sorted.tsv```

### Result

Importing both pizzly and oncofuse tsvs as dataframe in R to find common fusion gene pairs betwen both tools. The intersection between both results in quite low (35 pairs). Another point is the final filtered calls produced by pizzly are ~100 and oncofuse is ~500, which is an indication of huge amount of false positive predictions.

Thus, we need to find ways to shrink down the list of candidates from fusion detection tools, thus focusing on a reduced set of highly reliable fusions with a potential driver impact into the disease.

## Exploring fusion prioritization tool 

[FuGePrior](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5260008/) 












 

