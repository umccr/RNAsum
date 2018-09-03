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

For comparison purposes, I sorted the paircount column of this file using:

```sort -k5 -n -r MH17T001P013-oncofuse-test-flat-filtered.tsv > MH17T001P013-oncofuse-test-flat-filtered-sorted.tsv```

### Oncofuse

Oncofuse is a post-processing tool that takes fusions called by other programs (such as  STAR) as an input and predict oncogenic potential i.e. the probability of being 'driver' events to fusion sequences.



 

