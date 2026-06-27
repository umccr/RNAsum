# Read CNV Gene File

Reads the copy-number gene file (e.g. `purple.cnv.gene.tsv`), which
summarises copy number alterations of each gene in the HMF panel (see
https://github.com/hartwigmedical/hmftools/tree/master/purple#gene-copy-number-file).

## Usage

``` r
ppl_cnv_som_gene_read(x)
```

## Arguments

- x:

  Path to `purple.cnv.gene.tsv` file.

## Value

The input file as a tibble.
