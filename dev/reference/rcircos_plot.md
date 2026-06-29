# RCircos Plot

RCircos Plot

## Usage

``` r
rcircos_plot(df_circos, df_circos_pairs, cyto.info)
```

## Arguments

- df_circos:

  Dataframe with 4 columns: Chromosome, chromStart, chromEnd, and Gene.

- df_circos_pairs:

  Dataframe with 3 columns for each fusion mate: Chromosome, chromStart
  and chromEnd for geneA and geneB.

- cyto.info:

  A dataframe with hg38 cytoband/ideogram information from RCircos.

## Value

An RCircos plot object.
