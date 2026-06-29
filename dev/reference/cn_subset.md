# Copy Number Subset Table

Copy Number Subset Table

## Usage

``` r
cn_subset(gene.mut = NULL, cn_data, expr_data.z, expr_data.perc)
```

## Arguments

- gene.mut:

  Tibble with gene and alterations, if PCGR was run.

- cn_data:

  Tibble with gene and MeanCopyNumber

- expr_data.z:

  Tibble with SYMBOL and z-score expr values.

- expr_data.perc:

  Tibble with SYMBOL and percentile expr values.

## Value

A tibble with Gene, CN, Diff_Perc, Diff_Z_score, and Alterations if PCGR
was run.
