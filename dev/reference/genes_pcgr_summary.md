# Get Mutated Genes from PCGR

Get Mutated Genes from PCGR

## Usage

``` r
genes_pcgr_summary(pcgr_tbl, tiers = c(1:4), splice_vars = TRUE)
```

## Arguments

- pcgr_tbl:

  Parsed PCGR TSV tibble.

- tiers:

  PCGR tiers to keep. Default: 1, 2, 3, 4.

- splice_vars:

  Include non-coding splice region variants reported in PCGR?

## Value

Character vector of genes.
