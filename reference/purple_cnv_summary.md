# Get CNV Summary

Get CNV Summary

## Usage

``` r
purple_cnv_summary(tbl, cancer_genes_symbol, cn_bottom, cn_top)
```

## Arguments

- tbl:

  Gene copy-number tibble

- cancer_genes_symbol:

  Character vector of gene symbols to filter on.

- cn_bottom:

  CN threshold value to classify genes within lost regions.

- cn_top:

  CN threshold value to classify genes within gained regions.

## Value

List with genes, top and bottom copy number thresholds.
