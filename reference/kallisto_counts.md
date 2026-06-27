# Get Kallisto Counts

Get kallisto counts via tximport.

## Usage

``` r
kallisto_counts(x, tx2gene = NULL)
```

## Arguments

- x:

  Path to `abundance.tsv` file with abundances. See
  [tximport::tximport](https://rdrr.io/pkg/tximport/man/tximport.html).

- tx2gene:

  data.frame with tx_name and gene_id columns. See
  [tximport::tximport](https://rdrr.io/pkg/tximport/man/tximport.html).

## Value

Tibble with the counts per gene transcript, or NULL if any of the input
params are NULL.

## Examples

``` r
x <- system.file("rawdata/test_data/kallisto/abundance.tsv", package = "RNAsum")
tx2gene <- NULL
(kc <- kallisto_counts(x, tx2gene)) # NULL since no tx2gene specified
#> NULL
```
