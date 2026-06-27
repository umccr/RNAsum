# Read DRAGEN Mapping Metrics

Reads the `mapping_metrics.csv` file output by DRAGEN.

## Usage

``` r
dragen_mapping_metrics_read(x = NULL)
```

## Arguments

- x:

  Path to the `mapping_metrics.csv` file.

## Value

A tibble with the contents of the input TSV file, or NULL if input is
NULL or TSV has no data rows.

## Examples

``` r
x <- system.file("rawdata/test_data/dragen",
  "test.mapping_metrics.csv",
  package = "RNAsum"
)
(d1 <- dragen_mapping_metrics_read(x))
#> # A tibble: 113 x 6
#>    category rg    var                                   var_abbrev   count   pct
#>    <chr>    <chr> <chr>                                 <chr>        <dbl> <dbl>
#>  1 summary  TOTAL Total input reads                     reads_tot~  2.40e8 100  
#>  2 summary  TOTAL Number of duplicate marked reads      reads_num~  0        0  
#>  3 summary  TOTAL Number of duplicate marked and mate ~ reads_num~ NA       NA  
#>  4 summary  TOTAL Number of unique reads (excl. duplic~ reads_num~  2.40e8 100  
#>  5 summary  TOTAL Reads with mate sequenced             reads_w_m~  2.40e8 100  
#>  6 summary  TOTAL Reads without mate sequenced          reads_wo_~  0        0  
#>  7 summary  TOTAL QC-failed reads                       reads_qcf~  0        0  
#>  8 summary  TOTAL Mapped reads                          reads_map~  1.69e8  70.3
#>  9 summary  TOTAL Mapped reads adjusted for filtered m~ reads_map~  1.70e8  70.7
#> 10 summary  TOTAL Mapped reads R1                       reads_map~  8.43e7  70.3
#> # i 103 more rows
```
