# Read DRAGEN Fusions

Reads the `fusion_candidates.final` file output by DRAGEN.

## Usage

``` r
dragen_fusions_read(x = NULL)
```

## Arguments

- x:

  Path to the `fusion_candidates.final` file.

## Value

A tibble with the contents of the input TSV file, or NULL if input is
NULL or TSV has no data rows.

## Examples

``` r
x <- system.file("rawdata/test_data/dragen",
  "test_sample_WTS.fusion_candidates.final",
  package = "RNAsum"
)
(a <- dragen_fusions_read(x))
#> # A tibble: 2 × 15
#>   gene1 gene2   Score LeftBreakpoint RightBreakpoint Gene1Location Gene2Location
#>   <chr> <chr>   <dbl> <chr>          <chr>           <chr>         <chr>        
#> 1 ATAD2 FBXO32  0.966 chr8:12338052… chr8:123534814… IntactExon    IntactExon   
#> 2 NRIP1 AF1275… 0.658 chr21:1506474… chr21:14857708… IntactExon    IntactExon;I…
#> # ℹ 8 more variables: Gene1Sense <chr>, Gene2Sense <chr>, Gene1Id <chr>,
#> #   Gene2Id <chr>, NumSplitReads <dbl>, NumSoftClippedReads <dbl>,
#> #   NumPairedReads <dbl>, ReadNames <chr>
```
