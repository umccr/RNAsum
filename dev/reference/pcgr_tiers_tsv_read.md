# Read PCGR TSV File

Reads the `snvs_indels.tiers.tsv` (or v2's `snv_indel_ann.tsv.gz`) file
output by PCGR.

## Usage

``` r
pcgr_tiers_tsv_read(x = NULL)
```

## Arguments

- x:

  Path to PCGR `snvs_indels.tiers.tsv`/`snv_indel_ann.tsv.gz` file.

## Value

A tibble with the contents of the input TSV file, or NULL if x is NULL.

## Examples

``` r
x <- system.file(
  "rawdata/test_data/small_variants",
  "TEST-snvs_indels.tiers.tsv",
  package = "RNAsum"
)
(ptsv <- pcgr_tiers_tsv_read(x))
#> # A tibble: 99 × 7
#>    GENOMIC_CHANGE VARIANT_CLASS SYMBOL PROTEIN_CHANGE CONSEQUENCE AF_TUMOR TIER 
#>    <chr>          <chr>         <chr>  <chr>          <chr>          <dbl> <chr>
#>  1 3:g.179218303… SNV           PIK3CA p.Glu545Lys    missense        0.22 2    
#>  2 7:g.7674894G>A SNV           TP53   p.Arg213Ter    stop gained     0.47 3    
#>  3 3:g.110715614… SNV           ING1   p.Phe57Leu     missense        0.12 3    
#>  4 3:g.110719444… SNV           ING1   p.Glu118Gln    missense        0.2  3    
#>  5 8:g.20250093C… SNV           LZTS1  p.Glu474Lys    missense        0.3  3    
#>  6 8:g.13086457C… SNV           DLC1   p.Trp1433Ter   stop gained     0.24 3    
#>  7 0:g.31521795G… SNV           ZEB1   p.Gln820His    missense        0.29 3    
#>  8 9:g.132927250… SNV           TSC1   p.Ser54Phe     missense        0.26 3    
#>  9 3:g.47121452G… SNV           SETD2  p.Pro1062Ser   missense        0.18 3    
#> 10 6:g.117365075… SNV           ROS1   p.Arg1035Gly   missense        0.24 3    
#> # ℹ 89 more rows
```
