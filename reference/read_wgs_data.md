# Read WGS Data

Reads WGS data, including PCGR `tiers.tsv`, copy-number `cnv.gene.tsv`,
and `sv.prioritised.tsv`. If the file path has been specified in the
RNAsum params and is valid, it is returned. As a fallback, if the
umccrise directory param has been specified, then there is an attempt to
detect the file pattern in there.

## Usage

``` r
read_wgs_data(p)
```

## Arguments

- p:

  RNAsum params list.

## Value

A list of the input sample data.

## Examples

``` r
p <- list(
  umccrise = system.file("rawdata/test_data", package = "RNAsum"),
  pcgr_tiers_tsv = system.file(
    "rawdata/test_data/small_variants", "TEST-snvs_indels.tiers.tsv",
    package = "RNAsum"
  ),
  cn_gene_tsv = system.file(
    "rawdata/test_data/copy_number", "TEST.cnv.gene.tsv",
    package = "RNAsum"
  ),
  sv_tsv = system.file(
    "rawdata/test_data/structural", "TEST-sv.tsv",
    package = "RNAsum"
  )
)
(res <- read_wgs_data(p))
#> Warning: The following named parsers don't match the column names: chromosome, start, end, unused, somaticRegions, germlineHomDeletionRegions, germlineHetToHomDeletionRegions, transcriptId, transcriptVersion, chromosomeBand, minRegions, minRegionStart, minRegionEnd, minRegionStartSupport, minRegionEndSupport, minRegionMethod, minMinorAlleleCopyNumber
#> Rows: 99 Columns: 1
#> -- Column specification --------------------------------------------------------
#> Delimiter: "\t"
#> chr (1): annotation
#> 
#> i Use `spec()` to retrieve the full column specification for this data.
#> i Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> $pcgr_tiers_tsv
#> # A tibble: 99 x 7
#>    GENOMIC_CHANGE VARIANT_CLASS SYMBOL PROTEIN_CHANGE CONSEQUENCE AF_TUMOR TIER 
#>    <chr>          <chr>         <chr>  <chr>          <chr>          <dbl> <chr>
#>  1 3:g.179218303~ SNV           PIK3CA p.Glu545Lys    missense        0.22 2    
#>  2 7:g.7674894G>A SNV           TP53   p.Arg213Ter    stop gained     0.47 3    
#>  3 3:g.110715614~ SNV           ING1   p.Phe57Leu     missense        0.12 3    
#>  4 3:g.110719444~ SNV           ING1   p.Glu118Gln    missense        0.2  3    
#>  5 8:g.20250093C~ SNV           LZTS1  p.Glu474Lys    missense        0.3  3    
#>  6 8:g.13086457C~ SNV           DLC1   p.Trp1433Ter   stop gained     0.24 3    
#>  7 0:g.31521795G~ SNV           ZEB1   p.Gln820His    missense        0.29 3    
#>  8 9:g.132927250~ SNV           TSC1   p.Ser54Phe     missense        0.26 3    
#>  9 3:g.47121452G~ SNV           SETD2  p.Pro1062Ser   missense        0.18 3    
#> 10 6:g.117365075~ SNV           ROS1   p.Arg1035Gly   missense        0.24 3    
#> # i 89 more rows
#> 
#> $cn_gene_tsv
#> # A tibble: 25,417 x 3
#>    gene        minCopyNumber maxCopyNumber
#>    <chr>               <dbl>         <dbl>
#>  1 DDX11L1              1.98          1.98
#>  2 WASH7P               1.98          1.98
#>  3 MIR6859-1            1.98          1.98
#>  4 MIR1302-2HG          1.98          1.98
#>  5 MIR1302-2            1.98          1.98
#>  6 FAM138A              1.98          1.98
#>  7 OR4F5                1.98          1.98
#>  8 FO538757.1           1.98          1.98
#>  9 MIR6859-2            1.98          1.98
#> 10 OR4F29               1.98          1.98
#> # i 25,407 more rows
#> 
#> $sv_tsv
#> $sv_tsv$melted
#> # A tibble: 57 x 1
#>    Genes        
#>    <chr>        
#>  1 ADAMTSL4     
#>  2 ANXA9        
#>  3 ARNT         
#>  4 ARNT-RPS27AP6
#>  5 BRDT         
#>  6 BTBD6P1      
#>  7 CERS2        
#>  8 CNR2         
#>  9 CNR2-MIR378F 
#> 10 CTSK         
#> # i 47 more rows
#> 
#> $sv_tsv$total_variants
#> [1] 99
#> 
#> 
```
