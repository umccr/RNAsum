# Read Arriba TSV File

Reads the TSV file output by Arriba.

## Usage

``` r
arriba_tsv_read(x = NULL)
```

## Arguments

- x:

  Path to Arriba TSV file.

## Value

A tibble with the contents of the input TSV file, or NULL if input is
NULL or TSV has no data rows.

## Examples

``` r
x <- system.file("rawdata/test_data/dragen/arriba/fusions.tsv", package = "RNAsum")
(a <- arriba_tsv_read(x))
#> # A tibble: 4 x 30
#>   gene1     gene2      `strand1(gene/fusion)` `strand2(gene/fusion)` breakpoint1
#>   <chr>     <chr>      <chr>                  <chr>                  <chr>      
#> 1 EP300-AS1 XRCC6      -/-                    +/+                    chr22:4119~
#> 2 LAPTM4A   AC098828.3 -/-                    -/-                    chr2:20051~
#> 3 FNIP1     AC112492.~ -/-                    ./+                    chr5:13164~
#> 4 ZFP64     MORC2      -/-                    -/-                    chr20:5209~
#> # i 25 more variables: breakpoint2 <chr>, site1 <chr>, site2 <chr>, type <chr>,
#> #   split_reads1 <dbl>, split_reads2 <dbl>, discordant_mates <dbl>,
#> #   coverage1 <dbl>, coverage2 <dbl>, confidence <chr>, reading_frame <chr>,
#> #   tags <chr>, retained_protein_domains <chr>,
#> #   closest_genomic_breakpoint1 <chr>, closest_genomic_breakpoint2 <chr>,
#> #   gene_id1 <chr>, gene_id2 <chr>, transcript_id1 <chr>, transcript_id2 <chr>,
#> #   direction1 <chr>, direction2 <chr>, filters <chr>, ...
```
