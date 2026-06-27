# Read Arriba PDF File

Reads the PDF output by Arriba containing one fusion plot per page and
converts each page to a separate PNG file.

## Usage

``` r
arriba_pdf_read(pdf = NULL, fusions = NULL, outdir = NULL)
```

## Arguments

- pdf:

  Output PDF file from Arriba.

- fusions:

  Tibble with fusions from Arriba.

- outdir:

  Directory to write output PNGs to.

## Value

Tibble with paths to the created PNG images and their clean names. If
the PDF (or any of the required params) is NULL, returns NULL.

## Examples

``` r
pdf <- system.file("rawdata/test_data/dragen/arriba/fusions.pdf", package = "RNAsum")
tsv <- system.file("rawdata/test_data/dragen/arriba/fusions.tsv", package = "RNAsum")
fusions <- arriba_tsv_read(tsv)
(pngs <- arriba_pdf_read(pdf, fusions, tempdir()))
#> # A tibble: 4 x 2
#>   nm                                                                 png        
#>   <chr>                                                              <chr>      
#> 1 EP300.AS1__XRCC6_chr22.41195482_chr22.41636113                     /var/folde~
#> 2 LAPTM4A__AC098828.3_chr2.20051410_chr2.20064511                    /var/folde~
#> 3 FNIP1__AC112492.1.83716..DCAF8L1.218._chr5.131642783_chrX.27977775 /var/folde~
#> 4 ZFP64__MORC2_chr20.52096811_chr22.30963497                         /var/folde~
```
