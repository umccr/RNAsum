# Get Known Translocations

Flag known fusions based on info from:

- FusionGDB (https://ccsm.uth.edu/FusionGDB)

- Cancer Biomarkers database (CGI)
  (https://www.cancergenomeinterpreter.org/biomarkers)

## Usage

``` r
known_trans(kt_fusiongdb, kt_cgi)
```

## Arguments

- kt_fusiongdb:

  FusionGDB parsed tibble.

- kt_cgi:

  CGI parsed tibble.

## Value

Tibble with known translocations.
