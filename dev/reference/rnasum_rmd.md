# Generate UMCCR RNAsum Report

Generates a UMCCR RNAsum HTML report. It does so with the following
steps:

1.  copy the rmd into tmp/rnasum.Rmd

2.  render the rmd inside tmp/

3.  return the path to the output HTML

## Usage

``` r
rnasum_rmd(out_file = NULL, quiet = FALSE, pars)
```

## Arguments

- out_file:

  Path to output HTML file (needs '.html' suffix).

- quiet:

  Suppress printing during rendering.

- pars:

  List of named parameters that override custom params specified within
  the Rmd's YAML front-matter.

## Value

Path to rendered HTML report.
