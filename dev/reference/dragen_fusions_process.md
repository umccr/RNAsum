# Process DRAGEN Fusions

Process DRAGEN Fusions

## Usage

``` r
dragen_fusions_process(dragen.fusions, known_translocations, genes_cancer)
```

## Arguments

- dragen.fusions:

  Parsed tibble with DRAGEN fusions.

- known_translocations:

  Tibble with known translocations.

- genes_cancer:

  Character vector of cancer genes.

## Value

Tibble with post-processed DRAGEN fusion calls.
