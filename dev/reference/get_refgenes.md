# Get Reference Genes

Get a list of reference genes.

## Usage

``` r
get_refgenes(p)
```

## Arguments

- p:

  RNAsum params list.

## Value

List of tibbles corresponding to cancer, oncokb, immune and hrd genes

## Examples

``` r
p <- list(
  genes_cancer = system.file("extdata/genes/umccr_cancer_genes.2019-03-20.tsv.gz",
    package = "RNAsum.data"
  ),
  genes_oncokb = system.file(
    "extdata/OncoKB/CancerGenesList.txt.gz",
    package = "RNAsum.data"
  ),
  civic_var_summaries = system.file(
    "extdata/CIViC/01-Oct-2018-VariantSummaries.tsv.gz",
    package = "RNAsum.data"
  )
)
x <- get_refgenes(p)
```
