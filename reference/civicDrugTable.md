# Generate table with drugs targeting selected set of genes

Generates table with drugs targeting selected set of genes using info
from CIViC database (https://civicdb.org/).

## Usage

``` r
civicDrugTable(
  genes,
  civic_var_summaries,
  civic_clin_evid,
  evid_type = "Predictive",
  var_type = NULL
)
```

## Arguments

- genes:

  Genes of interest.

- civic_var_summaries:

  CIViC variant summaries.

- civic_clin_evid:

  CIViC clinical evidence.

- evid_type:

  Evidence type.

- var_type:

  Variant type.

## Value

Table with drugs targeting selected set of genes.
