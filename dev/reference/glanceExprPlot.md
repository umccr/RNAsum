# Generate boxplot presenting expression profiles for selected set of genes

Generates boxplot presenting expression profiles for selected set of
genes.

## Usage

``` r
glanceExprPlot(
  genes,
  data,
  targets,
  sampleName,
  int_cancer,
  ext_cancer,
  comp_cancer,
  add_cancer = NULL,
  hexcode,
  type = "z",
  sort = "diff",
  scaling = "gene-wise",
  report_dir
)
```

## Arguments

- genes:

  Genes of interest.

- data:

  Input data.

- targets:

  Target groups.

- sampleName:

  Sample name.

- int_cancer:

  Internal cancer group.

- ext_cancer:

  External cancer group.

- comp_cancer:

  Complete cancer group.

- add_cancer:

  Addtional cancer type.

- hexcode:

  Hexcode.

- type:

  Type for expression values.

- sort:

  Sort order.

- scaling:

  Scaling type.

- report_dir:

  Reprot directory.

## Value

Boxplot presenting expression profiles for selected set of genes.
