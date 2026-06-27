# Generate bar-plot for selected genes

Generates bar-plot for selected genes, highlighting samples of interest.

## Usage

``` r
barPlot(
  gene,
  data,
  targets,
  y_title = "Counts",
  sampleName,
  ext_cancer,
  int_cancer,
  add_cancer = NULL
)
```

## Arguments

- gene:

  User defined gene.

- data:

  Input data.

- targets:

  Targets.

- y_title:

  Title for y-axis.

- sampleName:

  Sample name.

- ext_cancer:

  External cancer group.

- int_cancer:

  Internal cancer group.

- add_cancer:

  Used for reordering groups.

## Value

Bar-plot for selected genes, highlighting samples of interest.
