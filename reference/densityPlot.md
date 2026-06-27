# Generate density and expression distribution plots for selected gene

Generate density and expression distribution plots for selected gene,
highlighting samples of interest.

## Usage

``` r
densityPlot(
  gene,
  data,
  main_title,
  x_title,
  sampleName,
  distributions = NULL,
  scaling = "gene-wise"
)
```

## Arguments

- gene:

  Selected gene.

- data:

  Input data.

- main_title:

  Main title.

- x_title:

  X-axis title.

- sampleName:

  Sample name.

- distributions:

  Distributions of interest.

- scaling:

  Scaling.

## Value

Density and expression distribution plots for selected gene.
