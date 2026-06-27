# Generate CDF plot for selected gene

Generates cumulative distribution function (CDF) plot for selected gene.
If option "addBoxPlot" = TRUE, then generate additional boxplot below to
show the data variance for selected gene in individual groups

## Usage

``` r
cdfPlot(
  gene,
  data,
  targets,
  sampleName,
  int_cancer,
  ext_cancer,
  comp_cancer,
  add_cancer = NULL,
  addBoxPlot = FALSE,
  scaling = "gene-wise",
  report_dir
)
```

## Arguments

- gene:

  Gene of interest.

- data:

  Input data.

- targets:

  Target data set.

- sampleName:

  Sample name.

- int_cancer:

  Internal cancer group.

- ext_cancer:

  External cancer group.

- comp_cancer:

  Complete cancer group.

- add_cancer:

  Used for reordering groups.

- addBoxPlot:

  Add boc plot (boolean).

- scaling:

  Gene-wise or sample-wise.

- report_dir:

  Output directory.

## Value

Cumulative distribution function (CDF) plot for selected gene.
