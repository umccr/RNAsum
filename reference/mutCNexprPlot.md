# Generate scatterplot with per-gene expression values (y-axis), CN values (x-axis) and mutation status info (colours), if provided

Generates scatterplot with per-gene expression values (y-axis), CN
values (x-axis) and mutation status info (colours), if provided

## Usage

``` r
mutCNexprPlot(
  data,
  alt_data = FALSE,
  cn_bottom = cn_bottom,
  cn_top = cn_top,
  comp_cancer,
  type = "z",
  report_dir
)
```

## Arguments

- data:

  Input data.

- alt_data:

  Boolean (Generates scatterplot with per-gene expression values
  (y-axis)).

- cn_bottom:

  Bottom value for copy number.

- cn_top:

  Top value for copy number.

- comp_cancer:

  Complete cancer group.

- type:

  Type.

- report_dir:

  Report directory.

## Value

Scatterplot with per-gene expression values (y-axis), CN values (x-axis)
and mutation status info (colours), if provided.
