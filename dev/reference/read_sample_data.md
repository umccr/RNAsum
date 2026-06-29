# Read Sample Data

Reads sample data, including Arriba fusions, Arriba plots, Salmon
counts, DRAGEN fusions, DRAGEN mapping metrics, SVs, copy-number data,
and PCGR SNVs.

## Usage

``` r
read_sample_data(p, results_dir, tx2gene = NULL)
```

## Arguments

- p:

  RNAsum params list.

- results_dir:

  Directory to output extracted Arriba PNGs to (created automatically if
  it does not already exist).

- tx2gene:

  data.frame with tx_name and gene_id columns. Required for salmon. See
  [tximport::tximport](https://rdrr.io/pkg/tximport/man/tximport.html).

## Value

A list of the input sample data.

## Examples

``` r
p <- list(
  dragen_wts_dir = system.file("rawdata/test_data/dragen", package = "RNAsum"),
  arriba_dir = system.file("rawdata/test_data/dragen/arriba", package = "RNAsum"),
  umccrise = system.file("rawdata/test_data", package = "RNAsum")
)
res <- read_sample_data(p, tempdir())
#> Error: Exactly one of salmon or kallisto count files are required
```
