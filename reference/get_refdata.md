# Get Reference Data File Paths

Get a list of paths to internal and external reference data, based on a
selected dataset name.

## Usage

``` r
get_refdata(dataset, batch_rm = FALSE)
```

## Arguments

- dataset:

  Reference RNAsum dataset of interest.

- batch_rm:

  Remove batch-associated effects between datasets.

## Value

List with paths to internal and external reference datasets (Counts,
Targets and Name).

## Examples

``` r
x <- get_refdata(dataset = "TEST")
```
