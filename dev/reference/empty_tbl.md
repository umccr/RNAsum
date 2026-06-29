# Create Empty Tibble Given Column Names

From https://stackoverflow.com/a/62535671/2169986. Useful for handling
edge cases with empty data.

## Usage

``` r
empty_tbl(cnames, ctypes = readr::cols(.default = "c"))
```

## Arguments

- cnames:

  Character vector of column names to use.

- ctypes:

  Character vector of column types corresponding to `cnames`.

## Value

A tibble with 0 rows and the given column names.
