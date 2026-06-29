# Process CGI Known Translocations

Processes parsed CGI known translocations tibble, removing duplicated
translocations that have a different source.

## Usage

``` r
known_translocations_cgi_process(kt_cgi)
```

## Arguments

- kt_cgi:

  Parsed CGI tibble.

## Value

Tibble with dups removed.
