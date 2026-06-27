# Wrapper for saveWidget bug fix

A wrapper to saveWidget which compensates for bug where the temporary
directory with widget dependencies is not deleted - see
https://github.com/ramnathv/htmlwidgets/issues/296.

## Usage

``` r
saveWidgetFix(widget, file, selfcontained = TRUE, ...)
```

## Arguments

- widget:

  Widget to save.

- file:

  Full path to save HTML into.

- selfcontained:

  Whether to save the HTML as a single self-contained file (with
  external resources base64 encoded) or a file with external resources
  placed in an adjacent directory.

- ...:

  Additional arguments passed to
  [htmlwidgets::saveWidget](https://rdrr.io/pkg/htmlwidgets/man/saveWidget.html).

## Value

Change work dir to compensate for saveWidget bug
