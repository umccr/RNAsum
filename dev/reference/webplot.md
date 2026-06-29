# Generate spider web plots to present immunogram genes

Generate spider web plots to present immunogram genes (from
http://www.statisticstoproveanything.com/2013/11/spider-web-plots-in-r.html).

## Usage

``` r
webplot(
  data,
  data.row = NULL,
  y.cols = NULL,
  main = NULL,
  add = F,
  col = "red",
  lty = 1,
  scale = T
)
```

## Arguments

- data:

  Dataframe or matrix.

- data.row:

  Row of data to plot (if NULL uses row 1).

- y.cols:

  Columns of interest (if NULL it selects all numeric columns).

- main:

  Title of plot (if NULL then rowname of data).

- add:

  Whether the plot should be added to an existing plot.

- col:

  Color of the data line

- lty:

  Lty of the data line

- scale:

  Scale.

## Value

Spider web plots to present immunogram genes
