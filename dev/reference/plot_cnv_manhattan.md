# Manhattan Plot for Copynumber Data

Manhattan Plot for Copynumber Data

## Usage

``` r
plot_cnv_manhattan(
  cnv_data,
  genes_to_highlight = NULL,
  title = "",
  show_legend = TRUE,
  point_size = 3,
  highlight_size = 8,
  cn_top = NULL,
  cn_bottom = NULL
)
```

## Arguments

- cnv_data:

  CN data frame including Gene, P, Diff_Perc, Diff_Z_score, Alterations,
  GENEBIOTYPE, ENSEMBL, SEQNAME, GENESEQSTART, GENESEQEND,
  Immunic_Cycle_Role.

- genes_to_highlight:

  Atomic vector of gene names to highlight.

- title:

  Plot title.

- show_legend:

  Show legend.

- point_size:

  Point size.

- highlight_size:

  Highlight size.

- cn_top:

  Top CN threshold.

- cn_bottom:

  Bottom CN threshold.

## Value

Plotly object.

## Examples

``` r
if (FALSE) { # \dontrun{
cnv_data <- data.annot
genes_to_highlight <- unique(data.sub$Gene)
} # }
```
