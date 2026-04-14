#' Manhattan Plot for Copynumber Data
#'
#' @param cnv_data CN data frame including Gene, P, Diff_Perc, Diff_Z_score,
#' Alterations, GENEBIOTYPE, ENSEMBL, SEQNAME, GENESEQSTART, GENESEQEND,
#' Immunic_Cycle_Role.
#' @param genes_to_highlight Atomic vector of gene names to highlight.
#' @param title Plot title.
#' @param show_legend Show legend.
#' @param point_size Point size.
#' @param highlight_size Highlight size.
#' @param cn_top Top CN threshold.
#' @param cn_bottom Bottom CN threshold.
#'
#' @returns Plotly object.
#'
#' @examples
#' \dontrun{
#' cnv_data <- data.annot
#' genes_to_highlight <- unique(data.sub$Gene)
#' }
#' @export
plot_cnv_manhattan <- function(
  cnv_data,
  genes_to_highlight = NULL,
  title = "",
  show_legend = TRUE,
  point_size = 3,
  highlight_size = 8,
  cn_top = NULL,
  cn_bottom = NULL
) {
  chr_lengths <- cnv_data |>
    dplyr::group_by(.data$SEQNAME) |>
    dplyr::summarize(max_pos = max(.data$GENESEQSTART), .groups = "drop") |>
    dplyr::arrange(.data$SEQNAME) |>
    dplyr::mutate(offset = cumsum(c(0, utils::head(.data$max_pos, -1))))

  plot_data <- cnv_data |>
    dplyr::left_join(chr_lengths, by = "SEQNAME") |>
    dplyr::mutate(
      pos_cum = .data$GENESEQSTART + .data$offset,
      is_highlighted = .data$Gene %in% genes_to_highlight
    )

  data_normal <- plot_data |>
    dplyr::filter(!.data$is_highlighted)
  data_highlight <- plot_data |>
    dplyr::filter(.data$is_highlighted)

  chr_centers <- plot_data |>
    dplyr::group_by(.data$SEQNAME) |>
    dplyr::summarize(
      center = (min(.data$pos_cum) + max(.data$pos_cum)) / 2,
      .groups = "drop"
    ) |>
    dplyr::arrange(.data$SEQNAME)

  chr_bands <- plot_data |>
    dplyr::group_by(.data$SEQNAME) |>
    dplyr::summarize(
      xmin = min(.data$pos_cum),
      xmax = max(.data$pos_cum),
      .groups = "drop"
    ) |>
    dplyr::arrange(.data$SEQNAME) |>
    dplyr::mutate(
      color = ifelse(
        seq_len(dplyr::n()) %% 2 == 0,
        "rgba(220, 220, 220, 0.3)",
        "rgba(180, 180, 180, 0.3)"
      )
    )

  y_max <- max(plot_data$P, cn_top, na.rm = TRUE) * 1.05
  y_min <- min(plot_data$P, cn_bottom, na.rm = TRUE)
  x_min <- min(plot_data$pos_cum, na.rm = TRUE)
  x_max <- max(plot_data$pos_cum, na.rm = TRUE)

  chr_label <- function(chr) {
    dplyr::case_when(
      chr == 23L ~ "X",
      chr == 24L ~ "Y",
      chr == 25L ~ "MT",
      .default = as.character(chr)
    )
  }

  p <- plotly::plot_ly()

  for (i in seq_len(nrow(chr_bands))) {
    p <- p |>
      plotly::add_ribbons(
        x = c(chr_bands$xmin[i], chr_bands$xmax[i]),
        ymin = rep(y_min, 2),
        ymax = rep(y_max, 2),
        fillcolor = chr_bands$color[i],
        line = list(width = 0),
        showlegend = FALSE,
        hoverinfo = "skip"
      )
  }

  for (chr in sort(unique(data_normal$SEQNAME))) {
    chr_data <- data_normal |>
      dplyr::filter(.data$SEQNAME == chr)
    color <- ifelse(chr %% 2 == 0, "#404040", "#808080")

    p <- p |>
      plotly::add_trace(
        data = chr_data,
        x = ~pos_cum,
        y = ~P,
        type = "scatter",
        mode = "markers",
        marker = list(size = point_size, color = color, opacity = 0.6),
        text = ~ paste0(
          "<b>",
          Gene,
          "</b><br>",
          "(",
          format(GENESEQSTART, big.mark = ","),
          ", ",
          round(P, 2),
          ") chr",
          chr_label(SEQNAME),
          "<br>",
          "Gene: ",
          Gene,
          "<br>",
          "Diff_Z_score: ",
          round(Diff_Z_score, 2),
          "<br>",
          "Diff_Perc: ",
          round(Diff_Perc, 1),
          "<br>",
          "Alterations: ",
          Alterations
        ),
        hoverinfo = "text",
        showlegend = FALSE,
        name = paste("Chr", chr_label(chr))
      )
  }

  if (nrow(data_highlight) > 0) {
    unique_genes <- unique(data_highlight$Gene)
    for (gene_name in unique_genes) {
      gene_data <- data_highlight |>
        dplyr::filter(.data$Gene == gene_name)
      p <- p |>
        plotly::add_trace(
          data = gene_data,
          x = ~pos_cum,
          y = ~P,
          type = "scatter",
          mode = "markers",
          marker = list(
            size = highlight_size,
            color = "#E74C3C",
            line = list(color = "#C0392B", width = 1),
            opacity = 0.9
          ),
          text = ~ paste0(
            "<b>",
            Gene,
            "</b><br>",
            "(",
            format(GENESEQSTART, big.mark = ","),
            ", ",
            round(P, 2),
            ") chr",
            chr_label(SEQNAME),
            "<br>",
            "Gene: ",
            Gene,
            "<br>",
            "Diff_Z_score: ",
            round(Diff_Z_score, 2),
            "<br>",
            "Diff_Perc: ",
            round(Diff_Perc, 1),
            "<br>",
            "Alterations: ",
            Alterations
          ),
          hoverinfo = "text",
          name = gene_name,
          showlegend = show_legend
        )
    }
  }

  if (!is.null(cn_top)) {
    p <- p |>
      plotly::add_segments(
        x = x_min,
        xend = x_max,
        y = cn_top,
        yend = cn_top,
        line = list(color = "gray", width = 1, dash = "dash"),
        inherit = FALSE,
        showlegend = FALSE,
        hoverinfo = "skip"
      )
  }
  if (!is.null(cn_bottom)) {
    p <- p |>
      plotly::add_segments(
        x = x_min,
        xend = x_max,
        y = cn_bottom,
        yend = cn_bottom,
        line = list(color = "gray", width = 1, dash = "dash"),
        inherit = FALSE,
        showlegend = FALSE,
        hoverinfo = "skip"
      )
  }

  p1 <- p |>
    plotly::layout(
      title = list(text = title, x = 0.5, xanchor = "center"),
      xaxis = list(
        title = list(text = "Chromosome"),
        tickmode = "array",
        tickvals = chr_centers$center,
        ticktext = vapply(chr_centers$SEQNAME, chr_label, character(1)),
        showgrid = FALSE,
        zeroline = FALSE
      ),
      yaxis = list(
        title = list(text = "CN value"),
        showgrid = TRUE,
        gridcolor = "rgba(200, 200, 200, 0.5)",
        zeroline = TRUE
      ),
      hovermode = "closest",
      plot_bgcolor = "white",
      paper_bgcolor = "white",
      legend = list(
        orientation = "v",
        yanchor = "top",
        y = 1,
        xanchor = "left",
        x = 1.02,
        bgcolor = "rgba(255, 255, 255, 0.8)",
        bordercolor = "rgba(0, 0, 0, 0.2)",
        borderwidth = 1
      ),
      margin = list(r = 150)
    )
  p1
}
