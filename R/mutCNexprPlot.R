#' Generate scatterplot with per-gene expression values (y-axis), CN values (x-axis) and mutation status info (colours), if provided
#'
#' Generates scatterplot with per-gene expression values (y-axis), CN values (x-axis) and mutation status info (colours), if provided
#'
#' @param data Input data.
#' @param alt_data Boolean (Generates scatterplot with per-gene expression values (y-axis)).
#' @param cn_bottom Bottom value for copy number.
#' @param cn_top Top value for copy number.
#' @param comp_cancer Complete cancer group.
#' @param type Type.
#' @param report_dir Report directory.
#'
#' @return Scatterplot with per-gene expression values (y-axis), CN values (x-axis) and mutation status info (colours), if provided.
#' @export
mutCNexprPlot <- function(data, alt_data = FALSE, cn_bottom = cn_bottom, cn_top = cn_top, comp_cancer, type = "z", report_dir) {
  ##### Extract info for genes to be annotated on the plot
  genes2annot <- data |>
    dplyr::filter(.data$CN >= cn_top | .data$CN <= cn_bottom) |>
    dplyr::pull(.data$Gene)

  if (length(genes2annot) == 0) {
    genes2annot <- ""
  }

  if (type == "z") {
    names(data)[names(data) %in% "Diff_Z_score"] <- "Expr"
    y_title <- paste0("mRNA expression (Z-score [Patient vs ", comp_cancer, "])")
  } else if (type == "perc") {
    names(data)[names(data) %in% "Diff_Perc"] <- "Expr"
    y_title <- paste0("mRNA expression (percentile [Patient vs ", comp_cancer, "])")
  }

  ##### Generate scatterplot with per-gene expression values (y-axis) (difference between Patient's and [comp_cancer] data), CN values (x-axis) and mutation status info (colours)
  if (alt_data) {
    p <- plotly::plot_ly(
      type = "scatter", mode = "markers", width = 800, height = 600, showlegend = FALSE
    ) |>
      plotly::add_markers(
        data = data, y = ~Expr, x = ~CN,
        name = ~Gene,
        text = paste0("Gene: ", data$Gene, "\nAlterations: ", data$Alterations),
        mode = "markers",
        marker = list(size = 10, symbol = "circle"),
        color = ~Gene,
        showlegend = TRUE,
        # legendtitle=TRUE,
        inherit = FALSE
      ) |>
      plotly::add_annotations(
        data = data |> dplyr::filter(.data$CN >= cn_top | .data$CN <= cn_bottom),
        text = genes2annot,
        x = ~CN, xanchor = "left",
        y = ~Expr, yanchor = "top",
        font = list(color = "Grey", size = 10),
        legendtitle = TRUE, showarrow = FALSE
      ) |>
      plotly::layout(
        xaxis = list(title = "CN value"), yaxis = list(title = y_title),
        margin = list(l = 50, r = 50, b = 50, t = 50, pad = 4),
        autosize = FALSE, legend = list(orientation = "v", x = 1, y = 0.97, yanchor = "top"),
        showlegend = TRUE
      )
    ##### Generate scatterplot with per-gene expression values (y-axis) and CN values (x-axis)
  } else {
    p <- plotly::plot_ly(
      data,
      x = ~CN, y = ~Expr, text = ~Gene, color = ~Gene, type = "scatter",
      mode = "markers", marker = list(size = 10, symbol = "circle"),
      width = 800, height = 600
    ) |>
      plotly::add_annotations(
        data = data |> dplyr::filter(.data$CN >= cn_top | .data$CN <= cn_bottom),
        text = ~Gene, x = ~CN, xanchor = "left", y = ~Expr, yanchor = "top",
        font = list(color = "Grey", size = 10),
        legendtitle = TRUE, showarrow = FALSE
      ) |>
      plotly::layout(
        xaxis = list(title = "CN value"), yaxis = list(title = y_title),
        margin = list(l = 50, r = 50, b = 50, t = 50, pad = 4), autosize = FALSE,
        legend = list(orientation = "v", y = 0.8, yanchor = "top"), showlegend = TRUE
      )
  }

  ##### Create directory for the plots
  mutCNexprPlotDir <- file.path(report_dir, "cn_expr_plot")
  fs::dir_create(mutCNexprPlotDir)

  ##### Save interactive plot as html file
  RNAsum::saveWidgetFix(widget = p, file = file.path(mutCNexprPlotDir, paste0("cn_expr_plot.", type, ".html")))
  return(p)
}
