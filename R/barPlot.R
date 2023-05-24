#' Generate bar-plot for selected genes
#'
#' Generates bar-plot for selected genes, highlighting samples of interest.
#'
#' @param gene User defined gene.
#' @param data Input data.
#' @param targets Targets.
#' @param y_title Title for y-axis.
#' @param sampleName Sample name.
#' @param ext_cancer External cancer group.
#' @param int_cancer Internal cancer group.
#' @param add_cancer Used for reordering groups.
#'
#' @importFrom magrittr %>%
#' @return Bar-plot for selected genes, highlighting samples of interest.
#' @export
barPlot <- function(gene, data, targets, y_title = "Counts", sampleName, ext_cancer, int_cancer, add_cancer = NULL) {
  ##### Used data for user-defined genes
  # uses rownames implicitly to filter
  data <- data[gene, , drop = FALSE]

  ##### Prepare data frame
  targets$Target[targets$Target == sampleName] <- "Patient"
  rownames(targets)[rownames(targets) == sampleName] <- "Patient"
  data.df <- data.frame(targets$Target, rownames(targets), as.numeric(data))
  colnames(data.df) <- c("Group", "Sample", "Data")

  ##### Reorder groups and add colours
  if (!is.null(add_cancer)) {
    data.df$Group <- factor(data.df$Group, levels = c(add_cancer, ext_cancer, int_cancer, "Patient"))
    group.colours <- c("forestgreen", "cornflowerblue", "red", "black")
  } else {
    data.df$Group <- factor(data.df$Group, levels = c(ext_cancer, int_cancer, "Patient"))
    group.colours <- c("cornflowerblue", "red", "black")
  }

  ##### The default order will be alphabetized unless specified as below
  data.df$Sample <- factor(data.df$Sample, levels = data.df[["Sample"]])
  p <- plotly::plot_ly(
    data.df,
    x = ~Sample, y = ~Data, color = ~Group, type = "bar",
    colors = group.colours, width = 750, height = 200
  ) %>%
    plotly::layout(
      title = "",
      xaxis = list(title = "", showticklabels = FALSE),
      yaxis = list(title = y_title), autosize = FALSE,
      legend = list(orientation = "h", y = 1.2),
      showlegend = TRUE
    )
  return(p)
}
