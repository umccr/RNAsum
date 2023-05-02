#' Perform PCA
#'
#' Perform PCA. Outputs a list with dataframe and samples colouring info ready for plotting.
#'
#' @param data Input data.
#' @param targets Target group.
#' @param title Title.
#' @param report_dir Report directory.
#' @param suffix Suffix for output file.
#'
#' @importFrom magrittr %>%
#' @return Outputs a list with data frame and samples colouring info ready for plotting.
#' @export
pca <- function(data, targets, title = "", report_dir, suffix = "") {
  ##### Keep only genes with variance > 0 across all samples
  rsd <- apply(data, 1, stats::sd)
  data.subset <- data[rsd > 0, ]

  ##### Perform PCA
  data.subset_pca <- stats::prcomp(t(data.subset), scale = FALSE)

  ##### Get variance importance for all principal components
  importance_pca <- summary(data.subset_pca)$importance[2, ]
  importance_pca <- paste(round(100 * importance_pca, 2), "%", sep = "")
  names(importance_pca) <- names(summary(data.subset_pca)$importance[2, ])

  ##### Prepare data frame
  data.subset_pca.df <- data.frame(
    targets$Target, targets$Dataset,
    data.subset_pca$x[, "PC1"],
    data.subset_pca$x[, "PC2"],
    data.subset_pca$x[, "PC3"]
  )
  colnames(data.subset_pca.df) <- c("Target", "Dataset", "PC1", "PC2", "PC3")

  ##### Assign colours to targets and datasets
  targets.colour <- getColours(targets$Target)
  datasets.colour <- getColours(targets$Dataset)

  ##### Create a list with dataframe and samples colouring info
  pca.list <- list(data.subset_pca.df, importance_pca, targets.colour, datasets.colour)
  names(pca.list) <- c("pca.df", "importance_pca", "targets", "datasets")

  ##### Change the datasets levels order
  data.subset_pca.df$Target <- factor(data.subset_pca.df$Target, levels = unique(data.subset_pca.df$Target))

  ##### Generate PCA 2-D plot
  pca_plot <- plotly::plot_ly(
    data.subset_pca.df,
    x = ~PC1, y = ~PC2, color = ~Target,
    text = paste(targets$Target, rownames(data.subset_pca.df), sep = ": "),
    colors = targets.colour[[1]], type = "scatter", mode = "markers",
    marker = list(size = 10, opacity = 0.7), width = 800, height = 500
  ) |>
    plotly::layout(
      title = title,
      xaxis = list(title = paste("PC1", " (", importance_pca["PC1"], ")", sep = "")),
      yaxis = list(title = paste("PC2", " (", importance_pca["PC2"], ")", sep = "")),
      margin = list(l = 50, r = 50, b = 50, t = 30, pad = 4), autosize = FALSE,
      showlegend = TRUE, legend = list(orientation = "v", y = 0.9)
    )

  ##### Generate Scree-plot
  data.subset_scree.df <- data.frame(
    paste0("PC ", c(1:length(importance_pca))),
    as.numeric(gsub("%", "", importance_pca))
  )
  colnames(data.subset_scree.df) <- c("PC", "Variances")

  ##### The default order will be alphabetized unless specified as below
  data.subset_scree.df$PC <- factor(data.subset_scree.df$PC, levels = data.subset_scree.df[["PC"]])

  scree_plot <- plotly::plot_ly(
    data.subset_scree.df,
    x = ~PC, y = ~Variances, type = "bar",
    width = 800, height = 350
  ) |>
    plotly::layout(
      title = title, xaxis = list(title = ""),
      margin = list(l = 50, r = 50, b = 100, t = 30, pad = 4), autosize = F
    )

  ##### Create directory for the plots
  PCAplotDir <- file.path(report_dir, "InputDataPlots")
  fs::dir_create(PCAplotDir)

  ##### Save interactive plot as html file
  saveWidgetFix(widget = pca_plot, file = file.path(PCAplotDir, glue::glue("pca_plot{suffix}.html")))
  saveWidgetFix(widget = scree_plot, file = file.path(PCAplotDir, glue::glue("scree_plot{suffix}.html")))

  list(pca.list, pca_plot, scree_plot)
}
