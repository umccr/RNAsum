##### Generate box-plot for selected genes, highlighting samples of interest
barPlot <- function(gene, data, targets, y_title = "Counts", sampleName,  ext_cancer = ext_cancer_group, int_cancer = int_cancer_group, comp_cancer = comp_cancer_group, add_cancer = NULL ) {

  ##### Used data for user-defined genes
  data <- data[ gene, ,drop=FALSE]

  ##### Prepare data frame
  targets$Target[ targets$Target==sampleName ] <- "Patient"
  rownames(targets)[ rownames(targets)==sampleName ] <- "Patient"
  data.df <- data.frame(targets$Target, rownames(targets), as.numeric(data))
  colnames(data.df) <- c("Group","Sample", "Data")

  ##### Reorder groups and add colours
  if ( !is.null(add_cancer) ) {
    data.df$Group <- factor(data.df$Group, levels=c( add_cancer, ext_cancer, int_cancer, "Patient"))
    group.colours <- c("forestgreen", "cornflowerblue", "red", "black")
  } else {
    data.df$Group <- factor(data.df$Group, levels=c(ext_cancer, int_cancer, "Patient"))
    group.colours <- c("cornflowerblue", "red", "black")
  }

  ##### The default order will be alphabetized unless specified as below
  data.df$Sample <- factor(data.df$Sample, levels = data.df[["Sample"]])
  p <- plot_ly(data.df, x = ~Sample, y = ~Data, color = ~Group, type = 'bar', colors = group.colours, width = 750, height = 200) %>%
    layout(title = "", xaxis = list( title = "", showticklabels = FALSE), yaxis = list(title = y_title), autosize = F, legend = list(orientation = 'h', y = 1.2), showlegend=TRUE)

  return( p )

  ##### Clean the space
  rm(list = ls(pattern='^data*'))

  #### Clear plots to free up some memory
  if(!is.null(dev.list())) invisible(dev.off())
}
