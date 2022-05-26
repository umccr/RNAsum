##### Generate density and expression distribution plots for selected gene, highlighting samples of interest
densityPlot <- function(gene, data, main_title, x_title, sampleName, distributions = NULL, scaling = "gene-wise") {

  if ( scaling == "gene-wise" ) {
    data.z <- t(scale(t(data)))
  } else {
    data.z <- scale(data, scale = FALSE)
  }

  ##### Used data for user-defined genes
  data.z <- data.z[ gene, ,drop=FALSE]

  ##### Create data frame and fill it with expression and density values for each sample for selected gene
  data.df <- data.frame(gene = "Observed distribution", sample = colnames(data.z)[order(data.z)], expr = sort(data.z), dens = density2freq(density(data.z, n=ncol(data.z))$y))

  ##### Generate values to generate various distributions
  if ( !is.null(distributions) ) {

    ##### Use the density values obtained from the expression values
    expr.sorted <- sort(data.z)

    ##### Get min and max values based on the expression data
    data.x <- seq(min(expr.sorted), max(expr.sorted), length.out = ncol(data.z))

    ##### Create empty data frame
    data.df.dist <- data.frame(matrix(ncol = 4, nrow = 0))
    colnames(data.df.dist) <- c("gene", "sample", "expr", "dens")

    ##### Generate y-values to mirror distributions of interest
    ##### Generate y-values for normal distribution. Useful resource https://stats.idre.ucla.edu/r/modules/probabilities-and-distributions/
    if ( "normal" %in% tolower(distributions) ) {
      data.y <- dnorm(data.x, mean = mean(data.x), sd = (max(data.x)-mean(data.x))/5)
      data.df.dist <- rbind(data.df.dist, data.frame(gene="Normal distribution", sample = colnames(data.z)[order(data.z)], expr=data.x, dens=density2freq(data.y)))
    }

    ##### Generate x- and y-values for binomial distribution. Useful link https://stat.ethz.ch/R-manual/R-devel/library/stats/html/Binomial.html
    if ( "binomial" %in% tolower(distributions) ) {
      data.x <- 1:ncol(data.z)
      data.y <- dbinom(data.x, ncol(data.z), 0.25)
      data.x <- rescale(data.x, to = c(min(expr.sorted), max(expr.sorted)))
      data.df.dist <- rbind(data.df.dist, data.frame(gene="Binomial distribution (p=0.25)", sample = colnames(data.z)[order(data.z)], expr=data.x, dens=density2freq(data.y)))

      data.x <- 1:ncol(data.z)
      data.y <- dbinom(data.x, ncol(data.z), 0.75)
      data.x <- rescale(data.x, to = c(min(expr.sorted), max(expr.sorted)))
      data.df.dist <- rbind(data.df.dist, data.frame(gene="Binomial distribution (p=0.75)", sample = colnames(data.z)[order(data.z)], expr=data.x, dens=density2freq(data.y)))
    }

    ##### Draw n/2 samples from a normal distributions with one median and another n/2 samples from a second normal distribution with a different median. Useful link                  https://stats.stackexchange.com/questions/355344/simulating-a-bimodal-distribution-in-the-range-of-15-in-r
    if ( "bimodal" %in% tolower(distributions) ){
      data.x1 <- seq(min(expr.sorted), median(expr.sorted), length.out = ncol(data.z)/2)
      data.x2 <- seq(median(expr.sorted), max(expr.sorted), length.out = ncol(data.z)/2)

      ##### Combine both normal distributions to generate a bimodal distribution. Make sure the the length of this vector is equal to the number samples in the data
      data.x <- c(data.x1, data.x2)
      data.x <- data.x[1:ncol(data.z)]

      ##### Generate y-values for bimodal distribution
      data.y <- c(dnorm(data.x1, mean = mean(data.x1), sd = (max(data.x1)-mean(data.x1))/3), dnorm(data.x2, mean = mean(data.x2), sd = (max(data.x2)-mean(data.x2))/3))
      data.y <- data.y[1:ncol(data.z)]

      ##### Add bimodal dist values to the distribution dataframe
      data.df.dist <- rbind(data.df.dist, data.frame(gene = "Bimodal distribution", sample = colnames(data.z)[order(data.z)], expr = data.x, dens = density2freq(data.y)))
    }

    data.df <- rbind(data.df, data.df.dist)

    ##### Extract expression for selected sample in the distributions dataframe
    data.df.selected <- data.df[ sampleName == data.df$sample, ]
  }

  ##### Get min and max values based on the expression data
  den.x <- sort(data.df$expr)
  den.y <- sort(data.df$dens)

  ##### Assign colours to distributions
  genes.colour <- getColours(rev(unique(data.df$gene)))

  ##### Generate interactive density plot
  p <- plot_ly(data.df, x = ~expr, y = ~dens, type = 'scatter', mode = 'lines', color = ~gene, colors = genes.colour[[1]], width = 750, height = 200) %>%
    add_markers(y = data.df.selected$dens, x = data.df.selected$expr,
                name = "Patient",
                text = "Patient",
                mode = 'markers',
                marker = list(size = 8, colors = data.df.selected$sample, color = rep(I("black"), each = nrow(data.df.selected)), line = list(color = "grey", width = 2)),
                showlegend = TRUE,
                inherit = FALSE) %>%
    layout(title = main_title,
           xaxis = list(title = x_title, range = c(den.x[1],den.x[length(den.x)])),
           yaxis = list (title = 'Weight', range = c(den.y[1],den.y[length(den.y)]), side = "right"),
           legend = list(orientation = 'h', y = 1.3))

  return( p )

  ##### Clean the space
  rm(gene, expr.sorted)
  rm(list = ls(pattern='^data*'))

  #### Clear plots to free up some memory
  if(!is.null(dev.list())) invisible(dev.off())
}
