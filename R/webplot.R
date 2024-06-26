#' Generate spider web plots to present immunogram genes
#'
#' Generate spider web plots to present immunogram genes
#' (from http://www.statisticstoproveanything.com/2013/11/spider-web-plots-in-r.html).
#'
#' @param data Dataframe or matrix.
#' @param data.row Row of data to plot (if NULL uses row 1).
#' @param y.cols Columns of interest (if NULL it selects all numeric columns).
#' @param main Title of plot (if NULL then rowname of data).
#' @param add Whether the plot should be added to an existing plot.
#' @param col Color of the data line
#' @param lty Lty of the data line
#' @param scale Scale.
#'
#' @return Spider web plots to present immunogram genes
#' @export
webplot = function(data, data.row = NULL, y.cols = NULL, main = NULL, add = F,
                   col = "red", lty = 1, scale = T) {
  if (!is.matrix(data) & !is.data.frame(data))
    stop("Requires matrix or data.frame")
  if (is.null(y.cols))
    y.cols = colnames(data)[sapply(data, is.numeric)]
  if (base::sum(!sapply(data[, y.cols], is.numeric)) > 0) {
    out = paste0("\"", colnames(data)[!sapply(data, is.numeric)], "\"",
                 collapse = ", ")
    stop(paste0("All y.cols must be numeric\n", out, " are not numeric"))
  }
  if (is.null(data.row))
    data.row = 1
  if (is.character(data.row))
    if (data.row %in% rownames(data)) {
      data.row = which(rownames(data) == data.row)
    } else {
      stop("Invalid value for data.row:\nMust be a valid rownames(data) or row-index value")
    }
  if (is.null(main))
    main = rownames(data)[data.row]
  if (scale == T) {
    data = scale(data[, y.cols])
    data = apply(data, 2, function(x) x/max(abs(x)))
  }
  data = as.data.frame(data)
  n.y = length(y.cols)
  min.rad = 360/n.y
  polar.vals = (90 + seq(0, 360, length.out = n.y + 1)) * pi/180

  if (add == F) {
    graphics::plot(0, xlim = c(-2.2, 2.2), ylim = c(-2.2, 2.2), type = "n", axes = F,
         xlab = "", ylab = "")
    graphics::title(main)
    lapply(polar.vals, function(x) graphics::lines(c(0, 2 * cos(x)), c(0, 2 * sin(x))))
    lapply(1:n.y, function(x) graphics::text(2.15 * cos(polar.vals[x]), 2.15 * sin(polar.vals[x]),
                                   y.cols[x], cex = 0.8))

    lapply(seq(0.5, 2, 0.5), function(x) graphics::lines(x * cos(seq(0, 2 * pi, length.out = 100)),
                                               x * sin(seq(0, 2 * pi, length.out = 100)), lwd = 0.5, lty = 2, col = "gray60"))
    graphics::lines(cos(seq(0, 2 * pi, length.out = 100)), sin(seq(0, 2 * pi, length.out = 100)),
          lwd = 1.2, col = "gray50")
  }

  r = 1 + data[data.row, y.cols]
  xs = r * cos(polar.vals)
  ys = r * sin(polar.vals)
  xs = c(xs, xs[1])
  ys = c(ys, ys[1])
  graphics::lines(xs, ys, col = col, lwd = 2, lty = lty)

  #### Clear plots to free up some memory
  if(!is.null(grDevices::dev.list())) invisible(grDevices::dev.off())
}
