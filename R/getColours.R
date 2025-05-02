#' Assign colours to different elements
#'
#' Assigns colours to different elements.
#'
#' @param elements Elements of interest.
#'
#' @return Colour assignment to different elements
#' @export
getColours <- function(elements) {
  ##### Predefined selection of colours for elements
  if (length(unique(elements)) == 2) {
    elements.colours <- c("cornflowerblue", "black")
  } else if (length(unique(elements)) == 3) {
    elements.colours <- c("cornflowerblue", "red", "black")
  } else if (length(unique(elements)) == 4) {
    elements.colours <- c("cornflowerblue", "forestgreen", "red", "black")
  } else {
    elements.colours <- grDevices::rainbow(length(elements))
  }

  f.elements <- factor(elements, levels = unique(elements))
  vec.elements <- elements.colours[seq_len(length(levels(f.elements)))]
  elements.colour <- rep(0, length(f.elements))
  for (i in seq_len(length(f.elements))) {
    elements.colour[i] <- vec.elements[f.elements[i] == levels(f.elements)]
  }

  list(vec.elements, elements.colour)
}
