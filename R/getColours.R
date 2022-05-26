##### Assign colours to different elements
getColours <- function(elements) {

  ##### Predefined selection of colours for elements
  if ( length(unique(elements)) == 3 ) {
    elements.colours <- c("powderblue", "red", "gray50")
  } else if ( length(unique(elements)) == 4 ) {
    elements.colours <- c("powderblue", "forestgreen", "red", "gray50")
  } else {
    elements.colours <- rainbow(length(elements))
  }

  f.elements <- factor(elements, levels = unique(elements))
  vec.elements <- elements.colours[1:length(levels(f.elements))]
  elements.colour <- rep(0,length(f.elements))
  for (i in 1:length(f.elements))
    elements.colour[i] <- vec.elements[ f.elements[i]==levels(f.elements)]

  return( list(vec.elements, elements.colour) )
}
