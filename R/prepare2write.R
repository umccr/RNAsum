##### Prepare object to write into a file
prepare2write <- function (x) {

  x2write <- cbind(rownames(x), x)
  colnames(x2write) <- c("",colnames(x))

  ##### Clean the space and return output
  rm(x)
  return(x2write)
}
