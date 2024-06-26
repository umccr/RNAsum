% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/salmon.R
\name{salmon_counts}
\alias{salmon_counts}
\title{Get Salmon Counts}
\usage{
salmon_counts(x, tx2gene = NULL)
}
\arguments{
\item{x}{Path to \code{quant.sf} or \code{quant.genes.sf} file with abundances. See \link[tximport:tximport]{tximport::tximport}.}

\item{tx2gene}{data.frame with tx_name and gene_id columns. See \link[tximport:tximport]{tximport::tximport}.}
}
\value{
Tibble with the counts per gene transcript, or NULL if any of the
input params are NULL.
}
\description{
Get salmon counts via tximport.
}
\examples{
x <- system.file("rawdata/test_data/dragen/TEST.quant.sf", package = "RNAsum")
tx2gene <- NULL
(sc <- salmon_counts(x, tx2gene)) # NULL since no tx2gene specified
}
