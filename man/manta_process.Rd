% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sv.R
\name{manta_process}
\alias{manta_process}
\title{Process Manta Object}
\usage{
manta_process(manta_tsv_obj)
}
\arguments{
\item{manta_tsv_obj}{Manta list object read via \code{sv_prioritize_old}.}
}
\value{
List with:
\itemize{
\item total Manta (unmelted) variants
\item tibble with melted variants
\item genes involved in multi-gene events
}
}
\description{
Process Manta Object
}
