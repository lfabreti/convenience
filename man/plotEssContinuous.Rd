% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotEssContinuous.R
\name{plotEssContinuous}
\alias{plotEssContinuous}
\title{Histogram of the ESS values}
\usage{
plotEssContinuous(
  x,
  per_run = FALSE,
  precision = 0.01,
  breaks = NULL,
  fill_color = NULL,
  filename = NULL,
  ...
)
}
\arguments{
\item{x}{A list of convenience.diag type}

\item{per_run}{If the plot should combine all runs or plot each run separately. Default value is FALSE}

\item{precision}{The precision of the mean estimates. Default is 0.01}

\item{breaks}{The number of bins to determine the intervals of the histogram}

\item{fill_color}{The color to fill the histogram bars}

\item{filename}{The name of the file to save the plot}

\item{...}{(various) Additional arguments passed to plot().}
}
\value{
Histogram
}
\description{
Plots the histogram of the ESS values for the continuous parameters
}
