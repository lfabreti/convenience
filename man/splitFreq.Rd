% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/splitFreq.R
\name{splitFreq}
\alias{splitFreq}
\title{Split Frequencies}
\usage{
splitFreq(runs, windows = FALSE)
}
\arguments{
\item{runs}{A list of rwty.chain type}

\item{windows}{A boolean to set if the calculation is within runs or between runs}
}
\value{
A list with name of the split, its frequency and the difference in frequency for the comparison in place
}
\description{
Calculate the split frequencies and difference in split frequencies between runs or between windows of the same run
}
