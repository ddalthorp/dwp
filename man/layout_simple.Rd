% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{layout_simple}
\alias{layout_simple}
\title{Example Simple Geometry Data Format for Site Layout}
\format{
A data frame illustrating an import format for a simple description of a
 site layout by turbine. Each turbine (\code{turbine}) is classed according to
 the shape of its search plot, either \code{circular}, \code{square}, or
 \code{RP} (search on the roads and turbine pad only). For \code{circular}
 plots, all ground within \code{radius} meters of the turbine is searched.
 For \code{square} plots, the \code{radius} is half the width of the square
 along the x-axis, NOT the distance to the corner. For \code{RP} plots, the
 \code{radius} is the maximum distance searched on the roads. The geometry of
 the \code{RP} also includes a circular turbine pad with radius \code{padrad},
 a road width of \code{roadwidth} meters, and the number of roads
 (\code{n_road}) searched out to \code{radius} meters from the turbine.
}
\usage{
layout_simple
}
\description{
Example Simple Geometry Data Format for Site Layout
}
\keyword{datasets}
