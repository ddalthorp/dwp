% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fun.R
\name{subset.shapeCarcass}
\alias{subset.shapeCarcass}
\alias{subset.shapeLayout}
\title{Subset Data from an Imported and Formatted Shape File}
\usage{
\method{subset}{shapeCarcass}(x, subset, select, ...)

\method{subset}{shapeLayout}(x, subset, select, ...)
}
\arguments{
\item{x}{object to be subsetted}

\item{subset}{values to subset by. For example, to subset \code{x} to include
only turbines \code{"t1"} and \code{"t2"}, then \code{subset = c("t1", "t2")}.
The name of the column with turbine names is given in \code{select}.}

\item{select}{the name of the column with the values to subset by. For example,
to subset \code{x} by turbines names "t1" and "t2" as found in the "turbine"
column in the data, use \code{select = "turbine"} and \code{subset = c("t1", "t2")}.}

\item{...}{ignored}
}
\value{
object of the same class as \code{x}, subsetted to values of \code{select}
 equal to some element in \code{subset}.
}
\description{
Subset Data from an Imported and Formatted Shape File
}
