% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fun.R
\name{ddPrint}
\alias{ddPrint}
\alias{print.dd}
\alias{print.ddArray}
\alias{print.fmod}
\title{Print S3 Objects in \code{dwp} Package}
\usage{
\method{print}{dd}(x, ...)

\method{print}{ddArray}(x, ...)

\method{print}{fmod}(x, ...)
}
\arguments{
\item{x}{a \code{\link[=ddFit]{ddArray}} or \code{\link[=ddFit]{ddArray}} object}

\item{...}{ignored}
}
\value{
no return value; output printed to the console
}
\description{
\code{\link[=ddFit]{dd}}, \code{\link[=ddFit]{ddArray}}, and 
 \code{\link[=modelFilter]{fmod}} objects are lists consisting of a great
 amount of data. Only a few of the elements are printed automatically. Other 
 elements of object \code{x} can be viewed and extracted as with other lists 
 in R, namely, by using the \code{x$element} or \code{x[[element]]} operator, 
 where \code{element} is the name of one of the elements of \code{x}, all of
 which can be viewed via \code{names(x)}.
}
