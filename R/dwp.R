#' @importFrom graphics axis box legend lines mtext par points polygon rect text
#' @importFrom grDevices colors
#' @importFrom magrittr %>% %<>%
#' @importFrom methods as is
#' @importFrom stats aggregate approxfun dbinom dchisq dexp dgamma dlnorm dnorm
#'  formula glm integrate pbinom pchisq pexp pgamma plnorm pnorm quantile rnorm
#'  runif uniroot var
#' @importFrom utils flush.console read.csv write.csv
utils::globalVariables(c(".", "cof_name", "constraints", "mod_all", "mod_color",
 "mod_name", "mod_offset", "mod_standard", "mod_xy", "natural", "par_default",
 "parm_name", "sieve_default", "mod_lty", "xyr"))

#' @title Density-Weighted Proportion
#'
#' @description This package is designed to analyze carcass dispersion data and
#'  fit models of carcass density as function of distance from turbine.
#'
#' @name dwp
#' 
#' @section Data sets:
#' \code{\link{carcass_polygon}}\cr
#' \code{\link{carcass_simple}}\cr
#' \code{\link{layout_eagle}}\cr
#' \code{\link{layout_polygon}}\cr
#' \code{\link{layout_simple}}\cr
#' \code{\link{layout_xy}}\cr
#' \code{\link{xyr}}\cr
#' \code{\link{sieve_default}}\cr
#'
#' @section Main Command-Line Functions:
#'  \describe{
#'   \item{\code{\link{initLayout}, \link{prepRing}, \link{readCarcass}, \link{addCarcass}}}{import and format data}
#'   \item{\code{\link{ddFit}}}{fit carcass distribution models}
#'   \item{\code{\link{estpsi}, \link{estdwp}}}{estimate probability that carcass 
#'      will lie in the searched area (psi) and the fraction of carcasses
#'      lying in the searched area (`dwp`)}
#'  \item{\code{\link{formatGenEst}, \link{exportGenEst}}}{format and export
#'     `dwphat` objects for use with GenEst}
#'   \item{\code{\link{aic}, \link{modelFilter}, \link{stats}, \link{ddCI}}}{statistics 
#'      for fitted models}
#'   \item{\code{plot}}{S3 function for \code{\link[=Plot]{ddArray}},
#'     \code{\link[=Plot]{dd}}, \code{\link[=Plot]{fmod}},
#'     \code{\link[=Plot]{polygonLayout}}, 
#'     \code{\link[=Plot]{psiHat}},
#'     \code{\link[=Plot]{dwpHat}} objects.}
#'   \item{\code{\link{ddd}, \link[=ddd]{pdd}, \link[=ddd]{qdd}, \link[=ddd]{rdd}, \link[=ddd]{rcd}}}{probability 
#'      functions for distance distributions}
#'  }
#' @section Potentially Useful Calculation and Editing Functions:
#'  \describe{ 
#'   \item{\code{\link{ddSim}, \link{dd2ddSim}}}{functions for simulating \code{dd} models},
#'   \item{\code{\link{getncarc}}}{extract the number of carcasses per turbine from a data set; 
#'      method for many types of objects},
#'   \item{\code{\link{cof2parms}, \link{cofOK}, \link[=cofOK]{cofOK0}, \link[=cofOK]{cofOKInf}, 
#'      \link{constraints}}}{functions for manipulating and checking model coefficients}
#'   \item{\code{\link{Acins}}}{calculate the area of the intersection of a circle
#'      and square sharing a common center}
#'   \item{\code{\link{rmat}, \link{off}}}{functions for constructing functions
#'      out of distribution information}
#'   \item{\code{\link{exclude}}}{simple function for excluding items from a 
#'      superset}
#'  }
NULL
