#' Incomplete Gamma Function
#' @param x non-negative numeric vector
#' @param a positive numeric vector
#' @param lower boolean for calculating lower or upper incomplete gamma function
#' @details The upper incomplete gamma function, following Wolfram Alpha, namely,
#'  incGamma(a, x) = \eqn{\int_x^\infty e^{-t} * t^{a - 1}dt}{%
#'  integral(exp(-t) * t^(a - 1) dt from x to Inf)},
#'  calculated using pgamma. Within the \code{dwp} package, \code{incGamma} is 
#'  used in the calculation of the cumulative distribution function (CDF) of the
#'  xep02 distribution (\code{pxep02}). NOTE: The function \code{pracma::incgam} also
#'  calculates incomplete gamma with \code{pracma::incgam(x, a) = incGamma(a, x)},
#'  but \code{pracma::incgam} is not vectorized and not used here.
#' @return scalar or vector of length = max(length(x), length(a)), with values
#'  of the shorter recycled to match the length of the longer a la pnorm etc.
#' @export
incGamma <- function(a, x, lower = FALSE) pgamma(x, a, lower.tail = lower) * gamma(a)

#' @name Distributions
#' @title Probability Distributions for Carcasses Versus Distance from Turbine
#'
#' @description PDFs and CDFs that are required by \code{ddd}, \code{pdd} and 
#' \code{qdd} but are not included among the standard R distributions. Relying on
#' custom code and included here are the Maxwell-Boltzmann (\code{pmb} and 
#' \code{dmb}), xep0 (Pareto), xep1, xepi0 (inverse gamma), xep2 (Rayleigh), 
#' xep02, xep12, xep012, xep123, and xep0123. Not included here are the 
#' distributions that can be calculated using standard probability functions from
#' base R, namely the exponential, truncated normal, lognormal, gamma (xep01), and 
#' chisquared distributions and the inverse gaussian, which is calculated using 
#' \code{statmod::dinvgauss} and \code{statmod::dinvgauss}. The functions are 
#' designed for vector \code{x} or \code{q} and scalar parameters. 
#' 
#' @details An xep distribution is calculated by dividing its kernel (for the 
#' densities) or the integral of its kernel (for the cumulative distributions) 
#' by the normalizing constant, which is the integral of the kernel from 0 to Inf. 
#' The kernel of an xep distribution is defined as \eqn{x e^{P(x)}}{%
#' x * exp(P(x))}, 
#' where \eqn{P(x)} is a polynomial with terms defined by the suffix on xep. For 
#' example, the kernel of xep12 would be \eqn{x e^{b_1*x + b_2*x^2}}{%
#' x * exp(b1*x + b2*x^2)}. A \code{0} in the suffix indicates a \code{log(X)} 
#' term and an \code{i} indicates a \code{1/x} term. The parameters of the xep 
#' distributions are some combination of \eqn{b_0, b_1, b_2, b_3}{%
#' bi, b0, b1, b2, and b3}. The parameterizations of the
#' inverse gamma (xepi0), Rayleigh (xep2), and Pareto (xep0) follow the standard
#' conventions of \code{shape} and \code{scale} for the inverse gamma, \code{s2} = 
#' \eqn{s^2} for the Rayleigh, and \code{a} = \eqn{a} for the Pareto (with
#' a scale or location parameter of 1 and PDF = \eqn{a/x^(x + 1)} 
#' with support (1, Inf).
#' 
#' The Maxwell-Boltzmann is a one-parameter family with parameter \code{a} and PDF 
#' \eqn{f(a) = \sqrt{2/\pi}\frac{x^2 e^{-x^2/(2a^2)}}{a^3}}{%
#' f(a) = sqrt(2/pi) * x^2 * exp(-x^2/(2a^2))/(a^3)}. The kernel is
#' \eqn{f(a) = x^2 e^{-x^2}}{x^2 * exp(-x^2)}, which has a simple closed-form 
#' integral that involves the error function (\code{pracma::erf}).
#' 
#' @param x,q vector of distances
#' @param a,b0,b1,b2,b3,shape,scale,s2 parameters used in the respective
#'  distributions.
#' @param const (optional) scalar normalizing constant for distributions that are
#'  numerically integrated using \code{integrate}, namely. Providing a \code{const} is
#'  not necessary but will improve the speed of calculation under certain
#'  conditions.
#' @return vector of probability densities or cumulative probabilities
#' @export
#'
dmb <- function(x, a) sqrt(2/pi) * x^2/as.vector(a)^3 * exp(-x^2/(2 * as.vector(a)^2))

#' @rdname Distributions
#' @export
pmb <- function(q, a)
  pracma::erf(q/(sqrt(2) * as.vector(a))) - sqrt(2/pi) * (q * exp(-q^2/(2 * as.vector(a)^2)))/as.vector(a)

#' @rdname Distributions
#' @export
#'
dxep1 <- function(x, b1){
  ans <- as.vector(b1)^2 * x * exp(as.vector(b1) * x)
  ans[x <= 0]  <- 0
  ans
}

#' @rdname Distributions
#' @export
pxep1 <- function(q, b1){
  totl <- max(length(q), length(b1))
  b1 <- rep(b1, length.out = totl)
  q <- rep(q, length.out = totl)
  exp(b1 * q) * (b1 * q - 1) + 1
}

#' @rdname Distributions
#' @export
pxep02 <- function(q, b0, b2){
  const <- 1/((-b2)^(-b0/2) * (-gamma(b0/2 + 1))/(2 * b2))
  ans <- const * q^b0 * (-b2 * q^2)^(-b0/2) *
    (incGamma(1/2*(b0 + 2), -b2 * q^2) - gamma(b0/2 + 1))/(2*b2)
  ans[q <= 0] <- 0
  ans
}

#' @rdname Distributions
#' @export
dxep02 <- function(x, b0, b2){
  const <- 1/((-b2)^(-b0/2) * (-gamma(b0/2 + 1))/(2 * b2))
  ans <- const * x * exp(b0 * log(x) + b2 * x^2)
  ans[x <= 0] <- 0
  ans
}


#' @rdname Distributions
#' @export
#'
dxep12 <- function(x, b1, b2){
  const <- 4 * b2^2/(sqrt(-b2 * pi) * exp(-b1^2/(4 * b2)) * b1 *
    (pracma::erf(0.5 * b1/sqrt(-b2)) + 1) - 2 * b2)
  ans <- numeric(length(x))
  ans <- const * x * exp(b1*x + b2*x^2)
  ans[x <= 0] <- 0
  ans
}

#' @rdname Distributions
#' @export
pxep12 <- function(x, b1, b2){
  const <- 1/(sqrt(-b2 * pi) * exp(-b1^2/(4 * b2)) * b1 *
    (pracma::erf(0.5 * b1/sqrt(-b2)) + 1) - 2 * b2)
  ans <- const * (sqrt(-b2 * pi) * exp(-b1^2/(4 * b2)) * b1 *
    (pracma::erf((-b2 * x - 0.5*b1)/sqrt(-b2)) + pracma::erf(0.5 * b1/sqrt(-b2))) +
    2 * b2 * (exp(b1 * x + b2 * x^2) - 1))
  ans[ans <= 0] <- 0
  ans
}

#' @rdname Distributions
#' @export
#'
dxep123 <- function(x, b1, b2, b3, const = NULL){
  if (is.null(const)) const <- 1/integrate(
    f = function(x) x * exp(b1 * x + b2 * x^2 + b3 * x^3),
    lower = 0, upper = Inf, rel.tol = .Machine$double.eps^0.5)$val
  ans <- const * x * exp(b1 * x + b2 * x^2 + b3 * x^3)
  ans[x <= 0] <- 0
  ans
}

#' @rdname Distributions
#' @export
pxep123 <- function(x, b1, b2, b3, const = NULL){
  if (is.null(const)) const <- tryCatch(1/integrate(
    f = function(x) x * exp(b1 * x + b2 * x^2 + b3 * x^3),
      lower = 0, upper = Inf, rel.tol = .Machine$double.eps^0.5)$val,
    error = function(e) NA
  )
  if (is.na(const)) return(rep(NA, length(x)))
  ans <- numeric(length(x))
  for (xi in which(x > 0)){
    ans[xi] <- tryCatch(const * integrate(
      f = function(r) r * exp(b1*r + b2*r^2 + b3*r^3), lower = 0, upper = x[xi]
      )$val,
      error = function(e) NA
    )
  }
  ans
}

#' @rdname Distributions
#' @export
#'
dxepi0 <- function(x, shape, scale){
  ans <- numeric(length(x))
  xx <- x
  if (any(x <= 0)) xx[x <= 0] <- 1e-4
  ans <- invgamma::dinvgamma(xx,
    shape = shape,
    scale =  1/scale)
  ans[x <= 0] <- 0
  ans
}

#' @rdname Distributions
#' @export
pxepi0 <- function(x, shape, scale){
  ans <- numeric(length(x))
  xx <- x
  if (any(x <= 0)) xx[x <= 0] <- 1e-4
  ans <- invgamma::pinvgamma(xx,
    shape = shape,
    scale = 1/scale)
  ans[x <= 0] <- 0
  ans
}

#' @rdname Distributions
#' @export
#'
dxep0123 <- function(x, b0, b1, b2, b3, const = NULL){
  if (is.null(const)) const <- 1/integrate(
    f = function(x) x * exp(b0 * log(x) + b1 * x + b2 * x^2 + b3 * x^3),
    lower = 0, upper = Inf, rel.tol = .Machine$double.eps^0.5)$val
  ans <- const * x * exp(b0*log(x) + b1*x + b2*x^2 + b3*x^3)
  ans[x <= 0] <- 0
  ans
}

#' @rdname Distributions
#' @export
pxep0123 <- function(x, b0, b1, b2, b3, const = NULL){
  # works for vector x but not beta parameters
  if (is.null(const)) const <- tryCatch(1/integrate(
    f = function(x) x * exp(b0 * log(x) + b1 * x + b2 * x^2 + b3 * x^3),
    lower = 0, upper = Inf, rel.tol = .Machine$double.eps^0.5)$val,
    error = function(e) NA)
  if(is.na(const)) return(rep(NA, length(x)))
  ans <- numeric(length(x))
  for (xi in which(x > 0)){
    ans[xi] <- const * integrate(
      f = function(r) r * exp(b0*log(r) + b1*r + b2*r^2 + b3*r^3),
      lower = 0, upper = x[xi]
    )$val
  }
  ans[x <= 0] <- 0
  ans
}

#' @rdname Distributions
#' @export
#'
dxep012 <- function(x, b0, b1, b2, const = NULL){
  if (is.null(const)) const <- 1/integrate(
    f = function(x) x * exp(b0 * log(x) + b1 * x + b2 * x^2),
    lower = 0, upper = Inf, rel.tol = .Machine$double.eps^0.5)$val
  ans <- const * x * exp(b0*log(x) + b1*x + b2*x^2)
  ans[x <= 0] <- 0
  ans
}

#' @rdname Distributions
#' @export
pxep012 <- function(x, b0, b1, b2, const = NULL){
  # works for vector x but not beta parameters
  if (is.null(const)) const <- tryCatch(1/integrate(
    f = function(x) x * exp(b0 * log(x) + b1 * x + b2 * x^2),
    lower = 0, upper = Inf, rel.tol = .Machine$double.eps^0.5)$val,
    error = function(e) NA)
  if(is.na(const)) return(rep(NA, length(x)))
  ans <- numeric(length(x))
  for (xi in which(x > 0)){
    ans[xi] <- const * integrate(
      f = function(r) r * exp(b0*log(r) + b1*r + b2*r^2),
      lower = 0, upper = x[xi]
    )$val
  }
  ans[x <= 0] <- 0
  ans
}

#' @rdname Distributions
#' @export
#'
dxep2 <- function(x, s2) x/as.vector(s2) * exp(-x^2/(2 * as.vector(s2)))

#' @rdname Distributions
#' @export
pxep2 <- function(x, s2)  1 - exp(-x^2/(2 * as.vector(s2)))

#' @rdname Distributions
#' @export
#'
dxep0 <- function(x, a){
  ans <- numeric(length(x))
  ans[x > 1] <- as.vector(a)/x^(as.vector(a) + 1)
  ans
}

#' @rdname Distributions
#' @export
pxep0 <- function(x, a){
  ans <- numeric(length(x))
  ans[x > 1] <- 1 - (1/x)^as.vector(a)
  ans
}

#' Simple Utility Function Used in Optimizing the GLM
#'
#' A simple utility function that is used in fitting a GLM, creating a matrix
#'  of "x" values for use in the polynomial part of a xep-type model.
#'
#' @param r vector of distances (>=0)
#' @param distr name of the distribution
#' @return array with \code{length(r)} rows and p columns, where p is the number 
#'  of parameters in the glm (including the intercept). The first column is all 
#'  1s, and the remaining columns are functions of r, specifically, log(r), r, 
#'  r^2, r^3, or 1/r, depending on what the distribution requires.
#' @export
rmat <- function(r, distr){
  switch(distr,
    xep01 = cbind(1, log(r), r),
    xep012 = cbind(1, log(r), r, r^2),
    xep02 = cbind(1, log(r), r^2),
    lognormal = cbind(1, log(r), log(r)^2),
    xep1 = cbind(1, r),
    xep12 = cbind(1, r, r^2),
    xep123 = cbind(1, r, r^2, r^3),
    xepi0 = cbind(1, 1/r, log(r)),
    xep0123 = cbind(1, log(r), r, r^2, r^3),
    xep2 = cbind(1, r^2),
    MaxwellBoltzmann = cbind(1, r^2),
    constant = matrix(1, nrow = length(r)) ,
    tnormal = cbind(1, r, r^2),
    exponential = cbind(1, r),
    xep0 = cbind(1, log(r)),
    chisq = cbind(1, log(r)),
    inverse_gaussian = cbind(1, 1/r, r)
  )
}

#' Utility Function for Constructing Offsets for GLMs
#'
#' This is a simple utility function for calculating offsets when exposure
#'  is assumed to be 100% at a given distance \code{r}. This is useful for
#'  calculating fitted distributions (\code{PDF} and \code{CDF}) but cannot be
#'  used in the fitting of the distributions themselves because it does not
#'  account for incomplete search coverages at given distances.
#'
#' @param r vector of distances
#' @param distr name of the distribution to calculate the offset for
#' @return A vector of offset values to use with distances \code{r} when fitting
#'  the \code{distr} model.
#' @export
off <- function(r, distr){
  if (natural[distr]) return(log(r))
  if (distr %in% c("tnormal", "exponential")) return(0)
  if (distr == "MaxwellBoltzmann") return(2 * log(r))
  if (distr == "chisq") return(log(r) - r/2)
  if (distr == "inverse_gaussian") return(-1.5 * log(r))
}

#' Check Whether GLM Coefficients Give Proper Distribution
#' 
#' In order for a fitted GLM to convert to a proper distance distribution, its
#' integral from 0 to Inf must be finite. As a rule, when the leading coefficient 
#' is positive, the integral diverges as the upper bound of integration approaches
#' infinity, and \code{cofOK} would return \code{FALSE}. Likewise, in some cases,
#' the GLM coefficients yield an integral that diverges as the lower bound 
#' approaches 0, in which case \code{cofOK} returns \code{FALSE} as well. 
#' \code{cofOK0} and \code{cofOKInf} check the left and right tails of the 
#' candidate distribution, repectively, for convergence.
#' 
#' @param cof vector or matrix of named glm parameters (with \code{"r"} as the
#'  distance variable)
#' @param distr name of the distribution
#' @return boolean vector (or scalar)
#' @export
cofOK <- function(cof, distr){
  if (missing("distr")) stop("cofOK requires distr")
  if (is.vector(cof))
    cof <- matrix(cof, nrow = 1, dimnames = list(NULL, names(cof)))
  output <- rep(TRUE, nrow(cof))
  lim <- constraints[[distr]]
  for (ci in rownames(lim)){
    output[cof[, ci] <= lim[ci, "lower"] | cof[, ci] >= lim[ci, "upper"]] <- FALSE
  }
  output
}

#' @rdname cofOK
#' @export
cofOK0 <- function(cof, distr){
  if (missing("distr")) stop("cofOK0 requires distr")
  if (is.vector(cof))
    cof <- matrix(cof, nrow = 1, dimnames = list(NULL, names(cof)))
  output <- rep(TRUE, nrow(cof))
  lim <- constraints[[distr]]
  ci <- grep("log(r)", rownames(lim))
  if (length(ci) > 0) output[cof[, ci] <= lim[ci, "lower"]] <- FALSE
  ci <- grep("I(1/r)", rownames(lim))
  if (length(ci) > 0) output[cof[, ci] <= lim[ci, "lower"]] <- FALSE
  output
}

#' @rdname cofOK
#' @export
cofOKInf <- function(cof, distr){
  if (missing("distr")) stop("cofOKInf requires distr")
  if (is.vector(cof))
    cof <- matrix(cof, nrow = 1, dimnames = list(NULL, names(cof)))
  output <- rep(TRUE, nrow(cof))
  lim <- constraints[[distr]]
  ci <- rownames(lim)[nrow(lim)]
  for (ci in rownames(lim)){
    output[cof[, ci] >= lim[ci, "upper"]] <- FALSE
  }
  output
}

#' Check Parameter Value Validity the Distribution
#' 
#' Performs a ouick check on whether the parameters given in \code{parms} are
#'  valid for the \code{distr}. 
#' 
#' @param parms vector or matrix of named glm parameters (with \code{"r"} as the
#'  distance variable)
#' @param distr name of the distribution
#' @return vector (or scalar) of 0s, 1s, and 2s to indicate whether the parameters
#'  are non-extensible (i.e., flat-out bogus), valid for the distribution, or
#'  valid for the distribution and give finite point densities.
#' @export
parOK <- function(parms, distr){
  if (is.vector(parms))
    parms <- matrix(parms, nrow = 1, dimnames = list(NULL, names(parms)))
  output <- rep(TRUE, nrow(parms))
  lim <- constraints_par[[distr]]
  for (ci in rownames(lim)){
    output[parms[, ci] <= lim[ci, "lower"] | parms[, ci] >= lim[ci, "upper"]] <- FALSE
  }
  
  output
}

#' Remove Particular Names from a Longer List 
#' 
#' Removes specific values (\code{what}) from a longer vector of values (\code{from}).
#' By default, \code{from = mod_standard}, and the intent is to simplify the
#' subsetting of \code{ddArray} objects created with the default standard models.
#' For example, \code{dmod2 <- dmod[exclude("lognormal")]} would subset the list
#' of models in \code{mod_standard} to exclude \code{"lognomal"}. The default can
#' be overridden by providing a specific vector for \code{from} (for example,
#' \code{dmod[exclude("lognormal", from = names(dmod)])}).
#' 
#' @param what vector of distribution names to exclude
#' @param from vector of distribution names to be excluded from
#' @return vector of names from "from" after excluding "what"
#' @export
exclude <- function(what, from = mod_standard) setdiff(from, what)
 

#' Convert Distribution Name + Parameters to \code{ddSim} Object
#'
#' Utility function to format a distribution name and a vector of its parameter values
#'  to a \code{\link{ddSim}} object for use in the p/d/r/q family of functions
#' 
#' @param distr character string giving the name of one of the models fit by
#'  \code{\link{ddFit}}
#' @param parms vector of parameters for the \code{distr}, or, alternatively, an array
#'  of parameter sets. Parameterization may follow  either the GLM format or the 
#'  distribution format. For example, the xep01  model (gamma distribution) has
#'  GLM parameters for \code{log(r)} and \code{r}, which are the coefficients of 
#'  the polynomial in the xep01 format (namely, x * exp(b0*log(r) + b1*r)), or 
#'  the gamma distribution parameters, \code{shape} and \code{rate}. The 
#'  elements of parameter vector must be named. For example, 
#'  \code{parms = c(log(r) = -0.373, r = -0.0147)} for the GLM format for an 
#'  xep01 model or, equivalently, \code{parms = c(shape = 1.63, rate = 0.0147)} 
#'  for the distribution format. If  both formats are given, the GLM parameters 
#'  are used and the distribution parameters ignored.
#' @return a ddSim object with \code{srad = NA}
#' @export
mpp2ddSim <- function(distr, parms){
  # error-checking:
  if (!is.character(distr) || !length(distr) == 1)
    stop("pdd: distr must be a single name of a distribution included in mod_all")
  if (!distr %in% mod_all){
    if (distr %in% names(alt_names)){
      distr <- alt_names[distr]
    } else {
      stop("distr must be name of a distribution in mod_all")
    }
  }
  if (is.null(parms)) stop("pdd: parameters must be provided as 'parms'",
      "when 'distr' is a character string giving the name of the distribution")
  if (!is.numeric(parms) || (!is.vector(parms) & !is.matrix(parms))) 
    stop("'parms' must be a numeric vector or matrix of parameters")
  if (is.null(names(parms))) 
    stop("parms must be a vector of named parameters for distr")
  if (is.vector(parms)){
   parms <- matrix(parms, nrow = 1, dimnames = list(NULL, names(parms)))
  }
  cnm <- cof_name[[distr]]
  pnm <- parm_name[[distr]]
  out <- matrix(ncol = length(cnm) + length(pnm) + 1, nrow = nrow(parms),
    dimnames = list(NULL, c(cnm, pnm, "extensible")))
  if (all(cnm[-1] %in% colnames(parms))){ # "(Intercept)" not necessary
    out[, cnm[-1]] <- parms[, cnm[-1]]
    if (cnm[1] %in% colnames(parms)) out[, 1] <- parms[, cnm[1]]
    out[, pnm] <- cof2parms(out[ , cnm], distr)
    out[, "extensible"] <- cofOK(out[, pnm], distr)
  } else if (all(pnm %in% colnames(parms))){
    out[, pnm] <- parms[, pnm]
    out[, cnm] <- parms2cof(parms, distr)
    out[, "extensible"] <- parOK(out[ , pnm], distr)
  } else {
    stop("incomplete parameter set in 'parms'")
  }

  attr(out, "distr") <- distr
  attr(out, "srad") <- NA
  class(out) <- "ddSim"
  out
}

parms2cof <- function(x, ...) UseMethod("parms2cof", x)
#' @rdname cof2parms
#' @export
parms2cof.matrix <- function(x, distr, ...){
  ans <- suppressWarnings(switch(distr,
    xep01 = {
      cof <- cbind(x[, "shape"] - 2, -x[, "rate"])
      colnames(cof) <- c("log(r)", "r")
      cof
    },
    lognormal = {
      log2r <- unname(-0.5/(x[, "sdlog"])^2)
      cof <- cbind(unname(-2 * (log2r*x[, "meanlog"] + 1)), log2r)
      colnames(cof) <- c("log(r)", "I(log(r)^2)")
      cof
    },
    xep1 = {
      cof <- unname(x[, "b1", drop = FALSE])
      colnames(cof) <- "r"
      cof
    },
    xep12 = {
      cof <- x[, c("b1", "b2"), drop = FALSE]
      colnames(cof) <- c("r", "I(r^2)")
      cof
    },
    xep02 = {
      cof <- x[, c("b0", "b2"), drop = FALSE]
      colnames(cof) <- c("log(r)", "I(r^2)")
      cof
    },
    xep123 ={
      cof <- x[, c("b1", "b2", "b3"), drop = FALSE]
      colnames(cof) <- c("r", "I(r^2)", "I(r^3)")
      cof
    },
    xepi0 = {
      cof <- -x[, c("shape", "scale")]
      colnames(cof) <- c("log(r)", "I(1/r)")
      cof[, "log(r)"] <- cof[, "log(r)"] - 2
      cof
    },
    xep0123 = {
      cof <- x[, c("b0", "b1", "b2", "b3"), drop = FALSE]
      colnames(cof) <- c("log(r)", "r", "I(r^2)", "I(r^3)")
      cof
    },
    xep012 = {
      cof <- x[, c("b0", "b1", "b2"), drop = FALSE]
      colnames(cof) <- c("log(r)", "r", "I(r^2)")
      cof
    },
    xep2 = {
      cof <- -0.5/x[, "s2", drop = FALSE]
      colnames(cof) <- "I(r^2)"
      cof
    },
    MaxwellBoltzmann = {
      cof <- -0.5 /x[, "a", drop = FALSE]^2
      colnames(cof) <- "I(r^2)"
      cof
    },
    xep0 = {
      cof <- -x[, "a", drop = FALSE] - 2
      colnames(cof) <- "log(r)"
      cof
    },
    chisq = {
      cof <- (x[, "df", drop = FALSE] - 4)/2
      colnames(cof) <- "log(r)"
      cof
    },
    inverse_gaussian = {
      cof <- cbind(-0.5/x[, "dispersion"], -0.5/(x[, "dispersion"] * x[, "mean"]))
      colnames(cof) <- c("I(1/r)", "r")
      cof
    },
    exponential = {
      cof <- -x[, "rate", drop = FALSE]
      colnames(cof) <- "r"
      cof
    },
    tnormal = {
      r2 <- -0.5/x[, "sd"]^2
      cof <- cbind(-2 * x[, "mean"] * r2, r2)
      colnames(cof) <- c("r", "I(r^2)")
      cof
    },
    constant = c(b1 = NA)
  ))
  if (is.null(ans)) stop("distr must be a distribution name")
  ans <- cbind(0, ans)
  colnames(ans)[1] <- "(Intercept)"
  if (is.matrix(ans) && nrow(ans) == 1) {
    nm <- colnames(ans)
    ans <- as.vector(ans)
    names(ans) <- nm
  }
  return(ans)
}

parms2cof.numeric <- function(x, distr, ...){
  output <- parms2cof(matrix(x, nrow = 1, dimnames = list(NULL, names(x))), distr)
  return(output)
}
