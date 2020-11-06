#' Incomplete Gamma Function
#' @param x non-negative numeric vector
#' @param a positive numeric vector
#' @param lower boolean for calculating lower or upper incomplete gamma function
#' @details The upper incomplete gamma function, following Wolfram Alpha, namely,
#'  incGam(x, a) = Gamma = integral(exp(-t) * t^(a - 1) dt from x to Inf),
#'  calculated using pgamma.
#' @return scalar or vector of length = max(length(x), length(a)), with values
#'  of the shorter recycled to match the length of the longer a la pnorm etc.
#' @export
incGamma <- function(a, x, lower = FALSE) pgamma(x, a, lower = lower) * gamma(a)

#' Error function
#' @description The error function is closely related to the standard normal CDF and arises frequently in probability calculations
#' @param x numeric scalar or array
#' @details \eqn{erf(x) = 2/\sqrt\pi * \int_0^x {exp(-x\^2} dx}
#' @return erf(x) with same dimensions as x
#' @export
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

#' Calculate Exposure for Corners of Square Plots

#' Probability Distributions for Carcasses vs. Distance
#' @param x, q vector of quantiles
#' @param a, b0, b1, b2, b3, shape, scale, s2 parameters used in the respective
#'  distributions. See Detailin the Maxwell-Boltzmann distribution
#' @section The distributions:
#'  Except where otherwise noted, the parameterizations match those given in
#'  Wikipedia (and the user guide and vignettes)...section under construction.
#'  \describe{
#'    \item{Maxwell-Boltzmann}{The pdf of the Maxwell-Boltzmann is defined as
#'      f(x; a) = sqrt(2/pi)/a^3 * x^2 * exp(-x^2/(2*a^2)), and the offset used in
#'      the glm for carcass distances is \code{log(exposure) + log(r)}.}
#'    \item{log-linear}{The pdf of the log-linear is f(x; b1) = b1^2 * x * exp(b1 * x),
#'      where b1 is the}
#'  }
#' @return \code{dmb} gives the density, \code{pmb} gives the distribution function.
#'  The length of the result is the maximum of the lengths of the numerical arguments
#'  for the other functions. The numerical arguments are recycled to the length of
#'  the result.
#' @export
#'
dmb <- function(x, a) sqrt(2/pi) * x^2/as.vector(a)^3 * exp(-x^2/(2 * as.vector(a)^2))

#' @rdname dmb
#' @export
pmb <- function(q, a)
  erf(q/(sqrt(2) * as.vector(a))) - sqrt(2/pi) * (q * exp(-q^2/(2 * as.vector(a)^2)))/as.vector(a)

#' Integral in xep12 distribution
#' @param x non-negative numeric vector
#' @param a positive numeric vector
#' @details The upper incomplete gamma function, following Wolfram Alpha, namely,
#'  incGam(x, a) = Gamma = integral(exp(-t) * t^(a - 1) dt from x to Inf),
#'  calculated using pgamma.
#' @return scalar or vector of length = max(length(x), length(a)), with values
#'  of the shorter recycled to match the length of the longer a la pnorm etc.
#' @export
lQint <- function(x, b1, b2){
  Re(sqrt(pi) * exp(-b1^2/(4*b2)) * b1 *
   (pracma::erfi(0.5*b1/sqrt(b2 + 0i)) - pracma::erfi((x*b2 + 0.5*b1)/sqrt(b2 + 0i)))/
   (4 * (b2 + 0i)^1.5) + (exp(x^2 * b2 + x*b1) - 1)/(2*b2))
}

#' @param x, q vector of quantiles
#' @param b1 \code{b1} parameter in the log-linear distribution
#' @details The pdf of the log-linear distribution is defined as
#'  f(x; b1) = b1^2 * x * exp(b1 * x)
#' @return \code{dlogL} gives the density, \code{plogL} gives the distribution function.
#'  The length of the result is the maximum of the lengths of the numerical arguments
#'  for the other functions. The numerical arguments are recycled to the length of
#'  the result.
#' @rdname dmb
#' @export
#'
dxep1 <- function(x, b1){
  ans <- as.vector(b1)^2 * x * exp(as.vector(b1) * x)
  ans[x <= 0]  <- 0
  ans
}

#' @rdname dmb
#' @export
pxep1 <- function(q, b1){
  totl <- max(length(q), length(b1))
  b1 <- rep(b1, length.out = totl)
  q <- rep(q, length.out = totl)
  exp(b1 * q) * (b1 * q - 1) + 1
}

#' @rdname dmb
#' @export
pxep02 <- function(q, b0, b2){
  const <- 1/((-b2)^(-b0/2) * (-gamma(b0/2 + 1))/(2 * b2))
  ans <- const * q^b0 * (-b2 * q^2)^(-b0/2) *
    (incGamma(1/2*(b0 + 2), -b2 * q^2) - gamma(b0/2 + 1))/(2*b2)
  ans[q <= 0] <- 0
  ans
}

#' @rdname dmb
#' @export
dxep02 <- function(x, b0, b2){
  const <- 1/((-b2)^(-b0/2) * (-gamma(b0/2 + 1))/(2 * b2))
  ans <- const * x * exp(b0 * log(x) + b2 * x^2)
  ans[x <= 0] <- 0
  ans
}


#' @param x, q vector of quantiles
#' @param b1, b2 \code{b1, b2} parameters in the log-quadratic distribution
#' @details The pdf of the log-quadratic distribution is defined as
#'  f(x; b1, b2) = const * x * exp(b1 * x + b2 * x^2). The normalization constant
#'  (\code{const}) is a complicated expression but can be expressed using only
#'  elementary functions and \code{erf}.
#' @return \code{dlogQ} gives the density, \code{plogQ} gives the distribution function.
#'  The length of the result is the maximum of the lengths of the numerical arguments
#'  for the other functions. The numerical arguments are recycled to the length of
#'  the result.
#' @rdname dmb
#' @export
#'
dxep12 <- function(x, b1, b2){
  const <- 4 * b2^2/(sqrt(-b2 * pi) * exp(-b1^2/(4 * b2)) * b1 *
    (erf(0.5 * b1/sqrt(-b2)) + 1) - 2 * b2)
  ans <- numeric(length(x))
  ans <- const * x * exp(b1*x + b2*x^2)
  ans[x <= 0] <- 0
  ans
}

#' @rdname dmb
#' @export
pxep12 <- function(x, b1, b2){
  const <- 1/(sqrt(-b2 * pi) * exp(-b1^2/(4 * b2)) * b1 *
    (erf(0.5 * b1/sqrt(-b2)) + 1) - 2 * b2)
  ans <- const * (sqrt(-b2 * pi) * exp(-b1^2/(4 * b2)) * b1 *
    (erf((-b2 * x - 0.5*b1)/sqrt(-b2)) + erf(0.5 * b1/sqrt(-b2))) +
    2 * b2 * (exp(b1 * x + b2 * x^2) - 1))
  ans[ans <= 0] <- 0
  ans
}

#' @param x, q vector of quantiles
#' @param b1, b2 \code{b1, b2} parameters in the log-quadratic distribution
#' @details The pdf of the log-cubic distribution is defined as
#'  f(x; b1, b2, b3) = const * x * exp(b1 * x + b2 * x^2 + b3 * x^3). The
#'  normalization constant (\code{const}) is a complicated expression that
#'  cannot be expressed using only elementary functions and readily approximated
#'  functions in R. Calculation of constant is cumbersome, so it may be calculated
#'  externally and provided. NOTE: The log-cubic distribution is calculated via
#'  numerical integration. It is vectorized for x but not for b1, b2, b3.
#' @return \code{dlogC} gives the density, \code{plogC} gives the distribution
#'  function. The length of the result is the maximum of the lengths of the
#'  numerical arguments for the other functions. The numerical arguments are
#'  recycled to the length of the result.
#' @rdname dmb
#' @export
#'
dxep123 <- function(x, b1, b2, b3, const = NULL){
  if (is.null(const)) const <- 1/integrate(
    f = function(x) x * exp(b1 * x + b2 * x^2 + b3 * x^3),
    lower = 0, upper = Inf, rel.tol = .Machine$double.eps^0.5)$val
  ans <- const * x * exp(b1 * x + b2 * x^2 + b3 * x^3)
  ans[x <= 0] <- 0
}

#' @rdname dmb
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

#' @param x, q vector of quantiles
#' @param shape, scale \code{a > 0 (shape), b > 0 (scale)} parameters in the
#'  inverse gamma distribution
#' @details The pdf of the inverse gamma distribution is defined as
#'  f(x; a, b) = b^a/gamma(a) * x^(-a - 1) * exp(-b/x). Calculation is via
#'  \code{invgamma::dinvgamma} and \code{invgamma::pinvgamma} with modification
#'  to return 0 instead of NaN for x <= 0. NOTE: The parameters
#' @return \code{digam} gives the density, \code{pigam} gives the distribution function.
#'  The length of the result is the maximum of the lengths of the numerical arguments
#'  for the otherp functions. The numerical arguments are recycled to the length of
#'  the result.
#' @rdname dmb
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

#' @rdname dmb
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

#' @param x, q vector of quantiles
#' @param b0, b1, b2, b3, const parameters in the paranormal-gamma distribution
#' @details The pdf of the paranormal-gamma distribution is defined as
#'  f(x; b0, b1, b2, b3) = const * x * exp(b0 * log(x) + b1 * x + b2 * x^2 + b3 * x^3).
#'  Calculation is cumbersome and done via numerical integration and is not
#'  vectorized yet.
#' @return \code{dxep0123} gives the density, \code{ppng} gives the distribution function.
#'  The length of the result is the maximum of the lengths of the numerical arguments
#'  for the other functions. The numerical arguments are recycled to the length of
#'  the result.
#' @rdname dmb
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

#' @rdname dmb
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

#' @param x, q vector of quantiles
#' @param b0, b1, b2, const parameters in the xep012 distribution
#' @details The pdf of the xep012 distribution is defined as
#'  f(x; b0, b1, b2) = const * x * exp(b0 * log(x) + b1 * x + b2 * x^2).
#'  Calculation is cumbersome and done via numerical integration and is not
#'  vectorized (yet).
#' @return \code{dxep012} gives the density, \code{pxep012} gives the distribution function.
#'  The length of the result is the maximum of the lengths of the numerical arguments
#'  for the other functions. The numerical arguments are recycled to the length of
#'  the result.
#' @rdname dmb
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

#' @rdname dmb
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

#' @param x, q vector of quantiles
#' @param s2 \code{s2} parameter in the xep2 distribution
#' @details The pdf of the xep2 distribution is defined as
#'  f(x; s2) = x/s2 * exp(-x^2/(2 * s2))
#' @return \code{dRay} gives the density, \code{pRay} gives the distribution
#'  function. The length of the result is the maximum of the lengths of the
#'  numerical arguments for the other functions. The numerical arguments are
#'  recycled to the length of the result.
#' @rdname dmb
#' @export
#'
dxep2 <- function(x, s2) x/as.vector(s2) * exp(-x^2/(2 * as.vector(s2)))

#' @rdname dmb
#' @export
pxep2 <- function(x, s2)  1 - exp(-x^2/(2 * as.vector(s2)))

#' @param x, q vector of quantiles
#' @param a \code{a} parameter in the xep0 distribution with \code{scale = x_m}
#'  fixed at 1.
#' @details The pdf of the Pareto distribution with scale = 1 is defined as
#'  f(x; a) = a/x^(a + 1)
#' @return \code{dPare1} gives the density, \code{pPare1} gives the distribution function.
#'  The length of the result is the maximum of the lengths of the numerical arguments
#'  for the other functions. The numerical arguments are recycled to the length of
#'  the result.
#' @rdname dmb
#' @export
#'
dxep0 <- function(x, a){
  ans <- numeric(length(x))
  ans[x > 1] <- as.vector(a)/x^(as.vector(a) + 1)
  ans
}

#' @rdname dmb
#' @export
pxep0 <- function(x, a){
  ans <- numeric(length(x))
  ans[x > 1] <- 1 - (1/x)^as.vector(a)
  ans
}

#' simple utility function used in optimizing the truncated glm
#' @param r vector of distances (>=0)
#' @param distr name of the distribution
#' @return array with length(r) rows and p columns, where p is the number of
#'  parameters in the glm (including the intercept). The first column is all 1s,
#'  and the remaining columns are functions of r, specifically, log(r), r, r^2,
#'  r^3, or 1/r, depending on what the distribution requires.
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

#' utility function for calculating offset for integrals for glm to dd
#'
#' This is a simple utility function for calculating offsets when exposure
#'  is assumed to be 100% at a given distance \code{r}. This is useful for
#'  calculating fitted distributions (\code{PDF} and \code{CDF}) but cannot be
#'  used in the fitting of the distributions themselves because it does not
#'  account for incomplete search coverages at given distances.
#'
#' @param r vector of distances
#' @param distr name of the distribution to calculate the offset for
#' @export
off <- function(r, distr){
  if (natural[distr]) return(log(r))
  if (distr %in% c("tnormal", "exponential")) return(0)
  if (distr == "MaxwellBoltzmann") return(2 * log(r))
  if (distr == "chisq") return(log(r) - r/2)
  if (distr == "inverse_gaussian") return(-1.5 * log(r))
}

#' reset par to default as in par_default for the \code{dwp} package
#' @export
par_reset <- function() do.call(par, par_default)

#' check whether glm coefficients give proper distribution
#' @param cof vector or matrix of named glm parameters (with \code{"r"} as the
#'  distance variable.
#' @param distr name of the distribution
#' @return boolean vector (or scalar)
#' @export
cofOK <- function(cof, distr){
  if (is.vector(cof))
    cof <- matrix(cof, nrow = 1, dimnames = list(NULL, names(cof)))
  output <- rep(TRUE, nrow(cof))
  lim <- constraints[[distr]]
  for (ci in rownames(lim)){
    output[cof[, ci] <= lim[ci, "lower"] | cof[, ci] >= lim[ci, "upper"]] <- FALSE
  }
  output
}

#' remove particular distribution names from a longer list (for subsetting \code{ddArray})
#' @param what vector of distribution names to exclude
#' @param from vector of distribution names to be excluded from
#' @export
exclude <- function(what, from = mod_name) setdiff(from, what)
