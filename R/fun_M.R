#' @title Calculate posterior distribution of M and extract statistics (M* and CI)
#'
#' @description Calculation of the posterior distribution of total mortality
#'  (\code{M}) given the carcass count, overall detection probability (\code{g}),
#'  and prior distribtion; calculation of summary statistics from the
#'  posterior distribution of \code{M}, including \code{M*} and credibility
#'  intervals.
#'
#' @details  The functions \code{postM} and \code{postM.ab} return the posterior
#'  distributions of \eqn{M|(X, g)} and \eqn{M|(X, Ba, Bb)}, respectively, where
#'  \code{Ba} and \code{Bb} are beta distribution parameters for the estimated
#'  detection probability. \code{postM} and \code{postM.ab} include options to
#'  to specify a prior distribution for \eqn{M} and a limit for truncating the
#'  prior to disregard implausibly large values of \eqn{M} and make the
#'  calculations tractable in certain cases where they otherwise might not be.
#'  Use \code{postM} when \eqn{g} is fixed and known; otherwise, use \code{postM.ab}
#'  when uncertainty in \eqn{g} is characterized in a beta distribution with
#'  parameters \eqn{Ba} and \eqn{Bb}. The non-informative, integrated reference
#'  prior for binomial random variables is the default (\code{prior = "IbinRef"}).
#'  Other options include "binRef", "IbetabinRef", and "betabinRef", which are
#'  the non-integrated and integrated forms of the binomial and betabinomial
#'  reference priors (Berger et al., 2012). For \eqn{X > 2}, the integrated and
#'  non-integrated reference priors give virtually identical posteriors. However,
#'  the non-integrated priors assign infinite weight to \eqn{m = 0} and return a
#'  posterior of \eqn{Pr(M = 0| X = 0, \hat{g}) = 1}, implying absolute certainty
#'  that the total number of fatalities was 0 if no carcasses were observed. In
#'  addition, a uniform prior may be specified by prior = "uniform". Alternatively,
#'  a custom prior may be given as a 2-dimensional array with columns for \eqn{m}
#'  and \eqn{Pr(M = m)}, respectively. The first column (\code{m}) must be
#'  sequential integers starting at \eqn{m = 0}. The second column gives the
#'  probabilities associated with \eqn{m}, which must be non-negative and sum to 1.
#'  The named priors (\code{"IbinRef"}, \code{"binRef"}, \code{"IbetabinRef"},
#'  and \code{"betabinRef"}) are functions of \eqn{m} and defined on \eqn{m=0,1,2,...}
#'  without upper bound. However, the posteriors can only be calculated for a
#'  finite number of \eqn{m}'s up to a maximum of \code{mmax}, which is set by
#'  default to the smallest value of \eqn{m} such that
#'  \eqn{Pr(X \leq x | m, \hat{g}) < 0.0001}, where \eqn{x} is the observed
#'  carcass count, or, alternatively, \code{mmax} may be specified by the user.
#'
#' @param x carcass count
#'
#' @param g overall carcass detection probability
#'
#' @param mmax cutoff for prior of M (large max requires large computing resources
#'  but does not help in the estimation)
#'
#' @param Ba,Bb parameters for beta distribution characterizing estimated \eqn{g}
#'
#' @param prior prior distribution of \eqn{M}
#'
#' @param pMgX posterior distribution of \eqn{M}
#'
#' @param crlev,alpha credibility level (\eqn{1-\alpha}) and its complement (\eqn{\alpha})
#'
#' @return The functions \code{postM} and \code{postM.ab} return the posterior
#'  distributions of \eqn{M | (X, g)} and \eqn{M | (X, Ba, Bb)}, respectively.
#'  The functions \code{calcMstar} and \code{MCI} return \eqn{M^*} value and
#'  credibility interval for the given posterior distribution, \code{pMgX}
#'  (which may be the return value of \code{postM} or \code{postM.ab}) and
#'  \eqn{\alpha} value or credibility level.
#'
#' @export
#'
postM <- function(x, g, prior = "IbinRef", mmax = NA){
 ### repaired version of postM from eoa
 ###  in eoa 2.0.7, the function doesn't work with custom priors
 # choices for prior are:
 # IbinRef = integrated reference prior for binomial
 # binRef = reference prior for binomial (gives Inf for P(M = 0 | X = 0)
 # uniform
 # custom prior = numeric array with columns m and p(M = m)
  if (length(x) * length(g) != 1){
    warning("error in data: x and g must be scalars")
    return(NA)
  }
  if (is.na(g)){
    warning("error in data: x and g must be numeric")
    return(NA)
  }
  if (!is.numeric(x) | !is.numeric(g)){
    warning("error in data: x and g must be numeric")
    return(NA)
  }
  if (x < -0.000001 || abs(round(x) - x) > 0.000001){
    warning("error in data: x must be a non-negative integer")
    return(NA)
  }
  if (g < 0.00001){
    warning("error in data: g must be strictly greater than 0")
    return(NA)
  }
  if (g > 0.99999) return (x)
  if (!is.na(mmax) &&
    (!is.numeric(mmax) || abs(round(mmax)-x) < 0.000001 || round(mmax) < 1)){
    print("problem is here")
    warning("error in mmax")
    return(NA)
  }
  x <- round(x)
  if (!is.numeric(prior) && !is.character(prior)){ # custom prior
    warning("error in prior")
    return(NA)
  } else {
    if (is.numeric(prior)){
    # numeric => must be two-dimensional array
      if (MpriorOK(prior)){
        if (is.na(mmax)){
          mmax <- max(prior[, 1])
        } else {
          if (round(mmax) != round(max(prior[, 1]))) {
            warning("mmax != max of custom prior. Aborting calculation.")
            return(NA)
          } else {
            mmax <- round(max(prior[, 1]))
          }
        }
        M <- x:mmax
        prior_M <- prior[,2]/sum(prior[,2])
      } else { # custom prior is not correctly formatted
        return(NA)
      }
    } else { # named priors (following different philosophies of uninformed)
      if (prior %in% c("IbinRef","binRef", "uniform")) {
        mmax <- round(ifelse(is.na(mmax), fmmax(x,g), mmax))
        M <- x:mmax
        prior_M <- switch(prior,
          IbinRef = diff(sqrt(0:(mmax+1)))/sum(diff(sqrt(0:(1+mmax)))),
          binRef = 1/sqrt(0:mmax),
          uniform = rep(1, mmax + 1)
        )
      } else {
          warning("error in prior")
          return(NA)
      }
    }
  }
  if (prior_M[1] == Inf & x == 0) return (1)
  pXgM <- dbinom(x, M, g)
  pM <- prior_M[x:mmax+1]
  pMgX <- pXgM*pM; pMgX <- pMgX/sum(pMgX) # posterior distribution for M
                                      #(ignoring M < X, which has probability 0)
  c(rep(0, x), pMgX) # full posterior, counting from m = 0 to mmax
}

#' @rdname postM
#' @export
postM.ab <- function(x, Ba, Bb, prior = "IbinRef", mmax = NULL){
  # error-checking
  suppressWarnings({
    if (length(x) * length(Ba) * length(Bb) != 1){
      print("error in data")
      return(NA)
    }
    if (is.na(as.numeric(x)) * is.na(as.numeric(Ba)) * is.na(as.numeric(Bb)) != 0){
      print("error in data: all values must be numeric")
      return(NA)
    }
    if (x < 0 || abs(round(x)-x) > 0.000001){
      print("error in data: x must be a non-negative integer")
      return(NA)
    }
    if (Ba < 0.000001){
      print("error in data: Ba must be positive")
      return(NA)
    }
    if (Bb < 0.000001){
      print("error in data: Bb must be positive")
      return(NA)
    }
  })
  x <- round(x)
  # check whether the mmax provided is valid
  if (!missing(mmax)){ # then mmax must be a non-negative integer
    if (!is.numeric(mmax)){
      warning("mmax must be a non-negative integer")
      return(NA)
    }
    if (abs(round(mmax) - mmax) > 0.000001){
      warning("mmax must be an integer")
      return(NA)
    }
    if (mmax < -0.000001){
      warning("mmax must be non-negative")
      return(NA)
    }
  }
  # check prior
  # if a custom prior is entered, check the format; if string is entered, calculated prior
  if (!is.numeric(prior) && !is.character(prior)){ # format is neither custom (numeric) nor named (string)
    warning("error in prior")
    return(NA)
  } else {
    if (is.numeric(prior)){ 
    # numeric => custom prior; must be two-dimensional array with probabilities starting at m = 0
      if (MpriorOK(prior)){
        if (missing(mmax)){
          mmax <- max(prior[, 1])
        } else {
          if (round(mmax) != round(max(prior[, 1]))) {
            warning("mmax != max of custom prior. Aborting calculation.")
            return(NA)
          } else {
            mmax <- round(max(prior[, 1]))
          }
        }
        M <- x:mmax
        prior_M <- prior[, 2]
      } else { # custom prior is not correctly formatted (numeric but not two columns etc.)
        return(NA)
      }
    } else { # named priors (following different philosophies of uninformed)
      if (prior %in% c("IbinRef", "binRef", "IbetabinRef", "betabinRef", "uniform")) {
        mmax <- round(ifelse(missing(mmax), fmmax.ab(x, Ba, Bb), mmax))
        M <- x:mmax
        prior_M <- switch(prior,
          IbinRef = diff(sqrt(0:(mmax + 1)))/sum(diff(sqrt(0:(1 + mmax)))),
          binRef = 1/sqrt(0:mmax),
          IbetabinRef = diff(log(sqrt(0:(mmax + 1)+ Ba+ Bb)+sqrt(0:(mmax + 1)))),
          betabinRef = 1/sqrt(0:mmax*(0:mmax + Ba + Bb)),
          uniform = rep(1, mmax + 1)
        )
      } else {
          warning("error in prior")
          return(NA)
      }
    }
  }

  if (prior_M[1] == Inf & x == 0) return (1)
  pXgM <- VGAM::dbetabinom.ab(x, M, Ba, Bb)
  pM <- prior_M[x:mmax + 1]
  pMgX <- pXgM*pM; pMgX <- pMgX/sum(pMgX) # posterior distribution for M (ignoring M < X, which has probability = zero)
  pMgX <- c(rep(0, x), pMgX)
  pMgX
}

#' @rdname postM
#' @export
calcMstar <- function(pMgX, alpha){
  min(which(cumsum(pMgX) >= 1 - alpha)) - 1
}

#' @rdname postM
#' @export
MCI <- function(pMgX, crlev = 0.95){
  cs <- cumsum(pMgX)
  aM <- 1 - crlev
  lwrbnd <- min(which(cs > aM/2)) - 1
  lwrArea <- ifelse(lwrbnd == 0, 0, cs[lwrbnd])
  uprbnd <- min(which(cs > 1 - aM + lwrArea)) - 1
  c(lwrbnd, uprbnd)
}

#' @title Find suitable mmax for clipping improper priors for M
#'
#' @description Improper priors need to be clipped in order to be usable.
#'  \code{fmmax} and \code{fmmax.ab} find values of \eqn{m} that are large enough
#'  that the probability of exceeding is less than 0.0001 (depends on \eqn{g} and
#'  \eqn{X}).
#'
#' @param x carcass count
#'
#' @param g overall carcass detection probability
#'
#' @param pBa,pBb parameters for beta distribution characterizing estimated \eqn{g}
#'
#' @return integer \eqn{m} such that \eqn{Pr(M >= m) < 0.0001}
#'
#' @export
#'
fmmax <- function(x, g){ # find the maximum m to sum over in calculating posterior
  if (pbinom(x, size = 1e5, prob = g) > 0.0001) {
    # if too huge of M's are required, return a large one and give warning
    warning(paste0("P(X <= ", x," | g = ", g, ", m = 100000) = ",
      signif(pbinom(x, size = 100000, prob = g), 6), ". Using mmax = 100000..."))
    return(1e5)
  }
  mmax <- x
  while (1){
    m <- mmax:(mmax + 100)
    mmax <- mmax + 100
    if (pbinom(x, size = mmax, prob = g) < 0.0001){
      # pick mmax large enough so that greater m's are very unlikely
      mmax <- m[min(which(pbinom(x, m, g) < 0.0001))]
      break
    }
  }
  mmax
}
#' @rdname fmmax
#' @export
fmmax.ab <- function(x, pBa, pBb){
  # find the maximum m to sum over in calculating posterior (for beta-binomial)
  if (VGAM::pbetabinom.ab(x, 1e5, pBa, pBb) > 0.0001) {
    # if too huge of M's are required, return a large one and give warning
    g <- pBa/(pBa + pBb)
    warning(paste0("P(X <= ", x," | g = ", g, ", m = 10000) = ",
      signif(pbinom(x, 10000, g), 6), ". Taking mmax = 10000..."
  ))
    return(1e5)
  }
  mmax <- x
  while (1){
    m <- mmax:(mmax + 100)
    mmax <- mmax + 100
    if (VGAM::pbetabinom.ab(x, size = mmax, shape1 = pBa, shape2 = pBb) < 0.0001){
      mmax <- m[min(which(
        VGAM::pbetabinom.ab(x, size = m, shape1 = pBa, shape2 = pBb) < 0.0001
      ))]
      break
    }
  }
  mmax
}

#' @title Check validity of format of custom prior for M
#'
#' @param prior a custom prior for M must be a matrix with columns for M and
#'  and associated probabalities P(M = m). The M column must begin at 0 and the
#'  probabilities must sum to 1.
#'
#' @return boolean. Is the prior formatted properly?
#'
#' @export
#'
MpriorOK <- function(prior){
  if (is.numeric(prior)){ # numeric => must be two-dimensional array with probabilities starting at m = 0
    if (length(dim(prior)) != 2){
      warning("error in prior")
      return(F)
    }
    if (abs(sum(diff(prior[, 1]) - rep(1, length(prior[, 1]) - 1))) > 0.00000001){
      # differences in M values not uniformly = 1
      warning("error in prior")
      return(F)
    }
    if (prior[1, 1] != 0){
      warning("error in prior: m[1] must be zero")
      return(F)
    }
    if (abs(sum(prior[, 2]) - 1) > 0.00001){
      warning("prior: probabilities must sum to 1")
      return(F)
    }
  } else {
    return(F)
  }
  return(T)
}