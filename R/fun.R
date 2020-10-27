#'  Create a siteLayout
#' @description Read shape files for search area polygons, turbine locations,
#'  and carcass discovery locations. If there is more than 1 turbine, there must
#'  be a turbine ID column (\code{unitCol}) in each of the files. The individual
#'  turbine ID's must be syntactically valid R names, i.e., contain combinations
#'  of letters, numbers, dot ( . ), and underscores ( _ ) only and not begin with
#'  a number or a dot followed by a number. Other characters are not allowed: no
#'  spaces, hyphens, parentheses, etc. All three shape files (search area
#'  polygons, turbine locations, carcass discovery locations) must have their
#'  three standard component files (.shp, .shx, .dbf) stored in the same folder.
#' @param file_layout name of a shape file (polygons or multipolygons) that
#'  delineates areas searched at each turbine. Search areas may be defined for
#'  any number of search classes (associated with different detection
#'  probabilities) and turbines. Search classes are listed the \code{scVar}
#'  column and turbine IDs are specified in the (\code{unitCol}) column. If a
#'  \code{unitCol} is provided (required if the site consists of more than one
#'  turbine), each turbine should have a syntactically valid name. Search class
#'  polygons at each turbine should be non-intersecting; otherwise, the function
#'  exits with an error.
#' @param file_CO name of a shape file (points) for carcass observations
#'  containing carcass discovery locations. If a \code{unitCol} name is provided,
#'  turbine IDs must be provided for  carcass discoveries. Search classes are
#'  optional.
#' @param file_turbine a shape file (points) giving the locations of the turbine
#'  centers for turbines listed in the \code{unitCol} column in the \code{file_layout}
#' @param unitCol (optional) column name for turbine IDs. If a \code{unitCol} is
#'  provided, a \code{file_turbine} with the turbine centers must also be
#'  provided and the \code{unitCol} must be present in \code{file_turbine},
#'  \code{file_layout}, and \code{file_CO}. Turbine IDs in the \code{unitCol}
#'  must be syntactically valid R names (see Description above).
#' @return A \code{siteLayout} object, which is a list with spatial
#'  characteristics of the site, including
#'  \describe{
#'    \item{\code{$layout}}{turbine search area configurations (polygons and multipolygons)
#'      from \code{file_layout} shape file as an \code{sf} object.}
#'    \item{\code{$layoutAdj}}{polygons from \code{$layout} but all recentered at (0, 0)}
#'    \item{\code{$turbines}}{turbine centers (as \code{sf} object)}
#'    \item{\code{$carcasses}}{carcass discovery locations along with (optional)
#'      incidental information (turbine, visibility class, species,)}
#'    \item{\code{$unitCol, $X, $tset, $tcenter}}{name of the column with turbine IDs
#'      (\code{$unitCol}, character string), number of carcasses found at each turbine
#'      (\code{$X}, vector), turbine names (\code{$tset}, vector of character strings),
#'      locations of turbine centers (\code{$tcenter}, nturb by 2 array) UTMs of
#'      turbine centers, extracted and simplified from \code{$turbines}}
#'  }
#' @export
#'
readLayout <- function(file_layout, file_CO, file_turbine = NULL, unitCol = NULL){
  # read as shape files
  plotLayout <- sf::st_read(file_layout, stringsAsFactors = FALSE)
  carcasses <- sf::st_read(file_CO, stringsAsFactors = FALSE)
  plotLayout0 <- sf::st_drop_geometry(plotLayout) # easier quick subsetting
  carcasses0 <- sf::st_drop_geometry(carcasses)
  #### error checks and formatting
  ## shape files
  # format of unitCol
  if (!is.null(unitCol)){
    if (is.null(file_turbine))
      stop("if unitCol is not NULL, file_turbine must be provided")
    turbines <- sf::st_read(file_turbine, stringsAsFactors = FALSE)
    turbines0 <- sf::st_drop_geometry(turbines)
    if (!unitCol %in% names(turbines))
      stop(unitCol, " column not included in turbine data")
    if (!unitCol %in% names(plotLayout))
      stop(unitCol, " column not included in plot layout data")
    if (!unitCol %in% names(carcasses))
      stop(unitCol, " column not included in carcass data")
    if (!all(plotLayout0[, unitCol] %in% turbines0[, unitCol]))
      stop("a turbine that is included in ", file_layout, " is missing from ",
        file_turbine)
    if (!all(carcasses0[, unitCol] %in% plotLayout0[, unitCol]))
      stop("a turbine that is included in ", file_CO, " is missing from ",
        file_turbine)
    if (!all(make.names(carcasses0[, unitCol, ]) == carcasses0[, unitCol])){
      badind <- which(make.names(carcasses0[, unitCol]) != carcasses0[, unitCol])
      badnm <- unique(carcasses0[badind, unitCol])
      badnm <- badnm[1:min(length(badnm), 3)]
      message("\n\nNOTE: Not all turbine names are syntactically valid.\n",
        "These (", paste(badnm, collapse = ", "),
        "...) can be converted to syntactically valid names (",
        paste(make.names(badnm), collapse = ", "), "...) if desired.")
      tmp <- readline("convert [y]? or abort [enter anything besides y]? ")
      if (!identical(tolower(tmp), "y")) stop("aborting...")
      carcasses[, unitCol] <- make.names(carcasses0[, unitCol])
      carcasses0[, unitCol] <- make.names(carcasses0[, unitCol])
      turbines[, unitCol] <- make.names(turbines0[, unitCol])
      plotLayout[, unitCol] <- make.names(plotLayout0[, unitCol])
    }
    tset <- as.character(as.data.frame(turbines)[, unitCol]) # turbine names
    tcenter <- sf::st_coordinates(turbines) # turbine centers (matrix with "X", "Y")
    rownames(tcenter) <- tset
  } else {
    tset <- "t"
    tcenter <- matrix(c(0, 0), nrow = 1, dimnames = list(tset, c("X", "Y")))
  }
  layoutAdj <- plotLayout # Adj => coordinates relative to turbine centers
  for (ai in 1:nrow(plotLayout)){
    ti <- plotLayout[ai, unitCol, drop = TRUE]
    layoutAdj[ai, ] <- sf::st_set_geometry(layoutAdj[ai, ],
      sf::st_geometry(plotLayout[ai, ]) - tcenter[ti,]) # recentering
  }
  # maximum distance from turbine to a searched point (radius of interpolation)
  n <- numeric(length(unique(plotLayout[, unitCol, drop = TRUE])))
  names(n) <- gtools::mixedsort(unique(plotLayout[, unitCol, drop = TRUE]))
  tabn <- table(carcasses0[, unitCol])
  n[names(tabn)] <- as.vector(tabn)
  ans <- list(layout = plotLayout, layoutAdj = layoutAdj, turbines = turbines,
    carcasses = carcasses, unitCol = unitCol, X = n,
    tset = tset, tcenter = tcenter)
  class(ans) <- "siteLayout"
  return(ans)
}

#'  Create ringData structure for fitting distance models
#' @description Create necessary data structures from \code{siteLayout} object
#'  for modeling carcass density as a function of distance from turbine and for
#'  modeling DWP.
#' @param site_layout a \code{siteLayout} object, which summarizes the spatial
#'  characteristics of the site as listed in the associated GIS shape files
#'  (search area polygons, turbine centers, and carcass discovery locations)
#' @param scVar (optional) name of a column with search class delineator in
#'  \code{site_layout$layout}.
#' @param notSearched name(s) of the search classes that are not searched
#'  (vector of character strings)
#' @return List containing:
#'  \describe{
#'    \item{\code{$ringData}}{data frame with a turbine-aggregated layout with distance
#'      (outer radius of 1m rings), carcass counts, exposure (area seached in
#'      each ring), and an optional search class column;}
#'    \item{code{$dwpArea}}{list of turbines with the area searched in at each
#'      distance from 1 to maximim distance searched at any turbine;}
#'    \item{\code{$searchMax}}{maximum distance searched at any turbine (numeric
#'      scalar);}
#'    \item{\code{$scVar}}{name of the search class column (if any);}
#'    \item{\code{$X}}{number of carcasses found at each turbine (vector}
#'  }
#' @export
prepRing <- function(siteLayout, scVar = NULL, notSearched = NULL){
  if (!is.null(scVar)){
    if (!scVar %in% names(siteLayout$layout))
      stop(scVar, " not in siteLayout$layout")
    if (!is.null(notSearched)){ # remove unsearched areas from data
      if (!is.character(notSearched))
        stop("notSearched must be a character vector (or NULL)")
      siteLayout$layout <-  siteLayout$layout[
        !siteLayout$layout[, scVar, drop = TRUE] %in% notSearched, ]
      siteLayout$layoutAdj <- siteLayout$layoutAdj[
        !siteLayout$layoutAdj[, scVar, drop = TRUE] %in% notSearched, ]
      if (scVar %in% names(siteLayout$carcasses))
        siteLayout$carcasses <-  siteLayout$carcasses[
          !siteLayout$carcasses[, scVar, drop = TRUE] %in% notSearched, ]
      if (length(unlist(sf::st_within(siteLayout$carcasses, siteLayout$layout))) <
          nrow(siteLayout$carcasses)){
        stop("carcass(es) found in area not searched.")
      }
    }
  }
  searchMax <- sqrt(max(rowSums((
    sf::st_coordinates(siteLayout$layoutAdj)[, c("X", "Y")])^2)))
  theta <- seq(0, 2 * pi, length = 300)
  trig <- cbind(cos(theta), sin(theta)); trig[dim(trig)[1], ] <- trig[1, ]
  radi <- 1:searchMax
  rings <- sf::st_sfc(sf::st_polygon(list(trig)))
  for (i in 2:length(radi))
    rings <- c(rings, sf::st_sfc(sf::st_polygon(list(radi[i] * trig))))
  rings <- sf::st_sf(rings)
  for (i in length(radi):2) rings[i, ] <- sf::st_difference(rings[i, ], rings[i-1, ])
  rings$r <- radi
  sf::st_crs(rings) <- sf::st_crs(siteLayout$layoutAdj)
  ## collapse search classes into one, i.e., searched area = union(search classes)
  layoutUnion <- list()
  for (ti in siteLayout$tset){
    layoutUnion[[ti]] <- sf::st_union(siteLayout$layoutAdj[
      as.data.frame(siteLayout$layoutAdj)[, unitCol] == ti, ])
  }
  dwpPoly <- list()
  for (ti in siteLayout$tset){# for (each tubine)
  #   create a unionized search polygon and calculate the area of intersection by ring
    dwpPoly[[ti]] <- sf::st_intersection(sf::st_geometry(rings),
      sf::st_union(siteLayout$layoutAdj[
        as.data.frame(siteLayout$layoutAdj)[, unitCol] == ti, ])
    )
  }
  # area searched in each ring at each turbine (without regard to search class)
  dwpArea <- lapply(dwpPoly, sf::st_area)
  # for each ring, what is the area sampled?
  int <- unlist(lapply(layoutUnion, length))
  for (ti in names(layoutUnion)){
    if (int[ti] > 0){
      dwpArea[[ti]] <- cbind(
        rings$r[sf::st_intersects(layoutUnion[[ti]], rings, sparse = FALSE)],
        dwpArea[[ti]]
      )
    }
  }
  for (ti in names(dwpArea))
    if (int[ti]) colnames(dwpArea[[ti]]) <- c("r", "area")
  rAll <- sort(unique(unlist(lapply(dwpArea,
    FUN = function(x) if (is.null(dim(x))) return (x) else return (x[, 1])))))
  dwpArea <- lapply(dwpArea, FUN = function(x){
    if (length(x) == 0){
      tmp <- matrix(c(rAll, numeric(length(rAll))), ncol = 2,
        dimnames = list(NULL, c("r", "area")))
    } else {
      ind <- which(!rAll %in% x[, "r"])
      if (length(ind) > 0){
        tmp <- rbind(x, cbind(rAll[ind], 0))
      } else {
        tmp <- x
      }
    }
    return(tmp[order(tmp[, 1]),])
  })
  ### create a data frame with area in each search class at each distance and
  ## tally the number of carcasses
  trscArea <- list() # nturbines x nrings arrays of areas in each search class
  scset <- gtools::mixedsort(unique(as.data.frame(siteLayout$layoutAdj)[, scVar]))
  for (sci in scset){
    trscArea[[sci]] <- matrix(0,
      nrow = length(siteLayout$tset), ncol = dim(rings)[1])
    rownames(trscArea[[sci]]) <- siteLayout$tset
    for (ti in siteLayout$tset){
      # indices for turbine = ti and scVar = sci
      ind <- which(siteLayout$layout[, unitCol, drop = TRUE] == ti &
        siteLayout$layout[, scVar, drop = TRUE] == sci)
      # polygons for ind
      for (ii in ind){
        tring <- sf::st_geometry(siteLayout$layoutAdj[ii,])
        jj <- unlist(sf::st_intersects(siteLayout$layoutAdj[ii, ], rings))
        trscArea[[sci]][ti, jj] <- trscArea[[sci]][ti, jj] +
          as.vector(sf::st_area(sf::st_intersection(tring, rings)))
      }
    }
  }
  trscArea <- lapply(trscArea, colSums)
  ringData <- data.frame(array(0, dim = c(dim(rings)[1] * length(scset), 4)))
  names(ringData) <- c("r", "scVar", "exposure", "ncarc")
  ringData$r <- rep(radi, length(scset))
  ringData[, "scVar"] <- rep(scset, each = length(radi))
  ringData$exposure <- unlist(trscArea)
  tmp <- suppressWarnings(sf::st_intersection(
    siteLayout$layout, siteLayout$carcasses))
  carcassSumry <- tmp[, c("Turbine", "Class")]
  carcr <- sqrt(rowSums((sf::st_coordinates(siteLayout$carcasses) -
    siteLayout$tcenter[siteLayout$carcasses[, unitCol, drop = TRUE], ])^2))
  for (ci in 1:nrow(carcassSumry)){
    ind <- which(ringData$r == ceiling(carcr[ci]) &
      ringData$scVar == unlist(as.matrix(carcassSumry)[ci, scVar]))
    ringData[ind, "ncarc"] <- ringData[ind, "ncarc"] + 1
  }
  rd <- ringData[ringData$exposure > 1e-3 & !ringData$scVar %in% notSearched,]
  return(list(ringData = rd, dwpArea = dwpArea, searchMax = searchMax,
    scVar = scVar, X = siteLayout$X))
}


#' Fit Distance Distribution Model(s)
#'
#' @description Fit glm's for distance distribution models corresponding to
#'  gamma, lognormal, log-linear, log-quadratic, log-cubic, truncated normal,
#'  exponential, Rayleigh, Pareto, inverse gamma, paranormal gamma, Maxwell
#'  Boltzmann, chi-squared, and/or inverse Gaussian distributions. By default,
#'  \code{ddFit} fits all the aforementioned models, but any subset of the models
#'  may be fit as well.
#'
#'  The glm is converted to a probability distribution by dividing by a
#'  normalizing constant, namely the integral of the glm evaluated from 0 to
#'  infinity. In some cases (most notably when the leading coefficient of the
#'  glm is positive so the fitted curve does not converge to zero as x increases),
#   the integral does not converge to a finite value and the glm cannot be
#'  converted to a probability distribution. In these cases, the distribution
#'  parameters are given as \code{NA}.
#'
#' @param ringData Data frame with a tally of carcasses within 1 m rings, along
#'  with outer radii of the rings, (optional) search class column, and the amount
#'  of area searched in each search class and ring.
#' @param model names (vector of character strings) of glm distribution templates
#'  to fit. Default is \code{"all"}, which fits \code{gamma}, \code{lognormal},
#'  \code{logLinear}, \code{logQuadratic}, \code{logCubic}, \code{inverse_gamma},
#'  \code{paranormal_gamma}, \code{Rayleigh}, \code{MaxwellBoltzmann}, \code{constant},
#'  \code{tnormal}, \code{exponential}, \code{Pareto}, \code{chisq}, and
#'  \code{inverse_gaussian}. Any subset of these may also be fit with a single
#'  call to \code{ddFit}.
#' @param scVar name of column in \code{ringData} with (optional) subdivisions
#'  of rings into search classes (character strings)
#' @param rCol name of column in \code{ringData} with outer radius for 1m rings
#' @param expoCol name of "exposure" column of total area in each ring belonging
#'  to each search class.

#' @return A list of fitted glm models as \code{dd} objects in a \code{ddArray}
#'  object if a vector of distributions is fit or a single \code{dd} object if a
#'  single model is fit. The \code{dd} objects are lists with all the standard
#'  \code{\link[stats]{glm}} components, in addition to the following elements:
#'  \describe{
#'    \item{\code{$distr}}{name of the distribution ("gamma", etc.)}
#'    \item{\code{$parms}}{vector of distribution parameter estimates (or \code{NA}
#'      if the model based on the MLE is not extensible)}
#'    \item{\code{$varbeta}}{the variance-covariance matrix of the glm parameter
#'      estimates. NOTE: This is identical to the covariance matrix from the glm,
#'      which can be extracted via \code{summary(x)$cov.unscaled}}
#'    \item{\code{$scVar}}{name of the (optional) search class variable (or \code{NULL})}
#'    \item{\code{$ncarc}}{number of carcasses}
#'    \item{\code{$n}}{number of rings}
#'    \item{\code{$k}}{number of parameters}
#'    \item{{\code{$srad}}{search radius}
#'  }
#'  When a \code{dd} object is printed, only a small subset of the elements are
#'  shown. To see a full list of the objects, use \code{names(x)}. The elements
#'  can be extracted in the usual R way via \code{$} or \code{[[x]]}.

#' @export
ddFit <- function(ringData, model = "all", scVar = NULL,
      rCol = "r", expoCol = "exposure", ncarcCol = "ncarc"){
  dat <- ringData[, c(rCol, scVar, expoCol, ncarcCol)]
  names(dat) <- c("r", if(!is.null(scVar)) "scVar", "exposure", "ncarc")
  scl <- ifelse(is.null(scVar), "", paste0(" + ", scVar))
  plu <- ifelse(!is.null(scVar), " + ", " ")
  sufx <- paste0(plu, scVar,  "+ offset(log(", expoCol, "))")
  sufx0 <- paste0(plu, scVar, "+ ")
  pref <- paste0("ncarc ~ ")
  output <- list()
  if (identical(model, "all")) model <- mod_name # all model names; see parameters.r
  for (distr in model){
    form <- formula(paste0("ncarc ~ ", paste(cof_name[[distr]][-1], collapse = " + "),
      scl, " + ", paste0("offset(", mod_offset[[distr]], ")")))
    output[[distr]] <- glm(
      formula = form,
      data = dat,
      family = "poisson"
    )
    output[[distr]]$distr <- distr
    output[[distr]]$form <- form
    output[[distr]]$scVar <- scVar
    output[[distr]]$parms <- cof2parms(output[[distr]]$coefficients, distr)
    output[[distr]]$varbeta <- summary(output[[distr]])$cov.unscaled
    output[[distr]]$ncarc <- sum(dat$ncarc)
    output[[distr]]$n <- nrow(dat)
    output[[distr]]$k <- output[[distr]]$n - output[[distr]]$df.residual
    class(output[[distr]]) <- c("dd", class(output[[distr]]))
  }
  if(length(output) == 1) return(output[[1]])
  class(output) <- "ddArray"
  return(output)
}

#' Calculate Akaike Information Criterion (AIC) for Distance Distributions
#'
#' @description functions for calculating AIC and AICc for carcass dispersion
#'  models.
#' @param x list of models (\code{ddArray}) or single model (\code{dd}) to calculate AICs for
#' @param extent Include only the extensible models (\code{extent = "full"}) or
#'  all models, regardless of whether they can be extended beyond the search
#'  radius (\code{extent = "win"}).
#' @param ... ignored
#' @return Data frame with AIC, AICc, deltaAICc for all models in \code{x}
#' @export
aic <- function(x, ...) UseMethod("aic", x)

#' @rdname aic
#' @export
aic.ddArray <- function(x, extent = "full", ...){
  # column names for data frame to be returned
  nm <- c("model", "k", "AIC", "AICc", "deltaAICc")
  if (extent == "full"){
    incmod <- NULL
    for(mod in x){
      if (cofOK(mod$coefficients, mod$distr)) incmod <- c(incmod, mod$distr)
    }
  } else if (extent == "win"){
    incmod <- names(x)
    nm <- c(nm, "extensible")
  }

  output <- data.frame(array(dim = c(length(incmod), length(nm)),
    dimnames = list(NULL, nm)), stringsAsFactors = FALSE)
  output$model <- incmod
  ci <- c("k", "AIC", "AICc")
  for (distr in incmod){
    output[output$model == distr, ci] <- aic(x[[distr]])[ci]
    if (extent == "win")
      output$extensible[output$model == distr] <- 1 * cofOK(x[distr]$coefficients, distr)
  }
  output$deltaAICc <- output$AICc - min(output$AICc, na.rm = TRUE)
  output <- output[order(output$deltaAICc), ]
  attr(output, "extent") <- extent
  attr(output, "n") <- attr(x, "n")
  return(output)
}

#' @rdname aic
#' @export
aic.dd <- function(x, ...){
  aic0 <- x$aic
  return(c(
    k = x$k,
    AIC = aic0,
    AICc = aic0 + 2 * x$k * (x$k + 1)/(x$n - x$k - 1)
  ))
}

#'  estDWP...obsolete?
#' @param data data
#' @param model model
#' @param nboot nboot
#' @return answer
#' @export
estDWP <- function(data, model, nboot = 1000){
# if ringData:
  radi <- sort(unique(data$ringData$r))
  radi <- radi - diff(c(0, radi))/2
  parms <- predict(
    object = model,
    newdata = data.frame(
      r = radi,
      scVar = data$ringData$scVar[1],
      exposure = 1),
    se.fit = TRUE)
  z <- rnorm(nboot)
  deno <- colSums(2 * pi * (radi - 0.5) * exp(parms$fit + outer(parms$se.fit, z)))
  numo <- lapply(data$dwpArea, FUN = function(x){
    radi <- x[, "r"] - diff(c(0, x[, "r"]))/2
    parms <- predict(
      object = model,
      newdata = data.frame(
        r = radi,
        scVar = data$ringData[, "scVar"][1],
        exposure = 1),
      se.fit = TRUE)
    colSums(x[, "area"] * exp(parms$fit + outer(parms$se.fit, z)))
  })
  dwp <- lapply(numo, FUN = function(x) x/deno)
  dwp_means <- list()
  for (ti in names(data$dwpArea)) dwp_means[[ti]] <- mean(dwp[[ti]])
  dwp_bco <- mapply(FUN = function(x, p)
    1/(1 + exp(-GenEst::logit(p) - var(GenEst::logit(p))/2)),
    x = as.list(data$X), p = dwp)
  # bco gives the p that results in unbiased estimate of 1/dwp and should be
  # good for estimating M, BUT it is based on the assumption that x is a constant
  # rather than a realization of a random process, so the binomial variation in
  # x must also be taken into account. This is the double bootstrap of the RP
  # paper.
  # the dwp_bco2 approach does xtilde on bco
  dwp_bco2 <- suppressWarnings(mapply(FUN = function(x, p)
    x/((rbinom(length(p), round(x/p), p) - (round(x/p)*p - x))/p),
    x = as.list(data$X), p = as.list(data.frame(dwp_bco))))
  dwp_bco2[dwp_bco2 < 0] <- 0

}

#' Tally carcasses by ring
#' @description Convert a vector of carcass distances into a data frame that
#'  gives a tally of the number of carcasses in 1m rings, appropriate for modeling
#'  carcass dispersion. Assumption is that all area within the search radius was
#'  searched and with uniform detection probability. Covariates associated with
#'  carcass counts are not allowed in this simplified data structure.
#' @param rdata data frame or array with a column giving the distances at which
#'  carcasses were found. Alternatively, \code{rdata} may be a simple vector of
#'  distances.
#' @param srad Search radius
#' @param rCol name of the column in \code{rdata} that includes the carcass
#'  distances.
#' @return Data frame with columns:
#'  \code{r} = the outer radius of 1m rings from 1 to \code{ceiling(srad)},
#'  \code{exposure} = area in ring = \code{(2 * pi * (r - 0.5))}
#'  \code{ncarc} = number of carcasses observed in each ring
#' @export
prepRing_r <- function(rdata, srad, rCol = "r"){
  if (is.vector(rdata)) {
    rdata <- data.frame(rdata)
    names(rdata) <- rCol
  }
  r <- 1:ceiling(srad)
  ans <- data.frame(r = r, exposure = 2 * pi * (r - 0.5), ncarc = 0)
  ctab <- table(ceiling(rdata[, rCol]))
  ans$ncarc[as.numeric(names(ctab))] <- as.vector(ctab)
  return(ans)
}

#'  plot dd and ddArray objects
#' @description Plot CDF, PDF, or rcd (relative carcass density) for a single
#'  carcass dispersion glm model (\code{dd} object) or a list of models
#'  (\code{ddArray} object).
#' @param x model(s) to plot
#' @param type Type or representation of carcass dispersion to plot:
#'  \code{"CDF"}, \code{"PDF"}, or \code{"rcd"}. The \code{"CDF"} gives the
#'  fraction of carcasses falling within \code{r} meters from a turbine and
#'  \code{"PDF"} is the associated probability density. The \code{"rcd"} gives the
#'  relative carcass density at a point \code{r} meters from a turbine and is
#'  PDF/(2 * \pi * r).
#' @param extent Plot dispersions as fraction of total carcasses (\code{"full"})
#'  or as fraction of carcasses within the searched area (\code{"win"}).
#' @param set To show lines for all fitted models, use \code{set = "all"}; or
#'  \code{set} = vector of model names to show lines for.
#' @param xmax maximum distance to show in the graph; if \code{xmax = NULL}, the
#'  maximum distance is taken as the max distance in the data set to which the
#'  models were fit.
#' @param resolution The number of line segments to break the curves into when
#'  plotting (i.e., \code{x = seq(0, xmax, length.out = resolution)}). Higher
#'  resolutions give smoother-looking curves.
#' @param mod_highlight Character string giving the name of the model to
#'  highlight by plotting it last and with \code{lwd = 2}. If \code{NULL}, the
#'  curve associated with the lowest (best) AICc score is highlighted.
#' @param CL confidence level to show in a \code{dd} plot (ignored for
#'  \code{ddArray} objects)
#' @param nsim Number of simulation reps to use for estimating confidence bounds
#'  for \code{dd} plot (ignored for \code{ddArray} objects)
#' @details \code{ddArray} objects are plotted with lines in order of decreasing
#'  AICc, so that the "better" models are closer to the top and more prominent.
#'  The model with the lowest AICc ("best" model) is plotted last with a heavier
#'  line than the others.
#'
#'  \code{dd} plots the curve for the MLE of the parameters, along with a
#'  100\code{CL}\% confidence bounds determined for \code{nsim} simulation reps
#'
#' @name plot
#' @export
plot.ddArray = function(x, type = "CDF", extent = "full", set = "all",
    xmax = NULL, resolution = 250, mod_highlight = NULL, par_reset = TRUE, ...){
  if (identical(set, "all")) set <- names(x)
  if (extent == "full")
    set <- set[sapply(x, function(tmp) !any(is.na(tmp$parm)))]
  if (type == "rcd")
    set <- set[which(!set %in% c("exponential", "tnormal", "constant"))]
  if (!any(set %in% mod_name)) stop("plot.ddArray: some model(s) undefined")
  dd <- x[set]
  aic_proper <- aic(dd)
  if (!type %in% c("CDF", "PDF", "rcd"))
    stop("type (", deparse(substitute(type)),
          ") must be \"PDF\", \"CDF\", or \"rcd\"")
  mod_best <- aic_proper$model[1]
  if (is.null(mod_highlight)) mod_highlight <- mod_best
  # unpacking ... args
  arglist <- list(...)
  if ("xlim" %in% names(arglist)){
    xlim <- arglist$xlim
    xmax <- max(xlim)
    xmin <- min(xlim)
  } else {
    if (is.null(xmax) || !is.finite(xmax)) xmax <- max(dd[[1]]$data[, "r"])
    xmin <- 0
    xlim = c(0, xmax)
  }
  xseq <- seq(xmin, xmax, length = resolution)
  ncol <- ifelse(nrow(aic_proper) >= 6, 2, 1)
  sz <- ifelse(round(nrow(aic_proper)/ncol) <= 3, 0.13,
          ifelse(round(nrow(aic_proper)/ncol) <= 5, 0.17, 0.22))
    # size of the main graph; vertical, relative to plot window size
  do.call(par, par_default)
  par(fig = c(0, 1, 0, sz), mar = c(0, 1, 0.5, 0), family = "mono")
  plot(0, type = "n", axes = F, xlab = "", ylab = "")
  leglab <- character(nrow(aic_proper))
  lwd <- numeric(nrow(aic_proper)) + 1
  if(mod_highlight %in% aic_proper$model) lwd[1] <- 2
  for (i in 1:nrow(aic_proper))
   leglab[i] <- sprintf("%18-s%6.2f", aic_proper$model[i], aic_proper$deltaAIC[i])
  mtext(side = 3, line = 0, adj = 0, sprintf("%19-s%s", "Distribution", "\u0394AICc"))
  lty <- 1 + !natural[aic_proper$model]
  names(lty) <- aic_proper$model
  legend(x = "topleft", legend = leglab, cex = 0.8, lwd = lwd, bty = "n",
    ncol = ncol, lty = lty, col = mod_color[aic_proper$model])
  if (type == "CDF"){
    par(fig = c(0, 1, sz, 1), mar = c(4, 4, 2.5, 0.5), family = "sans",
      new = TRUE) # for main graph
    xlab = ifelse("xlab" %in% names(arglist), arglist$xlab,
      "Distance from Turbine")
    ylab = ifelse("ylab" %in% names(arglist), arglist$ylab,
      "P(carcass falls within x meters from turbine)")
    pprm <- arglist
    pprm$x = 0
    pprm$type = "n"
    pprm$xlim = xlim
    pprm$ylim = 0:1
    pprm$xlab = xlab
    pprm$ylab = ylab
    do.call(plot, pprm)
#      plot(0, type = "n", xlim = xlim, ylim = 0:1, xlab = xlab, ylab = ylab, ...)
    for (fi in aic_proper$model[nrow(aic_proper):1])
      lines(xseq, pdd(xseq, model = dd[fi], extent = extent),
        lty = lty[fi], col = mod_color[fi])
    if(mod_highlight %in% names(dd)){
      lines(xseq, pdd(xseq, model = dd[mod_highlight], extent = extent),
        col = mod_color[mod_highlight], lwd = 2, lty = lty[mod_highlight])
    }
    if (extent == "win")
      mtext(side = 1, line = -1, adj = 1, "within search radius", family = "serif")
  } else if (type == "PDF"){
    par(fig = c(0, 1, sz, 1), mar = c(4, 4, 2.5, 0.5), family = "sans",
      new = TRUE) # for main graph
    ymax = -Inf
    for (fi in aic_proper$model){
      ymax <- max(ymax, max(ddd(xseq, dd[fi], extent = extent)))
    }
    plot(0, type = "n", xlim = xlim, ylim = c(0, ymax),
      xlab = "r", ylab = "PDF", ...)
    for (fi in aic_proper$model[nrow(aic_proper):1])
      lines(xseq, ddd(xseq, dd[fi], extent = extent),
        col = mod_color[fi], lty = lty[fi])
    lines(xseq, ddd(xseq, dd[mod_highlight], extent = extent),
      col = mod_color[mod_highlight], lwd = 2, lty = lty[mod_highlight])
    leglab <- character(nrow(aic_proper))
    if (extent == "win")
      mtext(side = 3, adj = 1, "within search radius", family = "serif", line = -1)
   } else if (type == "rcd"){
    ymax = -Inf
    for (fi in aic_proper$model){
      ymax <- max(ymax, max(rcd(x = xseq, model = dd[fi], extent = extent)))
    }
    par(fig = c(0, 1, sz, 1), mar = c(4, 4, 2.5, 0.5), family = "sans",
      new = TRUE) # for main graph
    plot(0, type = "n", xlim = xlim, ylim = c(0, ymax),
      xlab = "r", ylab = "Relative Carcass Density", ...)
    for (fi in aic_proper$model[nrow(aic_proper):1])
      lines(xseq, rcd(xseq, model = dd[fi], extent = extent),
        col = mod_color[fi], lty = lty[fi])
    lines(xseq, rcd(xseq, model = dd[mod_highlight], extent = extent),
      col = mod_color[mod_highlight], lwd = 2, lty = lty[mod_highlight])
    if (extent == "win")
      mtext(side = 3, adj = 1, "within search radius  ", family = "serif", line = -1)
  }
}

#' @name plot
#' @export
plot.dd <- function(x, type = "CDF", extent = "full",
    xmax = NULL, resolution = 250, nsim = 1000, CL = 0.9, ...){
  if (extent == "full" && anyNA(x$parms))
    stop("GLM non-extensible to distribution. Cannot plot")
  arglist <- list(...)
  if ("xlim" %in% names(arglist)){
    xmax <- max(arglist$xlim)
    xmin <- min(arglist$xlim)
  } else {
    xmin <- 0
  }
  if (is.null(xmax)) xmax <- max(x$data$r)
  xseq <- seq(xmin, xmax, length.out = resolution)
  if (!type %in% c("CDF", "PDF"))
    stop("type = ", deparse(substitute(type)), " not supported")
  CI <- ddCI(mod = x, x = xseq, type = type, CL = CL, nsim = nsim, extent = extent)
  do.call(par, par_default)
  if (type == "CDF"){
    ymax <- 1
    ylab <- "CDF = fraction of carcasses with x meters of turbine"
  } else if (type == "PDF"){
    ymax <- max(CI[, 3])
    ylab <- "PDF = probability density function for carcasses at distance x"
  }
  plot(0, xlim = range(xseq), ylim = c(0, ymax), type = "n",
    yaxs = "i", xaxs = "i",
    xlab = "x = Distance from Turbine (meters)",
    ylab = ylab,
    main = paste0("Distribution of Carcasses\n[", x$distr, " ", type, "]"))
  polygon(CI[c(1:nrow(CI), nrow(CI):1) , "r"] , c(CI[, 2], CI[nrow(CI):1, 3]),
    col = colors()[350], border = NA)
  msg <- switch(extent, win = "limited to carcasses within search radius  ")
  if (type == "PDF"){
    lines(xseq, ddd(xseq, x, extent = extent), col = mod_color[x$distr], lwd = 2)
    mtext(side = 3, line = -1, adj = 1, msg, family = "serif", lty = 1 + !natural[distr])
  } else {
    lines(xseq, pdd(xseq, x, extent = extent), col = mod_color[x$distr], lwd = 2)
    mtext(side = 1, line = -1, adj = 1, msg, family = "serif")
  }
  box()
}

#' Simulation of dispersion parameters
#' @param x object to simulate from
#' @param nsim number of simulation draws
#' @param ... ignored (but required format for CRAN)
#' @return array with simulated beta parameters from the glm model, and their
#'  conversion to distribution parameters
#' @export
ddSim <- function(x, ...) UseMethod("ddSim", x)

#' @rdname ddSim
#' @export
ddSim.dd <- function(x, nsim = 1000, extent = "full", ...){
  dd <- x
  bmean <- dd$coefficients
  bvar <- dd$varbeta
  if (any(is.na(bmean))) stop("ddSim: bad parameter estimates. cannot simulate")
  if (any(is.na(bvar))) stop("ddSim: bad var(beta_hat). cannot simulate")
  if (length(bmean) == 1){
    beta_sim <- rnorm(nsim, mean = bmean, sd = sqrt(bvar))
    beta_sim <- array(beta_sim, dim = c(length(beta_sim), 1))
    colnames(beta_sim) <- names(bmean)
  } else {
    beta_sim <- mvtnorm::rmvnorm(nsim, mean = bmean, sigma = bvar, method = "svd")
  }
  distr <- dd$distr
  attr(beta_sim, "distr") <- distr
  ans <- cbind(beta_sim, cof2parms(beta_sim, distr), extensible =
    ifelse(distr %in% c("exponential", "tnormal"), 1, 2) * cofOK(beta_sim, distr))
  attr(ans, "distr") <- distr
  attr(ans, "srad") <- max(dd$data$r, na.rm = TRUE)
  class(ans) <- "ddSim"
  ans
}

#' Subset simulated dispersion parameters while preserving attributes
#' @param x object to subset
#' @param ... ignored (but required format for CRAN)
#' @details Subset the ddSim object as if it were a simple matrix or array
#' @return array with simulated beta parameters from the glm model, their
#'  conversion to distribution parameters. NOTE: subsetting to a column or a row
#'  returns a matrix rather than a vector. This simplifies the coding and makes
#'  it easier to maintain integrity of data structures, but behavior differs from
#'  what is done when subsetting standard R matrices and arrays to a single column
#'  or row. Also unlike with standard R arrays and matrices, the class structure
#'  and attributes are preserved upon subsetting.
#' @export
"[.ddSim" <- function(x, i, j,...){
  y <- NextMethod("[")
  noclass <- FALSE
  if (!missing(j)){
    if (is.numeric(j)) j <- colnames(x)[j]
    if (!all(j %in% colnames(x))) stop("cannot subset ddSim object on given column")
    if (length(j) == 1) y <- matrix(y, ncol = 1, dimnames = list(NULL, j))
#    if (is.vector(y)) y <- matrix(y, nrow = 1, dimnames = list(NULL, names(y)))
    attr(y, "distr") <- attr(x, "distr")
    attr(y, "srad") <- attr(x, "srad")
    class(y) <- "ddSim"
    return(y)
  }
  if (nrow(x) == 1 || (!missing(i) && length(i) == 1))
    y <- matrix(y, nrow = 1, dimnames = list(NULL, names(y)))
  if (NCOL(y) == 1 & !noclass){
    nm <- names(y)
    y <- matrix(y, nrow = 1, dimnames = list(NULL, nm))
  }
  if (length(y) == 1) y <- as.vector(y)
  attr(y, "distr") <- attr(x, "distr")
  attr(y, "srad") <- attr(x, "srad")
  class(y) <- "ddSim"
  return(y)
}

#' Subset set of fitted dispersion models (ddArray)
#' @param x object to subset
#' @param distr vector of names of distributions to extract from fitted models
#' @param ... ignored (but required format for CRAN)
#' @details Subset the ddArray object as if it were a simple vector
#' @return list of selected models with crucial statistics
#' @export
"[.ddArray" <- function(x, distr){
  if (!is.vector(distr)) stop("subsetting ddArray objects must be by a vector ",
    "of integer array indices or distribution names")
  if (length(distr) == 1){
    output <- x[[distr]]
    class(output) <- "dd"
  } else {
    output <- list()
    if (is.numeric(distr)) distr <- names(x)[distr]
    for (di in distr) output[[di]] <- x[[di]]
    attr(output, "hidden") <- c("scVar", "ncarc", "ncarcCol", "rCol")
    attr(output, "expoCol") <- attr(x, "expoCol")
    attr(output, "rCol") <- attr(x, "rCol")
    attr(output, "scVar") <- attr(x, "scVar")
    attr(output, "ncarcCol") <- attr(x, "ncarcCol")
    if (length(output) == 1) {
      output <- output[[di]]
      class(output) <- "dd"
    } else {
      class(output) <- "ddArray"
    }
  }
  return(output)
}

#' Print dd and ddArray objects
#'
#' @description \code{dd} objects are lists consisting of all the elements of
#'  \code{glm} objects and several additional elements related to the distribution
#'  and its extensibility. Only a few of the elements are printed automatically.
#'  Others elements of object \code{x} can be viewed and extracted as with other
#'  lists in R, namely, by using the \code{x$element} or \code{x[[element]])
#'  operator, where \code{element} is the name of one of the elements of
#'  \code{x}which can be viewed via \code{names(x)}.
#'
#' @param x a \code{dd} or \code{ddArray} object
#' @param ... ignored
#'
#' @export
print.dd <- function(x, ...){
  cat(paste0("Distribution: ", x$distr, "\n"))
  cat(paste0("Formula: ", deparse(x$formula), "\n\n"))
  cat("Parameters:\n")
  print(x$parms)
  cat("\nCoefficients:\n")
  print(x$coefficients)
  cat("\nVariance:\n")
  print(x$varbeta)
  cat("\n")
  if(any(is.na(x$parms))){
  cat("Extensible", x$distr, "distribution not fittable.\n")
  cat("Cannot extrapolate beyond search radius.")
  cat("\nCoefficients for fitted glm are applicable within search radius only\n")
  }
}

#' @rdname dd
#' @export
print.ddArray <- function(x, ...){
  for (distr in names(x)) {
    print(x[[distr]])
    cat("****************************************\n")
  }
}

#'  Calculate CI for CDF, PDF, or quantile
#' @param mod a \code{dd} object
#' @param x distance from turbine (scalar or vector) or probability (for quantile)
#' @param type "CDF", "PDF", or "quantile"
#' @param CL confidence level for the confidence interval(s)
#' @param nsim number of simulation draws to base the estimate of CI on
#' @param na.tol maximum fraction of invalid parameter sets to discard when
#'  constructing CIs; abort if \code{mean(mod[, "extensible"]) > na.tol}
#' @return array (\code{ddCI} class) with columns for distance and the CI bounds
#' @export
ddCI <- function(mod, x, type = "CDF", CL = 0.9, nsim = 1000,
    extent = "full", zrad = 150, na.tol = 0.1){
  if (!"dd" %in% class(mod))
    stop("mod (", deparse(substitute(mod)), ") must be a dd object")
  p0 <- ddSim(mod, nsim = nsim, extent = extent)
  tmp <- array(dim = c(length(x), nsim))
  if (type == "CDF"){
    for (simi in 1:nsim)
      tmp[, simi] <- pdd(x, p0[simi, ], extent = extent, zrad = zrad)
  } else if (type == "PDF"){
    for (simi in 1:nsim)
      tmp[, simi] <- ddd(x, p0[simi, ], extent = extent, zrad = zrad)
  } else if (type == "quantile"){
    for (simi in 1:nsim)
      tmp[, simi] <- qdd(x, p0[simi, ], extent = extent, zrad = zrad)
  }
  CI <- cbind(x, matrixStats::rowQuantiles(tmp,
    probs = c((1 - CL)/2, (1 + CL)/2), na.rm = TRUE))
  if (anyNA(CI)){
    pNA <- mean(is.na(CI))
    if(pNA <= na.tol) {
      warning(paste0(100*round(pNA, 3), "% NAs in CIs"))
    } else {
      warning(paste0("cannot calculate meaningful CI for ", distr, ": ",
        100*round(pNA, 3), "% NAs in CIs"))
#      CI[, 2:3] <- NA
    }
  }
  class(CI) <- "ddCI"
  CI
}

#' Convert glm distance parameters into named distribution parameters
#' @param x object (vector or matrix of parameters, dmod, or glm)
#'  with named glm parameters ("r", "I(r^2)", "I(r^3)", "log(r)", or "I(1/r)")
#' @param distr name of the distribution
#' @param ... ignored [but required for CRAN]
#' @return matrix of parameters
#' @export
cof2parms <- function(x, ...) UseMethod("cof2parms", x)

#' @rdname cof2parms
#' @export
cof2parms.matrix <- function(x, distr, ...){
  ans <- suppressWarnings(switch(distr,
    gamma = {
      parms <- cbind(
        shape = x[, "log(r)"] + 2,
        rate = -x[, "r"])
      parms[parms[, "shape"] <= 0 | parms[, "rate"] <= 0, ]  <- NA
      parms
    },
    lognormal = {
      parms <- cbind(
        meanlog = -0.5 * (2 + x[, "log(r)"])/x[, "I(log(r)^2)"],
        sdlog = sqrt(-0.5/x[, "I(log(r)^2)"]))
      parms
    },
    logLinear = {
      parms <- x[, "r", drop = FALSE]
      parms[parms >= 0] <- NA
      colnames(parms) <- "b1"
      parms
    },
    logQuadratic = {
      parms <- x[, c("r", "I(r^2)"), drop = FALSE]
      colnames(parms) <- c("b1", "b2")
      parms[parms[, "b2"] >= 0, ] <- c(NA, NA)
      parms
    },
    logCubic ={
      parms <- x[, c("r", "I(r^2)", "I(r^3)"), drop = FALSE]
      colnames(parms) <- c("b1", "b2", "b3")
      parms[parms[, "b3"] >= 0, ] <- c(NA, NA, NA)
      parms
    },
    inverse_gamma = {
      parms <- cbind(
        shape = -x[, "log(r)"] - 2,
        scale = -x[, "I(1/r)"])
      parms[parms[, "shape"] <= 0 | parms[, "scale"] <= 0, ] <- c(NA, NA)
      parms
    },
    paranormal_gamma = {
      parms <- x[, c("log(r)", "r", "I(r^2)", "I(r^3)"), drop = FALSE]
      colnames(parms) <- c("b0", "b1", "b2", "b3")
      parms[parms[, "b3"] >= 0, ] <- c(NA, NA, NA, NA)
      parms
    },
    Rayleigh = {
      parms <- -0.5/x[, "I(r^2)", drop = FALSE]
      parms[parms <= 0] <- NA
      colnames(parms) <- "s2"
      parms
    },
    MaxwellBoltzmann = {
      parms <- sqrt(-0.5 /x[, "I(r^2)", drop = FALSE])
      colnames(parms) <- "a"
      parms
    },
    Pareto = {
      parms <- -x[, "log(r)", drop = FALSE] - 2
      parms[parms <= 0] <- NA
      colnames(parms) <- "a"
      parms
    },
    chisq = {
      parms <- 2 * x[, "log(r)", drop = FALSE] + 4
      parms[parms <= 0] <- NA
      colnames(parms) <- "df"
      parms
    },
    inverse_gaussian = {
      parms <- cbind(
        mean = sqrt(x[, "I(1/r)"]/x[, "r"]),
        dispersion = -0.5/x[, "I(1/r)"])
      parms[parms[, "mean"] <= 0 | parms[, "dispersion"] <= 0, ] <- NA
      parms
    },
    exponential = {
      parms <- -x[, "r", drop = FALSE]
      parms[parms <= 0] <- NA
      colnames(parms) <- "rate"
      parms
    },
    tnormal = {
      parms <- cbind(
        mean = -0.5 * x[, "r"]/x[, "I(r^2)"],
        sd = sqrt(-0.5/x[, "I(r^2)"]))
      parms
    },
    constant = c(b1 = NA)
  ))
  if (is.null(ans)) stop("distr must be a distribution name")
  if (is.matrix(ans) && nrow(ans) == 1) {
    nm <- colnames(ans)
    ans <- as.vector(ans)
    names(ans) <- nm
  }
  return(ans)
}

#' @rdname cof2parms
#' @export
cof2parms.numeric <- function(x, distr, ...){ # vector
  output <- cof2parms(matrix(x, nrow = 1, dimnames = list(NULL, names(x))), distr)
  return(output)
}

#' @rdname cof2parms
#' @export
cof2parms.dd <- function(x, ...){
  return(cof2parms(x$coef, x$distr))
}

#' Function to mimic dnorm, dexp, etc. for \code{dd} models
#' @param x numeric x >= 0
#' @param model either a \code{dd} object or the name of a distance model
#'  (e.g., "gamma").
#' @param parms scalar, vector, or array of parameters. Or if \code{x} is a
#'  \code{dd} object, then \code{parms} may be NULL, in which case, \code{ddd}
#'   uses the MLE parameters found in the \code{dd} object.
#' @param srad the radius to use for truncated distributions
#'  (\code{srad = NULL} to use full, non-truncated distribution)
#' @return vector of PDF values. If \code{srad = NULL} and a full distribution
#'  is used, the length of the output is equal to the maximum of the length of
#'  \code{x} and the number of parameter sets are included in \code{model} (i.e.
#'  1 if \code{model} is a \code{dd} object, or \code{nrow(model)} if
#'  \code{model} is a \code{ddSim} object. If a truncated distribution is to be
#'  used, then \code{length(x) = 1} or the number of parameter sets must be 1
#'  (which occurs if \code{model} is a \code{dd} object or \code{nrow(model) = 1}
#'  if \code{model} is a \code{ddSim} object.
#' @export
ddd <- function(x, model, extent = "full", zrad = 150){ # model is either ddSim or dd
  if("dd" %in% class(model)) model <- dd2ddSim(model)
  if ("ddSim" %in% class(model)){
    parms <- model
    distr <- attr(parms, which = "distr")
    srad <- attr(parms, which = "srad")
  } else {
    stop("ddd: model must be dd or ddSim object")
  }
  if (extent == "full"){ # intent to integrate to Inf; use zrad if parms improper
    # output has dimensions length(x) x nrow(parms) [or a vector]
    # calculations must be split into two parts:
    #  1) proper parms integrate to Inf
    #  2) improper parms are integrated to zrad
    output <- matrix(0, nrow = length(x), ncol = nrow(parms))
    i0 <- which(parms[, "extensible"] == 0)
    i1 <- which(parms[, "extensible"] > 0)
    if(length(i1) > 0){ # extensible parms
      xx <- rep(x, length(i1))
      ppi <- rep(i1, each = length(x))
      output[, i1] <- switch(distr,
        gamma = matrix(
          dgamma(xx, shape = parms[ppi, "shape"], rate = parms[ppi, "rate"]),
          nrow = length(x)
        ),
        lognormal = matrix(
          dlnorm(xx, meanlog = parms[ppi, "meanlog"], sdlog = parms[ppi, "sdlog"]),
          nrow = length(x)
        ),
        logLinear = matrix(dlogL(xx, b1 = parms[ppi, "b1"]),
          nrow = length(x)
        ),
        logQuadratic = matrix(dlogQ(xx, b1 = parms[ppi, "b1"], b2 = parms[ppi, "b2"]),
          nrow = length(x)
        ),
        logCubic = {
          const <- numeric(length(i1))
          tmp <- matrix(0, nrow = length(x), ncol = length(i1))
          for (i in i1){
            const <- tryCatch(1/integrate(
              f = function(r)
                x * exp(b1 * x + b2 * x^2 + b3 * x^3),
              lower = 0, upper = Inf, rel.tol = .Machine$double.eps^0.5)$val,
              error = function(e) NA
            )
            tmp[, which(i1 == i)] <- dlogC(x,
              b1 = c(parms[i, "b1"]), b2 = c(parms[i, "b2"]), b3 = c(parms[i, "b3"]),
              const = const
            )
          }
          tmp
        },
        inverse_gamma = matrix(
          digam(xx, shape = parms[ppi, "shape"], scale = parms[ppi, "scale"]),
          nrow = length(x)
        ),
        paranormal_gamma = { # assumes b0, b1, b2, b3 are scalars
          const <- numeric(length(i1))
          tmp <- matrix(0, nrow = length(x), ncol = length(i1))
          for (i in i1){
            b0 <- c(parms[i, "b0"])
            b1 <- c(parms[i, "b1"])
            b2 <- c(parms[i, "b2"])
            b3 <- c(parms[i, "b3"])
            const <- tryCatch(1/integrate(f = function(x)
              x * exp(b0 * log(x) + b1 * x + b2 * x^2 + b3 * x^3),
              lower = 0, upper = Inf, rel.tol = .Machine$double.eps^0.5)$val,
              error = function(e) NA
            )
            tmp[, which(i1 == i)] <- dpng(x,
              b0 = c(parms[i, "b0"]),
              b1 = c(parms[i, "b1"]),
              b2 = c(parms[i, "b2"]),
              b3 = c(parms[i, "b3"]),
              const = const
            )
          }
          tmp
        },
        Rayleigh = matrix(dRay(xx, s2 = parms[ppi, "s2"]), nrow = length(x)),
        MaxwellBoltzmann = matrix(dmb(xx, a = parms[ppi, "a"]), nrow = length(x)),
        constant = {
          tmp <- 2 * xx/zrad^2
          tmp[xx > zrad] <- 0
          matrix(tmp, nrow = length(x))
        },
        Pareto = matrix(dPare1(xx, a = parms), nrow = length(x)),
        tnormal = {
          mu <- parms[ppi, "mean"]
          sig  <- parms[ppi, "sd"]
          Fa <- pnorm(0, mean = mu, sd = sig)
          ans <- numeric(length(xx))
          ans[xx >= 0] <- dnorm(xx[xx >= 0], mean = mu, sd = sig)/(1 - Fa)
          matrix(ans, nrow = length(x))
        },
        exponential = matrix(dexp(xx, rate = parms[ppi, "rate"]), nrow = length(x)),
        inverse_gaussian = matrix(statmod::dinvgauss(xx,
            mean = parms[ppi, "mean"], dispersion = parms[ppi, "dispersion"]),
          nrow = length(x)
        ),
        chisq = matrix(dchisq(xx, df = parms[ppi, "df"]), nrow = length(x))
      )
    }
    if (length(i0) > 0){
      for(i in i0){
        if(distr == "constant") output[, i] <- 2 * x/zrad^2
        if(distr == "paranormal_gamma"){
          if(parms[i, "log(r)"] <= -2){
            output[, i] <- NA
            next
          }
        }
        if(distr == "gamma" && parms[i, "log(r)"] <= -2){
          output[, i] <- NA
          next
        }
        if(distr == "inverse_gamma" && parms[i, "log(r)"] > -2){
          output[, i] <- NA
          next
        }
        deno <- integrate(f = function(r)
          exp(rmat(r, distr) %*% parms[i, cof_name[[distr]]] + off(r, distr)),
          lower = 0,
          upper = zrad
        )$val
        output[x > 0, i] <-
          exp(rmat(x[x > 0], distr) %*% parms[i, cof_name[[distr]]] + off(x[x > 0], distr))/deno
      }
      output[x > zrad, i0] <- 0
    }
  } else if (extent == "win"){ # truncation radius is numeric scalar 0 < srad < Inf
    output <- matrix(0, nrow = length(x), ncol = nrow(parms))
    for (i in 1:nrow(parms)){
      deno <- tryCatch(integrate(f = function(r)
          exp(rmat(r, distr) %*% parms[i, cof_name[[distr]]] + off(r, distr)),
          lower = 0,
          upper = srad
        )$val,
        error = function(e) NA
      )
      if (!is.na(deno)){
        output[x > 0, i] <-
          exp(rmat(x[x > 0], distr) %*% parms[i, cof_name[[distr]]] + off(x[x > 0], distr))/deno
      } else {
        output[, i] <- NA
      }
    }
    output[x > srad, ] <- 0
  } else {
    stop("extent must be \"full\" or \"win\"")
  }
  if (nrow(parms) == 1 | length(x) == 1) output <- unname(as.vector(output))
  output
}

#' Function to mimic pnorm, pexp, etc. but
#' @param x numeric x >= 0
#' @param model either a \code{dd} object or the name of a distance model
#'  (e.g., "gamma").
#' @param parms scalar, vector, or array of parameters. Or if \code{x} is a
#'  \code{dd} object, then \code{parms} may be NULL, in which case, \code{ddd}
#'  uses the MLE parameters found in the \code{dd} object.
#' @param extent "full" or "win" to account for all carcasses or only those
#'  within the search radius.
#' @return answer
#' @rdname ddd
#' @export

pdd <- function(x, model, extent = "full", zrad = 150){ # model is either ddSim or dd
  if("dd" %in% class(model)) model <- dd2ddSim(model)
  if ("ddSim" %in% class(model)){
    distr <- attr(model, which = "distr")
  } else {
    stop("pdd: model must be dd or ddSim object")
  }
  parms <- model
  if (extent == "full"){ # intent to integrate to Inf; use zrad if parms improper
    # output has dimensions length(x) x nrow(parms) [or a vector]
    # calculations must be split into two parts:
    #  1) proper parms integrate to Inf
    #  2) improper parms are integrated to zrad
    output <- matrix(0, nrow = length(x), ncol = nrow(parms))
    i0 <- unname(which(parms[, "extensible"] == 0))
    i1 <- unname(which(parms[, "extensible"] > 0))
    if(length(i1) > 0){ # extensible parms
      xx <- rep(x, length(i1))
      ppi <- rep(i1, each = length(x))
      output[, i1] <- switch(distr,
        gamma = matrix(
          pgamma(xx, shape = parms[ppi, "shape"], rate = parms[ppi, "rate"]),
          nrow = length(x)
        ),
        lognormal = matrix(
          plnorm(xx, meanlog = parms[ppi, "meanlog"], sdlog = parms[ppi, "sdlog"]),
          nrow = length(x)
        ),
        logLinear = matrix(plogL(xx, b1 = parms[ppi, "b1"]),
          nrow = length(x)
        ),
        logQuadratic = matrix(plogQ(xx, b1 = parms[ppi, "b1"], b2 = parms[ppi, "b2"]),
          nrow = length(x)
        ),
        logCubic = {
          const <- numeric(length(i1))
          tmp <- matrix(0, nrow = length(x), ncol = length(i1))
          for (i in i1){
            const <- tryCatch(1/integrate(
              f = function(r)
                (exp(rmat(r, distr)[, -1] %*% c(parms[i, cof_name[[distr]]][-1]) + 
                  off(r, distr))),
              lower = 0, upper = Inf, rel.tol = .Machine$double.eps^0.5)$val,
              error = function(e) NA
            )
            tmp[, which(i1 == i)] <- plogC(x,
              b1 = c(parms[i, "b1"]), b2 = c(parms[i, "b2"]), b3 = c(parms[i, "b3"]),
              const = const
            )
          }
          tmp
        },
        inverse_gamma = matrix(
          pigam(xx, shape = parms[ppi, "shape"], scale = parms[ppi, "scale"]),
          nrow = length(x)
        ),
        paranormal_gamma = { # assumes b0, b1, b2, b3 are scalars
          const <- numeric(length(i1))
          tmp <- matrix(0, nrow = length(x), ncol = length(i1))
          for (i in i1){
            b0 <- c(parms[i, "b0"])
            b1 <- c(parms[i, "b1"])
            b2 <- c(parms[i, "b2"])
            b3 <- c(parms[i, "b3"])
            const <- tryCatch(1/integrate(f = function(x)
              x * exp(b0 * log(x) + b1 * x + b2 * x^2 + b3 * x^3),
              lower = 0, upper = Inf, rel.tol = .Machine$double.eps^0.5)$val,
              error = function(e) NA
            )
            tmp[, which(i1 == i)] <- ppng(x,
              b0 = b0,
              b1 = b1 ,
              b2 = b2,
              b3 = b3,
              const = const
            )
          }
          tmp
        },
        Rayleigh = matrix(pRay(xx, s2 = parms[ppi, "s2"]), nrow = length(x)),
        MaxwellBoltzmann = matrix(pmb(xx, a = parms[ppi, "a"]), nrow = length(x)),
        constant = {
          tmp <- (xx/zrad)^2
          tmp[xx > zrad] <- 1
          matrix(tmp, nrow = length(x))
        },
        Pareto = matrix(pPare1(xx, a = parms), nrow = length(x)),
        tnormal = {
          mu <- parms[ppi, "mean"]
          sig  <- parms[ppi, "sd"]
          Fa <- pnorm(0, mean = mu, sd = sig)
          ans <- numeric(length(xx))
          ans[xx >= 0] <- (pnorm(xx[xx >= 0], mean = mu, sd = sig) - Fa)/(1 - Fa)
          matrix(ans, nrow = length(x))
        },
        exponential = matrix(pexp(xx, rate = parms[ppi, "rate"]), nrow = length(x)),
        inverse_gaussian = matrix(statmod::pinvgauss(xx,
            mean = parms[ppi, "mean"], dispersion = parms[ppi, "dispersion"]),
          nrow = length(x)
        ),
        chisq = matrix(pchisq(xx, df = parms[ppi, "df"]), nrow = length(x))
      )
    }
    if (length(i0) > 0){
      for(i in i0){
        if(distr == "constant") output[, i] <- (x/zrad)^2
        if(distr == "paranormal_gamma"){
          if(parms[i, "log(r)"] <= -2){
            output[, i] <- NA
            next
          }
        }
        if(distr == "gamma" && parms[i, "log(r)"] <= -2){
          output[, i] <- NA
          next
        }
        if(distr == "inverse_gamma" && parms[i, "log(r)"] > -2){
          output[, i] <- NA
          next
        }
        deno <- tryCatch(
          integrate(f = function(r)
            exp(rmat(r, distr) %*% parms[i, cof_name[[distr]]] + off(r, distr)),
            lower = 0,
            upper = zrad, rel.tol = .Machine$double.eps^0.5
          )$val,
          error = function(e) NA
        )
        if (any(x > 0)){
          for (xi in which(x > 0 & x < zrad)){
            output[xi, i] <- tryCatch(
              integrate(function(r)
                exp(rmat(r, distr) %*% parms[i, cof_name[[distr]]] + off(r, distr)),
                lower = 0,
                upper = x[xi], rel.tol = .Machine$double.eps^0.5
              )$val/deno,
              error = function(e) NA
            )
          }
        }
      }
      output[x >= zrad, i0] <- 1
    }
  } else if (extent == "win"){ # truncation radius is numeric scalar 0 < srad < Inf
    srad <- attr(model, "srad")
    output <- matrix(0, nrow = length(x), ncol = nrow(parms))
    if (all(x <= 0)){
      if (nrow(output) == 1 | ncol(output) == 1) output <- unname(as.vector(output))
      return(output)
    }
    for (i in 1:nrow(parms)){
      deno <- tryCatch(integrate(f = function(r)
          exp(rmat(r, distr) %*% parms[i, cof_name[[distr]]] + off(r, distr)),
          lower = 0,
          upper = srad
        )$val,
        error = function(e) NA
      )
      if (!is.na(deno)){
        for (xi in which(x > 0 & x <= srad)){
          output[xi, i] <- integrate(function(r)
            exp(rmat(r, distr) %*% parms[i, cof_name[[distr]]] + off(r, distr)),
            lower = 0,
            upper = x[xi]
          )$val/deno
        }
      } else {
        output[, i] <- NA
      }
    }
    output[x > srad, ] <- 1
  } else {
    stop("extent must be \"full\" or \"win\"")
  }
  if (nrow(parms) == 1 | length(x) == 1) output <- unname(as.vector(output))
  output[output > 1] <- 1
  output[output < 0] <- 0
  output
}

#' Function to calculate relative carcass densities
#' @param x numeric x >= 0
#' @param model either a \code{dd} object or a set of simulated parameters
#'  (i.e., a \code{ddSim} object)
#' @param srad the radius to use for truncated distributions
#'  (\code{srad = NULL} to use full, non-truncated distribution)
#' @details \code{rcd} gives relative carcass densities at a point \code{x}
#'  meters from a turbine. In general, rcd is proportional to PDF(x)/x,
#'  normalized so that the surface of rotation of rcd(x) has total area of 1.
#'  There are more stringent contstraints on the allowable parameters in the
#'  fitted (or simulated) glm's because the integral of PDF(x)/x must converge.
#' @return vector or array of \code{rcd} values for every parameter set and
#'  distance. If there is only one parameter set (e.g., if \code{model} is a
#'  \code{dd} object), then the function returns a vector representing the
#'  \code{rcd} for each distance in \code{x}. If \code{model} is a \code{ddSim}
#'  object, then the function returns an array, with each column representing
#'  the vector \code{rcd(x)} values for all distances for a single parameter set,
#'  and the rows representing the \code{rcd} values for all parameter sets at a
#'  given distance. NOTE: \code{PDF(0)/0} is dengenerate in R, and \code{rcd}
#'  approximates x < 0.001 with x = 0.001 to streamline the calculations.
#' @rdname ddd
#' @export
rcd <- function(x, model, extent = "full", zrad = 150){
  x[x < 0.001] <- 0.001
  if("dd" %in% class(model)) model <- dd2ddSim(model)
  if ("ddSim" %in% class(model)){
    distr <- attr(model, which = "distr")
  } else {
    stop("pdd: model must be dd or ddSim object")
  }
  parms <- model
  if (distr == "constant" & extent == "full")
    warning("rcd cannot be calculated for extent = \"full\" and distr = ", distr)
  if (distr %in% c("exponential", "tnormal"))
    warning("rcd cannot be calculated for distr = ", distr)
  upr <- ifelse(extent == "full", Inf, srad)
  output <- array(dim = c(length(x), nrow(parms)))
  for (pri in 1:nrow(parms)){
    output[, pri] <- as.vector(ddd(x = x, model = parms[pri, ], extent = extent, zrad = zrad))
    const <- tryCatch(1/integrate(f = function(x)
      2 * pi * c(ddd(x, model = parms[pri,], extent = extent, zrad = zrad))/x,
      lower = 0, upper = upr)$val,
      error = function(e) NA)
    output[, pri] <- const * output[, pri]/x
  }
  if (length(x) == 1 | nrow(parms) == 1) output <- as.vector(output)
  return(output)
}

#' Display a Table of Summary Statistics
#'
#' @param x list of models (\code{ddArray}) or single model (\code{dd}) to calculate summary statistics for
#' @param extent distributions within searched area (\code{"win"}) or extended beyond (\code{"full"})
#' @param zrad maximum distance that carcasses can lie
#'  (only used when glm parameters not extensible to Inf)
#' @param ... ignored
#' @return list (or list of lists if \code{x} is \code{ddArray}) with \code{$model}
#'  giving the model parameters and \code{$stats} giving the median, and 75th,
#'  90th, and 95th quantiles of carcass distances
#' @export
stats <- function(x, ...) UseMethod("stats", x)

#' @rdname stats
#' @export
stats.dd <- function(x, extent = "full", zrad = 150, ...){
  distr <- x$distr
#  if (distr == "constant" & extent == "full") return(list(model = NA, stats = NA))
  cof <- numeric(length(cof_name[[distr]])) + NA
  names(cof) <- cof_name[[distr]]
  qtls <- c(0.5, 0.75, 0.9, 0.95)
  names(qtls) <- c("median", "75%", "90%", "95%")
  ans <- list(
    model = data.frame(array(dim = c(1, 2 * length(cof)))),
    stats = numeric(5) + NA
  )
  rownames(ans$model) <- distr
  colnames(ans$model) <- c(parm_name[[distr]], names(cof), "offset")
  ans[["model"]][1, "offset"]  <- mod_offset[distr]
  names(ans[["stats"]]) <- c(names(qtls), "mode")
  cof <- x$coefficients[cof_name[[distr]]]
  parms <- x$parms
  ans[["model"]][1, names(c(cof, parms))] <- c(cof, parms)
  if (any(is.na(x$parms)) & extent == "full") return(ans)
  ans[["stats"]][names(qtls)] <- qdd(qtls, model = x, extent = extent, zrad = zrad)
  xmax <- 1000
  xseq <-seq(0.1, xmax, by = 0.1)
  y <- ddd(x = xseq, model = x, extent = extent, zrad = zrad)
  mi <- which(y == max(y))
  ans[["stats"]]["mode"] <- ifelse(y[length(y)] %in% mi, NA, mean(xseq[mi]))
  if(ans[["stats"]]["mode"] == 0.1) ans[["stats"]]["mode"] <- 0
  ans[["stats"]] <- round(ans[["stats"]], 2)
  ans
}

#' @export
stats.ddArray <- function(x, extent = "full", zrad = 150, ...){
  aic0 <- aic(x, extent = extent)
  rnm <- aic0$model
  cnm<- c("median", "75%", "90%", "95%", "mode", "deltaAICc")
  ans <- data.frame(array(dim = c(length(rnm), length(cnm))))
  rownames(ans) <- rnm
  colnames(ans) <- cnm
  ans$deltaAICc <- aic0$deltaAICc
  for (distr in aic0$model)
    ans[distr, c("median", "75%", "90%", "95%", "mode")] <-
      stats(x[distr], extent = extent, zrad = zrad)$stats
  ans
}

#' @export
dd2ddSim <- function(dd){
  if (!"dd" %in% class(dd)) stop("dd2ddSim takes \"dd\" objects only")
  parms <- matrix(
    c(dd$coefficients, dd$parms, 1 * cofOK(dd$coefficients, dd$distr)),
    nrow = 1,
    dimnames = list(NULL,
      c(names(dd$coefficients), names(dd$parms), "extensible")
    )
  )
  attr(parms, "srad") <- max(dd$data$r, na.rm = TRUE)
  attr(parms, "distr") <- dd$distr
  class(parms) <- "ddSim"
  parms
}

#' @rdname ddd
#' @export
qdd <- function(p, model, extent = "full", zrad = 150, subdiv = 1000){ # model is ddSim or dd
# find r such that pdd(r) = p
  if("dd" %in% class(model)) model <- dd2ddSim(model)
  if (!"ddSim" %in% class(model)) stop("qdd: model must be dd or ddSim object")
  srad <- attr(model, "srad")
  distr <- attr(model, "distr")
  if (nrow(model) == 1){
    if ((distr %in% c("gamma", "paranormal_gamma") && model[, "log(r)"] <= -2) |
      (distr %in% c("inverse_gamma", "Pareto") && model[, "log(r)"] >= -2)){
      return(rep(NA, length(p)))
    }
    if (length(p) == 1){
      qval = tryCatch(
        uniroot(
          f = function(r) pdd(r, model = model, extent = extent, zrad = zrad) - p,
          lower = 0.1, upper = ifelse(extent == "full", zrad, srad)
        )$root,
        error = function(e) NA
      )
      if(is.na(qval)){
        r0 <- zrad
        while(1){
          r0 <- 2 * r0
          if (pdd(r0, model = model, extent = extent, zrad = zrad) > p){
            qval <- uniroot(
              f = function(r)
                pdd(r, model = model, extent = extent, zrad = zrad) - min(p),
              lower = r0/2, upper = r0
            )$root
            break
          }
        }
      }
      return(qval)
    } else if (length(p) > 1){
    # p is a vector, solve by linear interpolation between minr and maxr

      minr <- tryCatch(uniroot( # try to find the root
          f = function(r) pdd(r, model = model, extent = extent) - min(p),
          lower = 0,
          upper = ifelse(extent == "full", zrad, srad -1)
        )$root,
        error = function(e) NA
      )
      if (is.na(minr)){
        # if unsuccessful in finding the root, expand the search interval
        r0 <- zrad
        while(1){
          r0 <- 2 * r0
          if (pdd(r0, model = model, extent = extent) > min(p)){
            minr <- uniroot(
              f = function(r)
                pdd(r, model = model, extent = extent) - min(p),
              lower = r0/2, upper = r0
            )$root
            break
          }
        }
      }
      maxr <- tryCatch(
        uniroot(
          f = function(r)
            pdd(r, model = model, extent = extent, zrad = zrad) - max(p),
          lower = 0, upper = ifelse(extent == "full", zrad, srad)
          )$root,
        error = function(e) NA
      )
      if (is.na(maxr)){
        r0 <- zrad
        while(1){
          r0 <- 2 * r0
          if (pdd(r0, model = model, extent = extent) > max(p)){
            maxr <- uniroot(
              f = function(r)
                pdd(r, model = model, extent = extent, zrad = zrad) - max(p),
              lower = r0/2, upper = r0
            )$root
            break
          }
        }
      }
      if ("try-error" %in% class(maxr))
        stop(paste0("cannot invert CDF in qdd for p = ", max(p)))
      y <- seq(minr, maxr, length.out = subdiv)
      x <- pdd(x = y, model = model, extent = extent)
      y <- c(0, y, min(50000, srad))
      x <- c(0, x, 1)
      return(approxfun(x, y)(p))
    }
  } else if (nrow(model) > 1){
    pp <- rep(p, length.out = nrow(model))
    qval <- numeric(nrow(model))
    for (ppi in 1:length(pp)){
      qval[ppi] = tryCatch(
        uniroot(f = function(r)
          pdd(r, model = model[ppi, ], extent = extent, zrad = zrad) - pp[ppi],
          lower = 0, upper = ifelse(extent == "full", 50000, srad))$root,
        error = function(e) NA
      )
    }
    return(qval)
  } else {
    stop("dimension mismatch between p and simulated parameters (= model) in qdd")
  }
}

#' @export
rdd <- function(n, model, extent = "full", zrad = 150, subdiv = 1000){ # model is ddSim or dd
  if ("dd" %in% class(model)) parms <- dd2ddSim(model)
  if ("ddSim" %in% class(model)){
    if (nrow(model) > 1 & nrow(model) != n)
      stop("nrow(model) must be 1 or n in rdd")
  } else {
    stop ("class(model) in rdd must be dd or ddSim")
  }
  qdd(runif(n), model = model, extent = extent, zrad = zrad, subdiv = subdiv)
}
