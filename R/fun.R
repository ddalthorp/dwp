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
#'  probabilities) and turbines. Search classes are listed the \code{scCol}
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
#' @return \code{siteLayout} object = list with spatial characteristics of the site
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
#' @param scCol (optional) name of a column with search class delineator in
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
#'    \item{\code{$scCol}}{name of the search class column (if any);}
#'    \item{\code{$X}}{number of carcasses found at each turbine (vector}
#'  }
#' @export
prepRing <- function(siteLayout, scCol = NULL, notSearched = NULL){
  if (!is.null(scCol)){
    if (!scCol %in% names(siteLayout$layout))
      stop(scCol, " not in siteLayout$layout")
    if (!is.null(notSearched)){ # remove unsearched areas from data
      if (!is.character(notSearched))
        stop("notSearched must be a character vector (or NULL)")
      siteLayout$layout <-  siteLayout$layout[
        !siteLayout$layout[, scCol, drop = TRUE] %in% notSearched, ]
      siteLayout$layoutAdj <- siteLayout$layoutAdj[
        !siteLayout$layoutAdj[, scCol, drop = TRUE] %in% notSearched, ]
      if (scCol %in% names(siteLayout$carcasses))
        siteLayout$carcasses <-  siteLayout$carcasses[
          !siteLayout$carcasses[, scCol, drop = TRUE] %in% notSearched, ]
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
  scset <- gtools::mixedsort(unique(as.data.frame(siteLayout$layoutAdj)[, scCol]))
  for (sci in scset){
    trscArea[[sci]] <- matrix(0,
      nrow = length(siteLayout$tset), ncol = dim(rings)[1])
    rownames(trscArea[[sci]]) <- siteLayout$tset
    for (ti in siteLayout$tset){
      # indices for turbine = ti and svar = sci
      ind <- which(siteLayout$layout[, unitCol, drop = TRUE] == ti &
        siteLayout$layout[, scCol, drop = TRUE] == sci)
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
  names(ringData) <- c("r", "svar", "exposure", "ncarc")
  ringData$r <- rep(radi, length(scset))
  ringData[, "svar"] <- rep(scset, each = length(radi))
  ringData$exposure <- unlist(trscArea)
  tmp <- suppressWarnings(sf::st_intersection(
    siteLayout$layout, siteLayout$carcasses))
  carcassSumry <- tmp[, c("Turbine", "Class")]
  carcr <- sqrt(rowSums((sf::st_coordinates(siteLayout$carcasses) -
    siteLayout$tcenter[siteLayout$carcasses[, unitCol, drop = TRUE], ])^2))
  for (ci in 1:nrow(carcassSumry)){
    ind <- which(ringData$r == ceiling(carcr[ci]) &
      ringData$svar == unlist(as.matrix(carcassSumry)[ci, scCol]))
    ringData[ind, "ncarc"] <- ringData[ind, "ncarc"] + 1
  }
  rd <- ringData[ringData$exposure > 1e-3 & !ringData$svar %in% notSearched,]
  return(list(ringData = rd, dwpArea = dwpArea, searchMax = searchMax,
    scCol = scCol, X = siteLayout$X))
}

#' Fit an array of candidate distance models
#'
#' @description Fit glm's for canonical distance models.
#'
#' Models are fit for glm's corresponding to gamma, lognormal, log-linear,
#' log-quadratic, log-cubic, truncated normal, exponential, Rayleigh, Pareto,
#' inverse gamma, paranormal gamma, Maxwell Boltzmann, chi-squared, and inverse
#' Gaussian distributions. These models apply strictly to the carcasses within
#' the search radius. Extrapolation outside the search radius requires the
#' additional step of fitting truncated distributions (\code{tfit.dmod()})
#'
#' @param ringData Data frame with a tally of carcasses within 1m rings, along
#'  with outer radii of the rings, (optional) search class column, and the amount
#'  of area searched in each search class and ring.
#' @param rCol name of column in \code{ringData} with outer radius for 1m rings
#' @param scCol name of column in \code{ringData} with (optional) subdivisions
#'  of rings into search classes (character strings)
#' @param expoCol name of "exposure" column of total area in each ring belonging
#'  to each search class.
#' @param modelSet names (vector of character strings) of glm distribution templates
#'  to fit. Default is \code{"all"}, which fits \code{gamma}, \code{lognormal},
#'  \code{logLinear}, \code{logQuadratic}, \code{logCubic}, \code{inverse_gamma},
#'  \code{paranormal_gamma}, \code{Rayleigh}, \code{MaxwellBoltzmann}, \code{constant},
#'  \code{tnormal}, \code{exponential}, \code{Pareto}, \code{chisq}, and
#'  \code{inverse_gaussian}. Any subset of these may also be fit with a single
#'  call to \code{canFit}.
#' @return list of fitted glm models as a \code{dmod} objects, which are standard
#'  \code{glm} objects with the name of the fitted distribution template attached
#'  as \code{$distr} and the name of the (optional) search class column name
#'  attached as \code{$scCol}. In addition, the data in the \code{glm} are
#'  coded as "ncarc" for number of carcasses, "r" for distance from turbine,
#'  "exposure" for amount of area in the given search class at the given distance.
#'
#' @export
canFit <- function(ringData, rCol = "r", scCol = NULL, expoCol = "exposure",
  ncarcCol = "ncarc", modelSet = "all"){
  dat <- ringData[, c(rCol, scCol, expoCol, ncarcCol)]
  names(dat) <- c("r", if(!is.null(scCol)) "svar", "exposure", "ncarc")
  scl <- scCol
  plu <- ifelse(!is.null(scCol), " + ", " ")
  sufx <- paste0(plu, scCol,  "+ offset(log(", expoCol, "))")
  sufx0 <- paste0(plu, scCol, "+ ")
  pref <- paste0("ncarc ~ ")
  dmod <- list()
  if (modelSet == "all") modelSet <- mod_name
  for (distr in modelSet){
    dmod[[distr]] <- switch(distr,
      gamma = glm(formula(paste0(pref, "log(r) + r", sufx)),
        data = dat, family = "poisson"),
      lognormal = glm(formula(paste0(pref, "log(r) + I(log(r)^2)", sufx)),
        data = dat, family = "poisson"),
      logLinear = glm(formula(paste0(pref, "r", sufx)),
        data = dat, family = "poisson"),
      logQuadratic = glm(formula(paste0(pref, "r + I(r^2)", sufx)),
        data = dat, family = "poisson"),
      logCubic = glm(formula(paste0(pref, "r + I(r^2) + I(r^3)", sufx)),
        data = dat, family = "poisson"),
      inverse_gamma = glm(formula(paste0(pref, "I(1/r) + log(r)", sufx)),
        data = dat, family = "poisson"),
      paranormal_gamma = glm(formula(paste0(pref, "log(r) + r + I(r^2) + I(r^3)", sufx)),
        data = dat, family = "poisson"),
      Rayleigh = glm(formula(paste0(pref, "I(r^2)", sufx)),
        data = dat, family = "poisson"),
      Pareto = glm(formula(paste0(pref, "log(r)", sufx)),
        data = dat, family = "poisson"),
      MaxwellBoltzmann = glm(formula(
        paste0(pref, "I(r^2)", sufx0, "offset(I(log(exposure * r)))")),
        data = dat, family = "poisson"),
      constant = glm(formula =
        formula(ifelse(is.null(scCol),
          paste0(pref, "offset(log(", expoCol, "))"),
          paste0(pref, scCol, " + offset(log(", expoCol, "))"))),
        data = dat, family = "poisson"),
      tnormal = glm(formula(paste0(pref, "r + I(r^2)", sufx0,
        "offset(I(log(exposure/r)))")),
        data = dat, family = "poisson"),
      exponential = glm(formula(paste0(pref, "r", sufx0,
        "offset(I(log(exposure/r)))")),
        data = dat, family = "poisson"),
      inverse_gaussian = glm(formula(paste0(pref, "I(1/r) + r", sufx0,
        "offset(I(log(exposure) - 5/2 * log(r)))")),
        data = dat, family = "poisson"),
      chisq = glm(formula(paste0(pref, "log(r)", sufx0,
        "offset(I(log(exposure) - r/2))")),
        data = dat, family = "poisson")
    )
  }
  for (distr in names(dmod)){
    dmod[[distr]]$distr <- distr # distr solidly attached to its model
    dmod[[distr]]$scCol <- scCol
    class(dmod[[distr]]) <- c("dmod", class(dmod[[distr]]))# glm w/distr and dat
  }
  return(dmod)
}

#' Function for fitting glm's and distance distributions to carcass dispersions
#'
#' The function first fits the glm models for carcass dispersion within the
#' search radius. Then it uses truncated glm to fit carcass distributions that
#' may be extrapolated (with appropriate caution) beyond the search radius. The
#' fitted models can then be used to estimate dwp or to research carcass dispersion
#' patterns in general.
#'
#' @param ringData data frame with a simplified representation of the seach areas
#'  and carcass locations. Columns include the outer radii of concentric rings
#'  around turbine centers (\code{ncarcCol}), the number of carcasses in each
#'  ring (\code{ncarcCol}), an optional column for search class (\code{scCol}),
#'  and the amount of area or "exposure" for each search class in each ring
#'  (\code{expoCol}.
#' @param scCol name of the (optional) search class column. Search classes typically
#'  represent ground with different search characteristics resulting in different
#'  detection probabilities. Search classes will affect the number of carcasses
#'  observed (via differences in detection probabilities) but are assumed not to
#'  affect the distribution of distances at which carcasses land after getting
#'  struck by turbines.
#' @param rCol name of the column representing the outer radii of concentric rings
#'  for which carcass counts, search exposure, and search class are tallied.
#' @param expoCol name of the column representing the amount of area in each ring
#'  and (optional) search class.
#' @param ncarcCol name of the column with the number of carcasses observed in
#'  each ring and (optional) search class.
#' @param modelSet a vector of names (character strings) with the names of the
#'  distributions to be fit. Default is \code{"all"} which fits \code{gamma},
#'  \code{lognormal}, \code{logLinear}, \code{logQuadratic}, \code{logCubic},
#'  \code{inverse_gamma}, \code{paranormal_gamma}, \code{Rayleigh},
#'  \code{MaxwellBoltzmann}, \code{constant}, \code{tnormal}, \code{exponential},
#'  \code{Pareto}, \code{chisq}, and \code{inverse_gaussian}. Any subset of
#'  these may also be fit with a single call to \code{ddArray}.
#' @return a \code{ddArray} object, which is a list of \code{dd} objects, one for
#'  each fitted model. Each \code{dd} object consists of the following components:
#'  \describe{
#'    \item{\code{$distr}}{name of the fitted distribution}
#'    \item{code{$parms}}{MLE for parameter vector for the fitted distribution
#'      (in cases where the model can be extended beyond the search radius) or
#'      NA (in cases where the model cannot be extended). \code{$parms} is a
#'      with names corresponding to distribution parameter names (e.g., for a
#'      truncated normal distribution, the vector would be something like
#'       \code{c(mean = 28.3, sd = 12.3)}.}
#'    \item{\code{$beta}}{MLE for parameter vector for the fitted truncated glm
#'      that may be extended beyond the search radius.}
#'    \item{\code{$varbeta}}{estimated variance of estimated \code{$beta} vector.
#'      This variance matrix is used in simulating parameters and generating
#'      confidence intervals.}
#'    \item{\code{$ncarc}}{total number of carcasses observed}
#'    \item{\code{$n}}{number of rings}
#'    \item{\code{$k}}{number of parameters estimated}
#'    \item{\code{$llik}}{the log-likelihood for the fitted MLE parameters}
#'    \item{\code{$dmod}}{\code{dmod} object for the base, non-extendable model,
#'      applicable only to the carcasses within the search radius. The \code{dmod}
#'      object itself is a list, which is a \code{glm} object with the name of
#'      the model as \code{$distr} (e.g., \code{"gamma"} or \code{"lognormal"}).
#'      The \code{dmod} list has all the usual stuff from a fitted glm, including
#'      the raw data, covariate information, etc.}
#'  }
#' @export
ddArray <- function(ringData, scCol = NULL, rCol = "r",
      expoCol = "exposure", ncarcCol = "ncarc", modelSet = "all"){
  dmod0 <- canFit(ringData = ringData, modelSet = modelSet,
    rCol = rCol, scCol = scCol, expoCol = expoCol)
  dd <- list()
  for (distr in names(dmod0)) dd[[distr]] <- tfit(dmod0[[distr]])
  output <- dd
  attr(output, "hidden") <- c("scCol", "ncarc")
  attr(output, "expoCol") <- expoCol
  attr(output, "rCol") <- rCol
  attr(output, "scCol") <- scCol
  attr(output, "ncarcCol") <- ncarcCol
  class(output) <- "ddArray"
  return(output)
}

#'  Akaike Information Criterion
#' @description functions for calculating AIC and AICc for carcass dispersion
#'  models.
#' @param x list of models to calculate AICs for
#' @param rCol name of the "r" variable in model formulas
#' @param scCol name of the search class variable in model formulas
#' @param expoCol name of the exposure variable in model formulas
#' @param ... ignored (but required for S3 functions)
#' @return Data frame with AIC, AICc, deltaAICc for all models in \code{x}
#' @export
aic <- function(x, ...) UseMethod("aic", x)

#' @rdname aic
#' @export
aic.ddArray <- function(x, ...){
  nm <- c("model", "k", "AIC", "AICc", "deltaAICc")
  output <- data.frame(array(dim = c(length(x), length(nm))),
    stringsAsFactors = FALSE)
  names(output) <- nm
  output$model <- names(x)
  for (distr in names(x)){
    ci <- c("k", "AIC", "AICc")
    output[output$model == distr, ci] <- aic(x[[distr]])[ci]
  }
  output$deltaAICc <- output$AICc - min(output$AICc, na.rm = TRUE)
  output <- output[order(output$deltaAICc), ]
  attr(output, "n") <- attr(x, "n")
  return(output)
}

#' @rdname aic
#' @export
aic.dd <- function(x, ...){
  aic0 <- AIC(x$dmod)
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
      svar = data$ringData$svar[1],
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
        svar = data$ringData[, "svar"][1],
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
    xmax = NULL, resolution = 250, mod_highlight = NULL, ...){
  if (identical(set, "all")) set <- names(x)
  set <- set[sapply(x, function(tmp) !any(is.na(tmp$parm)))]
  if (type == "rcd")
    set <- set[which(!set %in% c("exponential", "tnormal", "constant"))]
  if (!any(set %in% mod_name)) stop("plot.ddArray: some model(s) undefined")
  dd <- x[set]
  if ("ddArray" %in% class(dd)){
    aic_proper <- aic(dd)
    if (extent == "win") {
      aic_proper <- aic_proper[order(aic_proper$AICw),]
      aic_proper$dAIC <- aic_proper$AICw - min(aic_proper$AICw)
    } else {
      aic_proper$dAIC <- aic_proper$deltaAICc
    }
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
      if (is.null(xmax)) xmax <- max(dd[[1]]$dmod$data[, "r"])
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
    lwd[1] <- 2
    for (i in 1:nrow(aic_proper))
     leglab[i] <- sprintf("%18-s%6.2f", aic_proper$model[i], aic_proper$dAIC[i])
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
      lines(xseq, pdd(xseq, model = dd[mod_highlight], extent = extent),
        col = mod_color[mod_highlight], lwd = 2, lty = lty[mod_highlight])
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
      for (fi in aic_proper$model)
        ymax <- max(ymax, max(rcd(x = xseq, model = dd[fi], extent = extent)))
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
  } else {
    plot(dd, ...)
  }
}

#' @name plot
#' @export
plot.dd <- function(x, type = "CDF", extent = "full",
    xmax = NULL, resolution = 250, nsim = 1000, CL = 0.9, ...){
  if (extent == "full" && is.na(x$beta[1]))
    stop("improper distribution. cannot plot")
  arglist <- list(...)
  if ("xlim" %in% names(arglist)){
    xmax <- max(arglist$xlim)
    xmin <- min(arglist$xlim)
  } else {
    xmin <- 0
  }
  if (is.null(xmax)) xmax <- max(x$dmod$data$r)
  xseq <- seq(xmin, xmax, length.out = resolution)
  CI <- ddCI(x, r = xseq, type = type, CL = CL, nsim = nsim, extent = extent)
  do.call(par, par_default)
  if (type == "CDF"){
    ymax <- 1
    ylab <- "CDF = fraction of carcasses with x meters of turbine"
  } else if (type == "PDF"){
    ymax <- max(CI[, 3])
    ylab <- "PDF = probability density function for carcasses at distance x"
  } else {
    stop("type = ", deparse(substitute(type)), " not supported")
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
  if (extent == "full"){
    bmean <- x$beta
    bvar <- x$varbeta
  } else if (extent == "win"){
    bmean <- x$dmod$coefficients
    bvar <- summary(x$dmod)$cov.unscaled
  } else {
    stop("ddSim: \"extent\" must be \"full\" or \"win\"")
  }
  if (any(is.na(bmean))) stop("ddSim: bad parameter estimates. cannot simulate")
  if (any(is.na(bvar))) stop("ddSim: bad var(beta_hat). cannot simulate")
  if (length(bmean) == 1){
    beta_sim <- rnorm(nsim, mean = bmean, sd = sqrt(bvar))
    beta_sim <- array(beta_sim, dim = c(length(beta_sim), 1))
    colnames(beta_sim) <- names(bmean)
  } else {
    beta_sim <- mvtnorm::rmvnorm(nsim, mean = bmean, sigma = bvar, method = "svd")
  }
  distr <- x$distr
  attr(beta_sim, "distr") <- distr
#  ans <- cof2parms(beta_sim, distr)
  ans <- cbind(beta_sim, cof2parms(beta_sim, distr), proper =
    ifelse(distr %in% c("exponential", "tnormal"), 1, 2) * cofOK(beta_sim, distr))
  attr(ans, "distr") <- distr
  attr(ans, "trad") <- max(x$dmod$data$r, na.rm = TRUE)
  class(ans) <- "ddSim"
  ans
}

#' @rdname ddSim
#' @export
ddSim.dmod <- function(x, nsim = 1000, extent = "full", ...){
  if(!extent %in% c("full", "win"))
    stop("ddSim: \"extent\" must be \"full\" or \"win\"")
  bmean <- x$coefficients
  bvar <- summary(x)$cov.unscaled
  if (any(is.na(bmean))) stop("ddSim: bad parameter estimates. cannot simulate")
  if (any(is.na(bvar))) stop("ddSim: bad var(beta_hat). cannot simulate")
  if (length(bmean) == 1){
    beta_sim <- rnorm(nsim, mean = bmean, sd = sqrt(bvar))
    beta_sim <- array(beta_sim, dim = c(length(beta_sim), 1))
    colnames(beta_sim) <- names(bmean)
  } else {
    beta_sim <- mvtnorm::rmvnorm(nsim, mean = bmean, sigma = bvar, method = "svd")
  }
  distr <- x$distr
  attr(beta_sim, "distr") <- distr
#  ans <- cof2parms(beta_sim, distr)
  ans <- cbind(beta_sim, cof2parms(beta_sim, distr), proper =
    ifelse(distr %in% c("exponential", "tnormal"), 1, 2) * cofOK(beta_sim, distr))
  attr(ans, "distr") <- distr
  attr(ans, "trad") <- max(x$data$r, na.rm = TRUE)
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
    if (!all(colnames(x) %in% j)) noclass <- TRUE
    if (length(j) == 1) y <- matrix(y, ncol = 1, dimnames = list(NULL, j))
    attr(y, "distr") <- attr(x, "distr")
    attr(y, "trad") <- attr(x, "trad")
    class(y) <- "ddSim"
    return(y)
  }
  if (nrow(x) == 1 || (!missing(i) && length(i) == 1)) y <- matrix(y, nrow = 1, dimnames = list(NULL, names(y)))
  if (NCOL(y) == 1 & !noclass){
    nm <- names(y)
    y <- matrix(y, nrow = 1, dimnames = list(NULL, nm))
  }
  attr(y, "distr") <- attr(x, "distr")
  attr(y, "trad") <- attr(x, "trad")
  if(!noclass) class(y) <- "ddSim"
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
  if (length(distr) == 1){
    output <- x[[distr]]
    class(output) <- "dd"
  } else {
    output <- list()
    for (di in distr) output[[di]] <- x[[di]]
    attr(output, "hidden") <- c("scCol", "ncarc", "ncarcCol", "rCol")
    attr(output, "expoCol") <- attr(x, "expoCol")
    attr(output, "rCol") <- attr(x, "rCol")
    attr(output, "scCol") <- attr(x, "scCol")
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

#' @export
print.dd <- function(x, ...){
  bad <- ifelse(any(is.na(x$parms)), TRUE, FALSE)
  if(!bad){
    cat(paste0("Distribution: ", x$distr, "\n"))
    cat(paste0("Formula: ", deparse(x$dmod$formula), "\n\n"))
    cat("Parameters:\n")
    print(x$parms)
    cat("\nCoefficients:\n")
    print(x$beta)
    cat("\nVariance:\n")
    print(x$varbeta)
    cat("\n")
  } else {
    cat(paste0("Distribution: ", x$distr, " (improper)\n"))
    cat(paste0("Formula: ", deparse(x$dmod$formula), "\n\n"))
    cat("Proper", x$distr, "distribution not fittable.\n")
    cat("Cannot extrapolate beyond search radius.")
    cat("\nCoefficients for fitted glm within search radius:\n")
    print(x$dmod$coefficients)
    cat("\n")
  }
}

#' @export
print.dmod <- function(x, ...){
  cat(paste0("Formula: ", deparse(x$formula), "\n\n"))
  cat(paste0("Distribution: ", x$distr))
  cat("\n\n")
  cat("Coefficients:\n")
  print(x$coefficients, quote = FALSE)
}

### change so that there is better coherence, less clutter
#   (group essential features of each model together
#' @export
print.ddArray <- function(x, ...){
  for (distr in names(x)) {
    if (distr == "constant") next
    print(x[[distr]])
    cat("****************************************\n")
  }
}

#'  convert glm parameters into distr parameters
#' @description There is a natural connection between glm and exponential
#'  families of probability distributions. \code{cof2parms} converts the glm
#'  parameters into probability distribution parameters in cases where the
#'  glm converges to 0 as x approaches infininty and converges to x0 < Inf as
#'  x approaches 0.
#' @param x structure containing the parameters values to be converted
#' @param distr name of the distribution
#' @param rvar name of the distance variable
#' @return array with columns for the glm parameters, the corresponding
#'  probability distribution parameters,
#' @export
cof2parms <- function(x, ...) UseMethod("cof2parms", x)

#' @export
cof2parms.numeric <- function(x, distr, rvar = "r", ...){
  # x is an array (or vector) of glm parameters, named to match the required
  #  parameters for the given function.
  #### NOTE: tricky because x can be a vector of parameters for a 1-parameter
  ####  distribution or a single set of parameters for a multi-parameter model
  if (length(dim(x)) == 0){
    # If named entries, assume vector is a single set of parameters.
    # If entries not named, assume vector of parameters for 1-parameter distribution.
    gp <- array(x, dim = c(1, length(x)))
    if (length(names(x)) == 0){
      gp <- array(x, dim = c(length(x), 1))
      colnames(gp) <- NULL # unknown name but can infer from distr (later)
    } else {
      gp <- array(x, dim = c(1, length(x))) # gp = glm parameters
      colnames(gp) <- names(x)
    }
  } else {
    gp <- x
  }
  if (distr == "gamma"){
    if (!any(grepl(paste0("^", rvar, "$"), colnames(gp))))
      stop("rvar = ", rvar, " not found in columns of x = ", deparse(substitute(x)))
    if (!any(grepl(paste0("^log\\(", rvar, "\\)$"), colnames(gp))))
      stop("log(", rvar, ") not found in columns of x = ", deparse(substitute(x)))
    gp <- gp[, c(paste0("log(", rvar, ")"), rvar), drop = FALSE]
    colnames(gp) <- c("log(r)", "r")
    parms <- array(dim = c(nrow(gp), 2 * ncol(gp) + 1))
    colnames(parms) <- c(colnames(gp), "shape", "rate", "proper")
    parms[, colnames(gp)] <- gp
    parms[, c("shape", "rate")] <-
      cbind(gp[, "log(r)"] + 2, -gp[, "r"])
    parms[, "proper"] <- 2 * (parms[, "shape"] > 0 & parms[, "rate"] > 0)
  } else if (distr == "lognormal"){
    if (!any(grepl(paste0("^log\\(", rvar, "\\)$"), colnames(gp))))
      stop("log(", rvar, ") not found among columns of x = ", deparse(substitute(x)))
    if (!any(grepl(paste0("^I\\(log\\(", rvar, "\\)\\^2\\)$"), colnames(gp))))
      stop("log(", rvar, "^2) not found in columns of x = ", deparse(substitute(x)))
    gp <- gp[, c(paste0("log(", rvar, ")"), paste0("I(log(", rvar, ")^2)")),
      drop = FALSE]
    colnames(gp) <- c(paste0("log(", rvar, ")"), paste0("I(log(", rvar, ")^2)"))
    parms <- array(dim = c(nrow(gp), 2 * ncol(gp) + 1))
    colnames(parms) <- c(colnames(gp), "meanlog", "sdlog", "proper")
    parms[, colnames(x)] <- gp
    parms[, c("meanlog", "sdlog")] <- cbind(
      -0.5 * (2 + gp[, "log(r)"])/gp[, "I(log(r)^2)"],
      sqrt(-0.5/gp[, "I(log(r)^2)"]))
    parms[, "proper"] <- 2 * (parms[, "sdlog"] > 0)
  } else if (distr == "logLinear"){
    if (ncol(gp) == 1){
      colnames(gp) <- "r"
    } else {
      if (!any(grepl(paste0("^", rvar, "$"), colnames(gp))))
        stop("rvar = ", rvar, " not found in colums of x = ", deparse(substitute(x)))
      gp <- array(gp[, rvar], dim = c(nrow(gp), 1))
      colnames(gp) <- "r"
    }
    parms <- array(dim = c(nrow(gp), 2 * ncol(gp) + 1))
    colnames(parms) <- c(colnames(gp), "b1", "proper")
    parms[, colnames(gp)] <- gp
    parms[, "b1"] <- parms[, "r"]
    parms[, "proper"] <- 2 * (parms[, "r"] < 0)
  } else if (distr == "logQuadratic"){
    if (!any(grepl(paste0("^", rvar, "$"), colnames(gp))))
      stop("rvar = ", rvar,  " not found among columns of x = ", deparse(substitute(x)))
    if (!any(grepl(paste0("^I\\(", rvar, "\\^2\\)$"), colnames(gp))))
      stop("I(", rvar, "^2) not found in columns of x = ", deparse(substitute(x)))
    gp <- gp[, c(rvar, paste0("I(", rvar, "^2)")), drop = FALSE]
    parms <- array(dim = c(nrow(gp), 2 * ncol(gp) + 1))
    colnames(parms) <- c(colnames(gp), "b1", "b2", "proper")
    parms[, colnames(gp)] <- gp
    parms[, c("b1", "b2")] <- parms[, colnames(gp)]
    parms[, "proper"] <- 2 * (parms[, colnames(gp)[2]] < 0)
  } else if (distr == "logCubic"){ # this is muddled!
    if (!any(grepl(paste0("^", rvar, "$"), colnames(gp))))
      stop("rvar = ", rvar,  " not found among columns of x = ", deparse(substitute(x)))
    if (!any(grepl(paste0("^I\\(", rvar, "\\^2\\)$"), colnames(gp))))
      stop("I(", rvar, "^2) not found in columns of x = ", deparse(substitute(x)))
    if (!any(grepl(paste0("^I\\(", rvar, "\\^3\\)$"), colnames(gp))))
      stop("I(", rvar, "^3) not found in columns of x = ", deparse(substitute(x)))
    gp <- gp[, c(rvar, paste0("I(", rvar, "^2)"), paste0("I(", rvar, "^3)")),
      drop = FALSE]
    parms <- array(dim = c(nrow(gp), 2 * ncol(gp) + 1))
    colnames(parms) <- c(colnames(gp), paste0("b", 1:3), "proper")
    parms[, c("r", "I(r^2)", "I(r^3)")] <- gp
    parms[, paste0("b", 1:3)] <- parms[, names(gp)]
    parms <- array(dim = c(nrow(gp), ncol(gp) + 1))
    colnames(parms) <- c(colnames(gp), "proper")
    parms[, names(x$beta)] <- gp
    parms[, "proper"] <- 2 * (parms[, colnames(gp)[3]] < 0)
  } else if (distr == "inverse_gamma"){
   if (!any(grepl(paste0("^I\\(1\\/", rvar, "\\)$"), colnames(gp))))
      stop(paste0("I(1/", rvar, ") not found in columns of x = "), deparse(substitute(x)))
    if (!any(grepl(paste0("^log\\(", rvar, "\\)$"), colnames(gp))))
      stop("log(", rvar, ") not found in columns of x = ", deparse(substitute(x)))
    gp <- gp[, c(paste0("I(1/", rvar, ")"), paste0("log(", rvar, ")")),
      drop = FALSE]
    colnames(gp) <- c("I(1/r)", "log(r)")  # not redundant if rvar != "r"
    parms <- array(dim = c(nrow(gp), 2 * ncol(gp) + 1))
    colnames(parms) <- c(colnames(gp), "shape", "scale", "proper")
    parms[, colnames(gp)] <- gp
    parms[, c("shape", "scale")] <- cbind(-gp[, "log(r)"] - 2, -gp[, "I(1/r)"])
    parms[, "proper"] <- 2 * (parms[, "shape"] > 0 & parms[, "scale"] > 0)
  } else if (distr == "paranormal_gamma"){
    if (!any(grepl(paste0("^", rvar, "$"), colnames(gp))))
      stop("rvar = ", rvar, " not found in columns of x = ", deparse(substitute(x)))
    if (!any(grepl(paste0("^log\\(", rvar, "\\)$"), colnames(gp))))
      stop("log(", rvar, ") not found in columns of x = ", deparse(substitute(x)))
    if (!any(grepl(paste0("^I\\(", rvar, "\\^2\\)$"), colnames(gp))))
      stop("I(", rvar, "^2) not found in columns of x = ", deparse(substitute(x)))
    if (!any(grepl(paste0("^I\\(", rvar, "\\^3\\)$"), colnames(gp))))
      stop("I(", rvar, "^3) not found in columns of x = ", deparse(substitute(x)))
    gp <- gp[, c(paste0("log(", rvar, ")"), rvar, paste0("I(", rvar, "^2)"),
      paste0("I(", rvar, "^3)")), drop = FALSE]
    parms <- array(dim = c(nrow(gp), 2 * ncol(gp) + 1))
    colnames(parms) <- c("log(r)", "r", "I(r^2)", "I(r^3)", paste0("b", 0:3), "proper")
    parms[, c("log(r)", "r", "I(r^2)", "I(r^3)")] <- gp
    parms[, paste0("b", 0:3)] <- gp
    parms[, "proper"] <- 2 * (parms[, colnames(gp)[length(colnames(gp))]] < 0)
  } else if (distr == "Rayleigh"){
    if (!any(grepl(paste0("^I\\(", rvar, "\\^2\\)$"), colnames(gp))))
      stop("I(", rvar, "^2) not found in columns of x = ", deparse(substitute(x)))
    gp <- gp[, paste0("I(", rvar, "^2)"), drop = FALSE]
    colnames(gp) <- "I(r^2)"
    parms <- array(dim = c(nrow(gp), 2 * ncol(gp) + 1))
    colnames(parms) <- c(colnames(gp), "s2", "proper")
    parms[, colnames(gp)] <- gp
    parms[, "s2"] <- -0.5 /gp[, "I(r^2)"]
    parms[, "proper"] <- 2 * (parms[, "s2"] > 0)
  } else if (distr == "MaxwellBoltzmann"){
    if (!any(grepl(paste0("^I\\(", rvar, "\\^2\\)$"), colnames(gp))))
      stop("I(", rvar, "^2) not found in columns of x = ", deparse(substitute(x)))
    gp <- gp[, paste0("I(", rvar, "^2)"), drop = FALSE]
    colnames(gp) <- "I(r^2)"
    parms <- array(dim = c(nrow(gp), 2 * ncol(gp) + 1))
    colnames(parms) <- c(colnames(gp), "a", "proper")
    parms[, colnames(gp)] <- gp
    parms[, "a"] <- sqrt(-0.5 /gp[, "I(r^2)"])
    parms[, "proper"] <- 2 * (parms[, "a"] > 0)
   } else if (distr == "custom"){
    # what about interactions?!
    # maybe a different pathway...much more difficult management;
    # could just do predict.glm? anyway, a different pathway
    if(!all(grepl(paste0("^", rvar, "$"), colnames(gp)) |
        grepl(paste0("[^a-zA-Z0-9_\\.]", rvar, "[^a-zA-Z0-9_\\.]"), colnames(gp))))
      stop("some columns in x = ", deparse(substitute(x)), " missing rvar = ", rvar)
    parms <- array(dim = c(nrow(gp), ncol(gp) + 1))
    nm <- gsub(paste0("^", rvar, "$"), "r", colnames(gp))
    nm <- gsub(paste0("\\(", rvar, "\\)"), "(r)", nm)
    nm <- gsub(paste0("\\(", rvar, "\\^"), "(r^", nm)
    nm <- gsub(paste0("\\/", rvar, "\\)"), "/r)", nm)
    nm <- gsub(paste0("\\/", rvar, "\\^"), "/r^", nm)
    colnames(parms) <- c(nm, "proper")
    parms[, nm] <- gp
    for (i in 1:nrow(gp)){
      bmat <- parms[i, nm]
      tmp <- try(integrate(f = function(r){
          rmat <- eval(parse(text = c("cbind(", paste(nm, collapse = ", "), ")")))
          exp(rmat %*% bmat)
        }, lower = 0, upper = Inf)$val, silent = TRUE)
      parms[i, "proper"] <- ifelse("try-error" %in% class(tmp), 0, 1)
    }
  } else if (distr == "exponential"){
    gp <- gp[, rvar, drop = FALSE]
    colnames(gp) <- "r"
    parms <- array(dim = c(nrow(gp), 2 * ncol(gp) + 1))
    colnames(parms) <- c(colnames(gp), "rate", "proper")
    parms[, colnames(gp)] <- gp
    parms[, "rate"] <- -parms[, "r"]
    parms[, "proper"] <- 1 * (parms[, "rate"] > 0)
  } else if (distr == "tnormal"){
    gp <- gp[, c(rvar, paste0("I(", rvar, "^2)")), drop = FALSE]
    colnames(gp) <- c("r", "I(r^2)")
    parms <- array(dim = c(nrow(gp), 2 * ncol(gp) + 1))
    colnames(parms) <- c(colnames(gp), "mean", "sd", "proper")
    parms[, colnames(gp)] <- gp
    parms[, "mean"] <- -0.5 * parms[, "r"]/parms[, "I(r^2)"]
    parms[, "sd"] <- sqrt(-0.5/parms[, "I(r^2)"])
    parms[, "proper"] <- 1 * (parms[, "sd"] > 0)
  } else if (distr == "inverse_gaussian"){
    gp <- gp[, c(paste0("I(1/", rvar, ")"), rvar), drop = FALSE]
    colnames(gp) <- c("I(1/r)", "r")
    parms <- array(dim = c(nrow(gp), 2 * ncol(gp) + 1))
    colnames(parms) <- c(colnames(gp), "mean", "dispersion", "proper")
    parms[, colnames(gp)] <- gp
    parms[, "mean"] <- sqrt(parms[, "I(1/r)"]/parms[, "r"])
    parms[, "dispersion"] <- -0.5/parms[, "r"]
    parms[, "proper"] <- 2 * (parms[, "mean"] > 0 & parms[, "dispersion"] > 0)
  } else if (distr == "chisq"){
    gp <- gp[, paste0("log(", rvar, ")"), drop = FALSE]
    colnames(gp) <- "log(r)"
    parms <- array(dim = c(nrow(gp), 2 * ncol(gp) + 1))
    colnames(parms) <- c(colnames(gp), "df", "proper")
    parms[, colnames(gp)] <- gp
    parms[, "df"] <- 2 * parms[, "log(r)"] + 4
    parms[, "proper"] <- 2 * (parms[, "df"] > 0)
  } else if (distr == "constant"){
    parms <- c(b1 = NA)
  }
  attr(parms, "distr") <- distr
  parms
} # cof2parms.numeric

#' @export
cof2parms.glm <- function(x, distr = NULL, rvar = "r", scCol = NULL, ...){
  # do not parse the distr; assume "custom" if not given
  # will change in future...
  if (is.null(distr)) stop("under contstruction...supply a distr")
  dmod <- x
  if (distr == "constant"){
    ans <- list(
      beta = c(r = 0),
      varbeta = c(r = 0),
      parms = c(b1 = NA),
      dmod = dmod,
      proper = 0
    )
    class(ans) <- "dd"
    return(ans)
  }
  if (!is.null(scCol)){
    ind <- !grepl("\\(Intercept\\)", names(dmod$coef)) & !grepl(scCol, names(dmod$coef))
  } else {
    ind <- which(!grepl("\\(Intercept\\)", names(dmod$coef)))
  }
  beta <- dmod$coef[ind]
  varbeta <- summary(dmod)$cov.scaled[ind, ind]
  # define default parms and proper...update later if distr is good
  parms <- rep(NA, length(beta))
  names(parms) <- parm_name[[distr]]
  proper <- 0
  if (distr == "gamma"){
    if (beta["log(r)"] <= -2){
      names(parms) <- c("shape", "rate")
    } else {
      proper <- 2
      parms <- c(shape = unname(beta["log(r)"]) + 2, rate = unname(-beta["r"]))
    }
  } else if (distr == "lognormal"){
    if (beta[critical_parameter[distr]] < 0){ # then calculate distr parms
      parms <- c(
        meanlog = unname(-0.5 * (2 + beta["log(r)"])/beta["I(log(r)^2)"]),
        sdlog = unname(sqrt(-0.5/beta["I(log(r)^2)"])))
      proper <- 2 # all the lognormals are proper if log(r)^2 term is negative
    }
  } else if (distr == "logLinear"){
    if (beta[critical_parameter[distr]] < 0){ # then calculate distr parms
      proper <- 2 # always proper = 2 if b1 < 0
      parms <- beta; names(parms) <- "b1"
    }
  } else if (distr == "logQuadratic"){
    if (beta[critical_parameter[distr]] < 0){
      proper <- 2
      parms <- beta; names(parms) = c("b1", "b2")
    }
  } else if (distr == "logCubic"){
    if (beta[critical_parameter[distr]] < 0){
      proper <- 2
      parms <- beta
      names(parms) = c("b1", "b2", "b3")
    }
  } else if (distr == "inverse_gamma"){
    if (-beta[critical_parameter[distr]] >= 2){
      proper <- 2
      parms <- c(-beta["log(r)"] - 2, -beta["I(1/r)"])
      names(parms) <- c("shape", "scale")
    }
  } else if (distr == "paranormal_gamma"){
    if (beta[critical_parameter[distr]] < 0){
      proper <- 2
      parms <- beta
      names(parms) <- c("b0", "b1", "b2", "b3")
    }
  } else if (distr == "Rayleigh"){
    if (beta[critical_parameter[distr]] < 0){
      proper <- 2
      parms <- -0.5 /beta["I(r^2)"]
      names(parms) <- "s2"
    }
  } else if (distr == "MaxwellBoltzmann"){
    if (beta[critical_parameter[distr]] < 0){
      proper <- 2
      parms <- sqrt(-0.5 /beta["I(r^2)"])
      names(parms) <- "a"
    }
  } else if (distr == "Pareto"){
    if (beta[critical_parameter[distr]] < -2){
      proper <- 2
      parms <- -beta["log(r)"] - 2
      names(parms) <- "a"
    }
  } else if (distr == "chisq"){
    if (beta[critical_parameter[distr]] > -2){
      proper <- 1
      parms <- 2 * beta + 4
      names(parms) <- "df"
      if (parms > 2) proper <- 2
    }
  } else if (distr == "inverse_gaussian"){
    if (beta["I(1/r)"] < 0 & beta["r"] < 0){
      proper <- 2
      parms <- c(sqrt(beta["I(1/r)"]/beta["r"]), -0.5/beta["I(1/r)"])
      names(parms) <- c("mean", "dispersion")
    }
  } else if (distr == "exponential"){
    if (beta[critical_parameter[distr]] < 0){
      proper <- 1 # never
      parms <- -beta["r"]
      names(parms) <- "rate"
    }
  } else if (distr == "tnormal"){
    if (beta[critical_parameter[distr]] < 0){
      proper <- 1 # never proper = 2
      parms <- c(-0.5 * beta["r"]/beta["I(r^2)"], sqrt(-0.5/beta["I(r^2)"]))
      names(parms) <- c("mean", "sd")
    }
  } else if (distr == "constant"){
    proper <- 0
    parms <- c(b1 = NA)
  }
  ans <- list(beta = beta, varbeta = varbeta, distr = distr,
    parms = parms, dmod = dmod, proper = proper)
  class(ans) <- "dd"
  return(ans)
}

#'  Calculate CI for CDF, PDF, or rcd for a \code{dd} model at distance \code{r}
#' @param mod a \code{dd} object
#' @param r distance from turbine (scalar or vector)
#' @param type "CDF", "PDF", or "rcd"
#' @param CL confidence level for the confidence interval(s)
#' @param nsim number of simulation draws to base the estimate of CI on
#' @param na.tol maximum fraction of invalid parameter sets to discard when
#'  constructing CIs; abort if \code{mean(mod[, "proper"]) > na.tol}
#' @return array (\code{ddCI} class) with columns for distance and the CI bounds
#' @export
ddCI <- function(mod, r, type = "CDF", CL = 0.9, nsim = 1000,
    extent = "full", na.tol = 0.1){
  if (!any(c("dd", "dmod") %in% class(mod)))
    stop("mod (", deparse(substitute(mod)), ") must be a dd object")
  p0 <- ddSim(mod, nsim = nsim, extent = extent)
  if (extent == "full" && mean(p0[, "proper"] == 0) > na.tol)
    stop("too many improper parameter sets to construct meaningful CI ",
      "[", round(100 * mean(p0[, "proper"] == 0), 1), "% > na.tol = ",
      na.tol, "]")
  parms <- p0[p0[, "proper"] > 0 , ] # this gives result with nrow < nsim
#  rx <- rep(r, each = nrow(parms))
  if (type == "CDF"){
    tmp <- array(dim = c(length(r), nrow(parms)))
    for (simi in 1:nrow(parms))
      tmp[, simi] <- pdd(r, parms[simi, ], extent = extent)
    CI <- cbind(r, matrixStats::rowQuantiles(tmp,
      probs = c((1 - CL)/2, (1 + CL)/2)
    ))
  } else if (type == "PDF"){
    tmp <- array(dim = c(length(r), nrow(parms)))
    for (simi in 1:nrow(parms))
      tmp[, simi] <- ddd(r, parms[simi, ], extent = extent)
    CI <- cbind(r, matrixStats::rowQuantiles(tmp,
      probs = c((1 - CL)/2, (1 + CL)/2)
    ))
  }
  class(CI) <- "ddCI"
  CI
}

#' Convert glm distance parameters into distribution parameters
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
cof2parms.dmod <- function(x, ...){
  return(cof2parms(x$coef, x$distr))
}

#'  Convert a glm model into a dd object (under construction...)
#' @param dmod a \code{glm} object
#' @param rvar name of distance variable
#' @param svar name of search class variable
#' @param expovar name of exposure variable
#' @return \code{dd} object
##' @export
cof2parms.glm <- function(x, rvar = "r", svar = NULL, expovar = "exposure",
  ...){
  # quick check on the format of the model: all terms are r, sc, or expos?
  form <- deparse(dmod$formula)
  trms <- names(attributes(terms(dmod))$dataClasses)[-1]
  terms_r <- (grepl(paste0("^", rvar, "$"), trms) |
    grepl(paste0("[^a-zA-Z0-9_\\.]", rvar, "[^a-zA-Z0-9_\\.]"), trms)) &
    !grepl("offset", trms)
  terms_sc <- grepl(paste0("^", svar, "$"), trms)
  terms_expo <- grepl(expovar, trms) & grepl("offset", trms)
  if (!all(terms_r | terms_sc | terms_expo)){
    bad <- which(!(terms_r | terms_sc | terms_expo))
    stop(form, ": \nmodel terms for ", paste(trms[bad], collapse = ", "),
      " not accounted for in arglist ",
      "rvar = ", rvar, ", svar = ", ifelse(is.null(svar), "NULL", svar),
      ", expovar = ", ifelse(is.null(svar), "NULL", expovar))
  }
  # need to parse a single dmod and convert to a distribution
  # decide what the form of the model is (if any specified template works)
  # look at deparse(dmod$formula)
  # check whether the model fits the template of the known distributions:

  # If any term in the formula is missing rvar, svar, and expovar, the model
  # is not properly formatted for analysis. In addition, the rvar terms must have
  # rvar in isolation, without terms that are missing rvar, svar, and expovar?
  ############
  ############
  logr <- paste0("log\\(", rvar, "\\)")
  Ilogr <- paste0("I\\(log\\(", rvar, "\\)\\)")
  logr2 <- paste0("I\\(log\\(", rvar, "\\)\\^2\\)")
  r2 <- paste0("I\\(", rvar, "\\^2\\)")
  r3 <- paste0("I\\(", rvar, "\\^3\\)")
  rinv <- paste0("I\\(1\\/", rvar, "\\)")
  offs <- paste0("offset\\(log\\(", expovar, "\\)\\)")
  if (is.null(svar)){
    sc <- ""
  } else {
    sc <- paste0(svar, " \\+ ")
  }
  distr <- list(
    gamma = paste0("\\~ ", logr, " \\+ ", rvar, " \\+ ", sc, offs, "$"),
    lognormal = paste0("\\~ ", logr, " \\+ ", logr2, " \\+ ", sc, offs, "$"),
    logLinear = paste0("\\~ r", " \\+ ", sc, offs, "$"),
    logQuadratic = paste0("\\~ r", " \\+ ", r2, " \\+ ", sc, offs, "$"),
    logCubic = paste0("\\~ r", " \\+ ", r2, " \\+ ", r3, " \\+ ", sc, offs, "$"),
    inverse_gamma = paste0("\\~ ", rinv, " \\+ ", logr, " \\+ ", sc, offs, "$"),
    parnormal_gamma = paste0("\\~ ", logr, " \\+ ", "r", " \\+ ", r2, " \\+ ",
      r3, " \\+ ", sc, offs, "$"),
    Rayleigh = paste0("\\~ ", r2, " \\+ ", sc, offs, "$"),
    MaxwellBoltzmann = paste0("\\~ ", r2, " \\+ ", sc,
      "offset\\(I\\(log\\(exposure \\* \\(r \\- 0\\.5\\)\\)\\)\\)", "$"),
    constant = paste0("\\~ ", sc, offs, "$")
  )
  nm <- names(distr)[sapply(distr, FUN = function(x) grepl(x, form))]
  nm <- ifelse(length(nm) == 0, "custom", nm)
  beta <- dmod$coef[which(!names(dmod$coef) %in% "(Intercept)" &
    !grepl(svar, names(dmod$coef)))]
  varbeta <- summary(dmod)$cov.scaled[names(beta), names(beta)]
  if (nm == "gamma" & beta["r"] < 0 & beta["log(r)"])
  ans <- list(beta = beta, varbeta = varbeta, distr = nm)
  class(ans) <- "distr"
  return(ans)
}

#' Function to mimic dnorm, dexp, etc. for \code{dd} models
#' @param x numeric x >= 0
#' @param model either a \code{dd} object or the name of a distance model
#'  (e.g., "gamma").
#' @param parms scalar, vector, or array of parameters. Or if \code{x} is a
#'  \code{dd} object, then \code{parms} may be NULL, in which case, \code{ddd}
#'   uses the MLE parameters found in the \code{dd} object.
#' @param trad the radius to use for truncated distributions
#'  (\code{trad = NULL} to use full, non-truncated distribution)
#' @return vector of PDF values. If \code{trad = NULL} and a full distribution
#'  is used, the length of the output is equal to the maximum of the length of
#'  \code{x} and the number of parameter sets are included in \code{model} (i.e.
#'  1 if \code{model} is a \code{dd} object, or \code{nrow(model)} if
#'  \code{model} is a \code{ddSim} object. If a truncated distribution is to be
#'  used, then \code{length(x) = 1} or the number of parameter sets must be 1
#'  (which occurs if \code{model} is a \code{dd} object or \code{nrow(model) = 1}
#'  if \code{model} is a \code{ddSim} object.
#' @export
ddd <- function(x, model, extent = "full"){
#model<-ddw[[distr]]
  if ("dd" %in% class(model)){
    if (extent == "full" && any(is.na(model$parms)))
      stop("improper distribution")
    if (extent == "full"){ # parameters for distribution
      parms <- matrix(model$parms, nrow = 1,
        dimnames = list(NULL, names(model$parms)))
    } else { # coefficients for glm
      cof <- matrix(model$dmod$coef, nrow = 1,
        dimnames = list(NULL, names(model$dmod$coef)))
      trad <- max(model$dmod$data$r, na.rm = TRUE)
    }
    distr <- model$distr
  } else if ("ddSim" %in% class(model)){
    if (extent == "full"){
      parms <- model[, parm_name[[attr(model, "distr")]], drop = FALSE]
    } else {
      cof <- model[, cof_name[[attr(model, "distr")]], drop = FALSE]
      trad <- attr(model, which = "trad")
    }
    distr <- attr(model, which = "distr")
  }
  if (extent == "full"){
    output <- switch(distr,
      gamma = dgamma(x, shape = c(parms[, "shape"]), rate = c(parms[, "rate"])),
      lognormal = dlnorm(x,
        meanlog = c(parms[, "meanlog"]),
        sdlog = c(parms[, "sdlog"])),
      logLinear = dlogL(x, b1 = c(parms[, "b1"])),
      logQuadratic = dlogQ(x, b1 = c(parms[, "b1"]), b2 = c(parms[, "b2"])),
      logCubic = {
        b1 <- c(parms[, "b1"]); b2 <- c(parms[, "b2"]); b3 <- c(parms[, "b3"])
        const <- try(1/integrate(f = function(x)
          x * exp(b1 * x + b2 * x^2 + b3 * x^3), lower = 0, upper = Inf,
          rel.tol = .Machine$double.eps^0.5)$val, silent = TRUE)
        if ("try-error" %in% class(const))
          stop("invalid parameter set for ", model$distr, ". Integral for ",
               "normalizing constant does not converge.")
        dlogC(x, b1 = b1, b2 = b2, b3 = b3)
      },
      inverse_gamma = digam(x, shape = parms[, "shape"], scale = parms[, "scale"]),
      paranormal_gamma = {
        b0 <- c(parms[, "b0"])
        b1 <- c(parms[, "b1"])
        b2 <- c(parms[, "b2"])
        b3 <- c(parms[, "b3"])
        const <- try(1/integrate(
          f = function(x) x * exp(b0 * log(x) + b1 * x + b2 * x^2 + b3 * x^3),
          lower = 0, upper = Inf, rel.tol = .Machine$double.eps^0.5)$val,
          silent = TRUE)
        if ("try-error" %in% class(const))
          stop("invalid parameter set for ", distr, ". Integral for ",
               "normalizing constant does not converge.")
        const * dpng(x, b0 = b0, b1 = b1, b2 = b2, b3 = b3, const = const)
      },
      Rayleigh = dRay(x, s2 = parms[, "s2"]),
      MaxwellBoltzmann = dmb(x, a = parms[, "a"]),
      constant = NA,
      Pareto = dPare1(x, a = parms[, "a"]),
      tnormal = {mu <- parms[, "mean"]; sig <- parms[, "sd"];
        Fa <- pnorm(0, mean = mu, sd = sig)
        ans <- numeric(length(x))
        ans[x >=0]  <- dnorm(x[x >= 0], mean = mu, sd = sig)/(1 - Fa)
        ans
      },
      exponential = dexp(x, rate = parms[, "rate"]),
      inverse_gaussian = statmod::dinvgauss(x,
        mean = parms[, "mean"], dispersion = parms[, "dispersion"]),
      chisq = dchisq(x, df = parms[, "df"])
    )
  } else if (extent == "win"){
    if (length(x) > 1 && nrow(cof) > 1 && length(x) != nrow(cof))
      stop("either length(x) or number of parameter sets must be 1 or equal to ",
           "each other")
    y <- numeric(max(nrow(cof), length(x)))
    xx <- rep(x, length.out = length(y))
    output <- switch(distr,
      gamma = {
        b1 <- c(cof[, "log(r)"])
        b2 <- c(cof[, "r"])
        y <- b2^2 * x * exp(b1*log(x) + b2*x)/(trad^b1 * (-b2*trad)^(-b1) *
          (gamma(b1 + 2) - expint::gammainc(b1 + 2, -b2*trad)))
        y[xx > trad] <- 0
        y[xx <= 0] <- 0
        y
      },
      lognormal = {
        cof <- cof[, c("(Intercept)", "log(r)", "I(log(r)^2)"), drop = FALSE]
        deno <- numeric(nrow(cof))
        for (i in length(deno)){
          deno[i] <- integrate(
            f = function(r) exp(rmat(r, distr) %*% cof[i, ] + off(r, distr)),
            lower = 0, upper = trad)$val
        }
        y <- exp(rmat(xx, distr) %*% t(cof) + off(xx, distr))/deno
        y[xx > trad] <- 0
        y[xx <= 0] <- 0
        y
      },
      logLinear = {
        b1 <- c(cof[, "r"])
        y <- b1^2 * x * exp(b1*x)/(exp(b1*trad) * (b1*trad - 1) + 1)
        y[xx >= trad] <- 0
        y[xx < 0] <- 0
        y
      },
      logQuadratic =  {
        b1 <- c(cof[, "r"])
        b2 <- c(cof[, "I(r^2)"])
        y <- x * exp(b1*x + b2*x^2)/lQint(trad, b1, b2)
        y[xx > trad] <- 0
        y[xx < 0] <- 0
        y
      },
      logCubic = {
        b1 <- c(cof[, "r"])
        b2 <- c(cof[, "I(r^2)"])
        b3 <- c(cof[, "I(r^3)"])
        if (length(x) == 1 | length(x) == length(b1)){
          x <- rep(x, length.out = length(b1))
          for (bi in 1:length(b1)){
            y[bi] <- x[bi] * exp(b1[bi]*x[bi] + b2[bi]*x[bi]^2 + b3[bi]*x[bi]^3)/
              integrate(f = function(r)
                r * exp(b1[bi]*r + b2[bi]*r^2 + b3[bi]*r^3),
                lower = 0, upper = trad)$val
          }
        } else if (length(b1) == 1){
          y <- x * exp(b1*x + b2*x^2 + b3*x^3)/
            integrate(f = function(r)
              r * exp(b1*r + b2*r^2 + b3*r^3), lower = 0, upper = trad)$val
        }
        y[x >= trad] <- 0
        y
      },
      inverse_gamma = {
        b1 <- c(cof[, "I(1/r)"])
        b2 <- c(cof[, "log(r)"])
        numo <- x * exp(b1/x + b2*log(x))
        deno <- numeric(length(b1))
        if (length(x) == 1 | length(x) == length(b1)){
          for (bi in 1:length(b1)){
            deno[bi] <- integrate(f = function(r)
              r * exp(b1[bi]/r + b2[bi]*log(r)), lower = 0, upper = trad)$val
          }
          y <- numo/deno
        } else if (length(b1) == 1){
          supi <- which(x < trad & x > 0)
          y[supi] <- x[supi] * exp(b1/x[supi] + b2*log(x[supi]))/
            integrate(f = function(r) r * exp(b1/r + b2*log(r)),
              lower = 0, upper = trad)$val
        }
        y[x >= trad] <- 0
        y[x <= 0] <- 0
        y
      },
      paranormal_gamma =  {
        b0 <- c(cof[, "log(r)"])
        b1 <- c(cof[, "r"])
        b2 <- c(cof[, "I(r^2)"])
        b3 <- c(cof[, "I(r^3)"])
        numo <- x * exp(b0*log(x) + b1*x + b2*x^2 + b3*x^3)
        deno <- numeric(max(length(x), length(b1)))
        if (length(x) == 1 | length(x) == length(b1)){
          x <- rep(x, length.out = length(b1))
          for (bi in 1:length(b1)){
            deno[bi] <- integrate(f = function(r)
              r * exp(b0[bi] * log(r) + b1[bi]*r + b2[bi]*r^2 + b3[bi]*r^3),
              lower = 0, upper = trad)$val
          }
        } else if (length(b1) == 1){
          supi <- which(xx < trad & xx >= 0)
          deno <- integrate(f = function(r)
            r * exp(b0*log(r) + b1*r + b2*r^2 + b3*r^3),
            lower = 0, upper = trad)$val
          y[supi] <- numo[supi]/deno
        }
        y[xx >= trad] <- 0
        y[xx <= 0] <- 0
        y
      },
      Rayleigh = {
        b1 <- c(cof[, "I(r^2)"])
        y <- 2 * b1 * x * exp(b1*x^2)/(exp(b1*trad^2) - 1)
        y[x >= trad] <- 0
        y[x <= 0] <- 0
        y
      },
      Pareto = {
        b1 <- c(cof[, "log(r)"])
        if (length(x) == 1) x <- rep(x, length.out = length(b1))
        supi <- which(x >= 1)
        y[supi] <- (b1 + 2) * x[supi]*exp(b1 * log(x[supi]))/(trad^(b1 + 2) - 1)
        y[x >= trad] <- 0
        y[x < 1] <- 0
        y
      },
      MaxwellBoltzmann =  {
        b1 <- c(cof[, "I(r^2)"])
        y <- -2 * b1 * x^2*exp(b1*x^2)/(trad * exp(b1 * trad^2) -
            0.5*sqrt(pi/b1 + 0i)*pracma::erfi(trad * sqrt(b1 + 0i)))
        y[x >= trad] <- 0
        y[x <= 0] <- 0
        Re(y)
      },
      constant = {
        y <- 2 * x/trad^2
        y[x < 0] <- 0
        y[x > trad] <- 0
        y
      },
      tnormal = {
        b1 <- c(cof[, "r"])
        b2 <- c(cof[, "I(r^2)"])
        deno <- -Re(0.5*sqrt(pi/b2  + 0i) * exp(-b1^2/(4*b2)) *
          (pracma::erfi((trad * b2 + 0.5 * b1)/sqrt(b2 + 0i)) -
          pracma::erfi(0.5*b1/sqrt(b2 + 0i))))
        y <- exp(b1*x + b2*x^2)/deno
        y[x >= trad] <- 0
        y[x < 0] <- 0
        Re(y)
      },
      exponential = {
        b1 <- c(cof[, "r"])
        y <- b1 * exp(x * b1)/(exp(trad * b1) - 1)
        y[x >= trad] <- 0
        y[x <= 0] <- 0
        y
      },
      inverse_gaussian =  {
        b1 <- c(cof[, "I(1/r)"])
        b2 <- c(cof[, "r"])
        numo <- numeric(max(length(b1), length(x)))
        deno <- numeric(max(length(b1), length(x)))
        if (length(x) == 1 | length(x) == length(b1)){
          x <- rep(x, length.out = length(b1))
          xp <- which(x > 0)
          numo[xp] <- exp(b1[xp]/x[xp] + b2[xp]*x[xp] - 3/2 * log(x[xp]))
          for (bi in 1:length(b1)){
            deno[bi] <- integrate(f = function(r)
              exp(b1[bi]/r + b2[bi]*r - 3/2 * log(r)),
              lower = 0, upper = trad)$val
          }
          y <- numo/deno
        } else if (length(b1) == 1){
          supi <- which(x <= trad & x > 0)
          y[supi] <- exp(b1/x[supi] + b2*x[supi] - 3/2 * log(x[supi]))/
            integrate(f = function(r)
              exp(b1/r + b2*r - 3/2 * log(r)), lower = 0, upper = trad)$val
          y[x > trad] <- 0
          y[x <= 0] <- 0
        }
        y
      },
      chisq =  {
        b1 = c(cof[, "log(r)"])
        if (length(b1) == 1){
          b1 <- rep(b1, length.out = length(x))
        }
        supi <- which(xx <= trad & xx > 0)
        y[supi] <- xx[supi] * exp(b1[supi]*log(x[supi]) - x[supi]/2)/
          ((gamma(b1[supi] + 2) - expint::gammainc(b1[supi] + 2, trad/2)) * 2^(b1[supi] + 2))
        y[x >= trad] <- 0
        y[x <= 0] <- 0
        y
      }
    )
  } else {
    stop("extent must be 'win' or 'full'")
  }
  return(matrix(output, nrow = ifelse(extent == "full", nrow(parms), nrow(cof))))
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

pdd <- function(x, model, extent = "full"){ # model is either ddSim or dd
  if("dd" %in% class(model)){
    distr <- model$distr
  } else if ("ddSim" %in% class(model)){
    distr <- attr(model, which = "distr")
  } else {
    stop("pdd: model must be dd or ddSim object")
  }
  if (extent == "full"){
    if ("dd" %in% class(model)){
      parms <- matrix(c(model$beta, model$parms), nrow = 1,
        dimnames = list(NULL, c(names(model$beta), names(model$parms))))
    } else if ("ddSim" %in% class(model)){
      parms <- model
    }    # this assume that all the parameter sets are proper
    output <- switch(distr,
      gamma = pgamma(x, shape = parms[, "shape"], rate = parms[, "rate"]),
      lognormal = plnorm(x, meanlog = parms[, "meanlog"], sdlog = parms[, "sdlog"]),
      logLinear = plogL(x, b1 = parms[, "b1"]),
      logQuadratic = plogQ(x, b1 = parms[, "b1"], b2 = parms[, "b2"]),
      logCubic = {
        b1 <- parms[, "b1"]; b2 <- parms[, "b2"]; b3 <- parms[, "b3"]
        const <- try(1/integrate(
          f = function(r) r * exp(b1 * r + b2 * r^2 + b3 * r^3),
          lower = 0, upper = Inf, rel.tol = .Machine$double.eps^0.5)$val,
         silent = TRUE)
        if ("try-error" %in% class(const)) stop("invalid parameter set for ",
          distr, ". Integral for normalizing constant does not converge.")
        plogC(x, b1 = b1, b2 = b2, b3 = b3, const = const)
      },
      inverse_gamma = pigam(x, shape = parms[, "shape"], scale = parms[, "scale"]),
      paranormal_gamma = { # assumes b0, b1, b2, b3 are scalars
        b0 <- parms[, "b0"]
        b1 <- parms[, "b1"]
        b2 <- parms[, "b2"]
        b3 <- parms[, "b3"]
        const <- try(1/integrate(
          f = function(r) r * exp(b0 * log(r) + b1 * r + b2 * r^2 + b3 * r^3),
          lower = 0, upper = Inf, rel.tol = .Machine$double.eps^0.5)$val,
         silent = TRUE)
        if ("try-error" %in% class(const)) stop("invalid parameter set for ",
          distr, ". Integral for normalizing constant does not converge.")
        ppng(x, b0 = b0, b1 = b1, b2 = b2, b3 = b3, const = const)
      },
      Rayleigh = pRay(x, s2 = parms[, "s2"]),
      MaxwellBoltzmann = pmb(x, a = parms[, "a"]),
      constant = NA,
      Pareto = pPare1(x, a = parms),
      tnormal = {
        mu <- parms[, "mean"]
        sig  <- parms[, "sd"]
        Fa <- pnorm(0, mean = mu, sd = sig)
        ans <- numeric(length(x))
        ans[x >= 0] <- (pnorm(x[x >= 0], mean = mu, sd = sig) - Fa)/(1 - Fa)
        ans
      },
      exponential = pexp(x, rate = parms[, "rate"]),
      inverse_gaussian = statmod::pinvgauss(x,
        mean = parms[, "mean"], dispersion = parms[, "dispersion"]),
      chisq = pchisq(x, df = parms[, "df"])
    )
  } else if (extent == "win"){ # truncation radius is numeric scalar 0 < trad < Inf
    if ("dd" %in% class(model)){
      cof <- matrix(model$dmod$coefficients[cof_name[[distr]]], nrow = 1,
        dimnames = list(NULL, cof_name[[distr]]))
      trad <- max(model$dmod$data$r) # truncation radius
    } else if ("ddSim" %in% class(model)){
      cof <- model[, cof_name[[distr]]]
      trad <- attr(model, "trad")
    }    # this assume that all the parameter sets are proper
    y <- numeric(max(length(x), nrow(cof)))
    if (length(x) == 1){
      xx <- rep(x, length.out = nrow(cof))
    } else {
      xx <- x
    }
    xx[xx == 0] <- 0.001 # reduce jargon with this simplification
    output <- switch(distr,
      gamma = {  # vectorized
        b1 <- cof[, "log(r)"]
        b2 <- cof[, "r"]
        y <- (xx^b1 * (-b2*xx)^(-b1) * (gamma(b1 + 2) - expint::gammainc(b1 + 2, -b2*xx)))/
          (trad^b1 * (-b2*trad)^(-b1) * (gamma(b1 + 2) - expint::gammainc(b1 + 2, -b2*trad)))
        y[xx >= trad] <- 1
        y[xx <= 0.001] <- 0
        y
      },
      lognormal = { # semi-vectorized
        if (nrow(cof) == 1){
          if (cofOK(cof, distr)){
            b <- cof2parms(cof, distr)
            y <- plnorm(x, meanlog = b["meanlog"], sdlog = b["sdlog"])/
              plnorm(trad, meanlog = b["meanlog"], sdlog = b["sdlog"])
          } else {
            deno <- integrate(f = function(r)
              exp(rmat(r, distr) %*% cof[,cof_name[[distr]]] + off(r, distr)),
              lower = 0, upper = trad)$val
            y[x > 0] <- rmutil::int(function(r)
              exp(rmat(r, distr) %*% cof[,cof_name[[distr]]] + off(r, distr)),
              a = 0, b = x[x > 0])
            y[x < 0] <- 0
            y[x >= trad] <- 1
          }
        } else if (length(x) == 1){
          if (x <= 0) {y <- y} else if (x >= trad) {y <- y + 1} else {
          ip <- which(cofOK(cof))
          b <- mod2distr(cof, distr)
          if (length(ip) > 0)
            y[ip] <- plnorm(x, meanlog = b[, "meanlog"], sdlog = b[, "sdlog"])/
                plnorm(trad, meanlog = b[, "meanlog"], sdlog = b[, "sdlog"])
          if (length(ip) < nrow(cof)){
            nip <- (1:nrow(cof))[-ip]
            for (i in nip){
              deno <- integrate(f = function(r)
                exp(rmat(r, distr) %*% cof[i, cof_name[[distr]]] + off(r, distr)),
                lower = 0, upper = trad)$val
              y[i] <- integrate(f = function(r)
                exp(rmat(r, distr) %*% cof[i, cof_name[[distr]]] + off(r, distr)),
                lower = 0, upper = x)$val
            }
          }}
        } else if (length(x) == nrow(cof)){
          ip <- which(cofOK(cof))
          b <- mod2distr(cof, distr)
          if (length(ip) > 0)
            y[ip] <- plnorm(x[ip], meanlog = b[, "meanlog"], sdlog = b[, "sdlog"])/
                plnorm(trad, meanlog = b[, "meanlog"], sdlog = b[, "sdlog"])
          if (length(ip) < nrow(cof)){
            nip <- (1:nrow(cof))[-ip]
            for (i in nip){
              deno <- integrate(f = function(r)
                exp(rmat(r, distr) %*% cof[i, cof_name[[distr]]] + off(r, distr)),
                lower = 0, upper = trad)$val
              y[i] <- integrate(f = function(r)
                exp(rmat(r, distr) %*% cof[i, cof_name[[distr]]] + off(r, distr)),
                lower = 0, upper = x[i])$val
            }
          }
        }
      },
      logLinear = { # vecorized
        b1 <- cof[, "r"]
        y <- (exp(b1*x) * (b1*x  - 1) + 1)/(exp(b1*trad) * (b1*trad - 1) + 1)
        y[xx >= trad] <- 1
        y[xx < 0] <- 0
        y
      },
      logQuadratic = { # vectorized
        b1 <- cof[, "r"]
        b2 <- cof[, "I(r^2)"]
        supi <- which(xx < trad & xx > 0)
        y[supi] <- lQint(xx[supi], b1, b2)/lQint(trad, b1, b2)
        y[xx >= trad] <- 1
        y[xx <= 0] <- 0
        y
      },
      logCubic = { # not vectorized; all must be numerically integrated
        b1 <- cof[, "r"]
        b2 <- cof[, "I(r^2)"]
        b3 <- cof[, "I(r^3)"]
        if (length(x) == 1 | length(x) == length(b1)){
          x <- rep(x, length.out = length(b1))
          for (bi in 1:length(b1)){
            y[bi] <-
              integrate(f = function(r)
                r * exp(b1[bi]*r + b2[bi]*r^2 + b3[bi]*r^3),
                lower = 0, upper = x)$val/
              integrate(f = function(r)
                r * exp(b1[bi]*r + b2[bi]*r^2 + b3[bi]*r^3),
                lower = 0, upper = trad)$val
          }
        } else if (length(b1) == 1){
          supi <- which(x < trad & x > 0)
          # rmutil::int is 50x faster than "integrate" on vectorized bounds
          y[supi] <- rmutil::int(f = function(r) # not vectorized wrt b1, b2, b3
              r * exp(b1*r + b2*r^2 + b3*r^3), a = 0, b = x[supi])/
            integrate(f = function(r)
              #integrate is 3x faster than rmutil::int on scalar bounds
              r * exp(b1*r + b2*r^2 + b3*r^3), lower = 0, upper = trad)$val
        }
        y[x >= trad] <- 1
        y[x <= 0] <- 0
        y
      },
      inverse_gamma = {
        b1 <- cof[, "I(1/r)"]
        b2 <- cof[, "log(r)"]
        deno <- numeric(length(b1))
        if (length(x) == 1 | length(x) == length(b1)){
          for (bi in 1:length(b1)){
            y[bi] <- integrate(f = function(r)
              r * exp(b1[bi]/r + b2[bi]*log(r)), lower = 0, upper = x)$val
            y[bi] <- y[bi]/integrate(f = function(r)
              r * exp(b1[bi]/r + b2[bi]*log(r)), lower = 0, upper = trad)$val
          }
        } else if (length(b1) == 1){
          supi <- which(x < trad & x > 0)
          y[supi] <- rmutil::int(f = function(r) # vectorized, 50x faster than "integrate"
            r * exp(b1/r + b2*log(r)), a = 0, b = x[supi])
          y[supi] <- y[supi]/integrate(f = function(r)
            r * exp(b1/r + b2*log(r)), lower = 0, upper = trad)$val
        }
        y[x >= trad] <- 1
        y[x <= 0] <- 0
        y
      },
      paranormal_gamma = {
        b0 <- cof[, "log(r)"]
        b1 <- cof[, "r"]
        b2 <- cof[, "I(r^2)"]
        b3 <- cof[, "I(r^3)"]
        if (length(x) == 1 | length(x) == length(b1)){
          x <- rep(x, length.out = length(b1))
          for (bi in 1:length(b1)){
            y[bi] <-
              integrate(f = function(r)
                r * exp(b0[bi] * log(r) + b1[bi]*r + b2[bi]*r^2 + b3[bi]*r^3),
                lower = 0, upper = x[bi])$val/
              integrate(f = function(r)
                r * exp(b0[bi] * log(r) + b1[bi]*r + b2[bi]*r^2 + b3[bi]*r^3),
                lower = 0, upper = trad)$val
          }
        } else if (length(b1) == 1){
          supi <- which(x < trad & x > 0)
          # rmutil::int is vectorized and 50x faster than "integrate"
          y[supi] <- rmutil::int(f = function(r) # not vectorized wrt b1, b2, b3
              r * exp(b0*log(r) + b1*r + b2*r^2 + b3*r^3), a = 0, b = x[supi])/
            integrate(f = function(r)
              r * exp(b0*log(r) + b1*r + b2*r^2 + b3*r^3),
              lower = 0, upper = trad)$val
        }
        y[x >= trad] <- 1
        y[x <= 0] <- 0
        y
      },
      Rayleigh = {
        b1 <- cof[, "I(r^2)"]
        y <- (exp(b1*x^2) - 1)/(exp(b1*trad^2) - 1)
        y[x >= trad] <- 1
        y[x <= 0] <- 0
        y
      },
      Pareto = {
        b1 <- cof[, "log(r)"]
        y <- (x^(b1 + 2) - 1)/(trad^(b1 + 2) - 1)
        y[x >= trad] <- 1
        y[x <= 1] <- 0
        y
      },
      MaxwellBoltzmann = {
        b1 <- cof[, "I(r^2)"]
        if (length(b1) == 1){
          b1 <- rep(b1, length.out = length(x))
        } else if (length(x) == 1){
          x <- rep(x, length.out = length(b1))
        }
        xp <- which(x > 0 & x <= trad)
        y[xp] <-
          (x[xp]*exp(b1[xp]*x[xp]^2)/(2*b1[xp]) - sqrt(pi)/(4*(b1[xp] + 0i)^1.5)*
            pracma::erfi(x[xp] * sqrt(b1[xp] + 0i)))/
          (trad *exp(b1[xp]*trad^2)/(2*b1[xp]) - sqrt(pi)/(4*(b1[xp] + 0i)^1.5)*
            pracma::erfi(trad * sqrt(b1[xp] + 0i)))
        y[x >= trad] <- 1
        y[x <= 0] <- 0
        y <- Re(y)
        y
      },
      constant = {
        y <- (x/trad)^2
        y[x >= trad] <- 1
        y[x <= 0] <- 0
        y
      },
      tnormal = {
        b1 <- cof[, "r"]
        b2 <- cof[, "I(r^2)"]
        erfib <- pracma::erfi(0.5*b1/sqrt(b2 + 0i))
        y <- Re((pracma::erfi((x * b2 + 0.5 * b1)/sqrt(b2 + 0i)) - erfib)/
            (pracma::erfi((trad * b2 + 0.5 * b1)/sqrt(b2 + 0i)) - erfib))
        y[x >= trad] <- 1
        y[x <= 0] <- 0
        y
      },
      exponential = {
        b1 <- cof[, "r"]
        y <- (exp(x * b1) - 1)/(exp(trad * b1) - 1)
        y[x >= trad] <- 1
        y[x <= 0] <- 0
        y
      },
      inverse_gaussian = {
        b1 <- cof[, "I(1/r)"]
        b2 <- cof[, "r"]
        if (length(x) == 1 | length(x) == length(b1)){
          x <- rep(x, length.out = length(b1))
          for (bi in 1:length(b1)){
            y[bi] <-
              integrate(f = function(r)
                exp(b1[bi]/r + b2[bi]*r - 3/2 * log(r)),
                lower = 0, upper = x[bi])$val/
              integrate(f = function(r)
                exp(b1[bi]/r + b2[bi]*r - 3/2 * log(r)),
                lower = 0, upper = trad)$val
          }
        } else if (length(b1) == 1){
          supi <- which(x < trad & x > 0)
          # rmutil::int is vectorized and 50x faster than "integrate"
          y[supi] <- rmutil::int(f = function(r) # not vectorized wrt b1, b2, b3
            exp(b1/r + b2*r - 3/2 * log(r)), a = 0, b = x[supi])/
            integrate(f = function(r)
              exp(b1/r + b2*r - 3/2 * log(r)), lower = 0, upper = trad)$val
          y[x >= trad] <- 1
          y[x <= 0] <- 0
        }
        y
      },
      chisq = {
        b1 = cof[, "log(r)"]
        supi <- which(x < trad & x > 0)
        y[supi] <- (gamma(b1 + 2) - expint::gammainc(b1 + 2, x[supi]/2))/
             (gamma(b1 + 2) - expint::gammainc(b1 + 2, trad/2))
        y[x >= trad] <- 1
        y[x <= 0] <- 0
        y
      }
    )
    parms <- cof
  } else {
    stop("extent must be \"full\" or \"win\"")
  }
  if (nrow(parms) > 1){
    output <- matrix(unname(output), nrow = nrow(parms))
  } else {
    output <- unname(as.vector(output))
  }
  output
}

#' Function to calculate relative carcass densities
#' @param x numeric x >= 0
#' @param model either a \code{dd} object or a set of simulated parameters
#'  (i.e., a \code{ddSim} object)
#' @param trad the radius to use for truncated distributions
#'  (\code{trad = NULL} to use full, non-truncated distribution)
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
rcd <- function(x, model, extent = "full"){
  x[x < 0.001] <- 0.001
  if ("ddSim" %in% class(model)){
    parms <- model
    distr <- attr(model, which = "distr")
    trad <- attr(model, which = "trad")
  } else if ("dd" %in% class(model)){
    if (extent == "full"){
      parms <- matrix(c(model$beta, model$parms), nrow = 1,
       dimnames = list(NULL, c(names(model$beta), names(model$parms))))
    } else if (extent == "win"){
      parms <- matrix(model$dmod$coefficients, nrow = 1,
       dimnames = list(NULL, names(model$dmod$coefficients)))
      trad <- max(model$dmod$data$r, na.rm = TRUE)
      attr(parms, "trad") <- trad
    }
    attr(parms, "distr") <- model$distr
    class(parms) <- "ddSim"
    distr <- model$distr
  }
  if (distr == "constant" & extent == "full")
    warning("rcd cannot be calculated for extent = \"full\" and distr = ", distr)
  if (distr %in% c("exponential", "tnormal"))
    warning("rcd cannot be calculated for distr = ", distr)
  upr <- ifelse(extent == "full", Inf, trad)
  output <- array(dim = c(length(x), nrow(parms)))
  for (pri in 1:nrow(parms)){
    output[, pri] <- as.vector(ddd(x = x, model = parms[pri, ], extent = extent))
    const <- 1/integrate(f = function(x)
      2 * pi * c(ddd(x, model = parms[pri,], extent = extent))/x,
      lower = 0, upper = upr)$val
    output[, pri] <- const * output[, pri]/x
  }
  if (length(x) == 1 | nrow(parms) == 1) output <- as.vector(output)
  return(output)
}

#' @export
summary.dd <- function(x, trad = NULL, ...){
  distr <- x$distr
  if (distr == "constant") return(NA)
  cof <- numeric(length(cof_name[[distr]][-1])) + NA
  names(cof) <- cof_name[[distr]][-1]
  qtls <- c(0.5, 0.75, 0.9, 0.95)
  names(qtls) <- c("median", "75%", "90%", "95%")
  ans <- list(
    model = data.frame(array(dim = c(1, 2 * length(cof) + 1))),
    stats = numeric(5) + NA
  )
  rownames(ans$model) <- distr
  colnames(ans$model) <- c(parm_name[[distr]], names(cof), "offset")
  ans[["model"]][1, "offset"]  <- mod_offset[distr]
  names(ans[["stats"]]) <- c(names(qtls), "mode")
  if (any(is.na(x$beta))) return(ans)

  cof <- x$beta[cof_name[[distr]][-1]]
  parms <- x$parms
  ans[["model"]][1, names(c(cof, parms))] <- c(cof, parms)
  if (cofOK(cof, distr)){
    for (i in names(qtls)){
      ans[["stats"]][i] <- tryCatch(uniroot(f = function(r) pdd(r, x) - qtls[i],
        lower = 1, upper = 5000)$root, error = function(e) NA, silent = TRUE)
      if ("try-error" %in% class(ans[["stats"]][i])){
        ans[["stats"]][i:4] <- NA
        break
      }
    }
    xmax <- 1000
    xseq <-seq(0.1, xmax, by = 0.1)
    y <- ddd(x = xseq, model = x)
    mi <- which(y == max(y))
    ans[["stats"]]["mode"] <- ifelse(y[length(y)] %in% mi, NA, mean(xseq[mi]))
    ans[["stats"]] <- round(ans[["stats"]], 2)
  } else {
    ans[["stats"]] <- NA
  }
  ans
}

#' @export
summary.ddArray <- function(x, trad = NULL, ...){
  aic0 <- aic(x)
  ans <- data.frame(array(dim = c(length(x), 6)))
  rownames(ans) <- aic0$model
  colnames(ans) <- c("median", "75%", "90%", "95%", "mode", "deltaAICc")
  ans$deltaAICc <- aic0$deltaAICc
  for (distr in names(x)){
    if (distr == "constant") next
    tmp <- summary(x[distr])
    ans[distr, c("median", "75%", "90%", "95%", "mode")] <- tmp$stats
  }
  ans
}

#' @export
dd2ddSim <- function(dd){
  parms <- matrix(c(dd$beta, dd$parms), nrow = 1,
    dimnames = list(NULL, c(names(dd$beta), names(dd$parms))))
  attr(parms, "distr") <- dd$distr
  class(parms) <- "ddSim"
  parms
}

#' @rdname ddd
#' @export
qdd <- function(p, model, extent = "full", subdiv = 1000){ # model is ddSim or dd
# find r such that pdd(r) = p
  if("dd" %in% class(model)){
    trad <- max(model$dmod$data$r, na.rm = TRUE)
  } else if ("ddSim" %in% class(model)){
    trad <- attr(model, "trad")
  }
  if ("dd" %in% class(model) || nrow(model) == 1){
    if (length(p) == 1){
      qval = try(uniroot(
        f = function(r) pdd(r, model = model, extent = extent) - p,
        lower = 0, upper = ifelse(extent == "full", 50000, trad))$root,
      silent = TRUE)
      if ("try-error" %in% class(qval))
        stop(paste0("cannot invert CDF in qdd for p = ", p))
      return(qval)
    }
    minr <- try(uniroot( # p is a vector, solve by linear interpolation
        f = function(r) pdd(r, model = model, extent = extent) - min(p),
        lower = 0,
        upper = ifelse(extent == "full", 50000, trad))$root,
      silent = TRUE)
    if ("try-error" %in% class(minr))
      stop(paste0("cannot invert CDF in qdd for p = ", min(p)))
    maxr <- try(uniroot(
        f = function(r) pdd(r, model = model, extent = extent) - max(p),
        lower = 0, upper = ifelse(extent == "full", 50000, trad))$root,
      silent = TRUE)
    if ("try-error" %in% class(maxr))
      stop(paste0("cannot invert CDF in qdd for p = ", max(p)))
    y <- seq(minr, maxr, length.out = subdiv)
    x <- pdd(x = y, model = model, extent = extent)
    y <- c(0, y, min(50000, trad))
    x <- c(0, x, 1)
    return(approxfun(x, y)(p))
  } else if ("ddSim" %in% class(model) && nrow(model) > 1){
    pp <- rep(p, length.out = nrow(model))
    qval <- numeric(nrow(model))
    for (ppi in 1:length(pp)){
      qval[ppi] = tryCatch(
        uniroot(f = function(r)
          pdd(r, model = model[ppi, ], extent = extent) - pp[ppi],
          lower = 0, upper = ifelse(extent == "full", 50000, trad))$root,
        error = function(e) NA
      )
    }
    return(qval)
  } else {
    stop("dimension mismatch between p and simulated parameters (= model) in qdd")
  }
}

#' @export
rdd <- function(n, model, extent = "full", subdiv = 1000){ # model is ddSim or dd
  if ("dd" %in% class(model)){
    parms <- dd2ddSim(model)
  } else if ("ddSim" %in% class(model)){
    if (nrow(model) > 1 & nrow(model) != n)
      stop("nrow(model) must be 1 or n in rdd")
  } else {
    stop ("class(model) in rdd must be dd or ddSim")
  }
  qdd(runif(n), model = model, extent = extent, subdiv = subdiv)
}

#' fit a truncated distribution from a \code{dmod} object
#' @param x \code{dmod} object, which is a \code{glm} object from a Poisson
#'  with the distance variable named "r", the search class named "svar", and
#'  the name of the distribution added to the glm list objects as a character
#'  string.
#' @return a fitted distribution for carcass v distance
#' @export
tfit <- function(x, ...) UseMethod("tfit", x)

#' @rdname tfit
#' @export
tfit.dmod <- function(x,...){
  dmod <- x
  distr <- dmod$distr
  # trad should not be an option? Just clip at end of data instead because the
  # data are already in ringData format
  if (distr == "constant") {
    output <- list(distr = distr, parms = NA, beta = NA, varbeta = NA,
      ncarc = sum(dmod$data$ncarc), n = nrow(dmod$data), k = length(dmod$coef),
      llik = NA, dmod = dmod)
    class(output) <- "dd"
    return(output)
  }
  trad <- max(dmod$data$r)
  # starting pt for optim
##  tmp <- ddInit(dmod)
#  sc0 <- tmp[grepl("svar", names(tmp))] # svar not reoptimized
#  par0 <- tmp[!grepl("svar", names(tmp))] # distance parameters
#  mmat <- cbind(model.matrix(dmod), exposure = dmod$data$exposure, ncarc = dmod$data$ncarc)
#  extra <- numeric(nrow(mmat))
#  if (distr %in% c("exponential", "tnormal")) extra <- -log(dmod$data$r)
#  if (distr == "chisq") extra <- -dmod$data$r/2
#  if (distr == "inverse_gaussian") extra <- -2.5*log(dmod$data$r)
#  if (distr == "MaxwellBoltzmann") extra <- log(dmod$data$r)
#  fit <- tryCatch(
#    optim(par0, fn = llik_fun,
#      mmat = mmat, trad = trad, sc0 = sc0, distr = distr, extra = extra,
#      method = "L-BFGS-B",
#      lower = constraints[[distr]][, "lower"],
#      upper = constraints[[distr]][, "upper"],
#      control = list(fnscale = -1, parscale = constraints[[distr]][, "parscale"])
#    ),
#    error = function(e) NA
#  )
#  if (length(fit) == 1 ||
#      !is.finite(tryCatch(integrate(f = function(r)
#          exp(rmat(r, distr) %*% fit$par + off(r, distr)),
#          lower = ifelse(distr == "Pareto", 1, 0),
#          upper = Inf,
#          rel.tol = sqrt(.Machine$double.eps)
#        )$val, error = function(e) NA)) ||
#      !cofOK(fit$par, distr)){
#    output <- list(
#      distr = distr,
#      parms = NA,
#      beta = NA,
#      varbeta = NA,
#      ncarc = sum(dmod$data$r),
#      n = nrow(dmod$data),
#      k = length(dmod$coef),
#      llik = NA,
#      dmod = dmod
#    )
#  } else {  #
#    varbeta <- tryCatch(
#      -solve(numDeriv::hessian(llik_fun, x = fit$par,
#        mmat = mmat, trad = trad, sc0 = sc0, distr = distr, extra = extra)),
#      error = function(e) array(dim = c(length(par0), length(par0)))
#    )
#    rownames(varbeta) <- names(fit$par)
#    colnames(varbeta) <- names(fit$par)
    output <- list(
      distr = distr,
      parms = cof2parms(x$coefficients, distr),
      beta = x$coefficients[cof_name[[distr]]],
      varbeta = summary(x)$cov.unscaled,
      ncarc = sum(x$data$r),
      n = nrow(x$data),
      k = length(x$coef),
      dmod = x
    )
#  }
  class(output) <- "dd"
  output
}

#' @export
llik_fun <- function(x, mmat, trad, sc0, distr, extra = 0){
  logmu <- mmat[, 1:(ncol(mmat) - 2)] %*% c(x, sc0) + # <= main form
    log(mmat[, "exposure"]) + extra # <= offset
  sum(
    mmat[, "ncarc"] * logmu - exp(logmu) - # likelihood
    log(integrate(f = # denominator in g(x)/F(b) for truncated distribution
      function(r) exp(rmat(r, distr) %*% x + off(r, distr)),
      lower = ifelse(distr == "Pareto", 1, 0),
      upper = trad,
      rel.tol = sqrt(.Machine$double.eps))$val
    )
  )
}
#' @export
ddInit <- function(x, ...) UseMethod("ddInit", x)

# x is either the data or a dd object
# if data, then fit a dmod
#   -- use the unedited fitted parameters if the critical parameter <0;
#   -- otherwise, use the parameters for the svar weights and a generic set of
#     parameters that yield something like a mean of 20 and sd of 20 for the
#     distance parameters
#' @export
ddInit.data.frame <- function(x, distr,
    rCol = "r", scCol = NULL, expoCol = "exposure", ncarcCol = "ncarc", ...){
  ringData <- x
  data <- ringData[, c(rCol, scCol, expoCol, ncarcCol)]
  names(data) <- c("r", if(!is.null(scCol)) "svar", "exposure", "ncarc")
  plu <- ifelse(!is.null(scCol), " + ", " ")
  sufx <- paste0(plu, scCol,  "+ offset(log(", expoCol, "))")
  sufx0 <- paste0(plu, scCol, "+ ")
  prfx <- paste0("ncarc ~ ")
  dmod <- switch(distr,
    gamma = glm(formula(paste0("ncarc ~ ", "log(r) + r", sufx)),
      data = ringData, family = "poisson"),
    lognormal = glm(formula(paste0("ncarc ~ ", "log(r) + I(log(r)^2)", sufx)),
      data = data, family = "poisson"),
    logLinear = glm(formula(paste0("ncarc ~ ", "r", sufx)),
      data = data, family = "poisson"),
    logQuadratic = glm(formula(paste0("ncarc ~", "r + I(r^2)", sufx)),
      data = data, family = "poisson"),
    logCubic = glm(formula(paste0("ncarc ~ ", "r + I(r^2) + I(r^3)", sufx)),
      data = data, family = "poisson"),
    inverse_gamma = glm(formula(paste0("ncarc ~ ", "I(1/r) + log(r)", sufx)),
      data = data, family = "poisson"),
    paranormal_gamma = glm(formula(paste0("ncarc ~ ", "log(r) + r + I(r^2) + I(r^3)", sufx)),
      data = data, family = "poisson"),
    Rayleigh =  lm(formula(paste0("ncarc ~ ", "I(r^2)", sufx)),
      data = data, family = "poisson"),
    Pareto = glm(formula(paste0("ncarc ~ ", "log(r)", sufx)),
      data = data, family = "poisson"),
    MaxwellBoltzmann = glm(formula(
      paste0("ncarc ~ ", "I(r^2)", sufx0, "offset(I(log(exposure * r)))"),),
      data = data, family = "poisson"),
    constant = glm(formula = formula(ifelse(is.null(scCol),
        paste0("ncarc ~ ", "offset(log(", expoCol, "),)"),
        paste0("ncarc ~ ", scCol, " + offset(log(", expoCol, "),)"),)),
      data = data, family = "poisson"),
    tnormal = glm(formula(paste0("ncarc ~ ", "r + I(r^2)", sufx0,
      "offset(I(log(exposure/r)))"),),
      data = data, family = "poisson"),
    exponential =  lm(formula(paste0("ncarc ~ ", "r", sufx0,
      "offset(I(log(exposure/r)))"),),
      data = data, family = "poisson"),
    inverse_gaussian = glm(formula(paste0("ncarc ~ ", "I(1/r) + r", sufx0,
      "offset(I(log(exposure) - 5/2 * log(r)))"),),
      data = data, family = "poisson"),
    chisq = glm(formula(paste0("ncarc ~ ", "log(r)", sufx0,
      "offset(I(log(exposure) - r/2))"),),
      data = data, family = "poisson")
  )
  par0 <- dmod$coef
  if (par0[critical_parameter[distr]] >= 0){
    par0[1:which(names(cof0) == critical_parameter[distr])] <- 0
    par0[1] <- ifelse(!grepl("inverse", distr), -5, 2)
    if (distr == "chisq") par0[1] <- -40
    par0[critical_parameter[distr]] <- icp[distr]
  }
  return(par0)
}

#' @export
ddInit.dmod <- function(x, ...){ #dmod is glm w/distr and data has standard names
  data <- x$data
  bb <- x$coef
  good <- TRUE
  pchk <- constraints[[x$distr]] # parameter check array
  for (bi in rownames(pchk)[-1])
    if (bb[bi] < pchk[bi, "lower"] || bb[bi] > pchk[bi, "upper"]) good <- FALSE
  if (!good){
    bb[1:which(names(bb) == critical_parameter[x$distr])] <- 0
    bb[1] <- ifelse(!grepl("inverse", x$distr), -5, 2)
    if (x$distr == "chisq") bb[1] <- -40
    if (x$distr %in% c("inverse_gamma", "inverse_gaussian")) bb["I(1/r)"] <- -10
    bb[critical_parameter[x$distr]] <- icp[x$distr]
  }
  return(bb)
}
