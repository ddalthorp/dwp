#' Create a shapeLayout
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
#' @param scCol name of column in file_layout with names of search classes. This
#'  is used for excluding unsearched areas from the grid data (x, y). It is used
#'  ONLY with \code{dataType = "xy"} and used to remove rows with
#'  \code{x[, scCol] == notSearched}, where \code{x} is the search grid data frame.
#' @param notSearched a companion to \code{scCol} for removing unsearched areas
#'  from the search grid.
#' @return A \code{shapeLayout} object, which is a list with spatial
#'  characteristics of the site, including
#'  \describe{
#'    \item{\code{$layout}}{turbine search area configurations (polygons and multipolygons)
#'      from \code{file_layout} shape file as an \code{sf} object.}
#'    \item{\code{$layoutAdj}}{polygons from \code{$layout} but all recentered at (0, 0)}
#'    \item{\code{$turbines}}{turbine centers (as \code{sf} object)}
#'    \item{\code{$carcasses}}{carcass discovery locations along with (optional)
#'      incidental information (turbine, visibility class, species,)}
#'    \item{\code{$unitCol, $ncarc, $tset, $tcenter}}{name of the column with turbine IDs
#'      (\code{$unitCol}, character string), number of carcasses found at each turbine
#'      (\code{$ncarc}, vector), turbine names (\code{$tset}, vector of character strings),
#'      locations of turbine centers (\code{$tcenter}, nturb by 2 array) UTMs of
#'      turbine centers, extracted and simplified from \code{$turbines}}
#'  }
#' @export
#'
readLayout <- function(file_layout, dataType = "simple",  unitCol = "turbine",
    file_turbine = NULL, radCol = "radius", shapeCol = "shape", padCol = "padrad",
    roadwidCol = "roadwidth", nRoadCol = "n_road", xCol = "x", yCol = "y",
    scCol = NULL, notSearched = NULL){
  # traffic directing
  dataType <- tolower(dataType)
  if (dataType == "gis" || (length(file_layout) == 1 &&
      is.character(file_layout) && grepl(".shp", file_layout))){
    plotLayout <- sf::st_read(file_layout, stringsAsFactors = FALSE)
    plotLayout0 <- sf::st_drop_geometry(plotLayout) # easier quick subsetting
    turbines <- sf::st_read(file_turbine, stringsAsFactors = FALSE)
    turbines0 <- sf::st_drop_geometry(turbines)
    if (!is.null(unitCol) && !is.na(unitCol)){
      if (is.null(file_turbine))
        stop("if unitCol is not NULL, file_turbine must be provided")
      if (!unitCol %in% names(turbines))
        stop(unitCol, " column not included in turbine data")
      if (!unitCol %in% names(plotLayout))
        stop(unitCol, " column not included in plot layout data")
      if (!all(plotLayout0[, unitCol] %in% turbines0[, unitCol])){
        stop("a turbine that is included in ", file_layout,
             " is missing from ", file_turbine
        )
      }
      # format of unitCol names: must be syntactically valid R names
      # make.names(x) = x if x is a valid name, else converts x to a valid name
      if (!all(make.names(turbines0[, unitCol]) == turbines0[, unitCol])){
        badind <- which(make.names(turbines0[, unitCol]) != turbines0[, unitCol])
        badnm <- unique(turbines0[badind, unitCol])
        badnm <- badnm[1:min(length(badnm), 3)]
        goodnm <- gsub("^X", "t", badnm)
        message(
          "\n\nNOTE: Not all the turbine names are syntactically valid.\n",
          "These (", paste(badnm, collapse = ", "),
          "...) can be converted to syntactically valid names (", goodnm, "...) ",
          "if desired."
        )
        tmp <- readline("convert [<enter>]? or not [n]? ")
        if (!identical(tolower(tmp), "y") & tmp != "") stop("aborting...")
        mknm <- make.names(turbines0[, unitCol])
        badind <- which(mknm != turbines0[, unitCol])
        badnm <- mknm[badind]
        turbines0[badind, unitCol] <- gsub("^X", "t", badnm)
        turbines[, unitCol] <- turbines0[, unitCol]
        mknm <- make.names(plotLayout0[, unitCol])
        badind <- which(mknm != plotLayout0[, unitCol])
        badnm <- mknm[badind]
        plotLayout[badind, unitCol] <- gsub("^X", "t", badnm)
      }
      tset <- as.character(turbines0[, unitCol]) # turbine names
      tcenter <- sf::st_coordinates(turbines) # turbine centers (matrix with "X", "Y")
      rownames(tcenter) <- tset
    } else {
      stop("readLayout with shape files must have file_turbine and unitCol.")
    }
    layoutAdj <- plotLayout # Adj => coordinates relative to turbine centers
    for (ai in 1:nrow(plotLayout)){
      ti <- plotLayout[ai, unitCol, drop = TRUE]
      layoutAdj[ai, ] <- sf::st_set_geometry(layoutAdj[ai, ],
        sf::st_geometry(plotLayout[ai, ]) - tcenter[ti,]) # recentering
    }
    # maximum distance from turbine to a searched point (radius of interpolation)
   output <- list(layout = plotLayout, layoutAdj = layoutAdj, turbines = turbines,
      unitCol = unitCol, tset = tset, tcenter = tcenter)
    class(output) <- "shapeLayout"
    return(output)
  } else if (dataType == "simple"){
    if (length(file_layout) == 1 && is.character(file_layout)){
      slayout <- read.csv(file_layout, stringsAsFactors = FALSE)
    } else if (is.data.frame(file_layout)){
      slayout <- file_layout
    }
    if ("simpleLayout" %in% class(file_layout)) return(file_layout)
    # rudimentary error-checking:
    if (!all(c(unitCol, radCol, shapeCol) %in% names(slayout))){
      stop("Simple layout must include columns for unit, radius, and shape ",
           "with column names specified in arg list for readLayout.")
    }
    if (any(slayout[, shapeCol] == "RP") &&
       !all(c(padCol, roadwidCol, nRoadCol) %in% names(slayout))){
      stop(
        "Simple layout that includes RP turbines must include columns for pad ",
        "radius, road width, and number of roads in the arg list for readLayout."
      )
    }
    if (!is.numeric(slayout[, radCol]) || any(slayout[, radCol] <= 0))
      stop("All search radii in simple layout must be positive numbers.")
    if (!all(slayout[, shapeCol] %in% c("circular", "square", "RP")))
      stop("Simple \"shape\" must be circular, square, or RP")
    srad <- slayout[ , radCol]
    ind <- which(slayout[, shapeCol] == "square")
    srad[ind] <- sqrt(2) * srad[ind]
    colnm <- c(unitCol = unitCol, radCol = radCol, shapeCol = shapeCol,
      padCol = padCol, roadwidCol = roadwidCol, nRoadcol = nRoadCol)
    slayout <- slayout[, colnm]
    names(slayout) <- c("turbine", "radius", "shape", "padrad", "roadwidth", "n_road")
    attr(slayout, "srad") <- ceiling(max(srad))
    attr(slayout, "colnm") <- colnm
    class(slayout) <- c("simpleLayout", "data.frame")
    return(slayout)
  } else if (dataType == "polygon"){
    if (is.data.frame(file_layout) || is.matrix(file_layout)){
      playout <- file_layout
    } else if (length(file_layout) == 1 && is.character(file_layout)){
      playout <- read.csv(file_layout, stringsAsFactors = FALSE)
    } else {
      stop("file_layout must be a name of a file; alternatively, ",
           "file_layout may be a properly formatted data frame or matrix"
      )
    }
    if (!is.matrix(playout) && !is.data.frame(playout))
      stop("Data in file_layout must be a properly formatted matrix or data frame")
    if (!xCol %in% colnames(playout) || !yCol %in% colnames(playout))
      stop("Data in file_layout must have columns for x and y coordinates")
    if (is.data.frame(playout)){
      rnm <- playout[, unitCol]
      playout <- as.matrix(playout[, c(xCol, yCol)], ncol = 2)
      rownames(playout) <- rnm
    } else {
      if (is.null(row.names(playout)))
        stop("If file_layout is a matrix, row names must be turbine names")
      if (min(table(rownames(playout))) < 3)
        stop("file_layout data must be polygons")
      playout <- playout[, c(xCol, yCol)]
      rnm <- rownames(playout)
    }
    pgon <- list()
    for (ti in unique(rnm)) pgon[[ti]] <- playout[which(rnm == ti), c("x", "y")]
    attr(pgon, "srad") <- ceiling(sqrt(max(rowSums(playout^2))))
    class(pgon) <- "polygonLayout"
    return(pgon)
  } else if (dataType == "xy"){
    if (is.character(file_layout)){
      xylayout <- read.csv(file_layout, stringsAsFactors = FALSE)
    } else if (is.data.frame(file_layout)){
      xylayout <- file_layout
    }
    if (!unitCol %in% names(xylayout))
      stop("unitCol must be included in file_layout")
    if (!identical(make.names(xylayout[, unitCol]), xylayout[, unitCol])){
      stop("Turbine names in file_layout must be sytactically valid names in R.",
           "Must not start with a number or dot (.) followed by a number, and ",
           "must contain combinations of letters, numbers, dots( . ), and ",
           "underscores ( _ )."
      )
    }
    if (!is.null(scCol) && scCol %in% names(xylayout) && !is.null(notSearched)){
      xylayout <- xylayout[!xylayout[, scCol] %in% notSearched, ]
    }
    if (!is.null(file_turbine)){
      xyturbine <- read.csv(file_turbine, stringsAsFactors = FALSE)
    } else {
      xyturbine <- data.frame(
        turbine = gtools::mixedsort(unique(xylayout$turbine)),
        x = 0,
        y = 0
      )
      names(xyturbine)[grepl("turbine", names(xyturbine))] <- unitCol
    }
    if (!unitCol %in% names(xyturbine))
      stop("unitCol must be included in file_turbine")
    if (!all(xylayout[, unitCol] %in% xyturbine[, unitCol]))
      stop("Some units in file_layout are missing their counterparts in ",
           "file_turbine")

    for (ti in unique(xylayout[, unitCol])){
      lind <- which(xylayout[, unitCol] == ti)
      tind <- which(xyturbine[, unitCol] == ti)
      xylayout[lind, xCol] %<>% `-`(., xyturbine[tind, xCol])
      xylayout[lind, yCol] %<>% `-`(., xyturbine[tind, yCol])
    }
    xylayout$r <- sqrt(rowSums(xylayout[, c(xCol, yCol)]^2))
    nms <- xyturbine[, unitCol]
    xyturbine <- as.matrix(xyturbine[, c(xCol, yCol)])
    tbl <- table(xylayout[ , c("ncarc", "turbine")])
    zilch <- which(rownames(tbl) == 0)
    tbl <- tbl[-zilch, ]
    if (is.vector(tbl)) ncarc <- tbl else ncarc <- colSums(tbl)
    ncarc <- ncarc[gtools::mixedsort(names(tbl))]
    ncarc <- c(ncarc, total = sum(ncarc))
    rownames(xyturbine) <- nms
    output <- list(xydat = xylayout, tcenter = xyturbine, ncarc = ncarc,
      unitCol = unitCol, tset = unique(xylayout[, unitCol]))
    class(output) <- "xyLayout"
    return(output)
  }
}


#' different types of data but all end up in the same place: list of arrays by
#'  (by turbine) with r, x, y, scVar
#' @param rCol column for carcass distances from turbine
#' @export
readCarcass <- function(file_cod, unitCol = "turbine"){
  if (grepl(".shp", file_cod)){
#    if (is.null(file_turbine))
#      stop("readCarcass with shape files must include a file_turbine")
    carcasses <- sf::st_read(file_cod, stringsAsFactors = FALSE)
    carcasses0 <- sf::st_drop_geometry(carcasses)
    if (!unitCol %in% names(carcasses))
      stop("readCarcass: unitCol = ", unitCol, " not included in carcass data")
    if (!all(make.names(carcasses0[, unitCol, ]) == carcasses0[, unitCol])){
      badind <- which(make.names(carcasses0[, unitCol]) != carcasses0[, unitCol])
      badnm <- unique(carcasses0[badind, unitCol])
      badnm <- badnm[1:min(length(badnm), 3)]
      goodnm <- gsub("^X", "t", badnm)
      message(
        "\n\nNOTE: Not all the turbine names are syntactically valid.\n",
        "These (", paste(badnm, collapse = ", "),
        "...) can be converted to syntactically valid names (", goodnm, "...) ",
        "if desired."
      )
      tmp <- readline("convert [<enter>]? or not [n]? ")
      if (!identical(tolower(tmp), "y") & tmp != "") stop("aborting...")
      mknm <- make.names(carcasses0[, unitCol])
      badind <- which(mknm != carcasses0[, unitCol])
      badnm <- mknm[badind]
      carcasses0[badind, unitCol] <- gsub("^X", "t", badnm)
      carcasses[, unitCol] <- carcasses0[, unitCol]
    }
    tabn <- table(carcasses0[, unitCol])
    n <- as.vector(tabn)
    names(n) <- names(tabn)
    n <- n[gtools::mixedsort(names(n))]
    output <- list(carcasses = carcasses, unitCol = unitCol, ncarc = n)
    class(output) <- "shapeCarcass"
    return(output)
  }
}

#' @export
prepRing<- function(x, ...) UseMethod("prepRing", x)
#' @rdname aic
#' @export
prepRing.shapeLayout <- function(x, scVar = NULL, notSearched = NULL,
  silent = FALSE, ...){
  # remove polygons of unsearched area and unsearched turbines
  unitCol <- x$unitCol
  if (!is.null(scVar)){
    if (!scVar %in% names(x$layout))
      stop(scVar, " not in x$layout")
    if (!is.null(notSearched)){ # remove unsearched areas from data
      if (!is.character(notSearched))
        stop("notSearched must be a vector of names of scVar levels (or NULL)")
      x$layout <-  x$layout[!x$layout[, scVar, drop = TRUE] %in% notSearched, ]
      x$layoutAdj <- x$layoutAdj[!x$layoutAdj[, scVar, drop = TRUE] %in% notSearched, ]
      x$tset <- unique(x$layout[, x$unitCol, drop = TRUE])
      x$turbines <- x$turbines[x$turbines[, unitCol, drop = TRUE] %in% x$tset, ]
      x$tcenter <- x$tcenter[x$tset, ]
    }
  }
  shapeLayout <- x
  # create sf rings to intersect with searched areas
  srad <- ceiling(sqrt(max(rowSums(( # maximum searched distance at any turbine
    sf::st_coordinates(shapeLayout$layoutAdj)[, c("X", "Y")])^2))))
  theta <- seq(0, 2 * pi, length = 1000)
  trig <- cbind(cos(theta), sin(theta)); trig[dim(trig)[1], ] <- trig[1, ]
  radi <- 1:srad
  rings <- sf::st_sfc(sf::st_polygon(list(trig)))
  for (i in 2:length(radi))
    rings <- c(rings, sf::st_sfc(sf::st_polygon(list(radi[i] * trig))))
  rings <- sf::st_sf(rings)
  for (i in length(radi):2) rings[i, ] <- sf::st_difference(rings[i, ], rings[i-1, ])
  rings$r <- radi
  sf::st_crs(rings) <- sf::st_crs(shapeLayout$layoutAdj)

  ## calculate the area in each serach class for each ring:
  trsca <- list() # nturbines x nrings arrays of areas in each search class
  scset <- gtools::mixedsort(unique(as.data.frame(shapeLayout$layoutAdj)[, scVar]))
  if (length(scset) == 0) scset <- "all"
  if (!silent) {cat("calculating ring areas...\n"); flush.console()}
  for (sci in scset){
    if (!silent & sci != "all"){
      print(paste0(substitute(scVar), ' = ', sci), quote = FALSE)
      flush.console()
    }
    trsca[[sci]] <- matrix(0, #trcsa: turbne, r, search class, area
      nrow = length(shapeLayout$tset), ncol = dim(rings)[1])
    rownames(trsca[[sci]]) <- shapeLayout$tset
    for (ti in shapeLayout$tset){
      if (!silent){
        cat(" ", ti)
        flush.console()
      }
      # indices for turbine = ti and scVar = sci
      ind <- shapeLayout$layout[, unitCol, drop = TRUE] == ti
      if (sci != "all")
        ind <- ind & (shapeLayout$layout[, scVar, drop = TRUE] == sci)
      ind <- which(ind)
      # polygons for ind
      for (ii in ind){
        tring <- sf::st_geometry(shapeLayout$layoutAdj[ii,])
        jj <- unlist(sf::st_intersects(shapeLayout$layoutAdj[ii, ], rings))
        trsca[[sci]][ti, jj] <- trsca[[sci]][ti, jj] +
          as.vector(sf::st_area(sf::st_intersection(tring, rings)))
      }
    }
    cat("\n")
  }
  #Q: was turbine searched? include only those that were searched
  trsca <- lapply(trsca, FUN = function(x){
    x[rowSums(sapply(trsca, rowSums)) > 0, ]
  })
  tset <- rownames(trsca[[1]])

  ## format the results:
  rdat <- list()
  for (ti in tset){
    rdat[[ti]] <- data.frame(
      r = rep(1:srad, length(names(trsca))),
      scVar = rep(names(trsca), each = srad),
      exposure = 0,
      ncarc = 0,
      stringsAsFactors = F
    )
    if (!is.null(scVar))
      names(rdat[[ti]])[which(names(rdat[[ti]]) == "scVar")] <- scVar
    if (is.null(scVar)){
      rdat[[ti]][, "exposure"] <- trsca[[sci]][ti, 1:srad]
    } else {
      for (sci in scset){
        ind <- which(rdat[[ti]][, scVar] == sci)
        rdat[[ti]][ind, "exposure"] <- trsca[[sci]][ti, 1:srad]
      }
    }
  }

  ## search coverage probability by ring and turbine
  rpA <- lapply(rdat, FUN = function(x) aggregate(exposure ~ r, x, FUN = sum))
  for (ti in names(rpA)){
    rpA[[ti]][, "pinc"] <-
      round(rpA[[ti]][, "exposure"]/(2*pi*(rpA[[ti]][, "r"] - 0.5)), 3)
    rpA[[ti]][, "exposure"] <- NULL
  }
  rpA[["total"]] <- data.frame(
    r = 1:srad,
    pinc = rowMeans(matrix(unlist(lapply(rpA, "[", "pinc")), ncol = length(tset)))
  )
  # calculate rdat summary for "total" (aggregate over all turbines)
  nm <- c("r", scVar, "exposure", "ncarc")
  rdat[["total"]] <- data.frame(array(0,
    dim = c(nrow(rings) * length(scset), length(nm))))
  colnames(rdat[["total"]]) <- nm
  rdat[["total"]]$r <- rep(radi, length(scset))
  rdat[["total"]][, scVar] <- rep(scset, each = length(radi))
  rdat[["total"]][, "exposure"] <- unname(unlist(lapply(trsca, colSums)))

  # erase scVar column if no scVar was provided
  if (is.null(scVar))
    rdat <- lapply(rdat, FUN = function(x) {x[, "scVar"] <- NULL; x})
  # cull rows where pinc or exposure < 0.001
  rdat <- lapply(rdat, FUN = function(x) x[x$exposure >= 0.001, ])
  rpA <- lapply(rpA, FUN = function(x) x[x$pinc >= 0.001, ])
  class(rpA) <- "rpA"
  ncarc <- numeric(length(rdat))
  names(ncarc) <- names(rdat)
  output <- list(rdat = rdat, rpA = rpA, srad = srad, ncarc = ncarc,
    scVar = scVar, tcenter = shapeLayout$tcenter)
  class(output) <- "rings"
  return(output)
} # prepRing.shapeLayout

#' Prepare Ring Data for Data Frames of Simple Turbine Layouts
#' @param x a data frame with simple characterizations of turbines
#' @param cod carcass observation distances. Data frame with columns for
#'  distances that carcasses were found at and which turbines. If \code{cod = NULL},
#'  an empty template is created for carcass with descriptions of the turbines
#'  (i.e., \code{r}, \code{exposure}, \code{ncarc}, \code{pinc}), and (optional)
#'  \code{scVar}. A separate function can then be used to add carcasses distances later. This setup allows for easier
#'  simulation of carcasses counts using the predefined template for turbines.
#'  A "total" is added for the site...sums the carcasses and exposures
#' @param (optional) data frame or vector with carcass distances and turbine names
#' @export
prepRing.simpleLayout <- function(x){
  rdat <- list()
  rpA <- list()
  rownames(x) <- x$turbine
  srad <- attr(x, "srad")
  totA <- data.frame(r = 1:srad, pinc = 0)
  totr <- data.frame(r = 1:srad, exposure = 0, ncarc = 0)
  for (ti in x$turbine){
    rmax <- ceiling(x[ti, "radius"] * ifelse(x[ti, "shape"] == "square", sqrt(2), 1))
    rdat[[ti]] <- data.frame(r = 1:rmax, exposure = 0, ncarc = 0)
    rdat[[ti]]$exposure = switch(x[ti, "shape"],
      circular = 2*pi*(rdat[[ti]]$r - 1/2) * (rdat[[ti]]$r <= x[ti, "radius"]),
      square = diff(Acins(r = 0:rmax, s = x[ti, "radius"])),
      RP = {
        expo <- numeric(rmax) # 1 exposure for each radius = 1, ..., rmax
        rout <- ceiling(x[ti, "padrad"]:rmax)
        expo[rout] <- x[ti, "roadwidth"] * x[ti, "n_road"]
        rin <- which(rdat[[ti]]$r < x[ti, "padrad"])
        expo[rin] <- 2*pi*(rin - 1/2)
        expo
      }
    )
    rpA[[ti]] <- rdat[[ti]][, c("r", "exposure")]
    names(rpA[[ti]]) <- c("r", "pinc")
    rpA[[ti]]$pinc <- pmin(rpA[[ti]]$pinc/(2*pi*(1:rmax - 1/2)), 1)
    totA$pinc[1:rmax] <- totA$pinc[1:rmax] + rpA[[ti]]$pinc
    totr$exposure[1:rmax] <- totr$exposure[1:rmax] + rdat[[ti]]$exposure
  }
  rpA <- lapply(rpA, FUN = function(xi){
    xi$pinc <- round(xi$pinc, 4)
    xi <- xi[xi$pinc >= 0.001, ]
    xi
  })
  rdat <- lapply(rdat, FUN = function(xi){
    xi$exposure <- signif(xi$exposure, 4)
    xi <- xi[xi$exposure >= 0.001, ]
    xi
  })
  totA$pinc <- pmin(totA$pinc/nrow(x), 1)
  rpA[["total"]] <- totA
  class(rpA) <- "rpA"
  rdat[["total"]] <- totr
  # cull the rows with radii beyond the edge
  ncarc <- numeric(nrow(x) + 1)
  names(ncarc) <- c(x[, unitCol], "total")
  output <- list(rdat = rdat, rpA = rpA, srad = srad, ncarc = ncarc, scVar = NULL)
  class(output) <- "rings"
  return(output)
}

prepRing.numeric <- function(x, srad, ...){
  if (!is.vector(dim(x))) stop("x must be a vector in prepRing.numeric")
  if (!is.numeric(srad) || length(srad) != 1 || srad <= 1)
    stop("srad must be numeric scalar > 1")
  srad <- ceiling(srad)
  rpA <- list()
  rdat <- list()
  rvec <- x[x<=srad]
  tbl <- table(ceiling(rvec))
  r <- 1:srad
  rdat[["total"]] <- data.frame(r = r, exposure = 2*pi*(r - 1/2), ncarc = 0)
  rdat[["total"]]$ncarc[as.numeric(names(tbl))] <- unname(tbl)
  rpA[["total"]] <- data.frame(r = r, pinc = 1)
  class(rpA) <- "rpA"
  output <- list(rdat = rdat, rpA = rpA, srad = srad, ncarc = length(rvec),
    scVar = NULL)
  class(output) <- "rings"
  return(output)
}

prepRing.polygonLayout <- function(x, ...){
  srad <- attr(x, "srad")
  rdat <- list()
  rpA <- list()
  theta <- seq(0, 2*pi, length = 1000)
  rr <- 1:srad
  cx <- outer(rr, cos(theta))
  cy <- outer(rr, sin(theta))
  iarea <- list()
  rtot <- data.frame(r = 1:srad, exposure = 0, ncarc = 0)
  Atot <- data.frame(r = 1:srad, pinc = 0)
  for (ti in names(x)){
    rmax <- ceiling(max(sqrt(rowSums(x[[ti]]^2))))
    iarea[[ti]] <- numeric(rmax)
    for (ri in 1:rmax){
      iarea[[ti]][ri] <- gpclib::area.poly(gpclib::intersect(
        as(cbind(cx[ri, ], cy[ri, ]), "gpc.poly"),
        as(x[[ti]], "gpc.poly")
      ))
    }
    iarea[[ti]] <- c(0, iarea[[ti]])
    rdat[[ti]] <- data.frame(r = 1:rmax, exposure = diff(iarea[[ti]]), ncarc = 0)
    rpA[[ti]] <- data.frame(
      r = 1:rmax,
      pinc = rdat[[ti]][, "exposure"]/((1:rmax - 1/2) * 2 * pi)
    )
    rtot[1:rmax, "exposure"] %<>% `+` (., rdat[[ti]][, "exposure"])
    Atot[1:rmax, "pinc"] %<>% `+` (., rpA[[ti]][, "pinc"])
  }
  rdat[["total"]] <- rtot
  rpA[["total"]] <- Atot
  rpA[["total"]][, "pinc"] %<>% `/` (., length(x))
  rpA <- lapply(rpA, round, 4)
  class(rpA) <- "rpA"
  ncarc <- numeric(length(rdat))
  names(ncarc) <- names(rdat)
  output <- list(rdat = lapply(rdat, round, 4), rpA = rpA,
    ncarc = ncarc, srad = srad, scVar = NULL)
  class(output) <- "rings"
  return(output)
}


#' Fit Distance Distribution Model(s)
#'
#' @description Fit glm's for distance distribution models corresponding to
#'  xep01 (gamma), xep012, lognormal, xep1, xep12, xep123, truncated normal,
#'  exponential, xep2 (Rayleigh), xep0 (Pareto), xepi0 (inverse gamma), xep0123
#'  (normal gamma with x = tau), Maxwell Boltzmann, chi-squared, and/or inverse
#'  Gaussian distributions. By default, \code{ddFit} fits all the aforementioned
#'  models, but any subset of the models may be fit as well by using, for example,
#'  \code{model = exclude(c("xep123", "constant"))} to fit all except the
#'  \code{"xep123"} and \code{"constant"} models or \code{model = c("xep01", "lognormal")}
#'  to fit just the \code{"xep01"} and \code{"lognormal"} models.
#'
#'  The glm is converted to a probability distribution by dividing by a
#'  normalizing constant, namely the integral of the glm evaluated from 0 to
#'  infinity. In some cases (most notably when the leading coefficient of the
#'  glm is positive so the fitted curve does not converge to zero as x increases),
#   the integral does not converge to a finite value and the glm cannot be
#'  converted to a probability distribution. In these cases, the distribution
#'  parameters are given as \code{NA}.
#'
#' @param rdat Data frame with a tally of carcasses within 1 m rings, along
#'  with outer radii of the rings, (optional) search class column, and the amount
#'  of area searched in each search class and ring.
#' @param model names (vector of character strings) of glm distribution templates
#'  to fit. Default is \code{"all"}, which fits \code{xep01}, \code{lognormal},
#'  \code{xep1}, \code{xep12}, \code{xep123}, \code{xepi0},
#'  \code{xep0123}, \code{xep012}, \code{xep2}, \code{MaxwellBoltzmann}, \code{constant},
#'  \code{tnormal}, \code{exponential}, \code{xep0}, \code{chisq}, and
#'  \code{inverse_gaussian}. Any subset of these may also be fit with a single
#'  call to \code{ddFit}.
#' @param scVar name of column in \code{rdat} with (optional) subdivisions
#'  of rings into search classes (character strings)
#' @param rCol name of column in \code{rdat} with outer radius for 1m rings
#' @param expoCol name of "exposure" column of total area in each ring belonging
#'  to each search class.
#' @param ... ignored
#'
#' @return A list of fitted glm models as \code{dd} objects in a \code{ddArray}
#'  object if a vector of distributions is fit or a single \code{dd} object if a
#'  single model is fit. The \code{dd} objects are lists with all the standard
#'  \code{\link[stats]{glm}} components, in addition to the following elements:
#'  \describe{
#'    \item{\code{$distr}}{name of the distribution ("xep01", etc.)}
#'    \item{\code{$parms}}{vector of distribution parameter estimates (or \code{NA}
#'      if the model based on the MLE is not extensible)}
#'    \item{\code{$varbeta}}{the variance-covariance matrix of the glm parameter
#'      estimates. NOTE: This is identical to the covariance matrix from the glm,
#'      which can be extracted via \code{summary(x)$cov.unscaled}}
#'    \item{\code{$scVar}}{name of the (optional) search class variable (or \code{NULL})}
#'    \item{\code{$ncarc}}{number of carcasses}
#'    \item{\code{$n}}{number of rings}
#'    \item{\code{$k}}{number of parameters}
#'    \item{\code{$srad}}{search radius}
#'  }
#'  When a \code{dd} object is printed, only a small subset of the elements are
#'  shown. To see a full list of the objects, use \code{names(x)}. The elements
#'  can be extracted in the usual R way via \code{$} or \code{[[x]]}.
#' @export

ddFit <- function(x, ...) UseMethod("ddFit", x)

#' @export
ddFit.data.frame <- function(x, distr = "standard", scVar = NULL, rCol = "r",
    expoCol = "exposure", ncarcCol = "ncarc", silent = FALSE, ...){
  rdat <- x
  dat <- rdat[, c(rCol, scVar, expoCol, ncarcCol)]
  names(dat) <- c("r", if(!is.null(scVar)) scVar, "exposure", "ncarc")
  dat$r <- dat$r - 0.5
  scl <- ifelse(is.null(scVar), "", paste0(" + ", scVar))
  plu <- ifelse(!is.null(scVar), " + ", " ")
  sufx <- paste0(plu, scVar,  "+ offset(log(", expoCol, "))")
  sufx0 <- paste0(plu, scVar, "+ ")
  pref <- paste0("ncarc ~ ")
  output <- list()
  if (identical(distr, "all")) distr <- mod_all # see parameters.r
  if (identical(distr, "standard")) distr <- mod_standard # see parameters.r
  for (di in distr){
    form <- formula(paste0("ncarc ~ ", paste(cof_name[[di]][-1], collapse = " + "),
      scl, " + ", paste0("offset(", mod_offset[[di]], ")")))
    output[[di]] <- suppressWarnings(glm(
      formula = form, data = dat, family = "poisson"
    ))
    output[[di]]$distr <- di
    output[[di]]$form <- form
    output[[di]]$scVar <- scVar
    output[[di]]$parms <- cof2parms(output[[di]]$coefficients, distr = di)
    output[[di]]$varbeta <- summary(output[[di]])$cov.unscaled
    output[[di]]$ncarc <- sum(dat$ncarc)
    output[[di]]$n <- nrow(dat)
    output[[di]]$k <- output[[di]]$n - output[[di]]$df.residual
    class(output[[di]]) <- c("dd", class(output[[di]]))
  }
  if (!silent){
    anybad <- FALSE
    for (di in distr){
      if (!silent && anyNA(output[[di]]$parms)){
        if (!anybad){
          anybad <- TRUE
          cat("Non-extensible models:\n")
        }
        cat(" ", di, "\n")
      }
    }
    if (anybad) flush.console()
  }
  if(length(output) == 1) return(output[[1]])
  class(output) <- "ddArray"
  return(output)
}

#' @export
ddFit.rings <- function(x, distr = "standard", scVar = NULL, rCol = "r",
    expoCol = "exposure", ncarcCol = "ncarc", silent = FALSE, ...){
  arglist <- as.list(match.call())
  arglist[[1]] <- NULL
  arglist$x <- x$rdat$total
  do.call(ddFit, arglist)
}

#' @export
ddFit.list <- function(x, distr = "standard", scVar = NULL, rCol = "r",
    expoCol = "exposure", ncarcCol = "ncarc", silent = FALSE, ...){
  arglist <- as.list(match.call())
  arglist[[1]] <- NULL
  arglist$x <- x[[1]]
  do.call(ddFit, arglist)
}

#' @export
ddFit.xyLayout <- function(x, distr = "standard", scVar = NULL, notSearched = NULL,
    rCol = "r", ncarcCol = "ncarc", unitCol = "turbine", silent = FALSE, ...){
  xydat <- x$xydat
  if (!is.null(scVar)){
    dat <- xydat[xydat[, scVar] != notSearched, c(rCol, scVar, ncarcCol)]
  } else {
    dat <- xydat[, c(rCol, ncarcCol)]
  }
  names(dat) <- c("r", if(!is.null(scVar)) scVar, "ncarc")
  scl <- ifelse(is.null(scVar), "", paste0(" + ", scVar, collapse = ""))
  pref <- paste0("ncarc ~ ")
  output <- list()
  if (identical(distr, "all") | identical(distr, "standard")) distr <- mod_xy
  if (!all(distr %in% mod_xy)){
    stop(paste0(distr[!distr %in% mod_xy], " not in list of models ",
        "accommodated for xy data."))
  }
  for (di in distr){
    form <- formula(paste0("ncarc ~ ", paste(cof_name[[di]][-1], collapse = " + "),
      scl))
    output[[di]] <- suppressWarnings(glm(
      formula = form, data = dat, family = "poisson"))
    output[[di]]$distr <- di
    output[[di]]$form <- form
    output[[di]]$scVar <- scVar
    output[[di]]$parms <- cof2parms(output[[di]]$coefficients, distr = di)
    output[[di]]$varbeta <- summary(output[[di]])$cov.unscaled
    output[[di]]$ncarc <- sum(dat$ncarc)
    output[[di]]$n <- nrow(dat)
    output[[di]]$k <- output[[di]]$n - output[[di]]$df.residual
    class(output[[di]]) <- c("dd", class(output[[di]]))
  }
  if (!silent){
    anybad <- FALSE
    for (di in distr){
      if (!silent && anyNA(output[[di]]$parms)){
        if (!anybad){
          anybad <- TRUE
          cat("Non-extensible models:\n")
        }
        cat(" ", di, "\n")
      }
    }
    if (anybad) flush.console()
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
    AIC = round(aic0, 2),
    AICc = round(aic0 + 2 * x$k * (x$k + 1)/(x$n - x$k - 1), 2),
    extensible = 1 * cofOK(x$coefficients, x$distr)
  ))
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
plot.ddArray = function(x, type = "CDF", extent = "full", distr = "all",
    xmax = NULL, resolution = 250, mod_highlight = NULL, par_reset = TRUE, ...){
  if (identical(distr, "all")) distr <- names(x)
  if (extent == "full")
    distr <- distr[sapply(x[distr], function(tmp) !any(is.na(tmp$parm)))]
  if (type == "rcd")
    distr <- distr[which(!distr %in% c("exponential", "tnormal", "constant"))]
  if (!any(distr %in% mod_all)) stop("plot.ddArray: some model(s) undefined")
  dd <- x[distr]
  aic_proper <- aic(dd, extent = extent)
  if (!type %in% c("CDF", "PDF", "rcd"))
    stop("type (", deparse(substitute(type)),
          ") must be \"PDF\", \"CDF\", or \"rcd\"")
  distr <- aic_proper$model
  mod_best <- distr[1]
  if (is.null(mod_highlight)) mod_highlight <- mod_best
  ### graph is in two parts, in succession with par(fig...) and par(new = T)
  ###  0) preliminaries
  ###  1) legend [common to all parts]
  ###  2) main graph [lines]
  # extracting and parsing parameters from ... args
  arglist <- list(...)
  if (!"col" %in% names(arglist)){
    arglist$col <- mod_color
  } else {
    if (is.null(names(arglist$col)) && length(arglist$col) != length(distr) &&
      length(arglist$col) != 1){
        stop(
          "'col' in plot.ddArray must be a scalar, a vector with one color for ",
          "each model plotted, or a vector with named elements matching the ",
          "names of the models to be plotted."
        )
    } else if (length(arglist$col) == 1){
      arglist$col <- rep(arglist$col, length.out = length(distr))
      names(arglist$col) <- distr
    } else if (!is.null(names(arglist$col))){
      if (!all(distr %in% names(arglist$col))){
        stop(
          "in plot.ddArray, 'col' be NULL or include a color for each plotted model"
        )
      }
    } else if (length(arglist$col) == length(distr)){
      names(arglist$col) <- distr
    }
  }
  if (!"lty" %in% names(arglist)){
    arglist$lty <- 1 + !natural[aic_proper$model]
    names(arglist$lty) <- aic_proper$model
  } else {
    if (!is.vector(arglist$lty))
      stop("'lty' in plot.ddArray must be a scalar or vector")
    if (is.null(names(arglist$lty)) && length(arglist$lty) >= length(distr)){
        stop(
          "'lty' in plot.ddArray must be a scalar or vector shorter than ",
          "number of plotted models. Values are recycled to fill out the vector ",
          "of line types to match the number of models to be plotted. Enter ",
          "?plot.ddArray for more info."
        )
    } else {
      arglist$lty <- rep(arglist$lty, length.out = length(distr))
      names(arglist$lty) <- distr
    }
  }
  if (!"xlim" %in% names(arglist)){
    if (is.null(xmax) || !is.finite(xmax)) xmax <- max(dd[[1]]$data[, "r"])
    xmin <- 0
    arglist$xlim = c(0, xmax)
  }
  arglist$x = 0
  arglist$type = "n"
  ## part 1: legend
  do.call(par, par_default)
  ncol <- ifelse(nrow(aic_proper) >= 6, 2, 1)
  sz <- ifelse(round(nrow(aic_proper)/ncol) <= 3, 0.13,
          ifelse(round(nrow(aic_proper)/ncol) <= 5, 0.17,
            ifelse(round(nrow(aic_proper)/ncol) <= 7, 0.22, 0.26)))
  par(fig = c(0, 1, 0, sz), mar = c(0, 1, 1.5, 0), family = "mono")
  plot(0, type = "n", axes = F, xlab = "", ylab = "")
  leglab <- character(nrow(aic_proper))
  lwd <- numeric(nrow(aic_proper)) + 1
  names(lwd) <- aic_proper$model
  if(mod_highlight %in% aic_proper$model) lwd[mod_highlight] <- 2
  for (i in 1:nrow(aic_proper))
   leglab[i] <- sprintf("%18-s%6.2f", aic_proper$model[i], aic_proper$deltaAIC[i])
  mtext(side = 3, line = 0, adj = 0, sprintf("%19-s%s", "Distribution", "\u0394AICc"))
  legend(x = "topleft", legend = leglab, cex = 0.8, lwd = lwd, bty = "n",
    ncol = ncol, lty = arglist$lty, col = arglist$col[aic_proper$model])
  ## part 2a: general plot layot
  xseq <- seq(min(arglist$xlim), max(arglist$xlim), length = resolution)
  # size of the main graph; vertical, relative to plot window size:
  par(fig = c(0, 1, sz, 1), mar = c(4, 4, 2, 0.5), family = "sans", new = TRUE)
  ## part 2b: lines
  if (type == "CDF"){
    if (!"xlab" %in% names(arglist))
      arglist$xlab = "Distance from Turbine"
    if (!"ylab" %in% names(arglist))
      arglist$ylab = "P(carcass falls within x meters from turbine)"
    if (!"ylim" %in% names(arglist)) arglist$ylim <- 0:1
    do.call(plot, arglist)
#      plot(0, type = "n", xlim = xlim, ylim = 0:1, xlab = xlab, ylab = ylab, ...)
    for (fi in aic_proper$model[nrow(aic_proper):1])
      lines(xseq, pdd(xseq, model = dd[fi], extent = extent),
        lty = arglist$lty[fi], col = arglist$col[fi])
    if (mod_highlight %in% distr){
      lines(xseq, pdd(xseq, model = dd[mod_highlight], extent = extent),
        col = arglist$col[mod_highlight], lwd = 2, lty = arglist$lty[mod_highlight])
    }
    if (extent == "win")
      mtext(side = 1, line = -1, adj = 1, "within search radius", family = "serif")
  } else if (type == "PDF"){
    if (!"ylim" %in% names(arglist)){
          ymax = -Inf
      for (fi in aic_proper$model){
        ymax <- max(ymax, max(ddd(xseq, dd[fi], extent = extent)))
      }
      arglist$ylim <- c(0, ymax)
    }
    arglist$xlab = ifelse("xlab" %in% names(arglist), arglist$xlab, "r")
    arglist$ylab = ifelse("ylab" %in% names(arglist), arglist$ylab, "PDF")
    do.call(plot, arglist)
    for (fi in aic_proper$model[nrow(aic_proper):1]){
      lines(xseq, ddd(xseq, dd[fi], extent = extent),
        col = arglist$col[fi], lty = arglist$lty[fi])
    }
    if (mod_highlight %in% distr){
      lines(xseq, ddd(xseq, dd[mod_highlight], extent = extent),
        col = arglist$col[mod_highlight], lwd = 2, lty = arglist$lty[mod_highlight])
    }
    leglab <- character(nrow(aic_proper))
    if (extent == "win")
      mtext(side = 3, adj = 1, "within search radius", family = "serif", line = -1)
   } else if (type == "rcd"){
    if (!"ylim" %in% names(arglist)){
          ymax = -Inf
      for (fi in aic_proper$model){
        ymax <- max(ymax, max(rcd(xseq, dd[fi], extent = extent)))
      }
      arglist$ylim <- c(0, ymax)
    }
    arglist$xlab = ifelse("xlab" %in% names(arglist), arglist$xlab, "r")
    arglist$ylab = ifelse("ylab" %in% names(arglist),
      arglist$ylab, "Relative Carcass Density")
    do.call(plot, arglist)
    for (fi in aic_proper$model[nrow(aic_proper):1]){
      lines(xseq, rcd(xseq, model = dd[fi], extent = extent),
        col = arglist$col[fi], lty = arglist$lty[fi])
    }
    if (mod_highlight %in% distr){
    lines(xseq, rcd(xseq, model = dd[mod_highlight], extent = extent),
      col = arglist$col[mod_highlight], lwd = 2, lty = arglist$lty[mod_highlight])
    }
    if (extent == "win") mtext(side = 3, adj = 1, "within search radius  ",
      family = "serif", line = -1)
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
  polygon(CI[c(1:nrow(CI), nrow(CI):1) , "x"] , c(CI[, 2], CI[nrow(CI):1, 3]),
    col = colors()[350], border = NA)
  msg <- switch(extent, win = "limited to carcasses within search radius  ")
  if (type == "PDF"){
    lines(xseq, ddd(xseq, x, extent = extent), col = mod_color[x$distr], lwd = 2)
    mtext(side = 3, line = -1, adj = 1, msg, family = "serif", lty = 1 + !natural[distr])
  } else {
    lines(xseq, pdd(x = xseq, model = x, extent = extent),
      col = mod_color[x$distr], lwd = 2)
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
#'  lists in R, namely, by using the \code{x$element} or \code{x[[element]]}
#'  operator, where \code{element} is the name of one of the elements of
#'  \code{x} which can be viewed via \code{names(x)}.
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
    xep01 = {
      parms <- cbind(
        shape = x[, "log(r)"] + 2,
        rate = -x[, "r"])
      parms[!cofOK(x, "xep01"), ]  <- NA
      parms
    },
    lognormal = {
      parms <- cbind(
        meanlog = -0.5 * (2 + x[, "log(r)"])/x[, "I(log(r)^2)"],
        sdlog = sqrt(-0.5/x[, "I(log(r)^2)"]))
      parms
    },
    xep1 = {
      parms <- x[, "r", drop = FALSE]
      parms[parms >= 0] <- NA
      colnames(parms) <- "b1"
      parms
    },
    xep12 = {
      parms <- x[, c("r", "I(r^2)"), drop = FALSE]
      colnames(parms) <- c("b1", "b2")
      parms[parms[, "b2"] >= 0, ] <- NA
      parms
    },
    xep02 = {
      parms <- x[, c("log(r)", "I(r^2)"), drop = FALSE]
      colnames(parms) <- c("b0", "b2")
      parms[parms[, "b2"] >= 0, ] <- NA

      parms
    },
    xep123 ={
      parms <- x[, c("r", "I(r^2)", "I(r^3)"), drop = FALSE]
      colnames(parms) <- c("b1", "b2", "b3")
      parms[parms[, "b3"] >= 0, ] <- NA
      parms
    },
    xepi0 = {
      parms <- cbind(
        shape = -x[, "log(r)"] - 2,
        scale = -x[, "I(1/r)"])
      parms[parms[, "shape"] <= 0 | parms[, "scale"] <= 0, ] <- NA
      parms
    },
    xep0123 = {
      parms <- x[, c("log(r)", "r", "I(r^2)", "I(r^3)"), drop = FALSE]
      colnames(parms) <- c("b0", "b1", "b2", "b3")
      parms[parms[, "b3"] >= 0, ] <- NA
      parms
    },
    xep012 = {
      parms <- x[, c("log(r)", "r", "I(r^2)"), drop = FALSE]
      colnames(parms) <- c("b0", "b1", "b2")
      parms[parms[, "b2"] >= 0, ] <- NA
      parms
    },
    xep2 = {
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
    xep0 = {
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
#'  (e.g., "xep01").
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
        xep01 = matrix(
          dgamma(xx, shape = parms[ppi, "shape"], rate = parms[ppi, "rate"]),
          nrow = length(x)
        ),
        lognormal = matrix(
          dlnorm(xx, meanlog = parms[ppi, "meanlog"], sdlog = parms[ppi, "sdlog"]),
          nrow = length(x)
        ),
        xep1 = matrix(
          dxep1(xx, b1 = parms[ppi, "b1"]),
          nrow = length(x)
        ),
        xep12 = matrix(
          dxep12(xx, b1 = parms[ppi, "b1"], b2 = parms[ppi, "b2"]),
          nrow = length(x)
        ),
        xep123 = {
          const <- numeric(length(i1))
          tmp <- matrix(0, nrow = length(x), ncol = length(i1))
          for (i in i1){
            const <- tryCatch(1/integrate(
              f = function(r)
                x * exp(b1 * x + b2 * x^2 + b3 * x^3),
              lower = 0, upper = Inf, rel.tol = .Machine$double.eps^0.5)$val,
              error = function(e) NA
            )
            tmp[, which(i1 == i)] <- dxep123(x,
              b1 = c(parms[i, "b1"]), b2 = c(parms[i, "b2"]), b3 = c(parms[i, "b3"]),
              const = const
            )
          }
          tmp
        },
        xepi0 = matrix(
          dxepi0(xx, shape = parms[ppi, "shape"], scale = parms[ppi, "scale"]),
          nrow = length(x)
        ),
        xep0123 = { # assumes b0, b1, b2, b3 are scalars
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
            tmp[, which(i1 == i)] <- dxep0123(x,
              b0 = c(parms[i, "b0"]),
              b1 = c(parms[i, "b1"]),
              b2 = c(parms[i, "b2"]),
              b3 = c(parms[i, "b3"]),
              const = const
            )
          }
          tmp
        },
        xep012 = {
          const <- numeric(length(i1))
          tmp <- matrix(0, nrow = length(x), ncol = length(i1))
          for (i in i1){
            b0 <- c(parms[i, "b0"])
            b1 <- c(parms[i, "b1"])
            b2 <- c(parms[i, "b2"])
            const <- tryCatch(1/integrate(f = function(x)
              x * exp(b0 * log(x) + b1 * x + b2 * x^2),
              lower = 0, upper = Inf, rel.tol = .Machine$double.eps^0.5)$val,
              error = function(e) NA
            )
            tmp[, which(i1 == i)] <- dxep012(x,
              b0 = c(parms[i, "b0"]),
              b1 = c(parms[i, "b1"]),
              b2 = c(parms[i, "b2"]),
              const = const
            )
          }
          tmp
        },
        xep02 = matrix(
          dxep02(xx, b0 = parms[ppi, "b0"], b2 = parms[ppi, "b2"]),
          nrow = length(x)
        ),
        xep2 = matrix(
          dxep2(xx, s2 = parms[ppi, "s2"]),
          nrow = length(x)
        ),
        MaxwellBoltzmann = matrix(dmb(xx, a = parms[ppi, "a"]), nrow = length(x)),
        constant = {
          tmp <- 2 * xx/zrad^2
          tmp[xx > zrad] <- 0
          matrix(tmp, nrow = length(x))
        },
        xep0 = matrix(dxep0(xx, a = parms), nrow = length(x)),
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
        if(distr == "xep0123"){
          if(parms[i, "log(r)"] <= -2){
            output[, i] <- NA
            next
          }
        }
        if(distr == "xep01" && parms[i, "log(r)"] <= -2){
          output[, i] <- NA
          next
        }
        if(distr == "xepi0" && parms[i, "log(r)"] > -2){
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
#'  (e.g., "xep01").
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
    srad <- attr(model, which = "srad")
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
        xep01 = matrix(
          pgamma(xx, shape = parms[ppi, "shape"], rate = parms[ppi, "rate"]),
          nrow = length(x)
        ),
        lognormal = matrix(
          plnorm(xx, meanlog = parms[ppi, "meanlog"], sdlog = parms[ppi, "sdlog"]),
          nrow = length(x)
        ),
        xep1 = matrix(pxep1(xx, b1 = parms[ppi, "b1"]),
          nrow = length(x)
        ),
        xep12 = matrix(pxep12(xx, b1 = parms[ppi, "b1"], b2 = parms[ppi, "b2"]),
          nrow = length(x)
        ),
        xep02 = matrix(pxep02(xx, b0 = parms[ppi, "b0"], b2 = parms[ppi, "b2"]),
          nrow = length(x)
        ),
        xep123 = {
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
            tmp[, which(i1 == i)] <- pxep123(x,
              b1 = c(parms[i, "b1"]), b2 = c(parms[i, "b2"]), b3 = c(parms[i, "b3"]),
              const = const
            )
          }
          tmp
        },
        xepi0 = matrix(
          pxepi0(xx, shape = parms[ppi, "shape"], scale = parms[ppi, "scale"]),
          nrow = length(x)
        ),
        xep012 = { # assumes b0, b1, b2, b3 are scalars
          const <- numeric(length(i1))
          tmp <- matrix(0, nrow = length(x), ncol = length(i1))
          for (i in i1){
            b0 <- c(parms[i, "b0"])
            b1 <- c(parms[i, "b1"])
            b2 <- c(parms[i, "b2"])
            const <- tryCatch(1/integrate(f = function(x)
              x * exp(b0 * log(x) + b1 * x + b2 * x^2),
              lower = 0, upper = Inf, rel.tol = .Machine$double.eps^0.5)$val,
              error = function(e) NA
            )
            tmp[, which(i1 == i)] <- pxep012(x,
              b0 = b0,
              b1 = b1 ,
              b2 = b2,
              const = const
            )
          }
          tmp
        },
        xep0123 = { # assumes b0, b1, b2, b3 are scalars
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
            tmp[, which(i1 == i)] <- pxep0123(x,
              b0 = b0,
              b1 = b1 ,
              b2 = b2,
              b3 = b3,
              const = const
            )
          }
          tmp
        },
        xep2 = matrix(pxep2(xx, s2 = parms[ppi, "s2"]), nrow = length(x)),
        MaxwellBoltzmann = matrix(pmb(xx, a = parms[ppi, "a"]), nrow = length(x)),
        constant = {
          tmp <- (xx/zrad)^2
          tmp[xx > zrad] <- 1
          matrix(tmp, nrow = length(x))
        },
        xep0 = matrix(pxep0(xx, a = parms), nrow = length(x)),
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
        if(distr %in% c("xep01", "xep02", "xep0123", "xep012")){
          if(parms[i, "log(r)"] <= -2){
            output[, i] <- NA
            next
          }
        }
        if(distr == "xepi0" && parms[i, "log(r)"] > -2){
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
      if (grepl("01", distr) && parms[i, "log(r)"] <= -2){
        output[, i] <- NA
        next
      }
      deno <- tryCatch(integrate(f = function(r)
          exp(rmat(r, distr) %*% parms[i, cof_name[[distr]]] + off(r, distr)),
          lower = ifelse(distr == "xep0", 1, 0),
          upper = srad
        )$val,
        error = function(e) NA
      )
      if (!is.na(deno)){
        for (xi in which(x > ifelse(distr == "xep0", 1, 0) & x <= srad)){
          output[xi, i] <- integrate(function(r)
            exp(rmat(r, distr) %*% parms[i, cof_name[[distr]]] + off(r, distr)),
            lower = ifelse(distr == "xep0", 1, 0),
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
    srad <- attr(model, which = "srad")
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
      lower = ifelse(distr == "xep0", 1, 0), upper = upr)$val,
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
  attr(parms, "srad") <- ceiling(max(dd$data$r, na.rm = TRUE))
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
    if ((distr %in% c("xep01", "xep012", "xep0123") && model[, "log(r)"] <= -2) |
      (distr == "xepi0" && model[, "log(r)"] >= -2)){
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

#' Estimate Probability Carcass lands in Searched Area
#' @description Estimated probability that carcass lands in searched area. This
#'  is an intermediate step in estimating dwp but is also interesting in its own
#'  right.
#' @param x data. Either \code{rdat} data frame for a single turbine, site,
#'  or other entity, or a list of such entities.
#' @param model A fitted \code{dd} model
#' @param extent calculate dwp within searched area only (\code{"win"}) or
#'  for fulss complement of carcasses (\code{"full"}).
#' @param nsim number of parametric bootstrap iterations for estimating CIs
#' @param ... ignored
#' @return \code{psiHat} object = list with
#' @export
estpsi <- function(x, ...) UseMethod("estpsi", x)

#' @rdname estpsi
#' @export
estpsi.data.frame <- function(x, model, extent = "full", nsim = 1000, zrad = 150, ...){
  if ("dd" %in% class(model)) model <- ddSim(model, nsim = nsim)
  if (!"ddSim" %in% class(model)) stop("estpsi: model must be dd or ddSim object")
  output <- list(
    psi = as.vector(t(ddd(x$r - 1/2, model, extent = extent, zrad = zrad)) %*% x$pinc),
    ncarc = sum(x$ncarc)
  )
  output$psi[output$psi > 1] <- 1
  attr(output, "extent") <- extent
  attr(output, "zrad") <- zrad
  class(output) <- "psiHat"
  output
}

#' @rdname estpsi
#' @export
estpsi.rpA <- function(x, model, extent = "full", nsim = 1000, zrad = 150, ...){
  if ("dd" %in% class(model)) parmsim <- ddSim(model, nsim = nsim)
  tmp <- lapply(x, FUN = estpsi, model = parmsim,
    extent = extent, nsim = nsim, zrad = zrad)
  output <- sapply(tmp, "[[", "psi")
  attr(output, "extent") <- extent
  attr(output, "zrad") <- zrad
  class(output) <- "psiHat"
  output
}

#' @rdname estpsi
#' @export
estpsi.xyLayout <- function(x, model, extent = "full", nsim = 1000, zrad = 150,
    ...){
  if ("dd" %in% class(model)) parmsim <- ddSim(model, nsim = nsim)
  di <- model$distr
  ep <- function(r, di, pco) c(exp(rmat(r, di) %*% pco[1, cof_name[[di]]]))
  tset <- x$tset
  output <- matrix(NA, nrow = nsim, ncol = length(tset) + 1)
  colnames(output) <- c(tset, "total")
  for (simi in 1:nsim){
    pco <- parmsim[simi, ]
    if (!cofOK0(pco, di)) next
    deno <- integrate(f = function(r, di, pco) 2 * pi * r * ep(r, di, pco),
      lower = 0, upper = ifelse(cofOKInf(pco, di), Inf, zrad), di = di, pco = pco)$val
    tmp <- aggregate(ep(x$xydat$r, di, pco) ~ x$xydat$turbine, FUN = "sum")
    output[simi, tmp[, 1]] <- tmp[, 2]/deno
  }
  output[output > 1] <- 1
  output[, "total"] <- rowMeans(output[, exclude("total", colnames(output))], na.rm = TRUE)
  attr(output, "extent") <- extent
  attr(output, "zrad") <- zrad
  class(output) <- "psiHat"
  return(output)
}

#' Estimate DWP
#' @param x \code{psiHat} object, which is a matrix \code{$psi} giving the
#'  probability of a given carcass landing in the searched area at each turbine.
#'  Column names = turbine IDs.
#' @param ncarc vector of total carcass count at each turbine represented in x.
#' @param nboot number of parametric bootstrap iterations for estimating CIs
#' @param GenEst format the results for importing into GenEst
#' @param ... ignored
#' @return list
#' @export
estdwp <- function(x, ...) UseMethod("estdwp", x)

#' @rdname estpsi
#' @export
estdwp.psiHat <- function(x, ncarc, nboot = NULL, forGenEst = FALSE,
    silent = TRUE, ...){
  # the bootstrap for estimated dwp is independent of the bootstrap for psi and
  # will potentially have different sample sizes. The default is to have them
  # both be the same (e.g., nboot = nsim = 1000). Another possibility would be
  # to have a fixed model for psi (nsim = 1) but still want to incorporate the
  # uncertainty in M|{x, psi}. A third possibility would be to have a different
  # number of iterations for estimating psi than for estimating dwp.
  psi <- x
  if (length(psi) == length(ncarc)){
    if (is.null(names(psi))){
      if (is.null(names(ncarc))) names(ncarc) <- paste0("t", 1:length(ncarc))
      names(psi) <- names(ncarc)
    }
    if (is.null(names(ncarc))) names(ncarc) <- names(psi)
    # CASE 1: single, "known" psi for each turbine
    # Need to incorporate uncertainty in M|[x, psi}
    # in this case the posterior distribution of M is known, so the implied dwp
    # distribution is also known, but GenEst requires a simulated column of
    # values of length nboot, so calculate nboot quantiles
    nboot <- ifelse(is.null(nboot), 1000, nboot)
    output <- list(psi = matrix(0, nrow = nboot, ncol = length(ncarc)))
    colnames(output$psi) <- names(psi)
    # calculate posterior of M and extract quantiles
    qtls <- seq(0.5/nboot, 1 - 0.5/nboot, length = nboot)
    if (length(ncarc) == 1){
      names(ncarc) <- "total"
      names(x$psi) <- "total"
    }
    for (ti in names(ncarc)){
      output$psi[, ti] <- findInterval(qtls, cumsum(eoa::postM(ncarc[ti], psi[ti])))
    }
    output$ncarc <- ncarc
  } else if (is.null(nboot) || length(x)/length(ncarc) == nboot){
    nboot <- ifelse(is.null(nboot), length(x)/length(ncarc), nboot)
    # CASE 2: nsim and nboot are equal
    # generate one random dwp for each psi
    if (length(ncarc) == 1) {
      names(ncarc) <- "total"
      psi <- matrix(psi, ncol = 1, dimnames = list(NULL, "total"))
    }
    tname <- names(ncarc)
    output <- array(0, dim = dim(psi), dimnames = list(NULL, colnames(psi)))
    for (ti in tname){
      # try fitting beta distribution to psi (increases speed by 100x or so)
      muB <- mean(psi[, ti]); sig2B <- var(psi[, ti])
      Ba <- muB^2/sig2B*(1 - muB) - muB; Bb <- Ba*(1/muB - 1)
      ab <- suppressWarnings(try(MASS::fitdistr(x = psi[, ti], densfun = "beta",
        start = list(shape1 = Ba, shape2 = Bb)), silent = T))
      if (!"try-error" %in% class(ab)){
        pm <- eoa::postM.ab(x = ncarc[ti], Ba = Ba, Bb = Bb, prior = "IbetabinRef")
        qtls <- seq(0.5/nboot, 1 - 0.5/nboot, length = nboot)
        output[, ti] <- ncarc[ti]/sample(findInterval(qtls, cumsum(pm)))
      } else {
        for (bi in 1:nboot){
          pm <- cumsum(eoa::postM(ncarc[ti], psi[bi, ti]))
          output[bi, ti] <- ncarc[ti]/sample(x = 1:length(pm) - 1, size = 1,
            prob = pm)
        }
      }
    }
  } else if (attr(x, "nsim") != nboot){
    tname <- names(ncarc)
    output <- array(0, dim = dim(psi), dimnames = list(NULL, colnames(psi)))
    for (ti in tname){
      for (bi in 1:nboot){
        pm <- cumsum(eoa::postM(ncarc[ti], x$psi[bi, ti]))
        output[bi, ti] <- ncarc[ti]/sample(x = 1:length(pm) - 1, size = 1,
          prob = pm)
      }
    }
  }
  if (ncol(output) == 1) { # single turbine, nsim reps
    if (!forGenEst){
      output <- c(output)
      names(output) <- names(ncarc)
    } else {
      output <- cbind(turbine = "t1", dwp = c(output))
    }
  }
  if (forGenEst){
    nms <- exclude("total", names(ncarc))
    output <- data.frame(
      turbine = rep(nms, each = nboot),
      dwp = c(output[, exclude("total", names(ncarc))])
    )
  }
  output
}

#' Read and Format Simple Polygons Data
#' @export
readLayout_polygons <- function(filename, unitCol = "turbine"){
  tmp <- read.csv(filename, stringsAsFactors = F)
  rnm <- tmp[, unitCol]
  tmp <- as.matrix(tmp[, c("x", "y")], ncol = 2)
  rownames(tmp) <- rnm
  pgon <- list()
  for (ti in unique(rnm)) pgon[[ti]] <- tmp[which(rnm == ti), c("x", "y")]
  attr(pgon, "srad") <- ceiling(sqrt(max(rowSums(tmp^2))))
  class(pgon) <- "polygonLayout"
  pgon
}

#' @export
plot.layoutSimple <- function(x, ...){
  par(mfrow = c(ceiling(sqrt(nrow(x))), round(sqrt(nrow(x)))),
    mar = c(0, 0, 0, 0), oma = c(2, 2, 1, 1), mgp = c(2, 0.7, 0), tck = -0.015)
  rmax <- max(x$radius)
  theta <- seq(0, 2* pi, length = 500)
  col <- colors()[647]
  for (i in 1:nrow(x)){
    plot(0, xlim = rmax * c(-1, 1) * 1.1, ylim = rmax * c(-1, 1) * 1.1,
      asp = 1, type = "n", axes = F)
    axis(1); axis(2)
    box()
    with(x[i, ], {
      switch(shape,
        circular = plotrix::draw.circle(0, 0, radius, nv = 500, col = col),
        square = rect(-radius, -radius, radius, radius, col = col),
        RP = {
          if (n_road == 2){
            if (rad_pad > width_road){
              theta0 <- asin(width_road/(2 * rad_pad))
              xx <- c(
                sqrt(rad_pad^2 - (width_road/2)^2),
                rad_pad * cos(theta[theta >= theta0 & theta <= pi - theta0]),
                -sqrt(rad_pad^2 - (width_road/2)^2),
                -radius, -radius,
                -sqrt(rad_pad^2 - (width_road/2)^2),
                rad_pad * cos(theta[theta >= pi + theta0 & theta <= 2 * pi - theta0]),
                sqrt(rad_pad^2 - (width_road/2)^2),
                radius, radius,
                sqrt(rad_pad^2 - (width_road/2)^2)
              )
              yy <- c(
                width_road/2,
                rad_pad * sin(theta[theta >= theta0 & theta <= pi - theta0]),
                width_road/2,
                width_road/2, -width_road/2,
                -width_road/2,
                rad_pad * sin(theta[theta >= pi + theta0 & theta <= 2 * pi - theta0]),
                -width_road/2,
                -width_road/2, width_road/2, width_road/2
              )
              polygon(xx, yy, col = col)
            } else {
              rect(-radius, -radius, radius, radius, col = col)
            }
          } else if (n_road == 1){
            theta0 <- asin(width_road/(2 * rad_pad))
            xx <- c(
              rad_pad * cos(theta[theta >= theta0 & theta <= 2*pi - theta0]),
              radius, radius
            )
            yy <- c(
              rad_pad * sin(theta[theta >= theta0 & theta <= 2*pi - theta0]),
              -width_road/2, width_road/2
            )
            polygon(xx, yy, col = col)
          }
        }
      )
      mtext(side = 1, line = -1.5, paste0(turbine, ", r = ", radius, " m"), adj = 0.1)
    })
  }
}




#' Prepare Ring Data for Data Frames of Simple Turbine Layouts
#' @param x a data frame with simple characterizations of turbines
#' @param cod carcass observation distances. Data frame with columns for
#'  distances that carcasses were found at and which turbines. If \code{cod = NULL},
#'  an empty template is created for carcass with descriptions of the turbines
#'  (i.e., \code{r}, \code{exposure}, \code{ncarc}, \code{pinc}), and (optional)
#'  \code{scVar}. A separate function can then be used to add carcasses distances later. This setup allows for easier
#'  simulation of carcasses counts using the predefined template for turbines.
#'  A "total" is added for the site...sums the carcasses and exposures
#' @param (optional) data frame or vector with carcass distances and turbine names
#' @export
prepRing_simple <- function(x, cod = NULL, unitCol = "turbine", rCol = "r"){
  rdat <- list()
  rpA <- list()
  rownames(x) <- x[, unitCol]
  srad <- ceiling(max(
    slayout[, "radius"] * ifelse(slayout[, "shape"] == "square", sqrt(2), 1)
  ))
  totA <- data.frame(r = 1:srad, pinc = 0)
  totr <- data.frame(r = 1:srad, exposure = 0, ncarc = 0)
  for (ti in x[, unitCol]){
    rmax <- ceiling(ifelse(x[ti, "shape"] %in% c("circular", "square"),
      x[ti, "radius"], sqrt(2) * x[ti, "radius"]))
    rdat[[ti]] <- data.frame(r = 1:rmax, exposure = 0, ncarc = 0)
    expo <- rdat[[ti]]$exposure
    rdat[[ti]]$exposure = switch(x[ti, "shape"],
      circular = 2*pi*(rdat[[ti]]$r - 1/2) * (rdat[[ti]]$r <= x[ti, "radius"]),
      square = diff(Acins(r = 0:rmax, s = x[ti, "radius"])),
      RP = {
        expo <- numeric(rmax) # 1 exposure for each radius = 1, ..., rmax
        rout <- ceiling(x[ti, "pad_rad"]:rmax)
        expo[rout] <- x[ti, "width_road"] * x[ti, "n_road"]
        rin <- which(rdat[[ti]]$r < x[ti, "rad_pad"])
        expo[rin] <- 2*pi*(rin - 1/2)
        expo
      }
    )
    rpA[[ti]] <- rdat[[ti]][, c("r", "exposure")]
    names(rpA[[ti]]) <- c("r", "pinc")
    rpA[[ti]][, "pinc"] <- pmin(rdat[[ti]][, "exposure"]/(2*pi*(1:rmax - 1/2)), 1)
    totA[1:rmax, "pinc"] <- totA[1:rmax, "pinc"] + rpA[[ti]][, "pinc"]
    totr[1:rmax, "exposure"] <- totr[1:rmax, "exposure"] + rpA[[ti]][, "exposure"]
  }
  totA[, "pinc"] <- pmin(totA[, "pinc"]/nrow(x), 1)
  totr[, "exposure"] <- totr[, "exposure"]/nrow(x)
  rpA[["total"]] <- totA
  rdat[["total"]] <- totr
  ncarc <- numeric(nrow(x))
  names(ncarc) <- x[, unitCol]
  output <- list(rdat = rd, rpA = rpA, srad = srad, ncarc = ncarc, scVar = NULL)
  class(output) <- "rings"
  return(output)
}

#' Calculate Area of Intersection of Circle and Square with Common Center
#' @param r radius of the circle (vector or scalar)
#' @param s half-width of the square (scalar)
#' @examples
#'  # calculate area in annulus intersecting square
#'  diff(Acins(r = c(inner = 11, outer = 12), s = 10.5))
#'  # calculate area in series of 1 m annuli extending to corner of square
#'  s <- 10.5 # radius of square (center to side)
#'  diff(Acins(r = 0:ceiling(sqrt(2) * s), s))
#' @export
Acins <- function(r, s){
  ans <- numeric(length(r))
  ans[r <= s] <- pi * r[r <= s]^2
  im <- which(r > s & r < sqrt(2) * s)
  ans[im] <- 4*(s*sqrt(r[im]^2 - s^2) + r[im]^2*(pi/2 - 2*acos(s/r[im]))/2)
  ans[r >= sqrt(2) * s] <- 4 * s^2
  ans
}

############
#' Add Carcasses to a Site Layout
#'
# addCarcasses:
#' @export
addCarcass<- function(x, ...) UseMethod("addCarcass", x)
#' @rdname addCarcass
#' @export
addCarcass.shapeCarcass <- function(x, data_ring, plotLayout = NULL,
    ncarcReset = TRUE, ...){
  unitCol <- x$unitCol
  scVar <- data_ring$scVar
  turbi <- x$carcasses[, unitCol, drop = T]
  if (any(!x$carcasses[, unitCol, drop = T] %in% names(data_ring$ncarc))){
    catturb <- x$carcasses[, unitCol, drop = T]
    turbnot <- catturb[!catturb %in% names(data_ring$ncarc)]
    warning("Carcasses found at unsearched turbines.")
    if (all(!x$carcasses[, unitCol, drop = T] %in% names(data_ring$ncarc)))
      stop("All carcasses were found at unsearched turbines. Mismatched ",
           "turbine names in carcass data (x) and data_ring?")
  }
  if (!is.null(scVar)){
    # need layout data (in addition to rings) OR carcasses in x to have scVar
    # levels that match up with distances and classes in data_ring
    if (!is.null(plotLayout)){
      if (!"shapeLayout" %in% class(plotLayout))
        stop("plotLayout must be a shapeLayout object in addCarcass.shapeCarcass")
      # check whether carcasses (at x, y) were found in areas searched, and
      #  attach scVar levels to them
      # indices for polygons that carcasses lie in:
      ind <- unlist(sf::st_within(x$carcasses, plotLayout$layout))
      if (length(ind) < nrow(x$carcasses)){
        stop("carcasses found in area not searched.")
      }
      rdat0 <- data.frame( # data for rdat [not properly tallied or formatted]
        turbine = x$carcasses[, unitCol, drop = T],
        scVar = plotLayout$layout[ind, scVar, drop = T],
        r = ceiling(unname(sqrt(rowSums(
          (plotLayout$tcenter[turbi, ] - sf::st_coordinates(x$carcasses))^2
        ))))
      )
      names(rdat0) %<>% gsub("scVar", scVar, .)
    } else {
      # check whether carcass search classes are represented at the distances
      # given in x
      if (!scVar %in% names(x$carcasses))
        stop("scVar not a column in carcass data in addCarcass.shapeCarcass")
      if (any(!x$carcasses[, scVar, drop = T] %in% unique(data_ring$rdat$total[, scVar])))
        stop("Search classes in x$carcasses do not match those of plotLayout")
      rdat0 <- data.frame(
        turbine = turbi,
        scVar = x$carcasses[, scVar, drop = T],
        r = ceiling(unname(sqrt(rowSums(
          (data_ring$tcenter[turbi, ] - sf::st_coordinates(x$carcasses))^2
        ))))
      )
      names(rdat0) %<>% gsub("scVar", scVar, .)
    }
  } else {
    rdat0 <- data.frame(
      turbine = x$carcasses[, unitCol, drop = T],
      r = ceiling(unname(sqrt(rowSums(
        (data_ring$tcenter[turbi,] - sf::st_coordinates(x$carcasses))^2
      ))))
    )
  }
  if (ncarcReset){ # remove previous carcasses from the layout
    data_ring$rdat <- lapply(data_ring$rdat, FUN = function(ri){
      ri$ncarc <- 0
      ri
    })
    data_ring$ncarc %<>% `*` (., 0)
  }
  # add carcass tallies to the ncarc column in ringData
  for(ti in exclude("total", names(data_ring[["rdat"]]))){
    tbl <- table(rdat0[rdat0[, "turbine"] == ti, c(scVar, "r")])
    if (is.null(scVar)){
      data_ring$rdat[[ti]][names(tbl), "ncarc"] %<>% `+` (., unname(tbl))
      data_ring$rdat[["total"]][names(tbl), "ncarc"] %<>% `+` (., unname(tbl))
    } else {
      for (sci in rownames(tbl)){
        ind <- which(data_ring$rdat[[ti]][, scVar] == sci)
        for (ri in colnames(tbl)){
          r <- as.numeric(ri)
          if (r %in% data_ring$rdat[[ti]][ind, "r"]){
            data_ring$rdat[[ti]][r, "ncarc"] %<>% `+` (., tbl[sci, ri])
            tmp <- data_ring$rdat[["total"]]
            tmp$ncarc[tmp$r == r & tmp[, scVar] == sci] %<>% `+`(., tbl[sci, ri])
            data_ring$rdat[["total"]] <- tmp
          }
        }
      }
    }
  }
  data_ring$ncarc <- sapply(data_ring$rdat, FUN = function(ti) sum(ti$ncarc))
  class(data_ring$rdat) <- "rings"
  data_ring
}

#' @param x data frame with carcass turbines and distances
#' @param data_ring ring data for adding carcasses to
#' @rdname addCarcass
#' @export
addCarcass.data.frame <- function(x, data_ring,
    ncarcReset = TRUE, unitCol = "turbine", rCol = "r", ...){
  if (ncarcReset){
    data_ring$rdat <- lapply(data_ring$rdat, FUN = function(ri){
      ri$ncarc <- 0
      ri
    })
    data_ring$ncarc %<>% `*` (., 0)
  }
  for (ti in unique(x[, unitCol])){
    tbl <- table(ceiling(x[x[, unitCol] == ti, rCol]))
    ind <- match(as.numeric(names(tbl)), data_ring$rdat[[ti]][, rCol])
    data_ring$rdat[[ti]][ind, "ncarc"] %<>% `+`(., unname(tbl))
    itot <- match(as.numeric(names(tbl)), data_ring$rdat[["total"]][, rCol])
    data_ring$rdat[["total"]][itot, "ncarc"] %<>% `+`(., unname(tbl))
  }
  data_ring$ncarc <- sapply(data_ring$rdat, FUN = function(ti) sum(ti[, "ncarc"]))
  class(data_ring$rdat) <- "rings"
  data_ring
}

