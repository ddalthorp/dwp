#' Create a Data Structure or Map for the Site Layout
#'
#' @description Read plot layout data and perform premliminary error-checking 
#'  and formatting. Search plot layout data can come in any of several different
#'  formats, including shape files for search area polygons and turbine locations,
#'  R polygons, (x, y) coordinates, or simple description of search plot type for
#'  each turbine (square, circular, road & pad). A vector of distances along with
#'  a search radius is also accommodated by \code{dwp}, but these can be directly
#'  processed in \code{\link{prepRing}} without preprocessing in \code{initLayout}.

#' @param data_layout Either the name of a shape file (polygons or multipolygons)
#'  that delineates areas searched at each turbine; a .csv file with R polygons,
#'  (x, y) coordinates, or simple descriptions of search parameters at each turbine;
#'  or a data frame with r polygons, (x, y) coordinates, or simple plot layout
#'  descriptions. See "Details" for details.
#' @param dataType An identifier for the type of data to be imported: \code{"shape"},
#'  \code{"polygon"}, \code{"xy"}, or \code{"simple"}. If \code{data_layout} is
#'  the name of a shape file, the \code{dataType = "shape"} identifier is optional.
#' @param unitCol Column name for turbine IDs. If \code{data_layout} is the name
#'  of a shape file, then \code{file_turbine} must also be provided, giving
#'  turbine locations. The \code{unitCol} must be present in \code{data_layout}
#'  and in \code{file_turbine} (if provided). Turbine IDs in the \code{unitCol}
#'  must be syntactically valid R names (see Details below).
#' @param file_turbine The name of a shape file (points) giving the turbine
#'  locations for turbines listed in the \code{unitCol} column in the
#'  \code{data_layout} if \code{data_layout} is a shape file. If 
#'  \code{dataType = "xy"} and the grids in \code{data_layout} are all centered
#'  at (0, 0) with their turbines at the center, then \code{file_turbine} is not
#'  necessary and is ignored. Otherwise, if the grid coordinates are UTMs, 
#'  \code{file_turbine} is either (1) a data frame with turbine names (in `unitCol`)
#'  and the location of turbines in `x` and `y`, or (2) the name of a .csv file 
#'  with turbine locations (`unitCol`, `x`, and `y`). \code{file_turbine} is not 
#'  required (and is ignored) for other data types.
#' @param radCol for \code{dataType = "simple"} layouts: the name of the column in
#'  \code{data_layout} that gives the search radius for each turbine
#' @param shapeCol for \code{dataType = "simple"} layouts: the name of the column in
#'  \code{data_layout} that gives the plot shape for each turbine.
#' @param padCol for \code{dataType = "simple"} layouts: the name of the column in
#'  \code{data_layout} that gives the radius of the turbine pad
#' @param roadwidCol for \code{dataType = "simple"} layouts: the name of the column in
#'  \code{data_layout} that gives the width of the turbine access road(s)
#' @param nRoadCol for \code{dataType = "simple"} layouts: the name of the column in
#'  \code{data_layout} that gives the number of turbine access roads at each turbine
#' @param xCol for \code{dataType = "xy"} or \code{dataType = "polygon"} layouts: 
#'  the name of the column in \code{data_layout} that gives \code{x} coordinates 
#'  on the grid (for \code{dataType = "xy"}) or \code{x} coordinates of search area
#'  polygon (for \code{dataType = "polygon"})
#' @param yCol for \code{dataType = "xy"} or \code{dataType = "polygon"} layouts: 
#'  the name of the column in \code{data_layout} that gives \code{y} coordinates 
#'  on the grid (for \code{dataType = "yy"}) or \code{y} coordinates of search area
#'  polygon (for \code{dataType = "polygon"})
#' @param ncarcCol for \code{dataType = "xy"} layouts: the name of the column with
#'  carcass counts in each grid cell. The column is required but may be all zeros
#'  with carcasses added from a matrix of carcass locations later
#' @param scCol for \code{dataType = "xy"} layouts: the name of column in data_layout
#'  with names of search classes. This is used for excluding unsearched areas
#'  from the grid data (x, y). It is used ONLY with \code{dataType = "xy"} and
#'  used to remove rows with \code{x[, scCol] == notSearched}, where \code{x}
#'  is the search grid data frame.
#' @param notSearched for \code{dataType = "xy"} layouts: the name(s) of search 
#'  class(es) in \code{scCol} that are not searched (optional). Ignored for data 
#'  types other than \code{xy}.
#' @param quiet boolean for controlling whether progress of calculations and other
#'  notes are printed to the console as the function runs
#'
#' @details All the layout types (except for vector, which is addressed elsewhere)
#'  can accommodate patterns of seached and not searched areas. If the searched
#'  areas are subdivided into different search classes with different detection
#'  probabilities, then search plot layout data must be input either from shape files
#'  with non-intersecting polgons delineating the search classes or from x-y grid
#'  data. If there is more than one search class variable (for example, ground
#'  cover and search schedule), then the covariates may be entered in separate
#'  columns if the layout files give grid coordinates or may be combined into
#'  one column in the shape files. For example, ground visibility may be easy or
#'  difficult and search schedule may be 1-day or 14-day. These can be combined
#'  into a single column with values of, say, \code{easy1}, \code{easy14},
#'  \code{difficult1}, and \code{difficult14}.
#'
#'  There must be a turbine ID column (\code{unitCol}) in each of the files. The 
#'  individual turbine ID's must be syntactically valid R names, i.e., contain 
#'  combinations of letters, numbers, dot ( . ), and underscores ( _ ) only and 
#'  not begin with a number or a dot followed by a number. Other characters are 
#'  not allowed: no spaces, hyphens, parentheses, etc.
#'
#'  If shape files are to be imported, both shape files (search area polygons and
#'  turbine locations) must have their three standard, mandatory component files
#'  (.shp, .shx, .dbf) stored in the same folder. Only the name of the .shp should
#'  be entered in the function call. Other components are automatically searched
#'  for and processed if available.
#'
#' @return A list or data frame with components appropriate to the type of data
#'  imported. The data structure is returned as an S3 class object, which other
#'  functions in \code{dwp} can recognize and properly process. There is minimal
#'  processing on the data after importing, but the structures are lightly error-checked
#'  and formatted for more thorough processing, depending on data type and analysis
#'  objectives. Typically, the layout data will be later processed by a call to
#'  \code{\link{prepRing}} to create a characterization of the searched area at 
#'  the site by "rings", with tallies of searched area, search classes, and 
#'  fraction of area searched in concentric, 1 meter rings around each turbine. 
#'  The format of the output depends on the format of the input. There are 
#'  several possibilities, including, each of which is an S3 object with 
#'  characteristics specific to the imported data:
#'  \describe{
#'    \item{\code{shapeLayout}}{List with elements:
#'      \itemize{
#'        \item \code{$layout} = turbine search area configurations (polygons 
#'          and multipolygons) from \code{data_layout} shape file as an 
#'          \code{\link[sf]{sf}} object.
#'        \item \code{$layoutAdj} = polygons from \code{$layout} but recentered at (0, 0)
#'        \item \code{$turbines} = turbine data (as \code{\link[sf]{sf}} object)
#'        \item \code{$unitCol} = name of the column with turbine IDs (character string)
#'        \item \code{$tset} = turbine names (vector of character strings)
#'        \item \code{$tcenter} = locations of turbine centers (nturb by 2 array)
#'          with UTMs of turbine locations, extracted and simplified from
#'          \code{$turbines}. Column names are \code{X} and \code{Y}, measuring
#'          meters from a reference point. Row names are the names of the turbines.
#'      }
#'    }
#'    \item{\code{simpleLayout}}{Data frame with columns:
#'      \itemize{
#'        \item \code{turbine} = turbine IDs (syntactically valid R names)
#'        \item \code{radius} = search radius. If \code{shape = "square"}, then
#'          radius is 1/2 the width of the square.
#'        \item \code{shape} = general descriptor of the shape of the search plot as
#'          \code{"square"}, \code{"circular"}, or \code{"RP"} (for roads and pads search).
#'        \item \code{padrad} = radius of the turbine pad (assumed circular)
#'        \item \code{roadwidth} = width of the access road(s)
#'        \item \code{n_road} = number of access roads
#'      }
#'    }
#'    \item{\code{polygonLayout}}{
#'      List of polygons, one for each turbine. The maximum search
#'      radius at any turbine is assigned as an attribute (\code{attr(, "rad")}).
#'    }
#'    \item{\code{xyLayout}}{List with elements:
#'      \itemize{
#'        \item \code{xydat} = data frame with columns for turbine names, x and y
#'          coordinates of 1m grid centers spanning the searched area, number of
#'          carcasses found in each grid cell, and optional covariates.
#'        \item \code{tcenter} = matrix giving turbine locations (\code{x, y}), with
#'          row names = turbine names.
#'        \item \code{ncarc} = vector giving the number of carcasses found at each
#'          turbine.
#'        \item \code{unitCol} = name of the column where turbine IDs are stored in
#'          \code{xydat}.
#'        \item \code{tset} = names of the searched turbines
#'      }
#'    }
#'  }
#'
#' @examples
#' data(layout_simple)
#' # converts properly formatted dataframe to 'simpleLayout' object
#' initLayout(layout_simple) 
#' 
#' data(layout_xy)
#' initLayout(layout_xy, dataType = "xy")
#'
#' data(layout_polygon)
#' initLayout(layout_polygon, dataType = "polygon", unitCol = "turbine")
#'
#' @export

initLayout <- function(data_layout, dataType = "simple",  unitCol = "turbine",
    file_turbine = NULL, radCol = "radius", shapeCol = "shape", padCol = "padrad",
    roadwidCol = "roadwidth", nRoadCol = "n_road", xCol = "x", yCol = "y",
    ncarcCol = "ncarc", scCol = NULL, notSearched = NULL, quiet = FALSE){
  # traffic directing
  dataType <- tolower(dataType)
  if (dataType == "shape" || (length(data_layout) == 1 &&
      is.character(data_layout) && grepl(".shp", data_layout))){
    suparg <- NULL
    if (!identical(radCol, "radius") && !is.null(radCol) && !is.na(radCol)) 
      suparg <- c(suparg, "radCol")
    if (!identical(shapeCol, "shape") && !is.null(shapeCol) && !is.na(shapeCol)) 
      suparg <- c(suparg, "shapeCol")
    if (!identical(padCol, "padrad") && !is.null(padCol) && !is.na(padCol)) 
      suparg <- c(suparg, "padCol")
    if (!identical(roadwidCol, "roadwidth") && !is.null(roadwidCol) && !is.na(roadwidCol)) 
      suparg <- c(suparg, "roadwidCol")
    if (!identical(nRoadCol, "n_road") && !is.null(nRoadCol) && !is.na(nRoadCol)) 
      suparg <- c(suparg, "nRoadCol")
    if (!identical(xCol, "x") && !is.null(xCol) && !is.na(xCol)) 
      suparg <- c(suparg, "xCol")
    if (!identical(yCol, "y") && !is.null(yCol) && !is.na(yCol)) 
      suparg <- c(suparg, "yCol")
    if (!identical(ncarcCol, "ncarc") && !is.null(ncarcCol) && !is.na(ncarcCol)) 
      suparg <- c(suparg, "ncarcCol")
    if (!identical(notSearched, "y") && !is.null(notSearched) && !is.na(notSearched)) 
      suparg <- c(suparg, "notSearched")
    if (!identical(scCol, "y") && !is.null(scCol) && !is.na(scCol)) 
      suparg <- c(suparg, "scCol")
 
    if (length(suparg) > 0){ 
      isare <- ifelse(length(suparg) == 1, "is", "are")
      message(c(paste0(suparg, collapse = ", "), " ", isare, 
        " superfluous for dataType = \"", dataType, "\" and ", isare, " ignored.\n\n"))
    }

    plotLayout <- sf::st_zm(sf::st_read(
      data_layout, stringsAsFactors = FALSE, quiet = quiet),
      drop = T, what = "ZM")

    plotLayout0 <- sf::st_drop_geometry(plotLayout) # easier quick subsetting
    turbines <- sf::st_zm(sf::st_read(
      file_turbine, stringsAsFactors = FALSE, quiet = quiet),
      drop = T, what = "ZM")
    turbines0 <- sf::st_drop_geometry(turbines)
    if (!is.null(unitCol) && !is.na(unitCol)){
      if (is.null(file_turbine))
        stop("if unitCol is not NULL, file_turbine must be provided")
      if (!unitCol %in% names(turbines))
        stop(unitCol, " column not included in turbine data")
      if (!unitCol %in% names(plotLayout))
        stop(unitCol, " column not included in plot layout data")
      if (!all(plotLayout0[, unitCol] %in% turbines0[, unitCol])){
        stop("a turbine that is included in ", data_layout,
             " is missing from ", file_turbine
        )
      }
      # format of unitCol names: must be syntactically valid R names
      # make.names(x) = x if x is a valid name, else converts x to a valid name
      if (!all(make.names(turbines0[, unitCol]) == turbines0[, unitCol])){
        badind <- which(make.names(turbines0[, unitCol]) != turbines0[, unitCol])
        badnm <- unique(turbines0[badind, unitCol])
        badnm <- badnm[1:min(length(badnm), 3)]
        goodnm <- gsub("^X", "t", make.names(badnm))
        if(!quiet){
          message(
            "\n\nNOTE: Not all the turbine names are syntactically valid.\n",
            "These (", paste(badnm, collapse = ", "),
            "...) will be converted to syntactically valid names (",
            paste(goodnm, collapse = ", "), "...). "
          )
        }
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
      stop("initLayout with shape files must have file_turbine and unitCol.")
    }
    layoutAdj <- plotLayout # Adj => coordinates relative to turbine centers
    for (ai in 1:nrow(plotLayout)){
      ti <- plotLayout[ai, unitCol, drop = TRUE]
      layoutAdj[ai, ] <- sf::st_set_geometry(layoutAdj[ai, ],
        sf::st_geometry(plotLayout[ai, ]) - tcenter[ti,]) # recentering
    }
    # maximum distance from turbine to a searched point (radius of interpolation)
    output <- list(layout = plotLayout, layoutAdj = layoutAdj,
      turbines = turbines, unitCol = unitCol, tset = tset, tcenter = tcenter)
    class(output) <- "shapeLayout"
    return(output)
  } else if (dataType == "simple"){
    if (length(data_layout) == 1 && is.character(data_layout)){
      slayout <- read.csv(data_layout, stringsAsFactors = FALSE)
    } else if (is.data.frame(data_layout)){
      slayout <- data_layout
      slayout[] <- lapply(slayout, 
        function(x) if(is.factor(x)) as.character(x) else x)
    }
    if ("simpleLayout" %in% class(data_layout)) return(data_layout)
    # rudimentary error-checking:
    suparg <- NULL

    if (!is.null(file_turbine) && !is.na(file_turbine)) 
      suparg <- c(suparg, "file_turbine")
    if (!identical(xCol, "x") && !is.null(xCol) && !is.na(xCol)) 
      suparg <- c(suparg, "xCol")
    if (!identical(yCol, "y") && !is.null(yCol) && !is.na(yCol)) 
      suparg <- c(suparg, "yCol")
    if (!identical(ncarcCol, "ncarc") && !is.null(ncarcCol) && !is.na(ncarcCol)) 
      suparg <- c(suparg, "ncarcCol")
    if (!is.null(scCol) && !is.na(scCol)) 
      suparg <- c(suparg, "scCol")
    if (!is.null(notSearched) && !is.na(notSearched)) 
      suparg <- c(suparg, "notSearched")
    
    if (length(suparg) > 0){ 
      isare <- ifelse(length(suparg) == 1, "is", "are")
      message(c(paste0(suparg, collapse = ", "), " ", isare, 
        " superfluous for dataType = \"", dataType, "\" and ", isare, " ignored.\n\n"))
    }

    if (!all(c(unitCol, radCol, shapeCol) %in% names(slayout))){
      stop("Simple layout must include columns for unit, radius, and shape ",
           "with column names specified in arg list for initLayout.")
    }
    if (any(slayout[, shapeCol] == "RP") &&
       !all(c(padCol, roadwidCol, nRoadCol) %in% names(slayout))){
      stop(
        "Simple layout that includes RP turbines must include columns for pad ",
        "radius, road width, and number of roads in the arg list for initLayout."
      )
    }
    if (!is.numeric(slayout[, radCol]) || any(slayout[, radCol] <= 0))
      stop("All search radii in simple layout must be positive numbers.")
    if (!all(slayout[, shapeCol] %in% c("circular", "square", "RP")))
      stop("Simple \"shape\" must be circular, square, or RP")
    srad <- slayout[, radCol]
    ind <- which(slayout[, shapeCol] == "square")
    srad[ind] <- sqrt(2) * srad[ind]
    # sterotype the names:
    names(slayout)[match(unitCol, names(slayout))] <- "turbine"
    names(slayout)[match(radCol, names(slayout))] <- "radius"
    names(slayout)[match(shapeCol, names(slayout))] <- "shape"
    names(slayout)[match(padCol, names(slayout))] <- "padrad"
    names(slayout)[match(roadwidCol, names(slayout))] <- "roadwidth"
    names(slayout)[match(nRoadCol, names(slayout))] <- "n_road"
    attr(slayout, "srad") <- ceiling(max(srad))
    attr(slayout, "unitCol") <- unitCol
    class(slayout) <- c("simpleLayout", "data.frame")
    return(slayout)
  } else if (dataType == "polygon"){
    if (is.data.frame(data_layout) || is.matrix(data_layout)){
      playout <- data_layout
      playout[] <- lapply(playout, 
        function(x) if(is.factor(x)) as.character(x) else x)
    } else if (length(data_layout) == 1 && is.character(data_layout)){
      playout <- read.csv(data_layout, stringsAsFactors = FALSE)
    } else {
      stop("data_layout must be a name of a file; alternatively, ",
           "data_layout may be a properly formatted data frame or matrix"
      )
    }
    suparg <- NULL
    if (!is.null(file_turbine) && !is.na(file_turbine)) 
      suparg <- c(suparg, "file_turbine")
    if (!identical(radCol, "radius") && !is.null(radCol) && !is.na(radCol)) 
      suparg <- c(suparg, "radCol")
    if (!identical(shapeCol, "shape") && !is.null(shapeCol) && !is.na(shapeCol)) 
      suparg <- c(suparg, "shapeCol")
    if (!identical(padCol, "padrad") && !is.null(padCol) && !is.na(padCol)) 
      suparg <- c(suparg, "padCol")
    if (!identical(roadwidCol, "roadwidth") && !is.null(roadwidCol) && !is.na(roadwidCol)) 
      suparg <- c(suparg, "roadwidCol")
    if (!identical(nRoadCol, "n_road") && !is.null(nRoadCol) && !is.na(nRoadCol)) 
      suparg <- c(suparg, "nRoadCol")
    if (!identical(ncarcCol, "ncarc") && !is.null(ncarcCol) && !is.na(ncarcCol)) 
      suparg <- c(suparg, "ncarcCol")
    if (!identical(scCol, "sc") && !is.null(scCol) && !is.na(scCol)) 
      suparg <- c(suparg, "scCol")
    if (!identical(notSearched, "sc") && !is.null(notSearched) && !is.na(notSearched)) 
      suparg <- c(suparg, "notSearched")
    if (length(suparg) > 0){ 
      isare <- ifelse(length(suparg) == 1, "is", "are")
      message(c(paste0(suparg, collapse = ", "), " ", isare, 
        " superfluous for dataType = \"", dataType, "\" and ", isare, " ignored.\n\n"))
    }
    if (!is.matrix(playout) && !is.data.frame(playout))
      stop("Data in data_layout must be a properly formatted matrix or data frame")
    if (!xCol %in% colnames(playout) || !yCol %in% colnames(playout))
      stop("Data in data_layout must have columns for x and y coordinates")
    if (is.data.frame(playout)){
      rnm <- playout[, unitCol]
      playout <- as.matrix(playout[, c(xCol, yCol)], ncol = 2)
      rownames(playout) <- rnm
    } else {
      if (is.null(row.names(playout)))
        stop("If data_layout is a matrix, row names must be turbine names")
      if (min(table(rownames(playout))) < 3)
        stop("data_layout data must be polygons")
      playout <- playout[, c(xCol, yCol)]
      rnm <- rownames(playout)
    }
    pgon <- list()
    for (ti in unique(rnm)) pgon[[ti]] <- playout[which(rnm == ti), c("x", "y")]
    attr(pgon, "srad") <- ceiling(sqrt(max(rowSums(playout^2))))
    class(pgon) <- "polygonLayout"
    return(pgon)
  } else if (dataType == "xy"){
    if (is.character(data_layout)){
      xylayout <- read.csv(data_layout, stringsAsFactors = FALSE)
    } else if (is.data.frame(data_layout)){
      xylayout <- data_layout
      xylayout[] <- lapply(xylayout, 
        function(x) if (is.factor(x)) as.character(x) else x)      
    }
    suparg <- NULL
    if (!identical(radCol, "radius") && !is.null(radCol) && !is.na(radCol)) 
      suparg <- c(suparg, "radCol")
    if (!identical(shapeCol, "shape") && !is.null(shapeCol) && !is.na(shapeCol)) 
      suparg <- c(suparg, "shapeCol")
    if (!identical(padCol, "padrad") && !is.null(padCol) && !is.na(padCol)) 
      suparg <- c(suparg, "padCol")
    if (!identical(roadwidCol, "roadwidth") && !is.null(roadwidCol) && !is.na(roadwidCol)) 
      suparg <- c(suparg, "roadwidCol")
    if (!identical(nRoadCol, "n_road") && !is.null(nRoadCol) && !is.na(nRoadCol)) 
      suparg <- c(suparg, "nRoadCol")
    if (length(suparg) > 0){ 
      isare <- ifelse(length(suparg) == 1, "is", "are")
      message(c(paste0(suparg, collapse = ", "), " ", isare, 
        " superfluous for dataType = \"", dataType, "\" and ", isare, " ignored.\n\n"))
    }
    if (!unitCol %in% names(xylayout))
      stop("unitCol must be included in data_layout")
    if (!identical(make.names(xylayout[, unitCol]), xylayout[, unitCol])){
      stop("Turbine names in data_layout must be sytactically valid names in R.",
           "Must not start with a number or dot (.) followed by a number, and ",
           "must contain combinations of letters, numbers, dots( . ), and ",
           "underscores ( _ )."
      )
    }
    if (!is.null(scCol) && scCol %in% names(xylayout) && !is.null(notSearched)){
      xylayout <- xylayout[!xylayout[, scCol] %in% notSearched, ]
    }
    if (!is.null(file_turbine)){
      if (is.character(file_turbine)){
        xyturbine <- read.csv(file_turbine, stringsAsFactors = FALSE)
      } else if (is.data.frame(file_turbine)){
        xyturbine[] <- lapply(xyturbine, 
          function(x) if (is.factor(x)) as.character(x) else x)
      }
    } else {
      tmp <- sqrt(xylayout$x^2 + xylayout$y^2)
      if (min(tmp) > 200 | max(tmp) > 2000){ 
        stop(
          "initLayout: coordinates in data_layout do not appear to be ",
          "relative to turbine locations, i.e., with each turbine located at ",
          "(0, 0) relative to its grid coordinates. Either revise the ",
          "coordinates or provide a file (or data frame) with the turbine ",
          "locations in the file_turbine argument."
        )
      }
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
      stop("Some units in data_layout are missing their counterparts in ",
           "file_turbine")
    # reset the xy coordinates to be relative to turbine locations
    for (ti in unique(xylayout[, unitCol])){
      lind <- which(xylayout[, unitCol] == ti)
      tind <- which(xyturbine[, unitCol] == ti)
      xylayout[lind, xCol] %<>% `-`(., xyturbine[tind, xCol])
      xylayout[lind, yCol] %<>% `-`(., xyturbine[tind, yCol])
    }
    xylayout$r <- sqrt(rowSums(xylayout[, c(xCol, yCol)]^2))
    nms <- xyturbine[, unitCol]
    xyturbine <- as.matrix(xyturbine[, c(xCol, yCol)])
    tbl <- aggregate(ncarc ~ turbine, xylayout, sum)
    ncarc <- tbl$ncarc
    names(ncarc) <- tbl$turbine
    ncarc <- ncarc[gtools::mixedsort(names(ncarc))]
    ncarc <- c(ncarc, total = sum(ncarc))
    rownames(xyturbine) <- nms
    output <- list(xydat = xylayout, tcenter = xyturbine, ncarc = ncarc,
      unitCol = unitCol, tset = unique(xylayout[, unitCol]))
    class(output) <- "xyLayout"
    return(output)
  }
}

#' Import Carcass Observations Locations from Shape Files
#'
#' Carcass coordiates (x, y) and turbine IDs are read using \code{\link[sf]{st_read}} 
#' and formatted for adding to \code{rings} data structures for analysis.
#'
#' @param file_cod name of the file with carcass observation data. Currently, the
#'  function requires a shape file, which gives the carcass locations on the same
#'  coordinate system that is used for the turbines. The geometry is a simple
#'  features points file, consisting of at least the three mandatory files standard
#'  components (.shp, .shx, .dbf) stored in the same directory. Only the name of the
#'  .shp is required (for example, \code{file_cod = "carcasses.shp"}). Other
#'  components are automatically searched for and processed if available.
#' @param unitCol name of column with turbine IDs. Column name and turbine IDs
#'  must match those of the \code{data_layout} and \code{file_turbine} used in
#'  the call to \code{\link{initLayout}}.
#' @param quiet boolean for directing the function to print calculation progress
#'  updates and messages to the console. This should be set to \code{FALSE} unless
#'  you know clearly why you want to turn off the messaging.
#'
#' @return a \code{shapeCarcass} object, which is a list with \code{$carcasses},
#'  which is a \code{sf} representation of the shape file \code{file_cod} data;
#'  \code{$unitCol}, which is the name of the unit column; and \code{$ncarc}, which
#'  is a vector of carcass counts at the turbines listed in \code{unitCol}. The
#'  elements of the \code{$ncarc} are named by turbines at which they were found.
#'
#' @export
readCarcass <- function(file_cod, unitCol = "turbine", quiet = FALSE){
  if (grepl(".shp", file_cod)){
    carcasses <- sf::st_read(file_cod, stringsAsFactors = FALSE, quiet = quiet)
    carcasses0 <- sf::st_drop_geometry(carcasses)
    if (!unitCol %in% names(carcasses))
      stop("readCarcass: unitCol = ", unitCol, " not included in carcass data")
    if (!all(make.names(carcasses0[, unitCol, ]) == carcasses0[, unitCol])){
      badind <- which(make.names(carcasses0[, unitCol]) != carcasses0[, unitCol])
      badnm <- unique(carcasses0[badind, unitCol])
      badnm <- badnm[1:min(length(badnm), 3)]
      goodnm <- gsub("^X", "t", make.names(badnm))
      if(!quiet){
        message(
          "\n\nNOTE: Not all the turbine names are syntactically valid.\n",
          "These (", paste(badnm, collapse = ", "),
          "...) will be converted to syntactically valid names (", 
          paste(goodnm, collapse = ", "), "...). "
        )
      }
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

#' Format a Search Layout into Rings for Analysis
#'
#' A function for creating a characterization of the search plot at each turbine by
#'  rings. The ground around each turbine is divided into 1 meter concentric rings
#'  out to the limit of the search plot. The amount of area searched in each ring
#'  and search class (if a search class column is present in the data) at each
#'  turbine is calculated, along with the fraction of area searched in each ring.
#'  In addition, sum totals of the area in each ring and the average fraction of
#'  the area searched in each ring across all turbines at the site are tallied as well.
#'  This is a convenient structure for the Poisson regressions that are used to
#'  estimate the carcass distributions with respect to distance from the turbines,
#'  the probabilities of carcasses landing in the searched areas, and the fraction
#'  of carcasses in the searched area.
#'
#' @param x a search plot layout as imported and processed by 
#'  \code{\link{initLayout}} into a \code{shapeLayout}, \code{polygonLayout}, or 
#'  \code{simpleLayout} object, or a bare vector of carcass distances if search 
#'  plots are all circular with the same radius and no unsearched area within 
#'  the search radius.
#' @param scVar name of the search class variable (optional), a column in the
#'  shape file for the search polygons. \code{scVar} is ignored if \code{x} is not
#'  a \code{shapeLayout} object.
#' @param notSearched name of the search class(es) in \code{scVar} that represent
#'  unsearched areas. Applicable only if \code{x} is a \code{shapeLayout} object
#'  and \code{scVar} is provided. Polygons associated with \code{scVar} values in
#'  \code{notSearched} are not included in the rings characterization of the site.
#'  Also, turbines with no polygons that are not \code{notSearched} are not included
#'  in the rings.
#' @param silent Processing shape files into rings may take several minutes. By
#'  default, \code{prepRing} prints periodic notice of the progress of the
#'  calculations for shape files. To suppress these notices, use 
#'  \code{silent = TRUE}.
#' @param srad search radius for data when \code{x} = bare vector of carcass
#'  observation distances.
#' @param ... ignored
#' @return an object of class \code{rings}, which is a list with components
#'  \describe{
#'    \item{\code{$rdat}}{list of data frames giving the area searched 
#'      (\code{"exposure"}), in a 1 meter ring with outer radius \code{"r"} and 
#'      the number of carcasses found \code{"ncarc"} in each ring, with search 
#'      class \code{scVar} optional. There is also a summary data frame 
#'      \code{$rdat[["total"]]} that sums the exposures and carcass counts for 
#'      all turbines across the site. The \code{$rdat[["total"]]} is the data 
#'      frame used in fitting the GLMs.}
#'    \item{\code{$rpA}}{list of data frames giving the proportion of area 
#'      included in the searches (\code{"pinc"}) in each ring (\code{"r"}). and 
#'      the number of carcasses found \code{"ncarc"} in each ring, with search 
#'      class \code{scVar} optional. There is also a summary data frame that 
#'      sums the exposures and carcass counts for all turbines across the site. 
#'      The \code{$rpA} data frames are used in estimating the probability of 
#'      carcasses falling in the searched area at each turbine, which, in turn 
#'      is used for calculating \code{dwp}}
#'    \item{\code{$srad}}{the maximum search radius at any of the turbines}
#'    \item{\code{$ncarc}}{vector of the number of carcasses at each turbine with
#'      names equal to the turbine names.}
#'    \item{\code{$scVar}}{name of the search class variable(s) or \code{NULL}}
#'    \item{\code{$tcenter}}{locations of turbine centers (nturb x 2 matrix) with
#'      UTMs of turbine locations. Column names are \code{X} and \code{Y}. Row
#'      names are the names of the turbines.}
#'  }
#'
#' @export
prepRing <- function(x, ...) UseMethod("prepRing", x)

#' @rdname prepRing
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
      cls <- unique(x$layout[, scVar, drop = TRUE])
      x <- subset(x, subset = exclude(notSearched, cls), select = scVar)
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
  for (i in length(radi):2) rings[i, ] <- sf::st_difference(rings[i, ], rings[i - 1, ])
  rings$r <- radi
  sf::st_crs(rings) <- sf::st_crs(shapeLayout$layoutAdj)

  ## calculate the area in each serach class for each ring:
  trsca <- list() # nturbines x nrings arrays of areas in each search class
  scset <- gtools::mixedsort(unique(as.data.frame(shapeLayout$layoutAdj)[, scVar]))
  if (length(scset) == 0) scset <- "all"
  if (!silent) {
    cat("calculating ring areas...\n")
    flush.console()
  }
  for (sci in scset){
    if (!silent & sci != "all"){
      print(paste0(substitute(scVar), " = ", sci), quote = FALSE)
      flush.console()
    }
    ctr <- NULL
    trsca[[sci]] <- matrix(0, #trcsa: turbne, r, search class, area
      nrow = length(shapeLayout$tset), ncol = dim(rings)[1])
    rownames(trsca[[sci]]) <- shapeLayout$tset
    for (ti in shapeLayout$tset){
      ctr <- paste0(ctr, "  ", ti)
      if (!silent){
        if(nchar(ctr) > options()$width){
          cat("\n", ti)
          ctr <- NULL
        }
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
        tring <- sf::st_geometry(shapeLayout$layoutAdj[ii, ])
        jj <- unlist(sf::st_intersects(shapeLayout$layoutAdj[ii, ], rings))
        trsca[[sci]][ti, jj] <- trsca[[sci]][ti, jj] +
          as.vector(sf::st_area(sf::st_intersection(tring, rings)))
      }
    }
    cat("\n")
  }
  #Q: was turbine searched? include only those that were searched
  trsca <- lapply(trsca, FUN = function(x){
    tmp <- sapply(trsca, rowSums)
    if (is.vector(tmp)) tmp <- matrix(tmp, nrow = 1)
    x[rowSums(tmp) > 0,, drop = FALSE]
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
  class(rdat) <- "rdat"
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

#' @rdname prepRing
#' @export
#'
prepRing.simpleLayout <- function(x, ...){
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
  class(rdat) <- "rdat"
  # cull the rows with radii beyond the edge
  ncarc <- numeric(nrow(x) + 1)
  names(ncarc) <- c(x[, "turbine"], "total")
  output <- list(rdat = rdat, rpA = rpA, srad = srad, ncarc = ncarc, scVar = NULL)
  class(output) <- "rings"
  return(output)
}

#' @rdname prepRing
#' @export
#'
prepRing.numeric <- function(x, srad, ...){
  if (!is.vector(x)) stop("x must be a vector in prepRing.numeric")
  if (!is.numeric(srad) || length(srad) != 1 || srad <= 1)
    stop("srad must be numeric scalar > 1")
  srad <- ceiling(srad)
  rpA <- list()
  rdat <- list()
  rvec <- x[x <= srad]
  tbl <- table(ceiling(rvec))
  r <- 1:srad
  rdat[["total"]] <- data.frame(r = r, exposure = 2*pi*(r - 1/2), ncarc = 0)
  rdat[["total"]]$ncarc[as.numeric(names(tbl))] <- unname(tbl)
  rpA[["total"]] <- data.frame(r = r, pinc = 1)
  class(rpA) <- "rpA"
  class(rdat) <- "rdat"
  output <- list(rdat = rdat, rpA = rpA, srad = srad, ncarc = length(rvec),
    scVar = NULL)
  class(output) <- "rings"
  return(output)
}

#' @rdname prepRing
#' @export
#'
prepRing.polygonLayout <- function(x, ...){
  srad <- attr(x, "srad")
  rdat <- list()
  rpA <- list()
  ndiv <- 1000
  theta <- seq(0, 2*pi, length = ndiv) # number of segments to divide the circle into for approx area
  rr <- 1:srad # outer radii of rings
  cx <- outer(rr, cos(theta)) # approx area in 1m ring with outer radius rr = 2 pi (rr - 0.5)
  cy <- outer(rr, sin(theta))
  iarea <- list()
  rtot <- data.frame(r = 1:srad, exposure = 0, ncarc = 0)
  Atot <- data.frame(r = 1:srad, pinc = 0)
  for(ti in names(x)){
    poly1 <- sf::st_polygon(list(poly1 = rbind(x[[ti]], x[[ti]][1,])))
    rmax <- ceiling(max(sqrt(rowSums(x[[ti]]^2))))
    iarea[[ti]] <- numeric(rmax)
    for(ri in 1:rmax){
      iarea[[ti]][ri] <- sf::st_length(cx[ri, ] %>% cbind(cy[ri, ]) %>% 
        sf::st_multipoint(., dim = "XY") %>%
        sf::st_cast(., "LINESTRING") %>%
        sf::st_intersection(., poly1))
      if(is.null(iarea[[ti]][ri])) next
    }
    rdat[[ti]] <- data.frame(
      r = 1:rmax, 
      exposure = (c(0, iarea[[ti]][1:(rmax - 1)]) + iarea[[ti]])/2, 
      ncarc = 0
    )
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
  rpA <- lapply(rpA, FUN = function(ti) ti[ti[, "pinc"] > 0, ])
  class(rpA) <- "rpA"
  class(rdat) <- "rdat"
  ncarc <- numeric(length(rdat))
  names(ncarc) <- names(rdat)
  output <- list(rdat = lapply(rdat, round, 4), rpA = rpA,
    ncarc = ncarc, srad = srad, scVar = NULL)
  class(output) <- "rings"
  return(output)
}

#' Fit Distance Distribution Model(s)
#'
#' @description Fit generalized linear models (glm) for distance distribution 
#'  models corresponding to standard forms [xep1, xep01 (gamma), xep2 (Rayleigh), 
#'  xep02, xep12, xep012, xep123, xep0123 (normal-gamma with x = tau), lognormal,  
#'  truncated normal, Maxwell Boltzmann, and constant] and supplentary forms 
#'  [exponential, chi-squared, inverse gamma, and inverse Gaussian].
#'
#'  The glm is converted to a probability distribution by dividing by a
#'  normalizing constant, namely the integral of the glm evaluated from 0 to
#'  infinity. In some cases (most notably when the leading coefficient of the
#'  glm is positive so the fitted curve does not converge to zero as x increases),
#   the integral does not converge to a finite value and the glm cannot be
#'  converted to a probability distribution. In these cases, the distribution
#'  parameters are given as \code{NA}, but the fitted model itself is saved.
#'
#' @param x a search plot layout object to fit carcass distribution models to. The
#'  layout may be a data frame with columns for ring radii, exposure (or searched
#'  area in each ring), search class variable (optional), and number of carcasses
#'  in each ring;
#'
#' @param distr names (vector of character strings) of glm distribution templates
#'  to fit. Default is \code{distr = "standard"} to fit the standard models listed in the
#'  description above. Setting \code{distr = "all"} will fit both the standard
#'  models and the supplementary models. Also, any subset of the models may be fit
#'  by using, for example, \code{distr = c("xep01", "lognormal")} to fit only
#'  the \code{"xep01"} and \code{"lognormal"} models, or
#'  \code{distr = exclude(c("xep123", "constant"))} to fit all standard models except
#'  \code{"xep123"} and \code{"constant"}, or \code{distr = exclude("lognormal", 
#'  mod_all)} to fit all the models except the lognormal.
#' @param scVar Search class variable to include in the model (optional). \code{scVar}
#'  is ignored if \code{x} is not a \code{shapeLayout} or \code{xyLayout} object.
#'  If \code{x} is a \code{shapeLayout} object, \code{scVar} may be either \code{NULL}
#'  or the name of a single column with search class data. If \code{x} is an \code{xyLayout}
#'  object, \code{scVar} may be either NULL or a vector of names of search class
#'  variables to include in the models.
#' @param notSearched the name of the level (if any) in \code{scVar} that
#'  indicates an unsearched area
#' @param rCol name of the distance column (which gives the outer radii of the rings).
#'  This will be correct by default for objects coming from \code{\link{prepRing}} 
#' and will rarely need to be explicitly specified.
#' @param expoCol name of the column with the exposure, which is the area in the ring
#'  with outer radius \code{rCol}. This will be correct by default for objects
#'  coming from \code{\link{prepRing}} and will rarely need to be 
#'  explicitly specified.
#' @param ncarcCol name of the column with tallies of carcasses by ring. This
#'  will be correct by default for objects coming from \code{\link{prepRing}} 
#'  and will rarely need to be explicitly specified.
#' @param silent set \code{silent = TRUE} to suppress information printed to the
#'  console as the calculations proceed, which may be useful when running 
#'  simulations.
#' @param unitCol name of the column with turbine IDs
#' @param ... ignored
#'
#' @return A list of fitted glm models as \code{\link[=ddFit]{dd}} objects in a 
#'  \code{\link[=ddFit]{ddArray}}
#'  object if a vector of distributions is fit, or a single \code{\link[=ddFit]{dd}} 
#'  object if a single model is fit. The \code{\link[=ddFit]{dd}} objects are 
#'  lists that include the following elements:
#'  \describe{
#'    \item{\code{\link[stats]{glm}}}{the fitted model}
#'    \item{\code{$distr}}{name of the distribution (\code{"xep01"}, etc.)}
#'    \item{\code{$parms}}{vector of distribution parameter estimates (or \code{NA}
#'      if the model based on the MLE is not extensible)}
#'    \item{\code{$varbeta}}{the variance-covariance matrix of the glm parameter
#'      estimates. NOTE: This is identical to the covariance matrix from the glm,
#'      which can be extracted via \code{summary(x)$cov.unscaled}}
#'    \item{\code{$scVar}}{name of the (optional) search class variable (or \code{NULL})}
#'    \item{\code{$ncarc}}{number of carcasses}
#'    \item{\code{$aicc}}{the AICc value of the fit}
#'    \item{\code{$n}}{number of rings}
#'    \item{\code{$k}}{number of parameters}
#'    \item{\code{$srad}}{search radius}
#'  }
#'  When a \code{dd} object is printed, only a small subset of the elements are
#'  shown. To see a full list of the objects, use \code{names(x)}. The elements
#'  can be extracted in the usual R way via \code{$} or \code{[[x]]}.
#'
#' @examples
#'  data(layout_simple) 
#'  data(carcass_simple)
#'  sitedata <- initLayout(layout_simple) # initialize
#'  ringdata <- prepRing(sitedata) # format site layout data for modeling
#'  ringsWithCarcasses <- addCarcass(carcass_simple, data_ring = ringdata) # add carcasses to site
#'  distanceModels <- ddFit(ringsWithCarcasses) # fit distance models

#' @export
#'
ddFit <- function(x, ...) UseMethod("ddFit", x)

#' @rdname ddFit
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
    output[[di]] <- list()
    form <- formula(paste0("ncarc ~ ", paste(cof_name[[di]][-1], collapse = " + "),
      scl, " + ", paste0("offset(", mod_offset[[di]], ")")))
    output[[di]]$glm <- suppressWarnings(glm(
      formula = form, data = dat, family = "poisson"
    ))
    output[[di]]$distr <- di
    output[[di]]$form <- form
    output[[di]]$scVar <- scVar
    output[[di]]$coefficients <- output[[di]]$glm$coefficients
    output[[di]]$parms <- cof2parms(output[[di]]$glm$coefficients, distr = di)
    output[[di]]$varbeta <- summary(output[[di]]$glm)$cov.unscaled
    output[[di]]$ncarc <- sum(dat$ncarc)
    output[[di]]$n <- nrow(dat)
    output[[di]]$k <- output[[di]]$n - output[[di]]$glm$df.residual
    output[[di]]$aicc<- output[[di]]$glm$aic + 2 *
      output[[di]]$k*(output[[di]]$k + 1)/(output[[di]]$n - output[[di]]$k - 1)
    output[[di]]$extensible <- cofOK(output[[di]]$coefficients, di)
    class(output[[di]]) <- "dd"
  }
  if (!silent){
    anybad <- FALSE
    cat("Extensible models:\n")
    for (di in distr){
      if (!anyNA(output[[di]]$parms)) cat(" ", di, "\n")
    }
    cat("\nNon-extensible models:\n")
    for (di in distr){
      if (anyNA(output[[di]]$parms)){
        anybad <- TRUE
        cat(" ", di, "\n")
      }
    }
    if (!anybad) cat(" none\n")
    flush.console()
  }

  if(length(output) == 1) return(output[[1]])
  class(output) <- "ddArray"
  return(output)
}

#' @rdname ddFit
#' @export
ddFit.rings <- function(x, distr = "standard", scVar = NULL, rCol = "r",
    expoCol = "exposure", ncarcCol = "ncarc", silent = FALSE, ...){
  ddFit(x$rdat$total, distr = distr, scVar = scVar, rCol = rCol,
    expoCol = expoCol, ncarcCol = ncarcCol, silent = silent)
}

#' @rdname ddFit
#' @export
ddFit.list <- function(x, distr = "standard", scVar = NULL, rCol = "r",
    expoCol = "exposure", ncarcCol = "ncarc", silent = FALSE, ...){
  arglist <- as.list(match.call())
  arglist[[1]] <- NULL
  arglist$x <- x[[1]]
  do.call(ddFit, arglist)
}

#' @rdname ddFit
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
  
  ################ new-->
  for (di in distr){
    output[[di]] <- list()
    form <- formula(paste0("ncarc ~ ", paste(cof_name[[di]][-1], collapse = " + "),
      scl))
    output[[di]]$glm <- suppressWarnings(glm(
      formula = form, data = dat, family = "poisson"
    ))
    output[[di]]$distr <- di
    output[[di]]$form <- form
    output[[di]]$scVar <- scVar
    output[[di]]$coefficients <- output[[di]]$glm$coefficients
    output[[di]]$parms <- cof2parms(output[[di]]$glm$coefficients, distr = di)
    output[[di]]$varbeta <- summary(output[[di]]$glm)$cov.unscaled
    output[[di]]$ncarc <- sum(dat$ncarc)
    output[[di]]$n <- nrow(dat)
    output[[di]]$k <- output[[di]]$n - output[[di]]$glm$df.residual
    output[[di]]$aicc<- output[[di]]$glm$aic + 2 *
      output[[di]]$k*(output[[di]]$k + 1)/(output[[di]]$n - output[[di]]$k - 1)
    output[[di]]$extensible <- cofOK(output[[di]]$coefficients, di)
    class(output[[di]]) <- "dd"
  }
  if (!silent){
    anybad <- FALSE
    cat("Extensible models:\n")
    for (di in distr){
      if (!anyNA(output[[di]]$parms)) cat(" ", di, "\n")
    }
    cat("\nNon-extensible models:\n")
    for (di in distr){
      if (anyNA(output[[di]]$parms)){
        anybad <- TRUE
        cat(" ", di, "\n")
      }
    }
    if (!anybad) cat(" none\n")
    flush.console()
  }
  if(length(output) == 1) return(output[[1]])
  class(output) <- "ddArray"
  return(output)
}

#' Calculate Akaike Information Criterion (AICc) for Distance Distributions
#'
#' @description functions for calculating AICc for carcass dispersion
#'  models.
#' @param x list of models (\code{\link[=ddFit]{ddArray}}) or single model 
#'  (\code{\link[=ddFit]{dd}}) to calculate AICs for
#' @param extent Include only the extensible models (\code{extent = "full"}) or
#'  all models (\code{extent = "win"}), whether or not they can be extended 
#'  beyond the search radius.
#' @param ... ignored
#' @return Data frame with AICc and deltaAICc for all models in \code{x}
#' @export
aic <- function(x, ...) UseMethod("aic", x)

#' @rdname aic
#' @export
#'
aic.ddArray <- function(x, extent = "full", ...){
  # column names for data frame to be returned
  nm <- c("k", "AICc", "deltaAICc")
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
    dimnames = list(incmod, nm)), stringsAsFactors = FALSE)
  rownames(output) <- incmod
  ci <- c("k", "AICc")
  for (di in incmod){
    output[di, ci] <- aic(x[[di]])[ci]
    if (extent == "win")
      output[di, "extensible"] <- 1 * cofOK(x[di][["coefficients"]], di)
  }
  output$deltaAICc <- output$AICc - min(output$AICc, na.rm = TRUE)
  output <- output[order(output$deltaAICc), ]
  attr(output, "extent") <- extent
  attr(output, "n") <- attr(x, "n")
  return(output)
}

#' @rdname aic
#' @export
#'
aic.ddArraycc <- function(x, extent = "full", ...){
  lapply(x, aic, extent = extent)
}

#' @rdname aic
#' @export
#'
aic.dd <- function(x, ...){
  return(data.frame(
    k = x$k, 
    AICc = round(x$aicc, 2), 
    extensible = 1 * cofOK(x$coefficients, x$distr), 
    row.names = x[["distr"]]
  ))
}

#' @name Plot
#' @title Plot dd and ddArray Objects
#'  
#' @description Plot CDF, PDF, or rcd (relative carcass density) for a single
#'  carcass dispersion glm model (\code{\link[=ddFit]{dd}} object) or a list of 
#'  models (\code{\link[=ddFit]{ddArray}} object).
#' @param x model(s) to plot
#' @param type Type or representation of carcass dispersion to plot:
#'  \code{"CDF"}, \code{"PDF"}, or \code{"rcd"}. The \code{"CDF"} gives the
#'  fraction of carcasses falling within \code{r} meters from a turbine and
#'  \code{"PDF"} is the associated probability density. The \code{"rcd"} gives the
#'  relative carcass density at a point \code{r} meters from a turbine and is
#'  PDF/(2 * pi * r).
#' @param extent Plot dispersions as fraction of total carcasses (\code{"full"})
#'  or as fraction of carcasses within the searched area (\code{"win"}).
#' @param distr vector of names of distributions to plot or \code{set = "all"}
#' @param xmax maximum distance to show in the graph; if \code{xmax = NULL}, the
#'  maximum distance is taken as the max distance in the data set to which the
#'  models were fit.
#' @param resolution The number of line segments to break the curves into when
#'  plotting (i.e., \code{x = seq(0, xmax, length.out = resolution)}). Higher
#'  resolutions give smoother-looking curves.
#' @param mod_highlight Character string giving the name of the model to
#'  highlight by plotting it last and with \code{lwd = 2}. If \code{NULL}, the
#'  curve associated with the lowest (best) AICc score is highlighted.
#' @param CL confidence level to show in a \code{\link[=ddFit]{dd}} plot (ignored 
#'  for \code{\link[=ddFit]{ddArray}} objects)
#' @param nsim Number of simulation reps to use for estimating confidence bounds
#'  for \code{\link[=ddFit]{dd}} plot (ignored for \code{\link[=ddFit]{ddArray}}
#'  objects)
#' @details \code{\link[=ddFit]{ddArray}} objects are plotted with lines in order 
#'  of decreasing AICc, so that the "better" models are closer to the top and 
#'  more prominent. The model with the lowest AICc ("best" model) is plotted 
#'  last with a heavier line than the others.
#'
#'  For \code{\link[=ddFit]{dd}} objects, the curve for the MLE of the parameters
#'  is plotted, along with a 100\code{CL}\% confidence bounds determined for 
#'  \code{nsim} simulation reps
#'
#'  The legend follows the ordering given by \code{\link{modelFilter}} with 
#'  the default sieve or, if \code{extent = "win"} by (1) delta AICc < 10, 
#'  (2) the absence of high-influence points, and (2) AICc. The best model 
#'  according to the filter is listed first, with a heavier line than the others; 
#'  the remaining distributions are listed in descending order, with the best 
#'  models in the leftmost column. 
#' @return Plot displayed; no return value.
#' @rdname Plot
#' @export
#'
plot.ddArray = function(x, type = "CDF", extent = "full", distr = "all",
    xmax = NULL, resolution = 250, mod_highlight = NULL, ...){
  if (identical(distr, "all")) distr <- names(x)
  if (extent == "full"){ # only plot the extensible functions (sieve = NULL)
    tmp <- modelFilter(x[distr], sieve = "default")    
    distr <- rownames(tmp$scores)[tmp$scores[, "extensible"] == TRUE]
  } else if (extent == "win"){
    tmp <- modelFilter(x, sieve = "win")
    distr <- rownames(tmp$scores)
  } else {
    stop("plot.ddArray: 'extent' must be specied as \"full\" or \"win\" in arg list")
  }
  if (type == "rcd")
    distr <- distr[which(!distr %in% c("exponential", "tnormal", "constant"))]
  if (!any(distr %in% mod_all)) stop("plot.ddArray: some model(s) undefined")
  dd <- x[distr]
  ndistr <- length(distr) # number of distributions to plot
#  aic_proper <- aic(dd, extent = extent) # no need for AIC? ordering by filter
  if (!type %in% c("CDF", "PDF", "rcd"))
    stop("type (", deparse(substitute(type)),
          ") must be \"PDF\", \"CDF\", or \"rcd\"")
  if (!is.null(mod_highlight) && !mod_highlight %in% distr){
    message("'mod_highlight' not among models to be graphed. Using modelFilter ",
            "to select highlighted model.")
    mod_highlight <- NULL
  }
  if (is.null(mod_highlight)){
    mod_best <- tmp$filtered$distr
    if (mod_best %in% mod_all){
      mod_highlight <- mod_best
    } else {
      mod_highlight <- distr[1]
    }
  }

    
  ### graph is in two parts, in succession with fig...) and par(new = T)
  ###  0) preliminaries
  ###  1) legend [common to all parts]
  ###  2) main graph [lines]
  # extracting and parsing parameters from ... args
  arglist <- list(...)
  if (!"col" %in% names(arglist)){
    arglist$col <- mod_color
  } else {
    if (is.null(names(arglist$col)) && length(arglist$col) != ndistr &&
      length(arglist$col) != 1){
        stop(
          "'col' in plot.ddArray must be a scalar, a vector with one color for ",
          "each model plotted, or a vector with named elements matching the ",
          "names of the models to be plotted."
        )
    } else if (length(arglist$col) == 1){
      arglist$col <- rep(arglist$col, length.out = ndistr)
      names(arglist$col) <- distr
    } else if (!is.null(names(arglist$col))){
      if (!all(distr %in% names(arglist$col))){
        stop(
          "in plot.ddArray, 'col' be NULL or include a color for each plotted model"
        )
      }
    } else if (length(arglist$col) == ndistr){
      names(arglist$col) <- distr
    }
  }
  if (!"lty" %in% names(arglist)){
    arglist$lty <- mod_lty[distr]
    names(arglist$lty) <- distr
  } else {
    if (!is.vector(arglist$lty))
      stop("'lty' in plot.ddArray must be a scalar or vector")
    if (is.null(names(arglist$lty)) && length(arglist$lty) >= ndistr){
        stop(
          "'lty' in plot.ddArray must be a scalar or vector shorter than ",
          "number of plotted models. Values are recycled to fill out the vector ",
          "of line types to match the number of models to be plotted. Enter ",
          "?plot.ddArray for more info."
        )
    } else {
      arglist$lty <- rep(arglist$lty, length.out = distr)
      names(arglist$lty) <- distr
    }
  }
  if (!"xlim" %in% names(arglist)){
    if (is.null(xmax) || !is.finite(xmax)) xmax <- max(dd[[1]]$glm$data[, "r"])
    xmin <- 0
    arglist$xlim = c(ifelse(type == "rcd", 1, 0), xmax)
  }
  arglist$x = 0
  arglist$type = "n"
  ## part 1: legend
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))            
  do.call(par, par_default) 
  ncol <- 1 + floor((ndistr - 1)/5)
  sz <- 0.17
  par(fig = c(0, 1, 0, sz), mar = c(0, 1, 1.2, 0), family = "sans")
  plot(0, type = "n", axes = F, xlab = "", ylab = "")
  lwd <- numeric(ndistr) + 1
  names(lwd) <- distr
  if(mod_highlight %in% distr) lwd[mod_highlight] <- 2

  legend(x = "topleft", legend = distr_names[distr], lwd = lwd, 
    lty = arglist$lty, col = arglist$col[distr], bty = "n", x.intersp = 0.65,
    inset = -0.02, ncol = ncol, cex = 0.8)

  ## part 2a: general plot layot
  xseq <- seq(min(arglist$xlim), max(arglist$xlim), length = resolution)
  # size of the main graph; vertical, relative to plot window size:
  par(fig = c(0, 1, sz, 1), mar = c(4, 4, 2, 0.5), family = "sans", new = TRUE)
  ## part 2b: lines
  if (type == "CDF"){
    if (!"xlab" %in% names(arglist))
      arglist$xlab = "Distance from turbine (meters)"
    if (!"ylab" %in% names(arglist))
      arglist$ylab = "Predicted fraction of carcasses within the given distance"
    if (!"ylim" %in% names(arglist)) arglist$ylim <- 0:1
    do.call(plot, arglist)
    for (fi in distr[ndistr:1]){
      lines(xseq, pdd(xseq, model = dd[fi], extent = extent),
        lty = arglist$lty[fi], col = arglist$col[fi])
    }
    if (mod_highlight %in% distr){
      lines(xseq, pdd(xseq, model = dd[mod_highlight], extent = extent),
        col = arglist$col[mod_highlight], lwd = 2, lty = arglist$lty[mod_highlight])
    }
    if (extent == "win")
      mtext(side = 1, line = -1, adj = 1, "within search radius", family = "serif")
  } else if (type == "PDF"){
    if (!"ylim" %in% names(arglist)){
          ymax = -Inf
      for (fi in distr){
        ymax <- max(ymax, max(ddd(xseq, dd[fi], extent = extent)))
      }
      arglist$ylim <- c(0, ymax)
    }
    arglist$xlab = ifelse("xlab" %in% names(arglist), arglist$xlab, "Distance from Turbine (meters)")
    arglist$ylab = ifelse("ylab" %in% names(arglist), arglist$ylab, "Probability Density (PDF)")
    do.call(plot, arglist)
    for (fi in distr[ndistr:1]){
      lines(xseq, ddd(xseq, dd[fi], extent = extent),
        col = arglist$col[fi], lty = arglist$lty[fi])
    }
    if (mod_highlight %in% distr){
      lines(xseq, ddd(xseq, dd[mod_highlight], extent = extent),
        col = arglist$col[mod_highlight], lwd = 2, lty = arglist$lty[mod_highlight])
    }
    if (extent == "win")
      mtext(side = 3, adj = 1, "within search radius", family = "serif", line = -1)
   } else if (type == "rcd"){
    if (!"ylim" %in% names(arglist)){
      ymax = -Inf
      for (fi in distr){
        ymax <- max(ymax, max(rcd(xseq, dd[fi], extent = extent)))
      }
      arglist$ylim <- c(0, ymax)
    }
    arglist$xlab = ifelse("xlab" %in% names(arglist), arglist$xlab, "Distance from Turbine (meters)")
    arglist$ylab = ifelse("ylab" %in% names(arglist),
      arglist$ylab, "Relative Carcass Density")
    do.call(plot, arglist)
    for (fi in distr[ndistr:1]){
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

#' @rdname Plot
#' @export
#'
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
  if (is.null(xmax)) xmax <- max(x$glm$data$r)
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
  polygon(CI[c(1:nrow(CI), nrow(CI):1), "x"], c(CI[, 2], CI[nrow(CI):1, 3]),
    col = colors()[350], border = NA)
  msg <- switch(extent, win = "limited to carcasses within search radius  ")
  if (type == "PDF"){
    lines(xseq, ddd(xseq, x, extent = extent), col = mod_color[x$distr], lwd = 2)
    mtext(side = 3, line = -1, adj = 1, msg, family = "serif", lty = mod_lty)
  } else {
    lines(xseq, pdd(xseq, model = x, extent = extent),
      col = mod_color[x$distr], lwd = 2)
    mtext(side = 1, line = -1, adj = 1, msg, family = "serif")
  }
  box()
}

#' @rdname Plot
#' @export
#' 
plot.fmod <- function(x, ...){
    plot(x[[1]], ...)
}

#' Simulation of Dispersion Parameters
#'
#' @param x object to simulate from
#' @param nsim number of simulation draws
#' @param extent simulate according to full distribution, including extrapolation
#'  beyond the search radius (\code{extent = "full"}); or restrict the
#'  distribution to the area within the search radius (\code{extent = "win"}).
#' @param ... ignored
#' @return array with simulated beta parameters from the glm model, and their
#'  conversion to distribution parameters
#' @export
#'
ddSim <- function(x, ...) UseMethod("ddSim", x)

#' @rdname ddSim
#' @export
#'
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
  attr(ans, "srad") <- max(dd$glm$data$r, na.rm = TRUE)
  class(ans) <- "ddSim"
  ans
}

#' Subset Simulated Dispersion Parameters while Preserving Attributes
#'
#' @param x object to subset
#' @param i,j row and column indices to subset
#' @param ... ignored
#' @details Subset the ddSim object as if it were a simple matrix or array
#' @return array with simulated beta parameters from the glm model, their
#'  conversion to distribution parameters. NOTE: subsetting to a column or a row
#'  returns a matrix rather than a vector. This simplifies the coding and makes
#'  it easier to maintain integrity of data structures, but behavior differs from
#'  what is done when subsetting standard R matrices and arrays to a single column
#'  or row. Also unlike with standard R arrays and matrices, the class structure
#'  and attributes are preserved upon subsetting.
#' @export
#'
"[.ddSim" <- function(x, i, j,...){
  y <- NextMethod("[")
  noclass <- FALSE
  if (!missing(j)){
    if (is.numeric(j)) j <- colnames(x)[j]
    if (!all(j %in% colnames(x))) stop("cannot subset ddSim object on given column")
    if (length(j) == 1) y <- matrix(y, ncol = 1, dimnames = list(NULL, j))
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

#' Subset a Set of Fitted Dispersion Models (\code{ddArray})
#'
#' @param x object to subset
#' @param distr vector of names or integer indices of distributions to extract 
#'  from fitted models
#' @details Subset the ddArray object as if it were a simple vector
#' @return list of selected models with crucial statistics
#' @export
#'
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

#' @name ddPrint
#' @title Print S3 Objects in \code{dwp} Package
#'
#' @description \code{\link[=ddFit]{dd}}, \code{\link[=ddFit]{ddArray}}, and 
#'  \code{\link[=modelFilter]{fmod}} objects are lists consisting of a great
#'  amount of data. Only a few of the elements are printed automatically. Other 
#'  elements of object \code{x} can be viewed and extracted as with other lists 
#'  in R, namely, by using the \code{x$element} or \code{x[[element]]} operator, 
#'  where \code{element} is the name of one of the elements of \code{x}, all of
#'  which can be viewed via \code{names(x)}.
#'
#' @param x a \code{\link[=ddFit]{ddArray}} or \code{\link[=ddFit]{ddArray}} object
#' @param ... ignored
#' @return no return value; output printed to the console
#' @rdname ddprint
#' @export
#'
print.dd <- function(x, ...){
  cat(paste0("Distribution: ", x$distr, "\n"))
  cat(paste0("Formula: ", deparse(x$glm$formula), "\n\n"))
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

#' @rdname ddprint
#' @export
print.ddArray <- function(x, ...){
  for (distr in names(x)) {
    print(x[[distr]])
    cat("****************************************\n")
  }
}

#' @rdname ddprint
#' @export
print.fmod <- function(x, ...){
  print(x$filtered)
  print(x$scores)
  if (!identical(x$note, "")) 
    print(x$note)
}
#' Calculate CI for CDF, PDF, or quantile
#'
#' Calculate a confidence interval for the CDF, PDF, or quantile of a carcass
#' distance distribution.
#'
#' @param mod a \code{\link[=ddFit]{dd}} object
#' @param x distance from turbine (scalar or vector) or probability (for quantile)
#' @param type \code{"CDF"}, \code{"PDF"}, or \code{"quantile"}
#' @param CL confidence level for the confidence interval(s)
#' @param nsim number of simulation draws to base the estimate of \code{CI} on
#' @param extent whether to calculate \code{CI} based on the full range of
#'  possible data and extrapolating beyond the search radius
#'  (\code{extent = "full"}) or restricting the distribution to the area within
#'  the search radius (\code{extent = "win"}).
#' @param zrad an ad hoc radius to integrate to when the (uncommon) simulated
#'  parameter estimates do not result in an extensible distribution. In effect,
#'  This replaces NAs with 1s in CDFs and with 0s in PDFs.
#' @param na.tol maximum fraction of invalid parameter sets to discard when
#'  constructing CIs; abort if \code{mean(mod[, "extensible"]) > na.tol}
#' @return array (\code{ddCI} class) with columns for distance and the CI bounds
#' @export
ddCI <- function(mod, x, type = "CDF", CL = 0.9, nsim = 1000,
    extent = "full", zrad = 200, na.tol = 0.1){
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
      warning(paste0("cannot calculate meaningful CI for ", mod$distr, ": ",
        100*round(pNA, 3), "% NAs in CIs"))
#      CI[, 2:3] <- NA
    }
  }
  class(CI) <- "ddCI"
  CI
}

#' Convert GLM Coefficients into Named Distribution Parameters
#'
#' @param x object (vector or matrix of parameters, dd, or glm)
#'  with named glm parameters (\code{"r"}, \code{"I(r^2)"}, \code{"I(r^3)"}, 
#'  \code{"log(r)"}, or \code{"I(1/r)"}). NOTE: This function has minimal 
#'  error-checking. 
#'
#' @param distr name of the distribution
#' @param ... ignored
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
        rate = -x[, "r"]
      )
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
cof2parms.numeric <- function(x, distr, ...){
  output <- cof2parms(matrix(x, nrow = 1, dimnames = list(NULL, names(x))), distr)
  return(output)
}

#' @rdname cof2parms
#' @export
cof2parms.dd <- function(x, ...){
  return(cof2parms(x$coef, x$distr))
}

#' Calculate Probability Functions for Distance Distributions
#'
#' Calculate the standard d/p/q/r family of R probability functions for distance
#' distributions (\code{\link[=ddFit]{dd}}) as well as the relative carcass 
#' density (\code{rcd}). Usage broadly parallels that of the d/p/q/r probability 
#' functions like \code{\link[=stats]{dnorm}}, \code{pnorm}, \code{qnorm}, and 
#' \code{rnorm}. 
#' 
#' The probability density function (PDF(\emph{x}) = \emph{f}(\emph{x}) = \code{ddd(x, ...)}) 
#' gives the probability that a carcass falls in a 1 meter ring centered at the 
#' turbine and with an outer radius of \code{x} meters. The cumulative distribution 
#' function [CDF(\emph{x}) = \emph{F}(\emph{x}) = \code{pdd(x, ...)}] gives the 
#' probability that a carcass falls within \code{x} meters from the turbine. For 
#' a given probability, \code{p}, the inverse CDF [\code{qdd(p,...)}] gives the 
#' \code{p} quantile of carcass distances. For example, \code{qdd(0.5,...)} 
#' gives the median carcass distance, and \code{qdd(0.9, ...)} gives the radius 
#' that 90\% of the carcasses are expected to fall in. Random carcass distances 
#' can be generated using \code{rdd}.
#' 
#' The relative carcass density function(\code{rcd}) gives relative carcass 
#' densities at a point \code{x} meters from a turbine. In general, rcd is 
#' proportional to PDF(x)/x, normalized so that the surface of rotation of rcd(x) 
#' has total volume of 1. There are more stringent contstraints on the allowable 
#' parameters in the fitted (or simulated) glm's because the integral of PDF(x)/x 
#' must converge.
#' 
#' Distributions may be extrapolated beyond the search radius to account for all
#' carcasses, including those that land beyond the search radius 
#' (\code{extent = "full"}), or may be restricted to carcasses falling within the
#' searched area (\code{extent = "win"}). Typically, in estimating \code{dwp} for
#' a fatality estimator like \code{eoa} or \code{GenEst}, the full distributions
#' would be used.
#'    
#' The probability functions have a number of purposes. A few of the more commonly
#' used are listed below.
#' \describe{
#'  \item{PDF and CDF (\code{ddd} and \code{pdd}):}{
#'   \itemize{
#'    \item{to calculate the probability that carcass lands at a distance
#'     \code{x} meters from the turbine (or, more precisely, within 0.5 meters of 
#'     \code{x}) or within \code{x} meters from the turbine, use a scalar value
#'     of \code{x} and a single model (\code{\link[=ddFit]{dd}} or \code{\link{ddSim}})
#'     with \code{ddd} or \code{pdd}, repspectively;}
#'    \item{to account for uncertainty in the probabilities at \code{x}, use 
#'     \code{ddd} or \code{pdd} for with scalar \code{x} and a simulated set of
#'     parameters from the fitted model (\code{\link{ddSim}} object). This would 
#'     be useful for calculating confidence intervals for the probabilities;}
#'    \item{to calculate probabilities for a range of \code{x} values according 
#'     to a single model, use a vector \code{x} with a \code{\link[=ddFit]{dd}} object or 
#'     a \code{\link{ddSim}} object with one row. This would be useful for 
#'     drawing graphs of PDFs or CDFs;}
#'    \item{to calculate simulated probabilites for a range of \code{x} values,
#'     use a vector \code{x} and a \code{ddSim} object of simulated parameter sets.
#'     This would be useful for drawing confidence regions around a fitted PDF or
#'     CDF.}
#'  }}
#'  \item{Inverse CDF (\code{qdd}):}{
#'   \itemize{
#'    \item{to calculate the distance that 100\code{p}\% of the carcasses are 
#'     expected to fall, use a scalar \code{p} in the interval (0, 1) and a 
#'     single model (\code{\link[=ddFit]{dd}}) or parameter set (\code{\link{ddSim}} with
#'     one row);}
#'    \item{to calculate account for the uncertainty in estimating the inverse 
#'     CDF for a given \code{p}, use a scalar \code{p} and a \code{\link{ddSim}}
#'     object. This would be useful for calculating a confidence interval for,
#'     say, the median or the expected 90th percentile of carcass distances;}
#'    \item{to calculate the inverse CDF for a range of probabilities for a single
#'     model, use a vector \code{p} and a single model (\code{\link[=ddFit]{dd}} or 
#'     \code{\link{ddSim}} object with one row.}
#'   }  
#'  }
#'  \item{Random Carcasses Distances (\code{rdd}):}{
#'    \itemize{
#'     \item{to generate \code{n} random carcass distances for a given (fixed)
#'      model, use a \code{\link[=ddFit]{dd}} object or a \code{\link{ddSim}} object with
#'      a single row;}
#'     \item{to generate \code{n} random carcass distances for a model and account
#'      for the uncertainty in estimating the model, use a \code{\link{ddSim}} 
#'      object with \code{n} rows, where \code{n} is also used as the \code{n}
#'      argument in the call to \code{rdd}.}
#'    }
#'  }
#'  \item{Relative Carcass Density (per m^2):}{
#'    \itemize{
#'     \item{to calculate the relative carcass density at a number of distances,
#'      use a vector \code{x}. This would be useful in generating maps of carcass
#'      density at a site.}
#'    }
#'  }
#' }
#' @param x,q,p,n numeric, \eqn{x \ge 0}
#' @param model either a \code{dd} object or a \code{ddSim} object
#' @param parms model parameters; required if model is specified as a character 
#'  string rather than a \code{dd} or \code{ddSim} object (otherwise optional and ignored)
#' @param extent for a full distribution extrapolated beyond the search radius
#'  to account for all carcasses, use \code{extent = "full"}; for a distribution
#'  restricted solely to carcasses falling within the search radius, use
#'  \code{extent = "win"}.
#' @param zrad the distance at which carcass density is assumed to be zero; to
#'  be used only in simulation reps in which simulated parameters do not yield
#'  extensible distributions, essentially returning 0 rather than NA for those
#'  pathological cases.
#' @param subdiv if the number of values to calculate with \code{rdd} or \code{qdd}
#'  is >1, the function uses breaks the PDF into \code{subdiv} subdivisions and
#'  interpolates to solve the inverse. More subdivisions gives greater accuracy
#'  but is slower.
#' @param silent If \code{TRUE}, then console messages are suppressed.
#' @return vector or matrix of values; a vector is returned unless \code{model} 
#'  is a \code{ddSim} object with more than one row and is to be calculated for
#'  more than one value (\code{x}, \code{q}, \code{p}), in which case an array 
#'  with dimensions \code{length(x)} by \code{nrow(model)} is returned (where 
#'  "\code{x}" is \code{x}, \code{q}, or \code{p}, depending on whether \code{ddd},
#'  \code{pdd}, or \code{qdd} is called).
#' 
#' @examples
#' data(layout_simple)
#' data(carcass_simple)
#' sitedata <- initLayout(layout_simple)
#' ringdata <- prepRing(sitedata)
#' ringsWithCarcasses <- addCarcass(carcass_simple, data_ring = ringdata)
#' distanceModels <- ddFit(ringsWithCarcasses)
#' modelEvaluations <- modelFilter(distanceModels)
#' bestModel <- modelEvaluations$filtered
#' pdd(100, model = bestModel) # estimated fraction of carcasses within 100m
#' ddd(1:150, model = bestModel) # estimated PDF of the carcass distances
#' qdd(0.9, model = bestModel) # estimated 0.9 quantile of carcass distances
#' rdd(1000, model = bestModel) # 1000 random draws from estimated carcass distribution

#' @export
ddd <- function(x, model, parms = NULL, extent = "full", zrad = 200){ # model is either ddSim or dd
  if ("dd" %in% class(model)) {
    model <- dd2ddSim(model)
  } else if (identical("character", class(model))) {
    model <- try(mpp2ddSim(model, parms))
    if (identical("try-error", class(model))) stop("ddd: improper model format")
  }
  if ("ddSim" %in% class(model)) {
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
            b1 = c(parms[i, "b1"])
            b2 = c(parms[i, "b2"])
            b3 = c(parms[i, "b3"])
            const <- tryCatch(1/integrate(
              f = function(r)
                r * exp(b1 * r + b2 * r^2 + b3 * r^3),
              lower = 0, upper = Inf, rel.tol = .Machine$double.eps^0.5)$val,
              error = function(e) NA
            )
            tmp[, which(i1 == i)] <- dxep123(x,
              b1 = b1, b2 = b2, b3 = b3,
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
        if(distr == "xepi0" && !cofOK0(parms, distr)){
          output[, i] <- NA
          next
        }
        deno <- try(integrate(f = function(r)
          exp(rmat(r, distr) %*% parms[i, cof_name[[distr]]] + off(r, distr)),
          lower = 0,
          upper = zrad
        )$val)
        if ("try-error" %in% class(deno)){
          output[x > 0, i] <- NA
        } else {
          output[x > 0, i] <- exp(rmat(x[x > 0], distr) %*%
            parms[i, cof_name[[distr]]] + off(x[x > 0], distr))/deno
        }
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

#' @rdname ddd
#' @export
pdd <- function(q, model,  parms = NULL, extent = "full", zrad = 200, silent = FALSE){ 
  if ("dd" %in% class(model)) {
    model <- dd2ddSim(model)
  } else if (identical("character", class(model))) {
    model <- try(mpp2ddSim(model, parms))
    if (identical("try-error", class(model))) stop("pdd: improper model format")
  }
  if ("ddSim" %in% class(model)) {
    distr <- attr(model, which = "distr")
    srad <- attr(model, which = "srad")
  } else {
    stop("pdd: model must be dd or ddSim object")
  }
  parms <- model
  x <- q
  if (extent == "full"){ # intent to integrate to Inf; use zrad if parms improper
    # output has dimensions length(q) x nrow(parms) [or a vector]
    # calculations must be split into two parts:
    #  1) proper parms integrate to Inf
    #  2) improper parms are integrated to zrad
    output <- matrix(0, nrow = length(q), ncol = nrow(parms))
    i0 <- unname(which(parms[, "extensible"] == 0))
    i1 <- unname(which(parms[, "extensible"] > 0))
    if(length(i1) > 0){ # extensible parms
      xx <- rep(q, length(i1))
      ppi <- rep(i1, each = length(q))
      output[, i1] <- switch(distr,
        xep01 = matrix(
          pgamma(xx, shape = parms[ppi, "shape"], rate = parms[ppi, "rate"]),
          nrow = length(q)
        ),
        lognormal = matrix(
          plnorm(xx, meanlog = parms[ppi, "meanlog"], sdlog = parms[ppi, "sdlog"]),
          nrow = length(q)
        ),
        xep1 = matrix(pxep1(xx, b1 = parms[ppi, "b1"]),
          nrow = length(q)
        ),
        xep12 = matrix(pxep12(xx, b1 = parms[ppi, "b1"], b2 = parms[ppi, "b2"]),
          nrow = length(q)
        ),
        xep02 = matrix(pxep02(xx, b0 = parms[ppi, "b0"], b2 = parms[ppi, "b2"]),
          nrow = length(q)
        ),
        xep123 = {
          const <- numeric(length(i1))
          tmp <- matrix(0, nrow = length(q), ncol = length(i1))
          for (i in i1){
            const <- tryCatch(1/integrate(
              f = function(r)
                (exp(rmat(r, distr)[, -1] %*% c(parms[i, cof_name[[distr]]][-1]) + 
                  off(r, distr))),
              lower = 0, upper = Inf, rel.tol = .Machine$double.eps^0.5)$val,
              error = function(e) NA
            )
            tmp[, which(i1 == i)] <- pxep123(q,
              b1 = c(parms[i, "b1"]), b2 = c(parms[i, "b2"]), b3 = c(parms[i, "b3"]),
              const = const
            )
          }
          tmp
        },
        xepi0 = matrix(
          pxepi0(xx, shape = parms[ppi, "shape"], scale = parms[ppi, "scale"]),
          nrow = length(q)
        ),
        xep012 = { # assumes b0, b1, b2, b3 are scalars
          const <- numeric(length(i1))
          tmp <- matrix(0, nrow = length(q), ncol = length(i1))
          for (i in i1){
            b0 <- c(parms[i, "b0"])
            b1 <- c(parms[i, "b1"])
            b2 <- c(parms[i, "b2"])
            const <- tryCatch(1/integrate(f = function(x)
              x * exp(b0 * log(x) + b1 * x + b2 * x^2),
              lower = 0, upper = Inf, rel.tol = .Machine$double.eps^0.5)$val,
              error = function(e) NA
            )
            tmp[, which(i1 == i)] <- pxep012(q,
              b0 = b0,
              b1 = b1,
              b2 = b2,
              const = const
            )
          }
          tmp
        },
        xep0123 = { # assumes b0, b1, b2, b3 are scalars
          const <- numeric(length(i1))
          tmp <- matrix(0, nrow = length(q), ncol = length(i1))
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
            tmp[, which(i1 == i)] <- pxep0123(q,
              b0 = b0,
              b1 = b1,
              b2 = b2,
              b3 = b3,
              const = const
            )
          }
          tmp
        },
        xep2 = matrix(pxep2(xx, s2 = parms[ppi, "s2"]), nrow = length(q)),
        MaxwellBoltzmann = matrix(pmb(xx, a = parms[ppi, "a"]), nrow = length(q)),
        constant = {
          tmp <- (xx/zrad)^2
          tmp[xx > zrad] <- 1
          matrix(tmp, nrow = length(q))
        },
        xep0 = matrix(pxep0(xx, a = parms), nrow = length(q)),
        tnormal = {
          mu <- parms[ppi, "mean"]
          sig  <- parms[ppi, "sd"]
          Fa <- pnorm(0, mean = mu, sd = sig)
          ans <- numeric(length(xx))
          ans[xx >= 0] <- (pnorm(xx[xx >= 0], mean = mu, sd = sig) - Fa)/(1 - Fa)
          matrix(ans, nrow = length(q))
        },
        exponential = matrix(pexp(xx, rate = parms[ppi, "rate"]), nrow = length(q)),
        inverse_gaussian = matrix(statmod::pinvgauss(xx,
            mean = parms[ppi, "mean"], dispersion = parms[ppi, "dispersion"]),
          nrow = length(q)
        ),
        chisq = matrix(pchisq(xx, df = parms[ppi, "df"]), nrow = length(q))
      )
    }
    if (length(i0) > 0){
      if (is.null(zrad) || is.na(zrad)){
        output[, i0] <- NA
      } else {
        for(i in i0){
          if(distr == "constant") output[, i] <- (q/zrad)^2
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
          if (any(q > 0) & !is.na(deno)){
            for (xi in which(q > 0 & q < zrad)){
              output[xi, i] <- tryCatch(
                integrate(function(r)
                  exp(rmat(r, distr) %*% parms[i, cof_name[[distr]]] + off(r, distr)),
                  lower = 0,
                  upper = q[xi], rel.tol = .Machine$double.eps^0.5
                )$val/deno,
                error = function(e) NA
              )
            }
          }
          
        }
        if(!silent){
          if (length(i0) == 1){
            msg <- paste0(
              "pdd: model in row ", i0, " is not extensible. \nResult for this model ",
              "is based on a truncated model with max radius of ", zrad, " meters.\n"
            )
          } else {
            msg <- paste0(
              "pdd: ", length(i0), " of the models are not extensible. \n",
              "results for these models are based on truncation of the models ",
              "at a maximum radius of ", zrad, " meters.\n"
            )
          }
          msg <- paste0(msg, "Use zrad = NULL to return NA rather than returning ",
            "results based on model truncation.")
          message(msg)
        }
      }
      output[q >= zrad, i0] <- 1
    }
  } else if (extent == "win"){ # truncation radius is numeric scalar 0 < srad < Inf
    srad <- attr(model, "srad")
    output <- matrix(0, nrow = length(q), ncol = nrow(parms))
    if (all(q <= 0)){
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
        for (xi in which(q > ifelse(distr == "xep0", 1, 0) & q <= srad)){
          output[xi, i] <- integrate(function(r)
            exp(rmat(r, distr) %*% parms[i, cof_name[[distr]]] + off(r, distr)),
            lower = ifelse(distr == "xep0", 1, 0),
            upper = q[xi]
          )$val/deno
        }
      } else {
        output[, i] <- NA
      }
    }
    output[q > srad, ] <- 1
  } else {
    stop("extent must be \"full\" or \"win\"")
  }
  if (nrow(parms) == 1 | length(q) == 1) output <- unname(as.vector(output))
  output[output > 1] <- 1
  output[output < 0] <- 0
  output
}

#' @rdname ddd
#' @export
qdd <- function(p, model, parms = NULL, extent = "full", zrad = 200, subdiv = 1000){ # model is ddSim or dd
  # find r such that pdd(r) = p
  if ("dd" %in% class(model)) {
    model <- dd2ddSim(model)
  } else if (identical("character", class(model))) {
    model <- try(mpp2ddSim(model, parms))
    if (identical("try-error", class(model))) stop("qdd: improper model format")
  }
  if ("ddSim" %in% class(model)) {
    distr <- attr(model, which = "distr")
    srad <- attr(model, which = "srad")
  } else {
    stop("pdd: model must be dd or ddSim object")
  }
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
      x <- pdd(y, model = model, extent = extent)
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

#' @rdname ddd
#' @export
rdd <- function(n, model, parms = NULL, extent = "full", zrad = 200, subdiv = 1000){ # model is ddSim or dd
  if ("dd" %in% class(model)) {
    model <- dd2ddSim(model)
  } else if (identical("character", class(model))) {
    model <- try(mpp2ddSim(model, parms))
    if (identical("try-error", class(model))) stop("qdd: improper model format")
  }
  if ("ddSim" %in% class(model)){
    if (nrow(model) > 1 & nrow(model) != n)
      stop("rdd: nrow(model) must be 1 or n")
  } else {
    stop ("class(model) in rdd must be dd or ddSim")
  }

  if (! "ddSim" %in% class(model)) stop("pdd: model must be dd or ddSim object")
  qdd(runif(n), model = model, extent = extent, zrad = zrad, subdiv = subdiv)
}

#' @rdname ddd
#' @export
rcd <- function(x, model, parms = NULL, extent = "full", zrad = 200){
  x[x < 0.001] <- 0.001
  if ("dd" %in% class(model)) {
    model <- dd2ddSim(model)
  } else if (identical("character", class(model))) {
    model <- try(mpp2ddSim(model, parms))
    if (identical("try-error", class(model))) stop("qdd: improper model format")
  }
  if ("ddSim" %in% class(model)){
    distr <- attr(model, which = "distr")
    srad <- attr(model, which = "srad")
  } else {
    stop("rcd: model must be dd or ddSim object")
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

#' Display a Tables of Summary Statistics for Distance Distributions
#' 
#' Calculate summary statistics for a single distance distribution 
#' (\code{\link[=ddFit]{dd}} object), an array of distance distributions all fit to the
#' same data set (\code{\link[=ddFit]{ddArray}}), or a list of arrays of distance 
#' distributions fit for different carcass classes but the same site layout
#' (\code{\link[=ddFit]{ddArraycc}}).
#'
#' @param x list of models (\code{ddArray}) or single model (\code{dd}) to calculate summary statistics for
#' @param extent distributions within searched area (\code{"win"}) or extended beyond (\code{"full"})
#' @param zrad maximum distance that carcasses can lie
#'  (only used when glm parameters not extensible to Inf)
#' @param ... ignored
#' @return list (or list of lists if \code{x} is \code{ddArray}) with \code{$model}
#'  giving the model parameters and \code{$stats} giving the median, and 75th,
#'  90th, and 95th quantiles of carcass distances and the estimated probability
#'  a carcass falls within the search area according to each model
#' @export
stats <- function(x, ...) UseMethod("stats", x)

#' @rdname stats
#' @export
stats.dd <- function(x, extent = "full", zrad = 200, ...){
  distr <- x$distr
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
  xseq <- seq(0.1, min(ans[["stats"]]["90%"], 1000), by = 0.1)
  y <- ddd(x = xseq, model = x, extent = extent, zrad = zrad)
  mi <- which(y == max(y))
  ans[["stats"]]["mode"] <- ifelse(y[length(y)] %in% mi, NA, mean(xseq[mi]))

  ans[["stats"]]["p_win"] <- pdd(max(x$glm$data$r), x) * 100
  if (ans[["stats"]]["mode"] == 0.1) ans[["stats"]]["mode"] <- 0
  if (ans[["stats"]]["mode"] > ans[["stats"]]["95%"])
    ans[["stats"]]["mode"] <- NA
  ans[["stats"]] <- round(ans[["stats"]], 1)
  ans[["stats"]]["p_win"] %<>% `*`(., 0.01)
  ans
}

#' @rdname stats
#' @export
stats.ddArray <- function(x, extent = "full", zrad = 200, ...){
  aic0 <- aic(x, extent = extent)
  rnm <- rownames(aic0)
  cnm<- c("median", "75%", "90%", "95%", "mode", "p_win", "deltaAICc")
  ans <- data.frame(array(dim = c(length(rnm), length(cnm))))
  rownames(ans) <- rnm
  colnames(ans) <- cnm
  ans$deltaAICc <- aic0$deltaAICc
  for (distr in rownames(aic0))
    ans[distr, c("median", "75%", "90%", "95%", "mode", "p_win")] <-
      stats(x[distr], extent = extent, zrad = zrad)$stats
  ans
}

#' @rdname stats
#' @export
stats.ddArraycc <- function(x, extent = "full", zrad = 200, ...){
  lapply(x, stats, extent = extent, zrad = zrad)
}

#' Extract Parameters from a Distance Model (\code{dd}) and Format as \code{ddSim} Object
#'
#' This is a utility function called internally by \code{dwp} functions
#' to extract parameters from a fitted \code{\link[=ddFit]{dd}} model and 
#' formats them for analysis and calculation.
#'
#' @param dd \code{dd} object
#' @return \code{ddSim} object with 1 row
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
  attr(parms, "srad") <- ceiling(max(dd$glm$data$r, na.rm = TRUE))
  attr(parms, "distr") <- dd$distr
  class(parms) <- "ddSim"
  parms
}


#' Estimate Probability Carcass lands in Searched Area
#'
#' @description Estimated probability that carcass lands in searched area. This
#'  is an intermediate step in estimating dwp but is also interesting in its own
#'  right. The estimation involves integrating the modeled carcass distribution 
#'  (\code{model}) over the  search plots at the turbines. Data for the search 
#'  plots is stored in the generic argument, \code{x}, which can take any of a
#'  number of different forms, as described in the \code{Arguments} section (below).
#' @param x data \describe{
#'  \item{\code{rings}}{a formatted site map created from raw data via function
#'    \code{\link{prepRing}} (or as a component of a list returned by 
#'    \code{\link{addCarcass}}).}
#'  \item{\code{ringscc}}{a list of \code{rings} objects, one for each carcass
#'    class; created from raw data via function \code{\link{prepRing}} (or as a 
#'    component of a list returned by \code{\link{addCarcass}}).}
#'  \item{\code{xyLayout}}{formatted site map data derived from (x, y) coordinates
#'    covering every square meter of searched areas at each turbine; derived from
#'    the function \code{\link{initLayout}}, when called with \code{xy} data.}
#'  \item{\code{rpA}}{(intended as an internal function that would rarely be 
#'    called directly by users) a list of data frames (one for each turbine) 
#'    giving the fraction of area searched (\code{pinc} at  each distance 
#'    (\code{r}). \code{rpA} data are embedded in \code{\link[=prepRing]{rings}} 
#'    objects that are created from site "maps" via \code{\link{prepRing}}.}
#'  \item{\code{rdat}}{(intended as an internal function that would rarely be 
#'    called directly by users) list of data frames giving the area searched 
#'    (\code{"exposure"}), in a 1 meter ring with outer radius \code{"r"} and the 
#'    number of carcasses found \code{"ncarc"} in each ring, with search class 
#'    \code{scVar} optional. There is also a summary data frame 
#'    \code{$rdat[["total"]]} that sums the exposures and carcass counts for all 
#'    turbines across the site. The \code{$rdat[["total"]]} is the data frame 
#'    used in fitting the GLMs. \code{rdat} objects are components of the return
#'    value of \code{\link{prepRing}}}
#'  \item{\code{data.frame}}{(intended as an internal function that would rarely be 
#'    called directly by users) a data frame giving the fraction of area searched 
#'    (\code{pinc} at  each distance (\code{r}).}
#' }
#' 
#' @param model A fitted \code{dd} model or an array of estimated parameters
#'  (\code{ddSim} object); or, if \code{x} is a \code{ringscc} object, a list
#'  of \code{dd} models (one for each carcass class), or a \code{ddArraycc}
#'  accompanied by a vector of model names to use (one for each carcass class).
#' @param modnames if \code{x} is a \code{ringscc} object, a vector of names of 
#' model to use for each carcass class; otherwise, \code{modnames} is ignored.
#' @param extent calculate dwp within searched radius only (\code{"win"}) or
#'  for full complement of carcasses (\code{"full"}), including those that fall
#'  outside the search radius.
#' @param nsim number of parametric bootstrap iterations for accounting for
#'  uncertainty in the estimator. Default is \code{nsim = 1000}. Use 
#'  \code{nsim = 0} for the estimate of \code{psi} based on the  MLE of
#'  the given model without accounting for uncertainty. 
#' @param zrad radius
#' @param ... ignored
#' @return A \code{psiHat} object, which is either 1) an array giving the 
#'  expected fraction  of carcasses lying in the searched area at each turbine 
#'  with \code{nsim} rows and one column for each turbine + one row for the 
#'  total; or 2) a list of such arrays, one for each carcass class if \code{x} 
#'  is a \code{ringscc} object. The uncertainty in the expected fractions is
#'  characterized by simulation and reflected in the variation in \code{psi}
#'  values within each column.
#' @export
estpsi <- function(x, ...) UseMethod("estpsi", x)

#' @rdname estpsi
#' @export
estpsi.rings <- function(x, model, extent = "full", nsim = 1000, zrad = 200, ...){
  estpsi(x$rpA, model = model, extent = extent, nsim = nsim, zrad = zrad)
}

#' @rdname estpsi
#' @export
estpsi.ringscc <- function(x, model, modnames = NULL, extent = "full",
                           nsim = 1000, zrad = 200, ...){
  output <- list()
  arglist <- as.list(match.call())
  arglist[[1]] <- NULL
  arglist[["modnames"]] <- NULL
  if ("ddArraycc" %in% class(model) && any(sapply(model, length) > 1)){
    if (is.null(modnames)){
      stop("estpsi: modnames must be provided if list of ",
           "models includes more than one model for any size class."
      )
    }
  }
  for (cc in names(x)){
    # if model is ddArraycc:
    if (!is.null(modnames)){ # then use modnames to select model
      arglist[["model"]] <- model[[cc]][[modnames[cc]]]
    } else { # then there is only one model for the given class
      arglist[["model"]] <- model[[cc]]
    }
    arglist[["x"]] <- x[[cc]][["rpA"]]
    output[[cc]] <- do.call(estpsi, arglist)
  }
  class(output) <- "psiHatcc"
  output
}

#' @rdname estpsi
#' @export
estpsi.xyLayout <- function(x, model, extent = "full", nsim = 1000, zrad = 200,
    ...){
  parmsim <- prepmod(model = model, nsim = nsim)
  di <- attr(parmsim, "distr")
  ep <- function(r, di, pco) c(exp(rmat(r, di) %*% pco[1, cof_name[[di]]]))
  tset <- x$tset
  output <- matrix(NA, nrow = nsim, ncol = length(tset) + 1)
  colnames(output) <- c(tset, "total")
  for (simi in 1:nsim){
    pco <- parmsim[simi, ]
    if (!cofOK0(pco, di)) next
    deno <- try(integrate(f = function(r, di, pco) 2 * pi * r * ep(r, di, pco),
      lower = 0, upper = ifelse(cofOKInf(pco, di), Inf, zrad),
      di = di, pco = pco)$val, silent = TRUE)
    if ("try-error" %in% class(deno)){
      deno <- try(integrate(f = function(r, di, pco) 2 * pi * r * ep(r, di, pco),
        lower = 0, upper = zrad, di = di, pco = pco)$val, silent = TRUE)
      if ("try-error" %in% class(deno)) {
        warning(
          "Difficulty integrating with simulated parameter set. ",
          "Discarding the offending rep."
        )
      }
    }
    tmp <- aggregate(ep(x$xydat$r, di, pco) ~ x$xydat$turbine, FUN = "sum")
    output[simi, tmp[, 1]] <- tmp[, 2]/deno
  }
  output[output > 1] <- 1
  output[, "total"] <- rowMeans(output[, exclude("total", colnames(output)), drop = F], na.rm = TRUE)
  attr(output, "extent") <- extent
  attr(output, "zrad") <- zrad
  class(output) <- c("psiHat", "matrix")
  return(output)
}

#' @rdname estpsi
#' @export
estpsi.rpA <- function(x, model, extent = "full", nsim = 1000, zrad = 200, ...){
  parmsim <- prepmod(model, nsim)
  tmp <- lapply(x, FUN = estpsi, model = parmsim,
    extent = extent, nsim = nsim, zrad = zrad)
  output <- sapply(tmp, "[[", "psi")
  attr(output, "extent") <- extent
  attr(output, "zrad") <- zrad
  class(output) <- c("psiHat", "matrix")
  output
}


#' @rdname estpsi
#' @export
estpsi.rdat <- function(x, model, extent = "full", nsim = 1000, zrad = 200, ...){
  rpA <- lapply(x, FUN = function(x){
    x[, "exposure"] <- x[, "exposure"]/((x[, "r"] - 0.5) * 2 * pi)
    x
  })
  rpA$total$exposure <- rpA$total$exposure/(length(x) - 1)
  rpA <- lapply(rpA, 
    FUN = function(x){
      colnames(x) <- gsub("exposure", "pinc", colnames(x))
      x[, "pinc"] <- round(pmin(x[, "pinc"], 1), 4)
      x[, "ncarc"] <- NULL
    }
  )
  class(rpA) <- "rpA"
  estpsi(rpA, model = model, extent = extent, nsim = nsim, zrad = zrad)
}

#' @rdname estpsi
#' @export
estpsi.data.frame <- function(x, model, extent = "full", nsim = 1000, zrad = 200,
    ...){
  parmsim <- prepmod(model, nsim)
  output <- list(
    psi = c(t(ddd(x$r - 1/2, parmsim, extent = extent, zrad = zrad)) %*% x$pinc),
    ncarc = sum(x$ncarc)
  )
  output$psi[output$psi > 1] <- 1
  attr(output, "extent") <- extent
  attr(output, "zrad") <- zrad
  class(output) <- c("psiHat", "matrix")
  output
}


#'  Internal Utility Function to Parse and Format Model for Calculating Psihat
#'
#' @param model \code{dd} or \code{ddSim} object
#' @param nsim number of simulation reps. If \code{nsim = 0}, return dd2ddSim
#' @return \code{\link{ddSim}} object of simulated distribution model parameters
#' @export
prepmod <- function(model, nsim){
  if ("dd" %in% class(model)){
    if (is.null(nsim) || is.na(nsim) || nsim == 0){
      parmsim <- dd2ddSim(model)
    } else {
      parmsim <- ddSim(model, nsim = nsim)
    }
  } else if ("ddSim" %in% class(model)){
    parmsim <- model
  } else {
    stop("estpsi: model must be dd or ddSim object")
  }
  parmsim
}

#' Estimate DWP
#'
#' Estimate the density-weighted proportion (DWP) of carcasses lying in the
#'  searched area at each turbine at a site. The calculation requires prior
#'  estimation of the expected proportion (\code{\link[=estpsi]{psi}}) and the
#'  number of carcasses found (\code{\link[=getncarc]{ncarc}}). NOTE: The 
#'  carcass counts affect the uncertainty in the estimate of the fraction
#'  of carcasses in the searched area (DWP), and \code{ncarc} is required for
#'  accounting for uncertainty in estimates of DWP. 
#'
#' @param x Either (1) \code{\link[=estpsi]{psiHat}} object, which is an 
#'  \code{nsim} by \code{nturbine} matrix that gives the estimated probability 
#'  of that a given carcass will land in the searched area at each turbine, with 
#'  turbine IDs as column names; or (2) a \code{\link[=estpsi]{psiHatcc}} 
#'  object, which is a list of \code{psiHat} objects, one for each carcass class.
#' @param ncarc vector of total carcass count at each turbine represented in x.
#' @param nboot number of parametric bootstrap iterations for estimating CIs
#' @param forGenEst format the results for importing into GenEst (boolean)
#' @param silent suppress messages from the fitting of a beta distribution in
#'  internal calculations that, if successful, increase the speed of the 
#'  calculations by 20-200x. The message would signal that this acceleration 
#'  cannot be applied.
#' @param ... ignored
#' @return list
#' @export
estdwp <- function(x, ...) UseMethod("estdwp", x)

#' @rdname estdwp
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
  if (length(psi) == length(ncarc)){ # nsim = 1
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
      output$psi[, ti] <- findInterval(qtls, cumsum(postM(ncarc[ti], psi[ti])))
    }
    output$ncarc <- ncarc
  } else if (is.null(nboot) || length(psi)/length(ncarc) == nboot){
    nboot <- ifelse(is.null(nboot), length(psi)/length(ncarc), nboot)
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
      dotog <- TRUE
      if (!any(is.na(psi[, ti]))){ # try doing them as betabinom
        muB <- mean(psi[, ti])
        sig2B <- var(psi[, ti])
        Ba <- muB^2/sig2B*(1 - muB) - muB
        Bb <- Ba*(1/muB - 1)
        ab <- suppressWarnings(try(
          MASS::fitdistr(x = psi[, ti], densfun = "beta",
            start = list(shape1 = Ba, shape2 = Bb)),
          silent = T
        ))
        if (!"try-error" %in% class(ab)){ # then betabinom works
          pm <- postM.ab(x = ncarc[ti], Ba = Ba, Bb = Bb, prior = "IbetabinRef")
          if (any(is.na(pm)) || mean(pm == 0) == 1){
            dotog <- FALSE
          } else {
            output[, ti] <- ncarc[ti]/sample(1:length(pm) - 1, size = nrow(output), 
              prob = pm, replace = TRUE)
          }
        } else { # betabinom doesn't work, so do each separately
          dotog <- FALSE
        }
      } else {
        dotog <- FALSE
      }
      if (!dotog){
        for (bi in 1:nboot){
          pm <- suppressWarnings(postM(ncarc[ti], psi[bi, ti]))
          if (!any(is.na(pm)) && sum(pm > 0) > 0){
            output[bi, ti] <- ncarc[ti]/sample(x = 1:length(pm) - 1, size = 1,
              prob = pm)
          }
        }
      }
    }
  } else if (nrow(psi) != nboot){
    tname <- names(ncarc)
    output <- array(0, dim = dim(psi), dimnames = list(NULL, colnames(psi)))
    for (ti in tname){
      for (bi in 1:nboot){
        pm <- postM(ncarc[ti], x$psi[bi, ti])
        if (!any(is.na(pm)) && sum(pm > 0) > 0){
          output[bi, ti] <- ncarc[ti]/sample(x = 1:length(pm) - 1, size = 1,
            prob = pm)
        }
      }
    }
  }
  output <- round(output, 3)
  output[output > 1] <- 1
  output[is.na(output)] <- psi[is.na(output)] # return psi when ncarc = 0 for a turbine
  if (ncol(output) == 1){ # single turbine, nsim reps
    if (forGenEst){
      output <- data.frame(turbine = "all", dwp = as.vector(output))
      class(output) <- c("dwphat", "GenEst", "data.frame")
      return(output)
    }
  }
  if (forGenEst){
   class(output) <- c("dwphat", "matrix")
   output <- formatGenEst(output)
  } else {
   class(output) <- c("dwphat", "notGenEst", "matrix")
  }
  return(output)
}

#' @rdname estdwp
#' @export
#x = psi_size; ncarc = getncarc(cod_free_size); forGenEst = T; nboot = NULL; silent = TRUE
estdwp.psiHatcc <- function(x, ncarc, nboot = NULL,
    forGenEst = FALSE, silent = TRUE, ...){
  if (!is.list(ncarc))
    stop("estdwp: ncarc must be list with element names matching those of x")
  if (!all(names(x) %in% names(ncarc)))
    stop("estdwp: not all x names found in ncarc")
  output <- list()
  arglist <- as.list(match.call())
  arglist[[1]] <- NULL
  for (sz in names(x)){
    arglist[["x"]] <- x[[sz]]
    arglist[["ncarc"]] <- ncarc[[sz]]
    output[[sz]] <- do.call(estdwp, arglist)
  }
  if (forGenEst){
    tmp <- data.frame(array(
      dim = c(nrow(output[[1]]), length(output) + 1),
      dimnames = list(NULL, c("turbine", names(output)))
     ))
    tmp$turbine <- output[[1]]$turbine
    for (sz in names(output)) tmp[, sz] <- output[[sz]][, "dwp"]
    output <- tmp
    class(output) <- c(class(output), "dwphat", "GenEst")
  } else {
    class(output) <- c(class(output), "dwphat", "notGenEst")
  }
  return(output)
}

#' Format DWP Estimate for Use in GenEst
#'
#' GenEst requires \code{dwp} data to be formatted as a data frame with columns
#'  for turbine ID and for estimated \code{dwp} for each carcass class. To
#'  incorporate uncertainty in the estimates, \code{nsim} simulated copies of
#'  the basic format are appended to the columns in the data set.
#'
#' @param dwphat a \code{dwphat} object
#' @return an \code{nsim*nturbine} by \code{nclass + 1} data frame, with columns
#'  for the turbine ID and for estimated \code{dwp} for each carcass class (e.g.,
#'  \code{large}, \code{medium}, \code{small}, \code{bat}).
#' @export
formatGenEst <- function(dwphat){
  if (!"dwphat" %in% class(dwphat))
    stop("formatGenEst: dwphat must be a dwphat object")
  if ("GenEst" %in% class(dwphat)) return(dwphat)
  if (is.matrix(dwphat)){
    if (ncol(dwphat) == 1){
      output <- data.frame(turbine = "all", dwp = as.vector(dwphat))
      class(output) <- c("dwphat", "GenEst", "data.frame")
      return(output)
    }
    nboot <- nrow(dwphat)
    nms <- exclude("total", colnames(dwphat))
    output <- data.frame(
      turbine = rep(nms, times = nboot),
      dwp = c(matrix(dwphat[, exclude("total", colnames(dwphat))], nrow = length(nms),
        byrow = T)),
      stringsAsFactors = FALSE
    )
  } else if (is.list(dwphat)){
    ccnm <- names(dwphat) # carcass classes
    nsim <- nrow(dwphat[[1]])
    nms <- c("turbine", ccnm)
    tnm <- colnames(dwphat[[1]])[-ncol(dwphat[[1]])]
    turbs <- rep(tnm, nsim)
    output <- data.frame(array(dim = c(length(turbs), length(nms)), 
      dimnames = list(NULL, nms)), stringsAsFactors = FALSE)
    output$turbine <- turbs
    for (sz in ccnm){
      output[, sz] <- c(t(dwphat[[sz]][, -ncol(dwphat[[sz]])]))
    }
  }
  class(output) <- c("dwphat", "GenEst", "data.frame")
  return(output)
}

#' @rdname Plot
#' @param ... arguments that may be passed to plotting functions
#' @export
plot.polygonLayout <- function(x, ...){
  dotlist <- list(...)
  arg_par <- list()
  arg_par[["mfrow"]] <- c(ceiling(sqrt(length(x))), round(sqrt(length(x))))
  arg_par[["oma"]] <- c(2, 2, 1, 1)
  arg_par[["mar"]] <- c(0, 0, 0, 0)
  for (nm in names(arg_par))
    if (nm %in% names(dotlist)) arg_par[[nm]] <- dotlist[[nm]]
  arg_par[["mgp"]] <- c(2, 0.7, 0)
  arg_par[["tck"]] <- -0.015
  do.call(par, arg_par)

  lmx <- max(sapply(x, FUN = function(ti) max(abs(ti))))
  arg_plot <- list()
  arg_plot[["xlim"]] <- lmx * c(-1, 1) * 1.1
  arg_plot[["ylim"]] <- lmx * c(-1, 1) * 1.1
  arg_plot[["bg"]] <- "white"
  arg_plot[["main"]] <- ""
  arg_plot[["sub"]] <- ""
  arg_plot[["asp"]] <- NA
  arg_plot[["asp"]] <- 1
  for (nm in names(arg_plot))
    if (nm %in% names(dotlist)) arg_plot[[nm]] <- dotlist[[nm]]
  arg_plot[["x"]] <- 0
  arg_plot[["type"]] <- "n"
  arg_plot[["axes"]] <- FALSE
  arg_polygon <- list()
  arg_polygon[["col"]] <- colors()[647]
  arg_polygon[["border"]] <- NA
  for (nm in names(arg_polygon))
    if (nm %in% names(dotlist)) arg_polygon[[nm]] <- dotlist[[nm]]

  arg_turb <- list()
  arg_turb[["pch"]] <- 1
  arg_turb[["cex"]] <- 1
  for (nm in names(arg_turb))
    if (nm %in% names(dotlist)) arg_turb[[nm]] <- dotlist[[nm]]
  arg_turb[["x"]] <- 0
  arg_turb[["y"]] <- 0
  theta <- seq(0, 2*pi, length = 500)
  r <- attr(x, "srad")
  for (ti in names(x)){
    do.call(plot, arg_plot)
    arg_polygon[["x"]] <- x[[ti]]
    do.call(polygon, arg_polygon)
    axis(1, at = pretty(par("usr")[1:2], n = 8))
    axis(2, at = pretty(par("usr")[3:4], n = 8), las = 2)
    do.call(points, arg_turb)
    lines(r * cos(theta), r * sin(theta), lty = 3)
    box()
    mtext(side = 1, line = -1.8, adj = 0.04, text = ti, cex = 1.4)
  }
}

#' @rdname Plot
#' @export
plot.layoutSimple <- function(x, ...){
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))            
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

#' Calculate Area of Intersection inside Circle and Square with Common Center
#'
#' @param r radius of the circle (vector or scalar)
#' @param s half-width of the square (scalar)
#' @return vector of intersections of interiors of circles with squaure
#' @examples
#'  # calculate area in annulus intersecting square
#' s <- 10 # radius or half-width of square
#' r <- c(11, 12) # inner and outer radii of circle
#' diff(Acins(r, s)) # intersection of square and annulus
#'# figure to illustrate the calculated area:
#' theta <- seq(0, 2 * pi, length = 1500)
#' plot(0, xlim = max(r) * c(-1, 1), ylim = max(r) * c(-1, 1),
#'   xlab = "x", ylab = "y", asp = 1, bty = "n", type = "n")
#' xi <- r[1] * cos(theta)
#' yi <- r[1] * sin(theta)
#' xo <- r[2] * cos(theta)
#' yo <- r[2] * sin(theta)
#' i1 <- which(abs(xi) <= s & abs(yi) <= s)
#' i2 <- which(abs(xo) <= s & abs(yo) <= s)
#' i2 <- sort(i2, decreasing = TRUE)
#' xi <- xi[i1]
#' yi <- yi[i1]
#' xo <- xo[i2]
#' yo <- yo[i2]
#' polygon(col = 8, border = NA,
#'   x = c(xi[xi >= 0 & yi >= 0], xo[xo >= 0 & yo >= 0]), 
#'   y = c(yi[xi >= 0 & yi >= 0], yo[xo >= 0 & yo >= 0]))
#' polygon(col = 8, border = NA, 
#'   x = c(xi[xi <= 0 & yi >= 0], xo[xo <= 0 & yo >= 0]), 
#'   y = c(yi[xi <= 0 & yi >= 0], yo[xo <= 0 & yo >= 0]))
#' polygon(col = 8, border = NA,
#'   x = c(xi[xi <= 0 & yi <= 0], xo[xo <= 0 & yo <= 0]), 
#'  y = c(yi[xi <= 0 & yi <= 0], yo[xo <= 0 & yo <= 0]))
#' polygon(col = 8, border = NA,
#'  x = c(xi[xi >= 0 & yi <= 0], xo[xo >= 0 & yo <= 0]), 
#'  y = c(yi[xi >= 0 & yi <= 0], yo[xo >= 0 & yo <= 0]))
#' lines(r[1] * cos(theta), r[1]* sin(theta))
#' lines(r[2]* cos(theta), r[2] * sin(theta))
#' rect(-s, -s, s, s)
  
#'  # calculate areas in series of 1 m annuli extending to corner of square
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

#' Add Carcasses to a Site Layout
#'
#' After the site layout is analyzed and structured by rings for analysis,
#'  carcass data may still need to be added to the site data. \code{addCarcass}
#'  grabs carcass location data from a shape file or data frame and formats it
#'  into ring data, with carcass tallies in every 1m ring from the turbine to
#'  the maximum search distance away from any turbine.
#'
#' @param x carcass data to insert into \code{data_ring}
#' @param data_ring ring data for receiving carcass data from \code{x}
#' @param plotLayout (optional) \code{shapeLayout} object to facilitate proper
#'  insertion of carcass data into ring structure.
#' @param ncarcReset boolean to direct the function to set the carcass counts
#'  in all the rings to 0 before adding the new carcasses (default) or to add the
#'  new carcasses to the old totals (\code{ncarcReset = FALSE}).
#' @param ccCol name of carcass class column (optional). Typically, the "carcass 
#'  class" would be for carcass characteristics that would be expected to affect
#'  distances that carcasses would fall from the turbine. For example, distances 
#'  would not be expected to be the same for large and small carcasses, and bats
#'  may have significantly different distance distributions than small birds. The ccCol
#'  could also be used for subsetting by any covariate that would be expected 
#'  to interact with carcass distance distributions, like season (if winds vary 
#'  by season) or turbine type (if the site has a diverse mix of turbines of 
#'  different sizes or types). Additionally, ccCol can be used to subset the data by area 
#'  (for example, NW, NE, SW, SE; or hilltop, river bank) or any other discrete 
#'  covariate that the user may be interested in.
#' @param unitCol name of unit column
#' @param rCol name of column with carcass distances
#' @param ... ignored
#' @return an object of class \code{rings} with a tally of the number of
#'  of carcasses discovered in each concentric 1m ring from the turbine to the 
#'  most distant point searched.
#' @examples
#'  data(layout_simple)
#'  data(carcass_simple)
#'  sitedata <- initLayout(layout_simple)
#'  ringdata <- prepRing(sitedata)
#'  ringsWithCarcasses <- addCarcass(carcass_simple, data_ring = ringdata)
#' @export
addCarcass <- function(x, ...) UseMethod("addCarcass", x)

#' @rdname addCarcass
#' @export
addCarcass.shapeCarcass <- function(x, data_ring, plotLayout = NULL,
    ncarcReset = TRUE, ccCol = NULL, ...){
  unitCol <- x$unitCol
  scVar <- data_ring$scVar
  turbi <- x$carcasses[, unitCol, drop = TRUE]
  if (any(!x$carcasses[, unitCol, drop = TRUE] %in% names(data_ring$ncarc))){
    catturb <- x$carcasses[, unitCol, drop = TRUE]
    warning("Carcasses found at unsearched turbines.")
    if (all(!x$carcasses[, unitCol, drop = TRUE] %in% names(data_ring$ncarc)))
      stop("All carcasses were found at unsearched turbines. Mismatched ",
           "turbine names in carcass data (x) and data_ring?")
  }
  if (!is.null(ccCol)){ # carcass class superstructure
    if (!ccCol %in% names(x$carcasses)) stop("ccCol not in carcass data frame")
    if (!is.character(x$carcasses[, ccCol, drop = TRUE]))
      stop("carcass classes must be class names (character)")
    szoutput <- list()
    for (sz in unique(x$carcasses[, ccCol, drop = TRUE])){
      # subset cod by ccCol
      szoutput[[sz]] <- addCarcass(subset(x, subset = sz, select = ccCol),
        data_ring = data_ring, plotLayout = plotLayout, ncarcReset = ncarcReset)
    }
    class(szoutput) <- "ringscc"
    return(szoutput)
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
        turbine = x$carcasses[, unitCol, drop = TRUE],
        scVar = plotLayout$layout[ind, scVar, drop = TRUE],
        r = ceiling(unname(sqrt(rowSums(
          (plotLayout$tcenter[turbi, ] - sf::st_coordinates(x$carcasses))^2
        ))))
      )
      names(rdat0) %<>% gsub("scVar", scVar, .)
    } else {
      # check whether carcass search classes are represented at the distances
      # given in x
      if (!scVar %in% names(x$carcasses))
        stop("scVar not a column in carcass data in x in addCarcass.shapeCarcass")
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
        (data_ring$tcenter[turbi, ] - sf::st_coordinates(x$carcasses))^2
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

#' @rdname addCarcass
#' @export
addCarcass.data.frame <- function(x, data_ring, ccCol = NULL,
    ncarcReset = TRUE, unitCol = "turbine", rCol = "r", ...){
  if (!is.null(ccCol)){
    if (!ccCol %in% names(x)) stop("ccCol not in carcass data frame")
    if (!is.character(x[, ccCol]))
      stop("carcass class must be a vector of names of carcass classes (char)")
    output <- list()
    for (sz in unique(x[, ccCol])){
      output[[sz]] <- addCarcass(x[x[, ccCol] == sz, ], data_ring = data_ring,
        ncarcReset = ncarcReset, unitCol = unitCol, rCol = rCol)
    }
    class(output) <- "ringscc"
    return(output)
  }
  if ("rings" %in% class(data_ring)){
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
  } else {
      stop("addCarcass: class(data_ring) must be rings")
  }
  return(data_ring)
}

#' Subset Data from an Imported and Formatted Shape File
#'
#' @param x object to be subsetted
#' @param subset values to subset by. For example, to subset \code{x} to include
#'  only turbines \code{"t1"} and \code{"t2"}, then \code{subset = c("t1", "t2")}.
#'  The name of the column with turbine names is given in \code{select}.
#' @param select the name of the column with the values to subset by. For example,
#'  to subset \code{x} by turbines names "t1" and "t2" as found in the "turbine"
#'  column in the data, use \code{select = "turbine"} and \code{subset = c("t1", "t2")}.
#' @param ... ignored
#' @return object of the same class as \code{x}, subsetted to values of \code{select}
#'  equal to some element in \code{subset}.
#' @rdname subset
#' @export
subset.shapeCarcass <- function(x, subset, select, ...){
  # like the previous subCarcass except S3 and substituting "subset" for "value"
  # and "select" for "ccCol"
  if (!select %in% names(x$carcasses)) stop("'select' not in carcass data")
  output <- list()
  output$carcasses <- x$carcasses[x$carcasses[, select, drop = TRUE] %in% subset, ]
  output$unitCol  <- x$unitCol
  output$ncarc <- table(output$carcasses[, x$unitCol, drop = TRUE])
  class(output) <- "shapeCarcass"
  return(output)
}

#' @rdname subset
#' @export
subset.shapeLayout <- function(x, subset, select, ...){
  if (length(select) != 1)
    stop("subset: 'select' must be the name of a single column in x")
  if (!select %in% names(x$layout))
    stop("subset: 'select' not in layout data")
  x$layout <- x$layout[x$layout[, select, drop = TRUE] %in% subset, ]
  x$layoutAdj <- x$layoutAdj[x$layoutAdj[, select, drop = TRUE] %in% subset, ]
  if (select %in% names(x$turbines)){
    x$turbines <- x$turbines[x$turbines[, select, drop = TRUE] %in% subset, ]
    x$tset <- x$turbines[, x$unitCol, drop = TRUE]
    x$tcenter <- x$tcenter[x$tset,, drop = FALSE]
  }
  class(x) <- "shapeLayout"
  return(x)
}

#' @rdname ddFit
#' @export
ddFit.ringscc <- function(x, distr = "standard", scVar = NULL, rCol = "r",
    expoCol = "exposure", ncarcCol = "ncarc", silent = FALSE, ...){
  output <- list()
  arglist <- as.list(match.call())
  arglist[[1]] <- NULL
  for (sz in names(x)){
    arglist$x <- x[[sz]]
    output[[sz]] <- do.call(ddFit, arglist)
  }
  class(output) <- "ddArraycc"
  return(output)
}

#' Simple Function to Extract Carcass Counts
#'
#' Carcass counts are easy to extract from any of the data structures, but it
#'  may be difficult to remember where to retrieve the data from for any particular
#'  structure. \code{getncarc} simplifies the task by having the same usage for
#'  all data types.
#'
#' @param x a data structure with \code{ncarc} buried in it somewhere
#' @param ... ignored
#'
#' @return \itemize{
#'    \item{scalar number of carcasses used in the fitted model 
#'     (\code{dd} and \code{ddArray} objects}
#'    \item{vector of numbers of carcasses of each size used in the fitted models 
#'     (\code{ddArraycc} objects)}
#'    \item{vector of carcass counts at each turbine and total at the site 
#'     (\code{xyLayout} and \code{rings} objects}
#'    \item{list of vectors of carcass counts at each turbine for each carcass class}
#'  }
#' @export
getncarc <- function(x, ...) UseMethod("getncarc", x)

#' @rdname getncarc
#' @export
getncarc.ringscc <- function(x, ...){
  lapply(x, "[[", "ncarc")
}

#' @rdname getncarc
#' @export
getncarc.rings <- function(x, ...){
  x$ncarc
}

#' @rdname getncarc
#' @export
getncarc.xyLayout <- function(x, ...){
  x$ncarc
}

#' @rdname getncarc
#' @export
getncarc.ddArray <- function(x, ...){
  x[[1]]$ncarc
}

#' @rdname getncarc
#' @export
getncarc.ddArraycc <- function(x, ...){
  ncarc <- numeric(length(x))
  names(ncarc) <- names(x)
  for (sz in names(ncarc)) ncarc[sz] <- getncarc(x[[sz]])
  return(ncarc)
}

#' @rdname getncarc
#' @export
getncarc.dd <- function(x, ...){
  return(x$ncarc)
}

#' @rdname Plot
#' @export
plot.psiHat <- function(x, ...){
  arglist <- list(...)
  arglist$x <- 0
  arglist$xlab <- "turbine"
  if ("xlab" %in% names(arglist)){
    xlab <- arglist$xlab
    arglist$xlab <- ""
  }
  if ("col" %in% names(arglist)){
    col <- arglist$col
    arglist$col <- 1
  } else {
    col <- 1
  }
  if (!"type" %in% names(arglist)) arglist$type <- "n"
  if (!"xlim" %in% names(arglist)) arglist$xlim <- c(0.5, ncol(x) + 0.5)
  if (!"ylim" %in% names(arglist)) arglist$ylim <- range(x)
  if (!"axes" %in% names(arglist)) arglist$axes <- FALSE
  if (!"ylab" %in% names(arglist)) arglist$ylab  <- expression(hat(psi))
  if (!"cex.lab" %in% names(arglist)) arglist$cex.lab <- 1.2
  bxwd <- 0.4
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))            
  par(mar = c(5.4, 5, 2, 1))
  do.call(plot, arglist)
#  plot(0, type = "n", xlim = c(1, ncol(x)), ylim = range(x), axes = FALSE,
#    xlab = "", ylab = expression(hat(psi)), cex.lab = 1.2)
  box()
  mtext(side = 1, line = 4, xlab, cex = 1.3)
  axis(1, at = 1:ncol(x), labels = FALSE)
  text(x = 1:ncol(x), y = par("usr")[3] - diff(par("usr")[3:4] * 0.04),
    labels = colnames(x), xpd = TRUE, srt = 40)
  axis(2)
  for (ti in 1:ncol(x)){
    qtls <- quantile(x[, ti], prob = c(0.005, 0.05, 0.25, 0.5, 0.75, 0.95, 0.995))
    rect(ti - bxwd, qtls[3], ti + bxwd, qtls[5], border = col)
    lines(ti + bxwd * c(-1, 1), rep(qtls[4], 2), lwd = 2, col = col)
    lines(rep(ti, 2), qtls[c(2, 3)], col = col)
    points(rep(ti, 2), qtls[c(2, 6)], pch = 3, col = col)
    lines(rep(ti, 2), qtls[c(5, 6)], col = col)
    lines(rep(ti, 2), qtls[c(2, 3) - 1], lty = 3, col = col)
    lines(rep(ti, 2), qtls[c(5, 6) + 1], lty = 3, col = col)
  }
}

#' @rdname Plot
#' @export
plot.dwphat <- function(x, ...){
  bxwd <- 0.4
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))            
  par(mar = c(5.4, 5, 2, 1))
  plot(0, type = "n", xlim = c(1, ncol(x)) + bxwd*c(-1, 1), ylim = range(x, na.rm = TRUE),
    axes = FALSE, xlab = "", ylab = expression(widehat(dwp)), cex.lab = 1.2)
  box()
  mtext(side = 1, line = 4, "turbine", cex = 1.3)
  axis(1, at = 1:ncol(x), labels = FALSE)
  text(x = 1:ncol(x), y = par("usr")[3] - diff(par("usr")[3:4] * 0.04),
    labels = colnames(x), xpd = TRUE, srt = 40)
  axis(2)
  for (ti in 1:ncol(x)){
    qtls <- quantile(x[, ti], prob = c(0.005, 0.05, 0.25, 0.5, 0.75, 0.95, 0.995),
      na.rm = TRUE)
    rect(ti - bxwd, qtls[3], ti + bxwd, qtls[5])
    lines(ti + bxwd * c(-1, 1), rep(qtls[4], 2), lwd = 2)
    points(rep(ti, 2), qtls[c(2, 6)], pch = 3)
    lines(rep(ti, 2), qtls[c(2, 3)])
    lines(rep(ti, 2), qtls[c(5, 6)])
    lines(rep(ti, 2), qtls[c(2, 3) - 1], lty = 3)
    lines(rep(ti, 2), qtls[c(5, 6) + 1], lty = 3)
  }
}

#' Export Estimated Density-Weighted Proportion to File in Proper GenEst Format
#'
#' GenEst imports DWP from files with comma-separated values (.csv), with a 
#'  column giving turbine ID and a column of DWP values for each carcass class.
#'  Column lengths are equal to the number of turbines times the number of 
#'  simulation reps (typically, \code{nsim = 1000}), giving \code{nsim} copies of 
#'  DWP values for each turbine in a single column. 
#'
#' NOTE: The \code{.csv} file uses the English convention (used in USA, UK, 
#'  Mexico, China, India, Australia, Japan, Korea, and others) with the comma ( , ) 
#'  and not the  semi-colon ( ; ) to separate values among different columns and 
#'  uses the period ( . ) as the decimal mark. Although GenEst can seamlessly
#'  accommodate either format, users in countries where the comma or other 
#'  character is used as a decimal mark may need to adjust their software settings
#'  or edit the data to be able to view it in a spreadsheet program (such 
#'  as Excel). 
#'
#' @param dwp a \code{dwphat} object
#' @param file name of file to export the \code{dwp} estimates to
#' @return The function writes the formatted data to a \code{.csv} file and 
#'  returns NULL.
#' @export
exportGenEst <- function(dwp, file){
  if (!"dwphat" %in% class(dwp)) stop("exportGenEst: 'dwp' must be a dwphat object")
  if (!"GenEst" %in% class(dwp)) dwp <- formatGenEst(dwp)
  write.csv(dwp, file = file, quote = FALSE, row.names = FALSE)
}

#' Run Models through a Sieve to Filter out Dubious Fits
#'
#' A set of fitted models (\code{\link[=ddFit]{ddArray}}) is filtered according 
#'  to a set of criteria that test for high AIC, high-influence points, and
#'  plausibility of the tail probabilities of each fitted distribution.
#'  \code{modelFilter} will either auto-select the best model according to a set
#'  of pre-defined, objective criteria or will will return all models that meet
#'  a set of user-defined, or default criteria. A table of how the models
#'  score according to each criterion is printed to the console.
#'
#' The criteria to test are entered in a list (\code{sieve}) with components:
#' \enumerate{
#'  \item \code{$rtail} = vector of probabilities that define a checkpoints on distributions
#'   to avoid situations where a model that may fit well within the range of data
#'   is nonetheless implausible because it predicts a significant or substantial 
#'   probability of carcasses falling great distances from the nearest turbine. 
#'   The default is to check whether or not a distribution predicts that less than 
#'   50\% of carcasses fall within 80 meters, 90\% within 120 meters, 95\% within 
#'   150 meters, or 99\% within 200 meters. Distributions that fall below any of 
#'   these points (for example predicting only 42\% within 80 meters or only 74\% 
#'   within 120 meters) fail the default \code{rtail} test. The format of the 
#'   default for  the test is \code{$rtail = c(p80 = 0.5, p120 = 0.90, 
#'   p150 = 0.95, p200 = 0.99)}. Users may override the default by using, for example, 
#'   \code{sieve = list(rtail = c(p80 = 0.8, p120 = 0.99, p150 = 0.99, p200 = 0.999))}
#'   in the argument list for a more stringent test or for a situation where
#'   turbines are small or winds are light. Alternatively, users may forego the
#'   test altogether by entering \code{sieve = list(rtail = FALSE)}. If specific
#'   probabilities are provided, they must be in a vector of length 4 with names
#'  "\code{p80}" etc. as in the examples above.
#'  \item \code{$ltail} = vector of probabilities that define checkpoints on distributions
#'   to avoid situations where the search radius is short and a distribution that
#'   fits the limited data set well but crashes to zero just outside the search
#'   radius. The default is to check whether or not a distribution predicts that
#'   greater than 50\% of carcasses fall with 20 meters or 90\% within 50 meters.
#'   Distributions that pass above either of these checkpoints (for example
#'   predicting 61\% of carcasses within 20 meters or 93\% within 50 meters)
#'   are eliminated by the default \code{ltail} test. The format of the default for
#'   the test is \code{$ltail = c(p20 = 0.5, p50 = 0.90)}. Users may override the
#'   default by using, for example, \code{sieve = list(rtail = c(p20 = 0.6, p50 = 0.8))}
#'   in the argument list for a situation where it is known that carcasses beyond
#'   50 meters are common.
#'  \item \code{$aic} = a numeric scalar cutoff value for model's delta AICc
#'   scores. Models with AICc scores exceeding the minimum AICc among all the 
#'   fitted models by \code{sieve$aic} or more fail the test. The default value 
#'   is 10. Users may override the default by using, for example, 
#'   \code{sieve = list(aic = 7)} in the argument list to use a delta AIC score 
#'   of 7 as the cutoff or may forego the test altogether by setting 
#'   \code{sieve = list(aic = FALSE)}
#'  \item \code{$hin} = \code{TRUE} or \code{FALSE} to test for high influence points,
#'   the presence of which cast doubt on the reliability of the model. The function
#'   defines "high influence" as models with high leverage points, namely, points
#'   with \eqn{\frac{h}{1 - h} >  \frac{2p}{n - 2p}}{h/(1 - h) >  2p/(n - 2p)} 
#'   (where \eqn{h} is leverage, \eqn{p} is the number of parameters in the model, 
#'   and \eqn{n} is the search radius) with Cook's distance \code{> 8/(n - 2*p)}. 
#'   The criteria for high influence points were adapted from Brian Ripley's GLM 
#'   diagnostics package \code{boot} (\code{\link[boot]{glm.diag}}). The test is 
#'   perhaps most valuable in identifying distributions with high probability of 
#'   carcasses landing well beyond what could reasonably be expected.
#' }
#' 
#' Several choices of pre-defined \code{sieve}s are available (or, as described 
#'  above, users may define their own criteria):
#' \describe{
#'  \item{\code{sieve = "default"}}{The models are ordered by the following
#'   criteria: 
#'   \enumerate{
#'     \item extensibility
#'    \item weight of right tail (discounting models that predict implausibly 
#'     high proportions of carcasses beyond the search radius)
#'    \item weight of the left tail (discounting models that predict implausibly 
#'     high proportions of carcasses near the turbines)
#'    \item AICc test (discounting models with delta AICc > 10)
#'    \item high influence points (discounting models in which one or more of the
#'     data points exert a high influence on the fitted model, according to 
#'     Ripley's GLM  diagnostics package \code{boot} (\code{\link[boot]{glm.diag}}))
#'    \item ranking by AICc
#'   }
#'   Precise definitions of the default sieve parameters are given in 
#'   \code{sieve_default}.}
#'  \item{\code{sieve = NULL}}{Returns a list of the extensible models without
#'   scoring them by other model selection criteria.}
#'  \item{\code{sieve = "win"}}{Sorts models by high-influence points and AICc}
#'  \item{\code{sieve = list(<custom>)}}{User provides a custom sieve, which may
#'   be a modification of the default sieve or de novo. To modify the default,
#'   use, for example, \code{sieve = list(hin = FALSE)} to disable the \code{hin}
#'   test but keep the other default tests, or \code{sieve = list(aic = 7)} to
#'   use 7 rather than 10 as the AIC cutoff, or 
#'   \code{sieve = list(ltail = c(p20 = 0.3, p50 = 0.8))} to use a more stringent 
#'   left tail test that requires CDF graphs to pass below the points (20, 0.3) 
#'   and (50, 0.8). Custom \code{ltail} and \code{rtail} parameters must match the 
#'   formats of the default tests, but their probabilities may vary. To turn off
#'   the \code{aic} filter, use \code{sieve = list(aic = Inf)}. To turn off the
#'   \code{ltail} filter, use  \code{sieve = list(ltail = c(p20 = 1, p50 = 1))}.
#'   To turn off the \code{rtail} filter, use 
#'   \code{sieve = list(rtail = c(p80 = 0, p120 = 0, p150 = 0, p200 = 0))}. These
#'   custom components may be mixed and matched as desired.}
#' }
#' @param dmod a \code{\link[=ddFit]{ddArray}} object
#' @param sieve a list of criteria for ordering models
#' @param quiet boolean to suppress (\code{quiet = TRUE}) or allow 
#'  (\code{quiet = FALSE}) messages from \code{modelFilter}
#' @return An \code{fmod} object, which is an unordered list of extensible models if
#'  \code{sieve = NULL}; otherwise, a list of class \code{fmod} with following 
#'  components:
#'  \describe{
#'    \item{\code{$filtered}}{the selected \code{dd} object or a \code{ddArray} list of
#'     models that passed the tests}
#'    \item{\code{$scores}}{a matrix with all models tested (rownames = model names) and 
#'      the results of each test (columns \code{aic_test}, \code{rtail}, 
#'      \code{ltail}, \code{hin}, \code{aic})}
#'    \item{\code{$sieve}}{the test criteria, stored in a list with
#'      \itemize{ 
#'        \item \code{$aic_test} = cutoff for AIC
#'        \item \code{$hin} = boolean to indicate whether high influence points were
#'         considered
#'        \item \code{$rtail} = numeric vector giving the probabilities that the
#'         right tail of the distribution must exceed at distances of 80, 120, 
#'         150, and 200 meters in order to pass
#'        \item \code{$ltail} = numeric vector giving the probabilities that the 
#'        left tail of the distribution must NOT exceed at distances of 20 and 
#'        50 meters in order to pass
#'      }
#'    }
#'    \item{\code{models}}{a list (\code{ddArray} object) of all models tested}
#'    \item{\code{note}}{notes on the tests}
#'  }
#'  When a \code{fmod} object is printed, only a small subset of the elements are
#'  shown. To see a full list of the objects, use \code{names(x)}, where \code{x}
#'  is the name of the \code{fmod} return value. The elements
#'  can be extracted in the usual R way via, for example, \code{x$sieve} or 
#'  \code{x[["sieve"]]}.
#'  
#' @examples
#'  data(layout_simple)
#'  data(carcass_simple)
#'  sitedata <- initLayout(layout_simple)
#'  ringdata <- prepRing(sitedata)
#'  ringsWithCarcasses <- addCarcass(carcass_simple, data_ring = ringdata)
#'  distanceModels <- ddFit(ringsWithCarcasses)
#'  stats(distanceModels)
#'  stats(distanceModels[["tnormal"]])
#'  stats(distanceModels[["lognormal"]])
#'  
#' @export
modelFilter <- function(dmod, sieve = "default", quiet = FALSE){
  if (!"ddArray" %in% class(dmod)) stop("modelFilter: dmod must be a ddArray")
  if (is.null(sieve))
    return(dmod[which(sapply(dmod, function(di) di[["extensible"]] == 1))])
  if (identical(sieve, "default") || identical(sieve, "auto_select")){
    fil <- sieve_default
  } else if (identical(sieve, "win")){
    fil <- sieve_win
  } else if (is.list(sieve)){
    fil <- sieve_default
    for (si in names(sieve)){
      # error-check
      if (si == "aic"){
        if (length(sieve[[si]]) == 1 && is(sieve[[si]], "logical")){
          if (sieve[[si]]) 
            stop("sieve['aic'] must be non-negative scalar or FALSE")
          fil[[si]] <- 1e6
          next
        }
        if (!is.numeric(sieve[[si]]) || length(sieve[[si]]) != 1 || sieve[[si]] < 0){
          stop("sieve[['aic']] must be non-negative scalar")
        } else {
          fil[[si]] <- sieve[[si]]
          next
        }
      }
      if (si == "hin"){
        if(!is(sieve[[si]], "logical") && length(sieve[[si]] == 1)){
          stop("sieve[['hin']] must be boolean")
        } else {
          fil[[si]] <- sieve[[si]]
          next
        }
      }
      if (si == "rtail"){
        if(length(sieve[[si]]) == 1 && is(sieve[[si]], "logical")){
          if (sieve[[si]]){
            fil[[si]] <- sieve_default[[si]]
          } else {
            fil[[si]] <- sieve_default[[si]] * 0
          }
          next
        }
        if(!is.numeric(sieve$rtail) || !length(sieve$rtail) == 4 ||
            any(sieve$rtail < 0) || any(sieve$rtail > 1)){
          stop("sieve[['rtail']] must be a numeric vector of length 4 with values in [0, 1]")
        }
        if (is.null(names(sieve$rtail)))
          names(sieve$rtail) <-  names(sieve_default$rtail)
        if(!identical(names(sieve$rtail), names(sieve_default$rtail))){
          stop(
            "names(sieve[['rtail']]) must be identical to ",
            "names(sieve_default[['rtail']])"
          )
        }
        fil[[si]] <- sieve[[si]]
        next
      }
      if (si == "ltail"){
        if(length(si) == 1 && is(si, "logical")){
          if (sieve[[si]]){
            fil[[si]] <- sieve_default[[si]]
          } else {
            fil[[si]] <- sieve_default[[si]]/sieve_default[[si]]
          }
          next
        }
        if(!is.numeric(sieve$ltail) || !length(sieve$ltail) == 2 ||
            any(sieve$ltail < 0) || any(sieve$ltail > 1))
          stop("sieve[['ltail']] must be a numeric vector of length 2 with values in [0, 1]")
        if (is.null(names(sieve$ltail)))
          names(sieve$ltail) <- names(sieve_default$ltail)
        if (!identical(names(sieve$ltail), names(sieve_default$ltail))){
          stop("names(sieve[['ltail']]) must be identical to ",
               "names(sieve_default[['ltail']])")
        }
        fil[[si]] <- sieve[[si]]
      }
    }
  } else {
    stop("modelFilter: sieve must be 'auto', 'default', NULL, or a custom list")
  }
  #extensible models
  critok <- c("extensible", "rtail", "ltail", "aicc", "hin", "deltaAICc")
  ptab <- array(1, dim = c(length(dmod), length(critok)), 
    dimnames = list(names(dmod), critok))
  if (identical(sieve, "win")){
    ptab[, "extensible"] <- 1
    ptab[, "rtail"] <- 1
    ptab[, "ltail"] <- 1
  } else {
    for (nm in rownames(ptab)){
      ptab[nm, "extensible"] <- dmod[nm]$extensible * 1
      if (ptab[nm, "extensible"] == 0){
        ptab[nm, "rtail"] <- 0
        ptab[nm, "ltail"] <- 1
      } else {
        ptab[nm, "rtail"] <- 1 * (
          pdd(q =  80, dmod[nm]) > fil[["rtail"]]["p80"] &
          pdd(q =  120, dmod[nm]) > fil[["rtail"]]["p120"] &
          pdd(q =  150, dmod[nm]) > fil[["rtail"]]["p150"] &
          pdd(q =  200, dmod[nm]) > fil[["rtail"]]["p200"]
        )
        ptab[nm, "ltail"] <- 1 * (
          pdd(q = 20, dmod[nm]) < fil[["ltail"]]["p20"] &
          pdd(q = 50, dmod[nm]) < fil[["ltail"]]["p50"]
        )
      }
    }
  }
  if (identical(fil$hin, TRUE) | identical(fil$hin, 1)){
   for (nm in rownames(ptab)){
      dnost <- boot::glm.diag(dmod[[nm]]$glm)
      h <- dnost$h
      cook <- dnost$cook
      n <- dmod[[nm]]$n
      p <- dmod[[nm]]$k
      ptab[nm, "hin"] <- 
        all(!((cook > 8/(n - 2*p)) & (h/(1 - h) > 2*p/(n - 2*p)))) * 1
    }
  }
  if (!identical(fil[["aic"]], FALSE)){
    ptab[, "aicc"] <- 
      sapply(dmod, "[[", "aicc") %>% "<="(., min(.) + fil[["aic"]]) * 1
    ptab[, "deltaAICc"] <- sapply(dmod, "[[", "aicc")
  }
  ptab[, "deltaAICc"] <- ptab[, "deltaAICc"] - min(ptab[, "deltaAICc"])
  ptab[, "aicc"] <- (ptab[, "deltaAICc"] <= 10) * 1
  if (!identical(fil[["aic"]], FALSE)){
    for (nm in rownames(ptab))
      ptab[nm, "aicc"] <- (ptab[nm, "deltaAICc"] <= fil[["aic"]]) * 1
  }
  ptab[is.na(ptab)] <- 1
  ptab <- ptab[order(
    ptab[, "extensible"],
    ptab[, "rtail"],
    ptab[, "ltail"],
    ptab[, "aicc"],
    ptab[, "hin"],
    -ptab[, "deltaAICc"],
    decreasing = TRUE
  ), ]
  if (identical(sieve, "win")) ptab[, c("extensible", "rtail", "ltail")] <- NA
  output <- list()
  output[["filtered"]] <- dmod[rownames(ptab)[1]]
  output[["scores"]] <- ptab
  output[["sieve"]] <- fil
  output[["models"]] <- dmod
  output[["note"]] <- ifelse(is.character(sieve), sieve, "")
  class(output) <- "fmod"
  if (!identical(sieve, "win") && sum(ptab[, "extensible"]) == 0){ 
    output[["note"]] <- c(output[["note"]], "None of the models are extensible")
    return(output)
  }
  if (!identical(sieve, "default") && !identical(sieve, "win")){
    if (all(rowSums(ptab[, 1:5]) < 5)){
      if (!quiet){
        message(
          "None of the models meet all the filter criteria in sieve. ",
          "Check data format, select less strict criteria, or set ",
          "sieve = 'default'."
        )
      }
    }
  }
  return(output)
}

#' Simple Extension of a \code{dd} Model beyond the Search Radius
#'
#' Extend a distance model beyond the search radius via multiplication by a fixed,
#' assumed constant rather than the default normalization used for extensible
#' models. \code{psi_extend} should not be used with \code{psiHat} objects that
#' were calculated with \code{extent = "full"}.
#'
#' @param psi \code{\link[=estpsi]{psiHat}} object
#' @param fwin fraction of carcasses assumed to lie within the search radius. If
#'  \code{psi} includes \code{psiHat} for multiple carcass classes, \code{fwin}
#'  should be either a vector with one value for each carcass class so that 
#'  \code{length(psi) = length(fwin)} or a scalar (which assumes all carcasses, 
#'  regardless of carcass class, have the same probability of landing outside the
#'  search radius).
#' @return \code{psiHat} object extended beyond the search radius
#' @export
psi_extend <- function(psi, fwin){
  if (!is.numeric(fwin)) stop("fwin must be numeric")
  if (any(fwin > 1) || any(fwin <= 0)) stop("fwin must be between 0 and 1")
  if (attr(psi, "extent") == "full")
    warning("psi has already been extended. psi_extend is extending psi again.")
  if ("psiHat" %in% class(psi)){
    if (is.list(psi)){
      if (length(fwin) == 1){
        fwin <- rep(fwin, length(psi))
      } else if (length(fwin) != length(psi)){
        stop(
          "psi_extend: fwin must be a scalar or a vector of values ",
          "corresponding to the carcass classes represented in psi."
        )
      }
      for (i in 1:length(psi)){
        attr(psi[[i]], "fwin") <- fwin[i]
        attr(psi[[i]], "extent") <- "full"
        psi[[i]] <- psi[[i]] * fwin[i]
      }
    } else {
      psi <- psi * fwin
      attr(psi, "fwin") <- fwin
      attr(psi, "extent") <- "full"
    }
  } 
  psi
}
