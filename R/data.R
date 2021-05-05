#' Names of the GLM Coefficients for Each Distribution
#'
#' @format A list with the names of the coefficients (vector of character strings)
#'  as they appear in the \code{$coefficients} value returned from \code{glm} for
#'  each model.
"cof_name"

#' Constraints on GLM Coefficients for Extensibility to a Distribution
#'
#' @format A list of matrices giving the upper and lower bounds that each
#'  model's coefficients must meet for the model to be extensible to a
#'  distribution. The \code{parscale} column may be used in \code{optim} for
#'  fitting a truncated weighted likelihood model.
"constraints"

#' Example Polygon Data for Site Layout
#'
#' @format A data frame illustrating the required raw data format for using
#'  standard R polygons to characterize a site layout. There must be three
#'  columns: one giving turbine IDs (\code{turbine}) and columns for the \code{x}
#'  and \code{y} coordinates that delineate the plot layouts. Turbine IDs must be
#'  syntactitally valid R names, that is, combinations of letters, numbers,
#'  underscores ( _ ) and periods ( . ) and no spaces, hyphens, or other special
#'  characters. Names must not begin with a number, so \code{1, 2, 3, ...},
#'  \code{3B1}, and \code{.72S} are NOT valid names. Coordinates shoud be in
#'  meters relative to a turbine at (0, 0).
"layout_polygon"

#' Example Bare Vector Format for Eagle Data
"layout_eagle"

#' Example Carcass Distances to Accompany the Polygon Layout Data Set
#'
#' @format A data frame illustrating an import format for carcass distances.
#'  There are columns with turbine IDs (\code{turbine}) and carcass distances
#'  from turbine (\code{r}). Distances (\code{r}) are in meters from the nearest
#'  turbine. The data set is used in the "polygon" example in the User Guide.
"carcass_polygon"

#' Example Simple Geometry Data Format for Site Layout
#'
#' @format A data frame illustrating an import format for a simple description of a
#'  site layout by turbine. Each turbine (\code{turbine}) is classed according to
#'  the shape of its search plot, either \code{circular}, \code{square}, or
#'  \code{RP} (search on the roads and turbine pad only). For \code{circular}
#'  plots, all ground within \code{radius} meters of the turbine is searched.
#'  For \code{square} plots, the \code{radius} is half the width of the square
#'  along the x-axis, NOT the distance to the corner. For \code{RP} plots, the
#'  \code{radius} is the maximum distance searched on the roads. The geometry of
#'  the \code{RP} also includes a circular turbine pad with radius \code{padrad},
#'  a road width of \code{roadwidth} meters, and the number of roads
#'  (\code{n_road}) searched out to \code{radius} meters from the turbine.
"layout_simple"

#' Example Data for Site Layout on an (x, y) Grid
#'
#' @format A data frame illustrating the required raw data format for using
#'  a grid format to characterize a site layout. There are five
#'  columns: one giving turbine IDs (\code{turbine}), columns for the \code{x}
#'  and \code{y} coordinates on 1 m. grids that overlay the search plot, a 
#'  column giving the carcass count in each grid cell, and a column giving the 
#'  distance of each cell from the turbine. There is only one turbine, searched
#'  on road and pad out. Coordinates are in meters relative to the turbine at 
#'  (0, 0).
"layout_xy"

#' Carcass Data to Accompany the Simple Geometry Data Format
"carcass_simple"

#' Names of All the Available Models
"mod_all"

#' Vector of Colors Used in Graphs of Fitted Models
"mod_color"

#' Vector of Line Types Used in Graphs of Fitted Models
"mod_lty"

#' Vector of GLM Offsets for Available Models
"mod_offset"

#' Vector of Names of Standard Models
"mod_standard"

#' Vector of Names of Models Available for Grid Layout
"mod_xy"

#' Vector of Names of Models with Natural Offset
#'
#' The natural offset is the area (m^2) searched at a given distance. Models
#'  that use the natural offset are referred to as "natural". Other models may
#'  use different offsets which alter the shape of the curve with distance in an
#'  a priori way that is unaffected by the data.
"natural"

#' List of Names of the Distribution Parameters Associated with Respective Models
"parm_name"

#' Default Graphics Parameters
"par_default"

#' Test Criteria for Model Selection
#'
#' @format A list containing the parameters used for test criteria in model 
#' selection an \code{ddArray} objects. The \code{sieve_default} values are used as
#' a default in \code{\link{modelFilter}}. If desired, users may create their 
#' own tests, using \code{sieve_default} as a template. The same list elements
#' must all be present and have the same structure as the defaults, namely:
#' \describe{
#'  \item{\code{$aic}}{the cutoff for DeltaAIC scores; models with higher scores 
#'    are removed from further consideration. Default is \code{$aic = 10}}
#'  \item{\code{$hin}}{a boolean to indicate whether or not to use high leverage 
#' points as a criterion for model selection. Default is \code{$hin = TRUE}}
#'  \item{\code{$rtail}}{a vector of probabilities that the fitted model must 
#'    exceed at 80, 120, 150, and 200 meters. Default is 
#'    \code{rtail = c(p80 = 0.50, p120 = 0.90, p150 = 0.95, p200 = 0.99)}.
#'    Custom test parameters must be a vector probabilities with "p80", "p120",
#'    "p150", and "p200" in the names.}
#'  \item{\code{ltail}}{a vector of probabilities that a fitted model must 
#'    not exceed at 20 and 50 meters. Default is 
#'    \code{ltail = c(p20 = 0.50, p50 = 0.90)}. Custom test parameters must be
#'    a vector of probabilities with "p20" and "p50" in \code{names}.}
#' }
#' 
"sieve_default"

#' Test Criteria for Model Selection within Search Area
#'
#' @format A list containing the parameters used for test criteria in model 
#' selection in \code{ddArray} objects. The \code{sieve_win} values are used 
#' when either \code{sieve = "win"} or \code{extent = "win"} in arg list of
#' \code{\link{modelFilter}}. The sieve parameters are:
#' \describe{
#'  \item{\code{$aic = 10}}{the cutoff for DeltaAIC scores; models with higher 
#'		scores are removed from further consideration.}
#'  \item{\code{$hin = T}}{a boolean to indicate whether or not to use high 
#'		leverage points as a criterion for model selection.}
#'  \item{\code{$rtail}}{Appropriate only for extrapolating beyond the search
#'		radius. Automatically disabled via\code{rtail = sieve_default$rtail * 0}.}
#'  \item{\code{ltail}}{Appropriate only for extrapolating beyond the search
#'		radius. Automaticall disabled via \code{ltail = sieve_default$ltail * 0 + 1}}
#' }
#' 
"sieve_win"

#' Locations of All Carcasses in Grid Data
#' @format matrix with columns x, y, r for all 100 carcasses in the simulation to
#'  generate the carcass data for the xy grid that was searched on road and pad.
"xyr"
