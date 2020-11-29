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

#' Example Carcass Distances to Accompany the Polygon Layout Data Set
#'
#' @format A data frame illustrating an import format for carcass distances.
#'  There are columns with turbine IDs (\code{turbine}) and carcass distances
#'  from turbine (\code{r}). Distances (\code{r}) are in meters from the nearest
#'  turbine. The data set is used in the "polygon" example in the User Guide.
"carcass_polygon"

#' Example Simple Data Format for Site Layout
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

#' Names of All the Available Models
"mod_all"

#' Vector of Colors Used in Graphs of Fitted Models
"mod_color"

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
