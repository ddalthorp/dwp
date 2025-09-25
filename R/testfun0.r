padrad = 10 # measuring from turbine center rather than edge
radius = 122 # from turbine center
roadwidth = 8
x <- tlay0 <- initLayout(data.frame(
  turbine = "t1",
  radius = radius, # "max search radius" is from edge of turbine, but glm assumes center
  shape = "RP",
  padrad = padrad,
  roadwidth = roadwidth,
  n_road = 1
))
