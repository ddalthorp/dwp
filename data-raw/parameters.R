layout_polygon <- read.csv(textConnection('
"turbine","x","y"
"t1",4.814164,69.209213
"t1",-12.278842,17.930197
"t1",-64.255532,15.139502
"t1",-54.836937,-54.627868
"t1",22.953680,-27.767431
"t1",70.395491,14.092991
"t1",20.162985,19.674381
"t2",-19.604416,-3.000014
"t2",-6.348615,-17.302325
"t2",14.581596,-2.302341
"t2",6.907185,10.255786
"t2",2.372306,25.953444
"t2",-1.116063,52.465045
"t2",-5.302105,58.744108
"t2",-17.860231,68.511540
"t2",-23.441621,81.767340
"t2",-35.650911,76.185950
"t2",-25.883479,33.627855
"t2",-35.999748,9.906949
'), as.is = T)

layout_simple <- read.csv(textConnection('
"turbine","radius","shape","padrad","roadwidth","n_road"
"t1",90,"circular",NA,NA,NA
"t2",65,"square",NA,NA,NA
"t3",120,"RP",15,5,2
"t4",100,"RP",20,4,1
'), as.is = T)

mod_color = c(
  xep1 = colors()[35], #brown3 colors()[585], #darkorchid2 #
  xep01 = 2,
  xep2 = colors()[257],
  xep02 = 1,
  xep12 = colors()[122],
  xep012 = 4,
  xep123 = colors()[148],
  xep0123 = colors()[148], #yellowgreen
  tnormal = colors()[26],
  MaxwellBoltzmann = colors()[257], #deepskyblue2
  lognormal = colors()[96],
  xepi0 = 5,
  xep0 = colors()[96],
  chisq = colors()[35],
  exponential = colors()[123], #brown3
  inverse_gaussian = 4,
  constant = 8
)
mod_lty = c(
  xep1 = 2, 
  xep01 = 1,
  xep2 = 1,
  xep02 = 1,
  xep12 = 2,
  xep012 = 1,
  xep123 = 1,
  xep0123 = 2,
  tnormal = 4,
  MaxwellBoltzmann = 2,
  lognormal = 1,
  xepi0 = 3,
  xep0 = 3,
  chisq = 3,
  exponential = 3, 
  inverse_gaussian = 3,
  constant = 3
)
degOrder <-c(
  "constant",
  "chisq",
  "xep0",
  "exponential",
  "xepi0",
  "inverse_gaussian",
  "lognormal",
  "xep1",
  "xep01",
  "xep2",
  "tnormal",
  "MaxwellBoltzmann",
  "xep02",
  "xep12",
  "xep012",
  "xep123",
  "xep0123"
)
mod_all = names(mod_color)
mod_standard = c(
  "xep1",
  "xep01",
  "xep2",
  "xep02",
  "xep12",
  "xep012",
  "xep123",
  "xep0123",
  "tnormal",
  "MaxwellBoltzmann",
  "lognormal",
  "constant"
)
mod_xy = c(
  "xep1",
  "xep01",
  "xep2",
  "xep02",
  "xep12",
  "xep012",
  "xep123",
  "xep0123"
)
parm_name <- list(
  xep1 = "b1",
  xep01 = c("shape", "rate"),
  xep2 = "s2",
  xep02 = c("b0", "b2"),
  xep12 = c("b1", "b2"),
  xep012 = c("b0", "b1", "b2"),
  xep123 = c("b1", "b2", "b3"),
  xep0123 = c("b0", "b1", "b2", "b3"),
  tnormal = c("mean", "sd"),
  MaxwellBoltzmann = "a",
  lognormal = c("meanlog", "sdlog"),
  xepi0 = c("shape", "scale"),
  xep0 = "a",
  exponential = "rate",
  chisq = "df",
  inverse_gaussian = c("mean", "dispersion"),
  constant = NULL
)
cof_name <- list(
  xep1 = c("(Intercept)", "r"),
  xep01 = c("(Intercept)", "log(r)", "r"),
  xep2 = c("(Intercept)", "I(r^2)"),
  xep02 = c("(Intercept)", "log(r)", "I(r^2)"),
  xep12 = c("(Intercept)", "r", "I(r^2)"),
  xep012 = c("(Intercept)", "log(r)", "r", "I(r^2)"),
  xep123 = c("(Intercept)","r", "I(r^2)", "I(r^3)"),
  xep0123 = c("(Intercept)", "log(r)", "r", "I(r^2)", "I(r^3)"),
  tnormal = c("(Intercept)", "r", "I(r^2)"),
  MaxwellBoltzmann = c("(Intercept)", "I(r^2)"),
  lognormal = c("(Intercept)", "log(r)", "I(log(r)^2)"),
  xep0 = c("(Intercept)", "log(r)"),
  xepi0 = c("(Intercept)", "I(1/r)", "log(r)"),
  exponential = c("(Intercept)", "r"),
  chisq = c("(Intercept)", "log(r)"),
  inverse_gaussian = c("(Intercept)", "I(1/r)", "r"),
  constant = c("(Intercept)")
)
distr_names <- c(
  xep1 = "xep1",
  xep01 = "xep01 (gamma)",
  xep2 = "xep2 (Rayleigh)",
  xep02 = "xep02",
  xep12 = "xep12",
  xep012 = "xep012",
  xep123 = "xep123",
  xep0123 = "xep0123 (normal-gamma)",
  tnormal = "truncated normal",
  MaxwellBoltzmann = "Maxwell-Boltzmann",
  lognormal = "lognormal",
  xep0 = "xep0 (Pareto)",
  xepi0 = "xepi0 (inverse gamma)",
  exponential = "exponential",
  chisq = "chi-squared",
  inverse_gaussian = "inverse Gaussian",
  constant = "constant"
)
alt_names <- c(
  "gamma" = "xep01",
  "Rayleigh" = "xep2",
  "normal-gamma" = "xep0123",
  "truncated normal" =  "tnormal",
  "Maxwell-Boltzmann" = "MaxwellBoltzmann",
  "Pareto" = "xep0",
  "inverse gamma" = "xepi0",
  "exponential" = "exponential",
  "chi-squared" = "chisq",
  "inverse Gaussian" = "inverse_gaussian"
)
natural = c( # distributions require offset log(exposure)
  xep1 = TRUE,
  xep01 = TRUE,
  xep2 = TRUE,
  xep02 = TRUE,
  xep12 = TRUE,
  xep012 = TRUE,
  xep123 = TRUE,
  xep0123 = TRUE,
  tnormal = FALSE,
  MaxwellBoltzmann = FALSE,
  lognormal = TRUE,
  xepi0 = TRUE,
  xep0 = TRUE,
  exponential = FALSE,
  chisq = FALSE,
  inverse_gaussian = FALSE,
  constant = TRUE
)
mod_offset <- c(
  xep1 = "log(exposure)",
  xep01 = "log(exposure)",
  xep2 = "log(exposure)",
  xep02 = "log(exposure)",
  xep12 = "log(exposure)",
  xep012 = "log(exposure)",
  xep123 = "log(exposure)",
  xep0123 = "log(exposure)",
  tnormal = "log(exposure) - log(r)",
  MaxwellBoltzmann = "log(exposure * r)",
  lognormal = "log(exposure)",
  xepi0 = "log(exposure)",
  xep0 = "log(exposure)",
  exponential = "log(exposure) - log(r)",
  chisq = "log(exposure) - r/2",
  inverse_gaussian = "log(exposure) - 2.5 * log(r)",
  constant = "log(exposure)"
)

constraints <- list(
  xep1 =  rbind(
    "(Intercept)" = c(lower = -Inf, upper = Inf, parscale = 5),
              "r" = c(lower = -Inf, upper = 0,   parscale = 0.1)),
  xep01 = rbind(
    "(Intercept)" = c(lower = -Inf, upper = Inf, parscale = 5),
         "log(r)" = c(lower = -2,   upper = Inf, parscale = 1),
              "r" = c(lower = -Inf, upper = 0,   parscale = 0.01)),
  xep2 = rbind(
    "(Intercept)" = c(lower = -Inf, upper = Inf, parscale = 5),
          "I(r^2)"= c(lower = -Inf, upper = 0,   parscale = 0.001)),
  xep02 = rbind(
    "(Intercept)" = c(lower = -Inf, upper = Inf, parscale = 5),
         "log(r)" = c(lower = -2,   upper = Inf, parscale = 1),
         "I(r^2)" = c(lower = -Inf, upper = 0,   parscale = 0.0001)),
  xep12 = rbind(
    "(Intercept)" = c(lower = -Inf, upper = Inf, parscale = 5),
              "r" = c(lower = -Inf, upper = Inf, parscale = 0.1),
          "I(r^2)"= c(lower = -Inf, upper = 0,   parscale = 0.001)),
  xep012 = rbind(
    "(Intercept)" = c(lower = -Inf, upper = Inf, parscale = 5),
         "log(r)" = c(lower = -2,   upper = Inf, parscale = 1),
              "r" = c(lower = -Inf, upper = Inf, parscale = 0.01),
         "I(r^2)" = c(lower = -Inf, upper = 0,   parscale = 0.0001)),
  xep123 = rbind(
    "(Intercept)" = c(lower = -Inf, upper = Inf, parscale = 5),
              "r" = c(lower = -Inf, upper = Inf, parscale = 0.01),
          "I(r^2)"= c(lower = -Inf, upper = Inf, parscale = 0.001),
          "I(r^3)"= c(lower = -Inf, upper = 0,   parscale = 0.00001)),
  xep0123 = rbind(
    "(Intercept)" = c(lower = -Inf, upper = Inf, parscale = 5),
         "log(r)" = c(lower = -2,   upper = Inf, parscale = 1),
              "r" = c(lower = -Inf, upper = Inf, parscale = 0.1),
          "I(r^2)"= c(lower = -Inf, upper = Inf, parscale = 0.001),
          "I(r^3)"= c(lower = -Inf, upper = 0,   parscale = 0.00001)),
  tnormal = rbind(
    "(Intercept)" = c(lower = -Inf, upper = Inf, parscale = 5),
              "r" = c(lower = -Inf, upper = Inf, parscale = 0.01),
          "I(r^2)"= c(lower = -Inf, upper = 0,   parscale = 0.0001)),
  MaxwellBoltzmann = rbind(
    "(Intercept)" = c(lower = -Inf, upper = Inf, parscale = 5),
          "I(r^2)"= c(lower = -Inf, upper = 0,   parscale = 0.001)),
  lognormal = rbind(
    "(Intercept)" = c(lower = -Inf, upper = Inf, parscale = 5),
         "log(r)" = c(lower = -Inf, upper = Inf, parscale = 1),
    "I(log(r)^2)" = c(lower = -Inf, upper = 0,   parscale = 0.1)),
  xepi0 = rbind(
    "(Intercept)" = c(lower = -Inf, upper = Inf, parscale = 5),
         "I(1/r)" = c(lower = -Inf, upper = 0,   parscale = 100),
         "log(r)" = c(lower = -Inf, upper = -2,  parscale = 1)),
  xep0 = rbind(
    "(Intercept)" = c(lower = -Inf, upper = Inf, parscale = 5),
          "log(r)"= c(lower = -Inf, upper = -2,   parscale = 1)),
  exponential = rbind(
    "(Intercept)" = c(lower = -Inf, upper = Inf, parscale = 5),
              "r" = c(lower = -Inf, upper = 0, parscale = 0.1)),
  chisq = rbind(
    "(Intercept)" = c(lower = -Inf, upper = Inf, parscale = 5),
          "log(r)"= c(lower =   -2, upper = Inf, parscale = 1)),
  inverse_gaussian =rbind(
    "(Intercept)" = c(lower = -Inf, upper = Inf, parscale = 5),
          "I(1/r)"= c(lower = -Inf, upper = 0,   parscale = 10),
              "r" = c(lower = -Inf, upper = 0,   parscale = 0.1)),
  constant = rbind(
    "(Intercept)" = c(lower =  Inf, upper = -Inf, parscale = 5))
)

constraints_par<- list( 
  xep1 = rbind(
    b1 = c(lower = -Inf, upper = 0)), 
  xep01 = rbind(
    shape = c(lower = 0, upper = Inf), 
    rate = c(0, Inf)), 
  xep2 = rbind(s2 = c(lower = 0, upper = Inf)),
  xep02 = rbind(
    b0 = c(lower = -2,   upper = Inf), 
    b2 =c(lower = -Inf, upper = 0)),
  xep12 = rbind(
    b1 = c(lower = -Inf, upper = Inf),
    b2 = c(lower = -Inf, upper = 0)),
  xep012 = rbind(
    b0 = c(lower = -2,   upper = Inf),
    b1 = c(lower = -Inf, upper = Inf),
    b2 = c(lower = -Inf, upper = 0)),
  xep123 = rbind(
    b1 = c(lower = -Inf, upper = Inf),
    b2 = c(lower = -Inf, upper = Inf),
    b3 = c(lower = -Inf, upper = 0)),
  xep0123 = rbind(
    b0 = c(lower = -2,   upper = Inf),
    b1 = c(lower = -Inf, upper = Inf),
    b2 = c(lower = -Inf, upper = Inf),
    b3 = c(lower = -Inf, upper = 0)),
  tnormal = rbind(
    mean = c(lower = -Inf, upper = Inf),
    sd = c(lower = 0, upper = Inf)),
  MaxwellBoltzmann = rbind(
    a = c(lower = 0, upper = Inf)),
  lognormal = rbind(
    meanlog = c(lower = -Inf, upper = Inf),
    sdlog = c(lower = 0, upper = Inf)),
  xepi0 = rbind(
    shape = c(lower = 0, upper = Inf),
    scale = c(lower = 0, upper = Inf)),
  xep0 = rbind(
    a = c(lower = 0, upper = Inf)),
  exponential = rbind(
    rate = c(lower = 0, upper = Inf)),
  chisq = rbind(
    df = c(lower = 0, upper = Inf)),
  inverse_gaussian =rbind(
    mean = c(lower = 0, upper = Inf),
    dispersion = c(lower = 0, upper = Inf)),
  constant = rbind(
    b1 = c(lower =  Inf, upper = -Inf))
)

par_default <- list(xlog = FALSE, ylog = FALSE, adj = 0.5, ann = TRUE, ask = FALSE,
  bg = "transparent", bty = "o", cex = 1, cex.axis = 1, cex.lab = 1,
  cex.main = 1.2, cex.sub = 1, col = "black", col.axis = "black",
  col.lab = "black", col.main = "black", col.sub = "black", crt = 0,
  err = 0, family = "", fg = "black", fin = c(6.999999, 6.999999),
  font = 1, font.axis = 1, font.lab = 1, font.main = 2, font.sub = 1,
  lab = c(5, 5, 7), las = 0, lend = "round", lheight = 1, ljoin = "round",
  lmitre = 10, lty = "solid", lwd = 1, mai = c(1.02, 0.82, 0.82, 0.42),
  mar = c(5.1, 4.1, 4.1, 2.1), mex = 1, mfcol = c(1, 1), mfg = c(1, 1, 1, 1),
  mfrow = c(1, 1), mgp = c(3, 1, 0), mkh = 0.001, new = FALSE, oma = c(0, 0, 0, 0),
  omd = c(0, 1, 0, 1), omi = c(0, 0, 0, 0), pch = 1,
  pin = c(5.759999, 5.159999), plt = c(0.1171429, 0.94, 0.1457143, 0.8828571),
  ps = 12, pty = "m", smo = 1, srt = 0, tck = NA, tcl = -0.5, usr = c(0, 1, 0, 1),
  xaxp = c(0, 1, 5), xaxs = "r", xaxt = "s", xpd = FALSE, yaxp = c(0, 1, 5),
  yaxs = "r", yaxt = "s", ylbias = 0.2)

# parameters for generating example carcass data sets
# these are not stored in the package database but are used to create data sets
# that are
dist_d = "gamma"
dparm = c(f50 = 0.6, f100 = 0.9)
dfbat <- suppressWarnings(optim(par = c(2, 50), fn = function(x){
   max(abs(pgamma(50, shape = x[1], scale=x[2]) - dparm["f50"]),
       abs(pgamma(100, shape = x[1], scale = x[2]) - dparm["f100"]))
})$par)
ncarc = 200
set.seed(20201111)   # 20201111 leads to awful data set
r <- rgamma(ncarc, shape = dfbat[1], scale = dfbat[2]) # distance
theta <- runif(ncarc) * 2 * pi # angle
# create carcass data frame for polygons:
turc <- sample(unique(layout_polygon$turbine), size = ncarc, replace = T)
# check whether the given distance was searched at the given turbine;
# if not, throw it out
tmpind <- numeric(ncarc)
for (ti in unique(turc)){
  i <- which(turc == ti)
  tmpind[i] <- splancs::inout(splancs::as.points(r[i] * cos(theta[i]), r[i] * sin(theta[i])),
    layout_polygon[layout_polygon$turbine == ti, ])
}
i <- which(tmpind == 1)
carcass_polygon <- data.frame(turbine = turc[i], r = round(r[i], 1),
  stringsAsFactors = FALSE)

# create carcass data for the simple geometry:
dist_d = "gamma"
dparm = c(f50 = 0.6, f100 = 0.9)
dfbat <- suppressWarnings(optim(par = c(2, 50), fn = function(x){
   max(abs(pgamma(50, shape = x[1], scale=x[2]) - dparm["f50"]),
       abs(pgamma(100, shape = x[1], scale = x[2]) - dparm["f100"]))
})$par)
ncarc = 100
set.seed(20200919)
r <- rgamma(ncarc, shape = dfbat[1], scale = dfbat[2]) # distance
theta <- runif(ncarc) * 2 * pi # angle
# create carcass data frame for polygons:
turc <- sample(unique(layout_simple$turbine), size = ncarc, replace = T)
carcass_simple0 <- data.frame(turbine = turc, r = r, theta = theta, found = 0)
rownames(layout_simple) <- layout_simple$turbine
carcass_simple0$x <- carcass_simple0$r * cos(carcass_simple0$theta)
carcass_simple0$y <- carcass_simple0$r * sin(carcass_simple0$theta)
t1found <- which(carcass_simple0$turbine == "t1" & 
  carcass_simple0$r <= layout_simple$radius[layout_simple$turbine == "t1"])
t2found <- which(carcass_simple0$turbine == "t2" & 
  abs(carcass_simple0$x) <= layout_simple$radius[layout_simple$turbine == "t2"] &
  abs(carcass_simple0$y) <= layout_simple$radius[layout_simple$turbine == "t2"]
)
tind <- which(layout_simple$turbine == "t3")
t3found <- which(carcass_simple0$turbine == "t3" & carcass_simple0$r <= layout_simple$radius[tind] & 
  (carcass_simple0$r <= layout_simple$padrad[tind] |
    (carcass_simple0$x >= 0 & abs(carcass_simple0$y) <= layout_simple$roadwidth[tind]/2) |
    (carcass_simple0$x <= 0 & abs(carcass_simple0$y) <= layout_simple$roadwidth[tind]/2)
  )
)

t4found <- which(carcass_simple0$turbine == "t4" &
  (carcass_simple0$r <= layout_simple$padrad[4] |
    (carcass_simple0$r <= layout_simple$radius[4] &
     carcass_simple0$x < 0 & abs(carcass_simple0$y) <= layout_simple$roadwidth[4]/2
    )
  )
)

carcass_simple0[c(t1found, t2found, t3found, t4found), "found"] <- 1
carcass_simple <- carcass_simple0[carcass_simple0$found == 1, c("turbine", "r")]
carcass_simple$r <- round(carcass_simple$r, 2)
carcass_simple <- carcass_simple[order(as.numeric(rownames(carcass_simple))), ]   
rownames(carcass_simple) <- 1:nrow(carcass_simple)

### xy data
# define carcass distribution
dist_d = "gamma"
dparm = c(f50 = 0.6, f100 = 0.9)
dfbat <- suppressWarnings(optim(par = c(2, 50), fn = function(x){
  max(abs(pgamma(50, shape = x[1], scale=x[2]) - dparm["f50"]),
      abs(pgamma(100, shape = x[1], scale = x[2]) - dparm["f100"]))
})$par)
ncarc = 100
set.seed(20201111)
# generate carcasses
r <- rgamma(ncarc, shape = dfbat[1], scale = dfbat[2])
theta <- runif(ncarc) * 2 * pi
xyr <- r * cbind(cos(theta), sin(theta), 1)
colnames(xyr) <- c("x", "y", "r")
# search plot parameters
radius = 75; padrad = 15; roadwidth = 5
# carcass searches
xyri <- xyr
xyri <- xyri[xyri[, "r"] <= radius + 0.5, ]
xyri <- xyri[xyri[, "r"] <= padrad | (xyri[, "x"] > 0 & abs(xyri[, "y"]) <= roadwidth/2), ]
# create layout
xygrid <- expand.grid(x = -radius:radius, y = -radius:radius)
xygrid <- xygrid[sqrt(rowSums(xygrid^2)) <= radius, ]
xygrid <- xygrid[sqrt(rowSums(xygrid^2)) <= padrad |
  (xygrid[, "x"] > 0 & abs(xygrid[, "y"]) <= roadwidth/2),]
xygrid <- as.matrix(xygrid)
xygrid <- cbind(xygrid, ncarc = 0)
xygrid <- cbind(xygrid, r = sqrt(rowSums(xygrid[, c("x", "y")]^2)))
xygrid <- xygrid[xygrid[, "r"] > 0.0001, ]
# assign carcasses to grid points
xyri[, "x"] <- round(xyri[, "x"])
xyri[, "y"] <- round(xyri[, "y"])
require(magrittr)
for (i in 1:nrow(xyri)){
  xygrid[which(xygrid[, "x"] == xyri[i, "x"] & xygrid[, "y"] == xyri[i, "y"]), "ncarc"] %<>% `+`(., 1)
}
xygrid[xygrid[, "x"] == -9 & xygrid[, "y"] == -12, "ncarc"] <- 1
xyri <- data.frame(xyri, stringsAsFactors = FALSE)
xyri$turbine <- "t1"
layout_xy <- as.data.frame(xygrid, stringsAsFactors = FALSE)
layout_xy$turbine <- "t1"
rownames(layout_xy) <- 1:nrow(layout_xy)

# eagle data
layout_eagle <- read.csv(textConnection('
"DateFound","turbine","r"
"03-Aug-05","t21",23
"10-Oct-05","t43",46
"31-Oct-05","t14",25
"09-Apr-06","t25",30
"28-Apr-06","t16",65
"03-May-06","t21",42
"04-May-06","t42",15
"05-May-06","t15",18
"01-Sep-06","t38",70
"02-May-06","t52",17
"09-May-07","t25",22
"13-Mar-08","t67",49
"03-Apr-08","t52",27
"16-Apr-08","t64",25
"16-Apr-08","t61",45
"22-Apr-08","t37",51
"30-Apr-08","t38",26
"25-Jun-08","t30",37
"12-Sep-08","t56",80
"01-Dec-08","t31",19
"04-Mar-09","t9",56
"11-Mar-09","t68",100
"08-Apr-09","t24",71
"28-Apr-09","t21",20
"01-May-09","t20",33
"07-Oct-09","t6",35
"26-Oct-09","t47",46
"29-Aug-12","t48",40
"28-Feb-14","t21",50
"25-Mar-14","t57",20
"09-Apr-14","t64",30
"28-Mar-15","t52",70
"06-Jun-15","t66",56
"24-Aug-15","t49",20
"08-Mar-16","t50",50
"22-Mar-16","t13",40
"05-Apr-16","t11",12
"19-Apr-16","t67",44
"23-Apr-16","t47",50
"01-May-16","t38",49
"02-May-16","t50",25
"06-May-16","t59",40
"21-May-16","t41",81
"05-Oct-16","t61",77
"16-Oct-16","t30",18
"26-Oct-16","t7",36
"10-Apr-18","t6",50
"19-Sep-17","t45",75
"08-Mar-18","t63",70
"16-Mar-18","t53",20
"28-Mar-18","t63",50
"08-Apr-18","t27",20
"29-Apr-18","t54",35
"15-Jun-18","t51",29
"15-Jun-18","t39",46
"15-Jun-18","t61",58
"15-Oct-18","t68",63
"27-Sep-18","t62",30
"16-Mar-19","t27",50
"23-Apr-19","t29",45
'), as.is = T)

sieve_default <- list(
  aic = 10,
  hin = TRUE,
  rtail = c(p80 = 0.50, p120 = 0.90, p150 = 0.95, p200 = 0.99),
  ltail = c(p20 = 0.50, p50 = 0.90)
)
sieve_win <- list(
  aic = 10,
  hin = TRUE,
  rtail = sieve_default$rtail * 0,
  ltail = sieve_default$ltail * 0 + 1
)

usethis::use_data(
  alt_names,
  carcass_polygon,
  carcass_simple,
  carcass_simple0,
  cof_name,
  constraints,
  constraints_par,
  degOrder,
  distr_names,
  layout_eagle,
  layout_polygon,
  layout_simple,
  layout_xy,
  mod_all,
  mod_color,
  mod_lty,
  mod_offset,
  mod_standard,
  mod_xy,
  natural,
  par_default,
  parm_name,
  sieve_default,
  sieve_win,
  xyr, 
  internal = FALSE, overwrite = TRUE)

