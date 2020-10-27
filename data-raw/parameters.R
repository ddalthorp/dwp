mod_color = c(
  gamma = 2,
  logLinear = colors()[97], #darkorchid2 #colors()[35], #brown3
  logQuadratic = 4,
  Rayleigh = colors()[123], #deepskyblue2 [dashed]
  logCubic = 3,
  paranormal_gamma = colors()[657], #yellowgreen
  lognormal = 1,
  inverse_gamma = colors()[148],  #goldenrod1
  Pareto = 7,
  chisq = 1,
  exponential = colors()[35], #brown3
  inverse_gaussian = 2,
  tnormal = 4,
  MaxwellBoltzmann = colors()[123], #deepskyblue2 [dashed]
  constant = 8
)
mod_name = names(mod_color)
parm_name <- list(
  gamma = c("shape", "rate"),
  lognormal = c("meanlog", "sdlog"),
  logLinear = "b1",
  logQuadratic = c("b1", "b2"),
  logCubic = c("b1", "b2", "b3"),
  inverse_gamma = c("shape", "scale"),
  paranormal_gamma = c("b0", "b1", "b2", "b3"),
  Rayleigh = "s2",
  MaxwellBoltzmann = "a",
  Pareto = "a",
  constant = NULL,
  tnormal = c("mean", "sd"),
  exponential = "rate",
  chisq = "df",
  inverse_gaussian = c("mean", "dispersion")
)
cof_name <- list(
  gamma = c("(Intercept)", "log(r)", "r"),
  lognormal = c("(Intercept)", "log(r)", "I(log(r)^2)"),
  logLinear = c("(Intercept)", "r"),
  logQuadratic = c("(Intercept)", "r", "I(r^2)"),
  logCubic = c("(Intercept)","r", "I(r^2)", "I(r^3)"),
  inverse_gamma = c("(Intercept)", "I(1/r)", "log(r)"),
  paranormal_gamma = c("(Intercept)", "log(r)", "r", "I(r^2)", "I(r^3)"),
  Rayleigh = c("(Intercept)", "I(r^2)"),
  MaxwellBoltzmann = c("(Intercept)", "I(r^2)"),
  Pareto = c("(Intercept)", "log(r)"),
  constant = c("(Intercept)"),
  tnormal = c("(Intercept)", "r", "I(r^2)"),
  exponential = c("(Intercept)", "r"),
  chisq = c("(Intercept)", "log(r)"),
  inverse_gaussian = c("(Intercept)", "I(1/r)", "r")
)
natural = c( # distributions that require offset other than log(exposure)
  gamma = TRUE,
  lognormal = TRUE,
  logLinear = TRUE,
  logQuadratic = TRUE,
  logCubic = TRUE,
  inverse_gamma = TRUE,
  paranormal_gamma = TRUE,
  Rayleigh = TRUE,
  Pareto = TRUE,
  MaxwellBoltzmann = FALSE,
  constant = TRUE,
  tnormal = FALSE,
  exponential = FALSE,
  chisq = FALSE,
  inverse_gaussian = FALSE
)
mod_offset <- c(
  gamma = "log(exposure)",
  lognormal = "log(exposure)",
  logLinear = "log(exposure)",
  logQuadratic = "log(exposure)",
  logCubic = "log(exposure)",
  inverse_gamma = "log(exposure)",
  paranormal_gamma = "log(exposure)",
  Rayleigh = "log(exposure)",
  MaxwellBoltzmann = "log(exposure * r)",
  Pareto = "log(exposure)",
  constant = "log(exposure)",
  tnormal = "log(exposure) - log(r)",
  exponential = "log(exposure) - log(r)",
  chisq = "log(exposure) - r/2",
  inverse_gaussian = "log(exposure) - 2.5 * log(r)"
)
constraints <- list(
  gamma = rbind(
    "(Intercept)" = c(lower = -Inf, upper = Inf, parscale = 5),
         "log(r)" = c(lower = -2,   upper = Inf, parscale = 1),
              "r" = c(lower = -Inf, upper = 0,   parscale = 0.01)),
  lognormal = rbind(
    "(Intercept)" = c(lower = -Inf, upper = Inf, parscale = 5),
         "log(r)" = c(lower = -Inf, upper = Inf, parscale = 1),
    "I(log(r)^2)" = c(lower = -Inf, upper = 0,   parscale = 0.1)),
  logLinear =  rbind(
    "(Intercept)" = c(lower = -Inf, upper = Inf, parscale = 5),
              "r" = c(lower = -Inf, upper = 0,   parscale = 0.1)),
  logQuadratic = rbind(
    "(Intercept)" = c(lower = -Inf, upper = Inf, parscale = 5),
              "r" = c(lower = -Inf, upper = Inf, parscale = 0.1),
          "I(r^2)"= c(lower = -Inf, upper = 0,   parscale = 0.001)),
  logCubic = rbind(
    "(Intercept)" = c(lower = -Inf, upper = Inf, parscale = 5),
              "r" = c(lower = -Inf, upper = Inf, parscale = 0.01),
          "I(r^2)"= c(lower = -Inf, upper = Inf, parscale = 0.001),
          "I(r^3)"= c(lower = -Inf, upper = 0,   parscale = 0.00001)),
  inverse_gamma = rbind(
    "(Intercept)" = c(lower = -Inf, upper = Inf, parscale = 5),
         "I(1/r)" = c(lower = -Inf, upper = 0,   parscale = 100),
         "log(r)" = c(lower = -Inf, upper = -2,  parscale = 1)),
  paranormal_gamma = rbind(
    "(Intercept)" = c(lower = -Inf, upper = Inf, parscale = 5),
         "log(r)" = c(lower = -2,   upper = Inf, parscale = 1),
              "r" = c(lower = -Inf, upper = Inf, parscale = 0.1),
          "I(r^2)"= c(lower = -Inf, upper = Inf, parscale = 0.001),
          "I(r^3)"= c(lower = -Inf, upper = 0,   parscale = 0.00001)),
  Rayleigh = rbind(
    "(Intercept)" = c(lower = -Inf, upper = Inf, parscale = 5),
          "I(r^2)"= c(lower = -Inf, upper = 0,   parscale = 0.001)),
  MaxwellBoltzmann = rbind(
    "(Intercept)" = c(lower = -Inf, upper = Inf, parscale = 5),
          "I(r^2)"= c(lower = -Inf, upper = 0,   parscale = 0.001)),
  Pareto = rbind(
    "(Intercept)" = c(lower = -Inf, upper = Inf, parscale = 5),
          "log(r)"= c(lower = -Inf, upper = -2,   parscale = 1)),
  constant = rbind(
    "(Intercept)" = c(lower =  Inf, upper = -Inf, parscale = 5)),
  tnormal = rbind(
    "(Intercept)" = c(lower = -Inf, upper = Inf, parscale = 5),
              "r" = c(lower = -Inf, upper = Inf, parscale = 0.01),
          "I(r^2)"= c(lower = -Inf, upper = 0,   parscale = 0.0001)),
  exponential = rbind(
    "(Intercept)" = c(lower = -Inf, upper = Inf, parscale = 5),
              "r" = c(lower = -Inf, upper = 0, parscale = 0.1)),
  chisq = rbind(
    "(Intercept)" = c(lower = -Inf, upper = Inf, parscale = 5),
          "log(r)"= c(lower =   -2, upper = Inf, parscale = 1)),
  inverse_gaussian =rbind(
    "(Intercept)" = c(lower = -Inf, upper = Inf, parscale = 5),
          "I(1/r)"= c(lower = -Inf, upper = 0,   parscale = 10),
              "r" = c(lower = -Inf, upper = 0,   parscale = 0.1))
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

usethis::use_data(
  cof_name,
  constraints,
  mod_color,
  mod_name,
  mod_offset,
  parm_name,
  natural,
  par_default,
  internal = FALSE, overwrite = TRUE
)
