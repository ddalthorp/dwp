mod_color = c(
  xep1 = colors()[35], #brown3 colors()[585], #darkorchid2 #
  xep01 = 2,
  xep2 = colors()[50],
  xep02 = 1,
  xep12 = colors()[122],
  xep012 = 4,
  xep123 = colors()[148],
  xep0123 = colors()[657], #yellowgreen
  tnormal = 4,
  MaxwellBoltzmann = colors()[257], #deepskyblue2 [dashed]
  lognormal = colors()[96],
  xepi0 = 5,
  xep0 = colors()[369],
  chisq = 1,
  exponential = colors()[123], #brown3
  inverse_gaussian = 2,
  constant = 8
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
  mod_all,
  mod_color,
  mod_offset,
  mod_standard,
  mod_xy,
  natural,
  par_default,
  parm_name,
  internal = FALSE, overwrite = TRUE
)
