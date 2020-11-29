#' @importFrom graphics axis box legend lines mtext par points polygon rect text
#' @importFrom grDevices colors
#' @importFrom magrittr %>% %<>%
#' @importFrom methods as is
#' @importFrom stats aggregate approxfun dbinom dchisq dexp dgamma dlnorm dnorm
#'  formula glm integrate pbinom pchisq pexp pgamma plnorm pnorm quantile rnorm
#'  runif uniroot var
#' @importFrom utils flush.console read.csv write.csv
utils::globalVariables(c(".", "cof_name", "constraints", "mod_all", "mod_color",
 "mod_name", "mod_offset", "mod_standard", "mod_xy", "natural", "par_default",
 "parm_name", "sieve_default"))

