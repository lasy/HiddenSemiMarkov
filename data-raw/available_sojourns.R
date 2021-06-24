available_sojourn = rbind(
  data.frame(distribution_type = "nonparametric",
             parameters = c("d"),
             required_parameter = TRUE,
             parameters_description = "non-parametric probability density",
             stringsAsFactors = FALSE),
  data.frame(distribution_type = "ksmoothed_nonparametric",
             parameters = c("d", "bw"),
             required_parameter = c(TRUE,FALSE),
             parameters_description =
               c(
                 "non-parametric probability density",
                 "smoothing bandwidth. See '?density' for default value."
               ),
             stringsAsFactors = FALSE),
  data.frame(distribution_type = "gamma",
             parameters = c("shape", "scale"),
             required_parameter =  c(TRUE,TRUE),
             parameters_description =
               c(
                 "shape of the gamma distribution",
                 "scale of the gamma distribution (rate = 1/scale)"
                 ),
             stringsAsFactors = FALSE),
  data.frame(distribution_type = "poisson",
             parameters = c("lambda", "shift"),
             required_parameter =  c(TRUE,FALSE),
             parameters_description =
               c("see '?dpois'",
                 "default = 0. Additional parameter allowing to shift the Poisson distribution."
                 ),
             stringsAsFactors = FALSE),
  data.frame(distribution_type = "lnorm",
             parameters = c("meanlog", "sdlog"),
             required_parameter =  c(TRUE,TRUE),
             parameters_description =
               c("mean of the log-normal distribution (see '?rlnorm')",
                 "standard deviation of the log-normal distribution (see '?rlnorm')"
               ),
             stringsAsFactors = FALSE),
  data.frame(distribution_type = "logarithmic",
             parameters = "shape",
             required_parameter =  c(TRUE),
             parameters_description = "decay parameter of the logarithmic distribution",
             stringsAsFactors = FALSE),
  data.frame(distribution_type = "nbinom",
             parameters = c("size", "prob", "shift"),
             required_parameter = c(TRUE, TRUE, FALSE),
             parameters_description = c(
               "size parameter of the negative binomial distribution (see ?dnbinom)",
               "probability of the negative binomial distribution (see ?dnbinom)",
               "default = 0. Additional parameter allowing to shift the distribution."
               ),
             stringsAsFactors = FALSE),
  data.frame(distribution_type = "geometric",
             parameters = "prob",
             required_parameter = TRUE,
             parameters_description = "see `rgeom` for the 'prob' parameter",
             stringsAsFactors = FALSE)
)
available_sojourn
