library(tidyverse)
library(magrittr)
library(HiddenSemiMarkov)

J = 5
state_names = c("iners_dominant", "iners_permissive", "limbo","non_lacto_dominant", "iners_recovering")
state_colors = c("green3","orange","hotpink","red", "steelblue")
init = rep(1/J, J)
transitions = matrix(0, J, J) %>% set_colnames(state_names) %>% set_rownames(state_names)
transitions["iners_dominant", "iners_permissive"] = 1
transitions["iners_permissive", "limbo"] = 1/2
transitions["iners_permissive", "non_lacto_dominant"] = 1/4
transitions["iners_permissive", "iners_recovering"] = 1/4
transitions["limbo", "non_lacto_dominant"] = 1/2
transitions["limbo", "iners_recovering"] = 1/2
transitions["non_lacto_dominant", "limbo"] = 1/2
transitions["non_lacto_dominant", "iners_recovering"] = 1/2
transitions["iners_recovering", "iners_dominant"] = 1
# sojourn
x = 1:(40*7)
d = matrix(0, nrow = length(x), ncol = J) %>% set_colnames(state_names)
d[,"iners_dominant"] = c(rep(0,4*7), rep(1, 36*7))
d[,"iners_permissive"] = dnorm(x, mean = 5*7, sd = 7)
d[,"limbo"] = c(rep(0,4*7), rep(1, 36*7))
d[,"non_lacto_dominant"] = c(rep(0,4*7), rep(1, 36*7))
d[,"iners_recovering"] = dnorm(x, mean = 5*7, sd = 7)
d = t(t(d)/colSums(d))
sojourn = list(type = "ksmoothed-nonparametric",
               d = d)

# emission par
marg_em_dist = list(
  iners.prop = list(
    type = "beta",
    params = list(
      shape1 = c(10, 1, 2,  1, 1),
      shape2 = c(1 , 1, 2, 10, 1)
    )
  ),
  d.iners.prop = list(
    type = "norm",
    params = list(
      mean = c(0,-0.2, 0, 0, 0.2),
      sd = c(0.1, 0.2, 0.2, 0.1, 0.2)
    )
  ),
  non.lactobacillus.prop = list(
    type = "beta",
    params = list(
      shape1 = c(1 , 1, 2, 10, 1),
      shape2 = c(10, 1, 2,  1, 1)
    )
  ),
  d.non.lactobacillus.prop = list(
    type = "norm",
    params = list(
      mean = c(0, 0.2, 0, 0, -0.2),
      sd = c(0.1, 0.2, 0.2, 0.1, 0.2)
    )
  )
)

sample_cat_hsmm = specify_hsmm(J = J,
                               state_names = state_names,
                               state_colors = state_colors,
                               init = init, transition = transitions,
                               sojourn = sojourn,
                               marg_em_probs = marg_em_dist)

X = simulate_hsmm(model = sample_cat_hsmm, n_state_transitions = 100)
plot_hsmm_seq(X = X, model = sample_cat_hsmm)


fitted_model_1 = fit_hsmm(model = sample_cat_hsmm, X = X, use_sojourn_prior = TRUE)

fitted_model_2 = fit_hsmm(model = sample_cat_hsmm, X = X, use_sojourn_prior = FALSE)

plot_hsmm_fit_param(fitted_model_1)
plot_hsmm_fit_param(fitted_model_2)



plot_hsmm_sojourn_dist(model = fitted_model_1$model, maxt = 300) + facet_grid(state ~ . , scale = "free")
plot_hsmm_sojourn_dist(model = fitted_model_2$model, maxt = 300) + facet_grid(state ~ . , scale = "free")



plot_hsmm_marg_dist(model = sample_cat_hsmm)


