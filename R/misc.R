

# generate random observations for all variables
# generate_random_obs = function(n = 10, parem, state){
#   r_obs = sapply(1:length(parem), function(i) generate_random_obs_par(n, parem[[i]],state))
#   r_obs = matrix(r_obs, nrow = n) # in case we generate for only 1 observation
#   colnames(r_obs) = names(parem)
#   r_obs = r_obs %>% as.data.frame()
#   r_obs = .check_types_and_values(data = r_obs, parem = parem)
#   return(r_obs)
# }


# generate random observations for one variable
# generate_random_obs_par = function(n, par, state){
#   if(par$type == "norm"){
#     x = rnorm(n, mean = par$param$mean[state], sd = par$param$sd[state])
#   }else if(par$type == "binom"){
#     x = rbinom(n, size = par$param$size[state], prob = par$param$prob[state])
#   }else if(par$type == "non-par"){
#     x = sample(rep(par$param$values, round(par$param$probs[,state]*5000)) , n, replace = TRUE)
#   }else{stop("This distributions hasn't been implemented yet")}
#   return(x)
# }


# impute
# impute_missing_data_at_random_given_state = function(model, X, state, seed = 1){
#   set.seed(seed)
#   n_missing_obs =  X %>% select(all_of(names(model$parms.emission))) %>% is.na() %>% colSums()
#   if(any(n_missing_obs > 0)){
#     for(var in names(model$parms.emission)[n_missing_obs > 0]){
#       j = which(is.na(X[, var]))
#       X[j, var] = generate_random_obs_par(n = n_missing_obs[var],state = state, par = model$parms.emission[[var]])
#     }
#   }
#   set.seed(Sys.time())
#   X
# }

# impute_missing_data_with_most_likely_given_state = function(model, X, state){
#   n_missing_obs =  X %>% select(all_of(names(model$parms.emission))) %>% is.na() %>% colSums()
#   if(any(n_missing_obs > 0)){
#     for(var in names(model$parms.emission)[n_missing_obs > 0]){
#       j = which(is.na(X[, var]))
#       val = most_probable_value(par = model$parms.emission[[var]], state = state)
#       X[j, var] = rep(val, n_missing_obs[var])
#     }
#   }
#   X
# }




# generate_missingness

# generate_missingness = function(n = 10, parem, state){
#   r_obs = sapply(1:length(parem), function(v) rbinom(n = n, size = 1, prob = parem[[v]]$missing_prob[state]))
#   r_obs = matrix(r_obs, nrow = n)
#   colnames(r_obs) = names(parem)
#   r_obs = r_obs %>% as.data.frame()
#   return(r_obs)
# }
#
#
# generate_censoring = function(n = 10, prob, state, parem){
#   mc = rbinom(n = n, size = 1, prob = prob[state])
#   mc = matrix(mc, nrow = n, ncol = length(parem))
#   colnames(mc) = names(parem)
#   mc = as.data.frame(mc)
#   mc
# }

# state_names = model$state_names
# state_cols = model$state_colors
#
# transition_mat = model$transition
# colnames(transition_mat) = state_names; rownames(transition_mat) = state_names
#
# g = graph_from_adjacency_matrix(transition_mat, weighted = TRUE, mode = "directed")
# if(is.null(layout)) layout = layout_nicely(g)
#
# #E(g)$width <- 1+E(g)$weight*5
# E(g)$width <- E(g)$weight
# plot(g,
#      vertex.size=50, vertex.frame.color="white",
#      vertex.label.color="white", vertex.label.family  = "sans",
#      vertex.color = state_cols,
#      layout = layout
# )




# if(model$marg_em_probs[[var]]$type == "non-par"){
#   X_var = X_var %>%  dplyr::mutate(y_num = y %>% as.numeric())
#   y_breaks = 1:length(model$marg_em_probs[[var]]$params$values)
#   y_limits = c(0,max(y_breaks))
#   y_labels = model$marg_em_probs[[var]]$params$values
#   y_ref = 0
#   add_line = TRUE
#   add_point = FALSE
#
#   rel_heights = c(rel_heights, 1+sqrt(length(model$marg_em_probs[[var]]$params$values)))
# }else if(model$marg_em_probs[[var]]$type == "norm"){
#   X_var = X_var %>%  dplyr::mutate(y_num = y)
#   y_limits = range(X_var$y_num); y_limits = y_limits+c(-1,1)*max(abs(y_limits))*0.01
#   y_breaks = seq(-10, 10, by = 1)
#   y_labels = y_breaks
#   y_ref = 0
#   add_line = TRUE
#   add_point = FALSE
#   rel_heights = c(rel_heights, 5)
# }else if(model$marg_em_probs[[var]]$type == "binom"){
#   X_var = X_var %>%  dplyr::mutate(y_num = y)
#   y_breaks = 0:max(model$marg_em_probs[[var]]$params$size)
#   y_limits = range(y_breaks)
#   y_labels = y_breaks
#   y_ref = 0
#   add_line = FALSE
#   add_point = TRUE
#   rel_heights = c(rel_heights, pmax(2,length(y_breaks)))
# }
#
# g_var = g_base
# if(add_line) g_var = g_var + geom_segment(data = X_var, aes(xend = t, y = y_num, yend = y_ref, col = y_num))
# if(add_point) g_var = g_var + geom_point(data = X_var, aes(y = y_num, col = y_num))
# g_var = g_var +
#   scale_y_continuous(breaks = y_breaks, labels = y_labels, minor_breaks = NULL, limits = y_limits)+
#   scale_color_gradient2(low = "blue", high = "red", mid = "gray90", midpoint = 0)+
#   scale_fill_identity()+
#   guides(col = FALSE)+
#   ggtitle(var)





# .re_estimate_emission_parameters = function(model, X, p, ground_truth = ground_truth, trust_in_ground_truth = 0.75, N0 = 200, verbose = FALSE){
#
#   # initialization of the weights matrix
#   weights = p
#
#   #### First, we need to combine the ground_truth with the results from the E-step
#   weights_ground_truth = weights
#   j = which(!is.na(ground_truth$state))
#   if(length(j)>0){
#     weights_ground_truth[j,] = 0
#     weights_ground_truth[cbind(j,ground_truth$state[j])] = 1
#   }
#
#   weights = (1-trust_in_ground_truth) * weights + trust_in_ground_truth * weights_ground_truth
#
#   state_seq = apply(weights,1,which.max)
#
#   # Check for un-visited states or rarely visited states
#   if(verbose){
#     if(any(table(state_seq)<(nrow(X)/model$J/10))){warning("Some states are rarely visited\n")}
#     if(any(table(state_seq)==0)){warning("Some states are never visited\n")}
#   }
#
#   #### Then we can update the model parameters
#   new_model = model
#
#   #### The global_censoring_probability needs to be updates
#   new_model$global_censoring_prob = .re_estimate_global_censoring_prob(model = model, X = X, w = weights, N0 = N0)
#
#   #### So do the marginal emission distributions and the marginal missingness probabilities
#   new_model = model
#   new_model$marg_em_probs = .re_estimate_emission_distributions(model = model, X = X, w = weights)
#
#   #### Finally, (and most importantly) we re-estimate the joint observation probabilities
#   new_model$b = .re_estimate_b(model = model, X = X, w = weights, N0 = N0, verbose = verbose)
#
#   ##### Return the new model with updated emission parameters
#   new_model
# }
#
# .re_estimate_global_censoring_prob = function(model = model, X = X, w = w, N0 = N0){
#   # first we get the observed censoring frequency
#   all_nas = X %>% dplyr::select(all_of(names(model$marg_em_probs))) %>% apply(.,1, function(x) all(is.na(x)))
#   observed_counts = colSums(all_nas * w)
#   updated_prob = (model$global_censoring_prob * N0 + observed_counts)/(N0 + nrow(X))
#   updated_prob
# }
#
# .re_estimate_emission_distributions = function(model = model, X = X, w = w){
#
#   # keeping only the observations (= variables for which we have marg_em_probs)
#   obs_names = names(model$marg_em_probs)
#   obs = X %>% dplyr::select(dplyr::all_of(obs_names)) %>% as.data.frame()
#
#   new_parem = .compute_new_em_parms(obs = obs, w = w, parem = model$marg_em_probs)
#   new_parem
# }
#
#
# #library(SDMTools)
#
# .compute_new_em_parms = function(obs, w, parem){
#   new.parem = parem
#   for(i in 1:length(parem)){
#     #cat(i, "\n")
#     if(new.parem[[i]]$type == "norm"){
#       # weighted mean
#       x_bar = sapply(1:ncol(w),function(state) weighted.mean(obs[,i],w = w[,state], na.rm = TRUE))
#       x_bar[is.na(x_bar)] = parem[[i]]$params$mean[is.na(x_bar)]
#       # weighted sd
#       # sd_bar = sapply(1:ncol(w),function(state) wt.sd(obs[,i],w = w[,state]))
#       sd_bar = sapply(1:ncol(w), function(state) sqrt( sum( w[,state]  * (obs[,i] - x_bar[state])^2)))
#       sd_bar[is.na(sd_bar)] = parem[[i]]$params$sd[is.na(sd_bar)]
#       # variables and hyper-parameters
#       n = apply(w,2,sum)
#       mu_0 = parem[[i]]$params$mu_0
#       n0 = parem[[i]]$params$n0
#       new_n0 = n + n0
#       new_alpha = parem[[i]]$params$alpha + n/2
#       new_beta = parem[[i]]$params$beta + (n-1)/2 * sd_bar^2 + n * n0 / (n + n0) * (x_bar - mu_0)^2 / 2
#
#       # posteriors
#       new_mean = (n0 * mu_0 + n*x_bar)/(n0 + n)
#       new_sd = sqrt( new_beta * (new_n0 + 1) / (new_n0 * new_alpha) )
#
#       # updating the model
#       new.parem[[i]]$params$mean = new_mean
#       new.parem[[i]]$params$sd = new_sd
#     }else if(new.parem[[i]]$type == "binom"){
#
#       # observed proportions
#       obs_prob = sapply(1:ncol(w),function(state) weighted.mean(obs[,i]/(parem[[i]]$param$size[state]),w = w[,state], na.rm = TRUE)) # would be better with an EM approach, but it's good enough for now
#       obs_prob[is.na(obs_prob)] = parem[[i]]$params$prob[is.na(obs_prob)]
#       n = apply(w,2,sum)
#       obs_successes = n*obs_prob
#       # hyper-parameters
#       new_alpha = parem[[i]]$params$alpha + obs_successes
#       new_beta = parem[[i]]$params$beta + n * parem[[i]]$param$size  - obs_successes
#
#       # posterior
#       new_prob = new_alpha / (new_alpha + new_beta)
#
#       # update the model
#       new.parem[[i]]$params$prob = new_prob
#
#     }else if(new.parem[[i]]$type == "non-par"){
#       # observed proportions
#       obs_probs = sapply(
#         1:ncol(w), # for each state
#         function(state) {
#           tt = parem[[i]]$params$probs[,state]
#           j = which(w[,state]>0)
#           if(length(j)>0) tt = table(sample(obs[j,i], prob = w[j,state], size = nrow(obs), replace = TRUE))
#           tt = tt/sum(tt)
#           return(t(tt))
#         }
#       )
#
#       # variables
#       n = apply(w,2,sum)
#       n0 = parem[[i]]$params$n0
#       probs_0 = parem[[i]]$params$probs_0
#       # posterior
#       new_probs = n0 * probs_0 + t(n * t(obs_probs))
#       new_probs = t(t(new_probs)/colSums(new_probs))
#
#       # updating the model
#       new.parem[[i]]$params$probs = new_probs
#     }else{ stop("This type of distribution has not been handled yet")}
#
#     # we also need to impute the missingness of this variable
#     # TODO
#     new_missing_prob = sapply(1:ncol(w), function(state) weighted.mean(is.na(obs[,i]), w = w[,state], na.rm = FALSE))
#     new_missing_prob[is.na(new_missing_prob)] = parem[[i]]$missing_prob[is.na(new_missing_prob)]
#     new.parem[[i]]$missing_prob = new_missing_prob # %>% pmin(.,0.9999) %>% pmax(.,0.0001)
#   }
#   return(new.parem)
# }
#
#
# .re_estimate_b = function(model, X, w, N0 = 200, verbose = FALSE){
#
#   if(verbose) cat("Re-estimating b\n")
#
#   all_levels = .get_all_possible_levels(model = model, with_missing = TRUE, continuous_var_binned = TRUE)
#   var_names = colnames(all_levels)
#   Xb = .cut_continuous_var(model = model, X = X)
#   Xb_with_w = Xb %>%
#     dplyr::select(seq_id, t, all_of(var_names)) %>% # we only keep the seq_id, the time-points and the observations
#     dplyr::full_join(., # we join with the weight, but we need to transform them to long format first
#                      w %>% as.data.frame() %>%
#                        magrittr::set_colnames(1:model$J) %>%
#                        dplyr::mutate(seq_id = X$seq_id, t = X$t) %>%
#                        tidyr::pivot_longer(col = c(-seq_id, -t), names_to = "state", values_to = "p") %>%
#                        dplyr::mutate(state = as.integer(state)),
#                      by = c("seq_id","t"))
#
#   # to re-estimate b, first we need to get the contingency tables for values and for missing.
#
#   # let's first do the missingness: update prob_missing
#   # we transform the observations as logical specifying if the observations were NA or not
#   observed_missing = Xb_with_w %>% dplyr::mutate(dplyr::across(tidyselect::all_of(var_names), is.na))
#
#   for(obs_comb in names(model$prob_missing)){
#     obs_names = stringr::str_split(obs_comb, "_") %>% unlist()
#     contingency_missing = observed_missing %>%
#       dplyr::select(tidyselect::all_of(obs_names), state, p) %>%
#       dplyr::group_by(.dots = c("state", obs_names)) %>%
#       dplyr::summarise(counts = sum(p), .groups = "drop") %>%
#       dplyr::group_by(state) %>%
#       dplyr::mutate(N = sum(counts))
#     updated_prob_missing =
#       dplyr::right_join(contingency_missing,
#                         model$prob_missing[[obs_comb]] %>% dplyr::select(state, tidyselect::all_of(obs_names), p0),
#                         by = c("state", obs_names)) %>%
#       dplyr::mutate(p = (p0*N0 + counts)/(N0+N)) %>%
#       dplyr::ungroup()
#     model$prob_missing[[obs_comb]] = updated_prob_missing
#   }
#
#   # now we can do the values: prob_values
#   for(obs_comb in names(model$prob_values)){
#     obs_names = stringr::str_split(obs_comb, "_") %>% unlist()
#     contingency_values = Xb_with_w %>%
#       dplyr::select(tidyselect::all_of(obs_names), state, p) %>%
#       dplyr::mutate(has_NA = observed_missing %>% dplyr::select(tidyselect::all_of(obs_names)) %>% apply(.,1,any)) %>%
#       dplyr::filter(!has_NA) %>% dplyr::select(-has_NA) %>%
#       dplyr::group_by(.dots = c("state", obs_names)) %>%
#       dplyr::summarise(counts = sum(p), .groups = "drop")
#
#     updated_prob_values =
#       dplyr::left_join(
#         model$prob_values[[obs_comb]] %>% dplyr::select(state, tidyselect::all_of(obs_names), p0),
#         contingency_values,
#         by = c("state", obs_names)) %>%
#       dplyr::mutate(counts = counts %>% tidyr::replace_na(0)) %>%
#       dplyr::group_by(state) %>%
#       dplyr::mutate(N = sum(counts)) %>%
#       dplyr::ungroup() %>%
#       dplyr::mutate(p = (p0*N0 + counts)/(N0+N))
#     model$prob_values[[obs_comb]] = updated_prob_values
#   }
#
#   # from these contingency tables, and updated emission and missing probabilities, we can update b
#   b = .compute_b_from_prob_values_and_prob_missing(model)
#   b
# }
#
#
#
# .get_value_prob = function(model, s, non_missing_vars, i, b){
#   p = 1
#   if(non_missing_vars != "")
#     values = b[i,stringr::str_split(non_missing_vars,"_") %>% unlist()]
#   p = dplyr::left_join(values,
#                        model$prob_values[[non_missing_vars]] %>% dplyr::filter(state == s),
#                        by = colnames(values)) %>% dplyr::select(p) %>% unlist()
#   p
# }
#
# .get_missing_prob = function(model, s, missing_vars){
#   p = 1
#   if(missing_vars != "")
#     p = model$prob_missing[[missing_vars]] %>% dplyr::filter(state == s) %>% dplyr::select(p) %>% unlist()
#   p
# }







#' ####### INITIALIZATION ---------------------
#'
#'
#' #' Initialize a hidden semi-Markov model.
#' #'
#' #' This function initializes the joint emission probabilities given the states,
#' #' i.e. determines the values of \eqn{Pr(X_i, E_i | S_j)}  where \eqn{X_i} are the observations at a given time-point \eqn{i}, \eqn{E_i} are the values of the augmented variables at time-point \eqn{i} and \eqn{S_j} is the state \eqn{j}.
#' #'
#' #' @param model a \code{hsmm_spec} object. A specified hidden semi-Markov model. Use the function \code{specify_hsmm()} to specify a hidden semi-Markov model.
#' #' @param nmin (optional) a strictly positive integer. XXXX
#' #' @param seq_sim_seed (optional) an integer. XXX
#' #' @param verbose (optional) a logical. If TRUE, the function prints the internal steps of the function.
#' #'
#' #' @return an object of class \code{hsmm}, which can be used to decode observation sequences with the function \code{predict_states_hsmm()} or which can be fitted to observations with the function \code{fit_hsmm()}.
#' #'
#' #' @importFrom magrittr %>%
#' #' @export
#' #' @examples
#' #' my_model_spec = simple_model_spec
#' #' class(my_model_spec)
#' #' my_model_init = initialize_hsmm(my_model_spec)
#' #' class(my_model_init)
#' #'
#' initialize_hsmm = function(model,
#'                            nmin = NULL,
#'                            verbose = FALSE,
#'                            seq_sim_seed = NULL){
#'
#'   # CHECKS
#'   if(!(class(model) %in% c("hsmm_spec","hsmm"))) stop("model must be of class 'hsmm_spec'. Use function 'specify_hsmm' to specify your model.")
#'   #TODO: check model
#'   # seeds
#'   if(is.null(seq_sim_seed)) seq_sim_seed = sample(1:10000,1)
#'   if(verbose) cat("Checks done\n")
#'
#'   # compute $prob_values
#'   if(verbose) cat("Initialize joint emission probabilities\n")
#'   model$prob_values = .initialize_prob_values(model = model, seq_sim_seed = seq_sim_seed, nmin = nmin)
#'
#'   # compute $prob_missing
#'   if(verbose) cat("Initialize joint missing probabilities\n")
#'   model$prob_missing = .initialize_prob_missing(model = model)
#'
#'   # compute b table from $prob_values and $prob_missing
#'   if(verbose) cat("Initialize observation probabilities\n")
#'   model$b = .compute_b_from_prob_values_and_prob_missing(model = model)
#'
#'   # upgrade the class of the model and return it
#'   class(model) = "hsmm"
#'   model
#' }
#'
#'
#'
#' #' Initialize the list of matrices with the joint emission probabilities
#' #'
#' #' This function returns a list of matrix. Each matrix gives the joint probabilities for a specific set of observed variables. These matrices provide the joint probabilities for all value combinations of a specific variable set.
#' .initialize_prob_values = function(model = model, seq_sim_seed = seq_sim_seed, nmin = nmin){
#'
#'   prob_values = list() # initialization
#'
#'   # we retrieve all possible values from all variables
#'   all_levels = .get_all_possible_levels(model = model, with_missing = TRUE, continuous_var_binned = TRUE)
#'   all_levels_no_missing = .get_all_possible_levels(model = model, with_missing = FALSE, continuous_var_binned = TRUE)
#'
#'   # and specify which variable can be missing
#'   var_names = colnames(all_levels)
#'   vars_with_missing_values = var_names[apply(all_levels, 2, function(x) any(is.na(x)))]
#'   vars_never_missing = var_names %>% setdiff(vars_with_missing_values)
#'
#'   # if there is a data augmentation function, we need to get their marginal emission probabilities by simulating sequences
#'   if(length(model$augment_data_fun(get_var_names_only = TRUE)) > 0){
#'     Xsim = simulate_hsmm(model = model, n_state_transitions = 100, all_states = TRUE, min_tp_per_state = nmin, seed = seq_sim_seed)
#'   }else{
#'     Xsim = data.frame()
#'   }
#'
#'   # we format the marginal probabilities so that the matrices are easy to build.
#'   marg_probs = purrr::map(var_names, .get_marginal_prob, model = model, Xsim = Xsim) # marg_probs is a list with one element for each variable. Each element is a data.table with the marginal probabilities in long format for the specific variable.
#'   names(marg_probs) = var_names
#'
#'   # we create a matrix for each combination of observable variable set.
#'   # we loop over the size of the set (l loop), then loop again over all possible combination of that size (k loop).
#'   start_l = ifelse(length(vars_never_missing) == 0, 1, 0)
#'   for(l in start_l:length(vars_with_missing_values)){ # the size of the set of observable variables is length(var_never_missing) + l
#'     var_combinations = combn(x = vars_with_missing_values, m = l)
#'     var_combinations = rbind(var_combinations,
#'                              matrix(vars_never_missing, ncol = ncol(var_combinations), nrow = length(vars_never_missing), byrow = FALSE))
#'     for(k in 1:ncol(var_combinations)){
#'       # first we retrieve all the marginal probabilities
#'       P = data.frame(state = 1:model$J)
#'       for(var in var_combinations[,k]) P = P %>% dplyr::full_join(., marg_probs[[var]] %>% magrittr::set_colnames(c("state",var,paste0("prob_",var))), by = c("state"))
#'       # then we compute the joint probabilities
#'       P$p = apply(P %>% dplyr::select(dplyr::starts_with("prob_")), 1, prod)
#'       P$p0 = P$p
#'       # we attach p to the list
#'       this_var_combination = stringr::str_c(var_names[var_names %in% var_combinations[,k]], collapse = "_")
#'       prob_values[[this_var_combination]] = P %>% dplyr::select(state, tidyselect::all_of(var_combinations[,k]), p0, p)
#'     }
#'   }
#'   prob_values
#' }



#' #' Initialize the list of matrices with the joint missing probabilities
#' #'
#' #' This function returns a list of matrix. Each matrix gives the joint missing probabilities for a specific set of observed variables.
#' .initialize_prob_missing = function(model = model){
#'
#'   prob_missing = list() # initialization
#'
#'   # we retrieve all possible values from all variables
#'   all_levels = .get_all_possible_levels(model = model, with_missing = TRUE, continuous_var_binned = TRUE)
#'
#'   # and specify which variable can be missing
#'   var_names = colnames(all_levels)
#'   vars_with_missing_values = var_names[apply(all_levels, 2, function(x) any(is.na(x)))]
#'   if(length(vars_with_missing_values) == 0) return(prob_missing)
#'
#'   # we create a model$J-row matrix for each combination of observable variable set.
#'   # we loop over the size of the set (l loop), then loop again over all possible combination of that size (k loop).
#'   for(l in 1:length(vars_with_missing_values)){ # the size of the set of observable variables is l
#'     var_combinations = combn(x = vars_with_missing_values, m = l)
#'     for(k in 1:ncol(var_combinations)){
#'       # first we retrieve all the marginal probabilities
#'       P = data.frame(state = 1:model$J)
#'       for(var in var_combinations[,k]) P = dplyr::bind_cols(P, data.frame(p_missing = model$marg_em_probs[[var]]$missing_prob) %>% magrittr::set_colnames(paste0("p_missing_",var)))
#'       P[,var_combinations[,k]] = TRUE
#'       # then we compute the joint probabilities
#'       P$p = (1-model$global_censoring_prob) * apply(P %>% dplyr::select(dplyr::starts_with("p_missing_")), 1, prod)
#'       if(l == length(vars_with_missing_values)) P$p = P$p + model$global_censoring_prob
#'       P$p0 = P$p
#'       # we attach p to the list
#'       this_var_combination = stringr::str_c(var_names[var_names %in% var_combinations[,k]], collapse = "_")
#'       prob_missing[[this_var_combination]] = P %>% dplyr::select(state, all_of(var_combinations[,k]), p0, p)
#'     }
#'   }
#'   prob_missing
#' }
#'
#'
#' .compute_b_from_prob_values_and_prob_missing = function(model){
#'
#'   all_levels = .get_all_possible_levels(model = model, with_missing = TRUE, continuous_var_binned = TRUE)
#'   var_names = colnames(all_levels)
#'
#'   # first we build the structure of b as a long data.frame
#'   b = tidyr::expand_grid(state = 1:model$J)
#'   for(var in var_names) b = b %>% tidyr::expand_grid(., all_levels %>% dplyr::select(dplyr::all_of(var)) %>% dplyr::distinct())
#'   b = b %>% dplyr::mutate(dplyr::across(tidyselect::all_of(var_names), is.na, .names = "{col}_is_missing"))
#'
#'
#'   # then we build a P_missing data.frame with the missing probabilities
#'   P_missing = data.frame()
#'   for(missing_comb in names(model$prob_missing)){
#'     tmp = model$prob_missing[[missing_comb]]
#'     missing_vars = missing_comb %>% stringr::str_split(., pattern = "_") %>% unlist()
#'     tmp[,setdiff(var_names, missing_vars)] = FALSE
#'     tmp = tmp %>% dplyr::select(state, dplyr::all_of(var_names), p)
#'     P_missing = rbind(P_missing, tmp)
#'   }
#'
#'   # similarly, we build a P_values data.frame with the observation probabilities
#'   P_values = data.frame()
#'   for(obs_comb in names(model$prob_values)){
#'     tmp = model$prob_values[[obs_comb]]
#'     obs_vars = obs_comb %>% stringr::str_split(., pattern = "_") %>% unlist()
#'     tmp[,setdiff(var_names, obs_vars)] = NA
#'     tmp = tmp %>% dplyr::select(state, dplyr::all_of(var_names), p) %>%
#'       dplyr::rename(prob_value = p)
#'     P_values = rbind(P_values, tmp)
#'   }
#'
#'   # Now, we join b and P_missing on whether observations are missing or not:
#'   if(nrow(P_missing) > 0){
#'     b = b %>%
#'       dplyr::left_join(.,
#'                        P_missing %>% magrittr::set_colnames(c("state",paste0(var_names,"_is_missing"), "prob_missing")),
#'                        by = c("state", stringr::str_c(var_names,"_is_missing")))
#'
#'   }else{
#'     b = b %>% dplyr::mutate(prob_missing = 1)
#'   }
#'
#'   if(nrow(P_values) > 0){
#'     b = b %>%
#'       dplyr::left_join(.,
#'                        P_values,
#'                        by = c("state", var_names))
#'   }else{
#'     b = b %>% dplyr::mutate(prob_value = 1)
#'   }
#'
#'   # combining missing and value probabilities
#'   b = b %>%
#'     dplyr::mutate(
#'       prob_missing = prob_missing %>% tidyr::replace_na(1),
#'       prob_value = prob_value %>% tidyr::replace_na(1),
#'       p = prob_missing * prob_value
#'     )
#'
#'   # and finally, we convert to wide format for faster computation when decoding
#'   b = b %>%
#'     dplyr::select(state, tidyselect::all_of(var_names), p) %>%
#'     tidyr::pivot_wider(names_from = state, values_from = p, names_prefix = "p_")
#'   b
#' }







# cat(length(fwbw_res$F1),"\n")
# cat("range(fwbw_res$F1) : ",range(fwbw_res$F1),"\n")
# cat("any(is.nan(fwbw_res$F1)) : ",any(is.nan(fwbw_res$F1)),"\n")
# cat("any(is.nan(fwbw_res$N)) : ",any(is.nan(fwbw_res$N)),"\n")
# cat("any(is.nan(fwbw_res$si)) : ",any(is.nan(fwbw_res$si)),"\n")
# cat("any(is.nan(fwbw_res$G))  : ",any(is.nan(fwbw_res$G)),"\n")
# cat("any(is.nan(fwbw_res$L1)) : ",any(is.nan(fwbw_res$L1)),"\n")
#
# cat("all(is.nan(fwbw_res$F1)) : ",all(is.nan(fwbw_res$F1)),"\n")
# F1 = matrix(fwbw_res$F1, ncol = J)
# matplot(F1, type = "l", lty = 1, col = model$state_colors)
# cat("all(is.nan(fwbw_res$F1)) : ",all(is.nan(fwbw_res$F1)),"\n")
# cat("apply(F1, 2, function(x) any(is.nan(x))) : ", apply(F1, 2, function(x) any(is.nan(x))), "\n")
# cat("apply(F1, 2, function(x) sum(is.nan(x))) : ", apply(F1, 2, function(x) sum(is.nan(x))), "\n")
# matplot(is.nan(F1), type = "l", lty = 1, col = model$state_colors)
# plot(apply(F1,1, sum), type = "l")
#
# plot(fwbw_res$N, type = "l", ylim = c(0,0.001))
# cat("range(fwbw_res$N)",range(fwbw_res$N), "\n")
#
# G = matrix(fwbw_res$G, ncol = J)
# cat("all(is.nan(G)) : ",all(is.nan(G)),"\n")
# cat("apply(G, 2, function(x) any(is.nan(x))) : ", apply(G, 2, function(x) any(is.nan(x))), "\n")
# cat("apply(G, 2, function(x) sum(is.nan(x))) : ", apply(G, 2, function(x) sum(is.nan(x))), "\n")
# matplot(is.nan(G), type = "l", lty = 1, col = model$state_colors)




# .compute_observed_probabilities = function(model, X, verbose){
#
#   all_levels = .get_all_possible_levels(model = model)
#   var_names = colnames(all_levels)
#
#   Xp = .cut_continuous_var(model = model, X = X)
#
#   P = Xp %>%
#     dplyr::select(dplyr::all_of(var_names), state, w) %>%
#     dplyr::group_by(.dots = dplyr::all_of(c(var_names,"state"))) %>%
#     dplyr::summarize(n = sum(w), .groups = "drop") %>%
#     dplyr::group_by(state) %>% dplyr::mutate(Tot = sum(n)) %>% dplyr::ungroup() %>%
#     dplyr::mutate(p = n / Tot) %>%
#     dplyr::group_by(state) %>%  dplyr::mutate(p_max = max(p)) %>% dplyr::ungroup() %>%
#     dplyr::mutate(p_n = p/p_max)
#
#   P
# }






#
# all_levels = .get_all_possible_levels(model = model, with_missing = TRUE, continuous_var_binned = TRUE)
# n_levels = apply(all_levels, 2, function(x) length(unique(x)))
# obs_probs_var = data.frame(state = 1:model$J)
# for(var in colnames(all_levels)){
#   cn = colnames(obs_probs_var)
#   obs_probs_var = tidyr::expand_grid(obs_probs_var, unique(all_levels[,var])) %>% magrittr::set_colnames(c(cn, var))
# }
#
# if(verbose) cat("Computing marginal emission distribution and censoring probabilities of model variables.\n")
# for(var in names(model$parms.emission)){
#   if(verbose) cat("\t ",var,"\t")
#   probs_this_var = tidyr::expand_grid(state = 1:model$J, x = unique(all_levels[,var])) %>%
#     dplyr::mutate(prob_missing = ifelse(is.na(x), model$parms.emission[[var]]$missing_prob[state], 1))
#   if(model$parms.emission[[var]]$type == "norm"){
#     x_continuous = seq(min(model$parms.emission[[var]]$params$mean - 5 * model$parms.emission[[var]]$params$sd),
#                        max(model$parms.emission[[var]]$params$mean + 5 * model$parms.emission[[var]]$params$sd),
#                        len = 10000)
#     marg_prob = tidyr::expand_grid(x_continuous = x_continuous,
#                                    state = 1:model$J)
#     marg_prob = marg_prob %>% dplyr::mutate(x = cut(x_continuous, breaks = model$parms.emission[[var]]$breaks),
#                                             d = dnorm(x = x_continuous,
#                                                       mean = model$parms.emission[[var]]$params$mean[state],
#                                                       sd = model$parms.emission[[var]]$params$sd[state])) %>%
#       dplyr::group_by(state, x) %>% dplyr::summarize(prob = sum(d), .groups = "drop") %>%
#       dplyr::group_by(state) %>% dplyr::mutate(tot = sum(prob), prob_value = prob/tot) %>% dplyr::select(-tot,-prob)
#     probs_this_var = probs_this_var %>%
#       dplyr::left_join(., marg_prob, by = c("state","x")) %>%
#       dplyr::mutate(prob_value = prob_value %>% tidyr::replace_na(1))
#   }else if(model$parms.emission[[var]]$type == "binom"){
#     probs_this_var = probs_this_var %>%
#       dplyr::mutate(prob_value = dbinom(x = x, size = model$parms.emission[[var]]$params$size[state], prob = model$parms.emission[[var]]$params$prob[state]))
#   }else if(model$parms.emission[[var]]$type == "non-par"){
#     probs_this_var = probs_this_var %>%
#       dplyr::mutate(prob_value = model$parms.emission[[var]]$params$probs[cbind(match(x,model$parms.emission[[var]]$params$values),state)])
#   }
#   probs_this_var = probs_this_var %>%
#     dplyr::mutate(prob_value = prob_value %>% tidyr::replace_na(1),
#                   prob = prob_missing * prob_value) %>%
#     magrittr::set_colnames(c("state",var,stringr::str_c(c("prob_missing_","prob_value_","p_"),var)))
#
#   obs_probs_var = obs_probs_var %>% dplyr::left_join(., probs_this_var, by = c("state",var))
# }
# if(verbose) cat("\n")
#
#
# if(length(model$augment_data_fun(get_var_names_only = TRUE)) > 0){
#   if(verbose) cat("Computing marginal emission distribution of the augmented variables by simulating sequences.\n")
#   nmin = 200
#   Xsim = simulate_hsmm(model = model, n_state_transitions = 500, all_states = TRUE, min_tp_per_state = nmin)
#   Xsima = .augment_data(model = model, X = Xsim, verbose = FALSE)
#
#   Evar_names = model$augment_data_fun(get_var_names_only = TRUE)
#   Evar_types = model$augment_data_fun(get_var_types_only = TRUE)
#   for(var in Evar_names){
#     if(verbose) cat("\t ",var,"\t")
#     probs_this_var = tidyr::expand_grid(state = 1:model$J, x = unique(all_levels[,var]))
#
#     contingency_this_var = Xsima %>% dplyr::select(state, dplyr::all_of(var)) %>%
#       dplyr::rename(x = matches(var))
#     if(Evar_types[[var]]$type != "cat") contingency_this_var = contingency_this_var %>% dplyr::mutate(x = cut(x, breaks = c(-Inf, seq(0.2,0.8,by = 0.2),Inf)))
#     contingency_this_var = contingency_this_var %>%
#       dplyr::group_by(state, x) %>%
#       dplyr::summarize(n = dplyr::n(), .groups = "drop") %>%
#       dplyr::group_by(state) %>% dplyr::mutate(tot = sum(n)) %>% dplyr::ungroup() %>%
#       dplyr::mutate(prob = n/tot) %>%
#       dplyr::select(state, x, prob)
#
#     probs_this_var = probs_this_var %>% dplyr::left_join(., contingency_this_var, by = c("state","x")) %>%
#       dplyr::mutate(prob = prob %>% tidyr::replace_na(0),
#                     prob_value = prob) %>%
#       magrittr::set_colnames(c("state",var,stringr::str_c(c("prob_value_","p_"), var)))
#
#     obs_probs_var = obs_probs_var %>% dplyr::left_join(., probs_this_var, by = c("state",var))
#   }
# }
# if(verbose) cat("\n")
#
#
#
#
# if(verbose) cat("Computing joint probabilities (b).\n")
#
# obs_probs_var$p = apply(obs_probs_var %>% dplyr::select(starts_with("p_")), 1, prod)
#
# model$b = obs_probs_var %>%
#   dplyr::select(state, dplyr::all_of(colnames(all_levels)), p) %>%
#   tidyr::pivot_wider(id_cols = dplyr::all_of(colnames(all_levels)), names_from = "state", values_from = c("p"), names_prefix = "p_")
#
# if(verbose) cat("Function to compute 'b' from the observations.\n")
#
# model$compute_obs_probs_fun = function(model, X){
#
#   # we need all variable names (included augmenting variables)
#   all_levels = .get_all_possible_levels(model = model, continuous_var_binned = TRUE)
#   var_names = colnames(all_levels)
#
#   # cut continuous variables
#   Xp = .cut_continuous_var(model = model, X = X)
#
#   # join X with model$b
#   Xb = dplyr::left_join(Xp, model$b, by = intersect(colnames(model$b), colnames(Xp)))
#   # select the prob columns and turn them into a matrix
#   p = Xb %>% dplyr::select(dplyr::matches("p_")) %>% as.matrix()
#   p
# }
#
# if(verbose) cat("Keeping the prior distributions.\n")
#
# # keeping the prior prob_value and prob_missing for the fit
# var_names = colnames(all_levels)
# vars_with_missing_values = var_names[apply(all_levels, 2, function(x) any(is.na(x)))]
# vars_never_missing = var_names %>% setdiff(vars_with_missing_values)
#
# # prob_missing : probability of variables being missing in each state.
# model$prob_missing = list()
# if(length(vars_with_missing_values)>0){
#   for(l in 1:length(vars_with_missing_values)){
#     var_combinations = combn(x = vars_with_missing_values, m = l)
#     for(k in 1:ncol(var_combinations)){
#       this_var_combination = stringr::str_c(var_names[var_names %in% var_combinations[,k]], collapse = "_")
#       P = tidyr::expand_grid(state = 1:model$J)
#       for(var in var_combinations[,k]){
#         cn = colnames(P)
#         P = tidyr::expand_grid(P, TRUE) %>% magrittr::set_colnames(c(cn, var))
#       }
#       this_comb_obs_probs = obs_probs_var %>%
#         dplyr::select(state, tidyselect::all_of(var_combinations[,k]), tidyselect::all_of(paste0("prob_missing_",var_combinations[,k]))) %>%
#         dplyr::mutate(dplyr::across(tidyselect::all_of(var_combinations[,k]), is.na)) %>%
#         dplyr::distinct()
#       P = P %>% dplyr::left_join(., this_comb_obs_probs, by = c("state",var_combinations[,k]))
#       prob = apply(P %>% dplyr::select(tidyselect::starts_with("prob_missing_")), 1, prod)
#       P$p0 = prob
#       model$prob_missing[[this_var_combination]] = P %>% dplyr::select(state, tidyselect::all_of(var_combinations[,k]), p0)
#     }
#   }
# }
#
# # prob_values : joint probability of variables
# model$prob_values = list()
# all_levels_no_missing = .get_all_possible_levels(model = model, with_missing = FALSE, continuous_var_binned = TRUE)
# for(l in 0:length(vars_with_missing_values)){
#   var_combinations = combn(x = vars_with_missing_values, m = l)
#   var_combinations = rbind(var_combinations,
#                            matrix(vars_never_missing, ncol = ncol(var_combinations), nrow = length(vars_never_missing), byrow = FALSE))
#   for(k in 1:ncol(var_combinations)){
#     this_var_combination = stringr::str_c(var_names[var_names %in% var_combinations[,k]], collapse = "_")
#     P = tidyr::expand_grid(state = 1:model$J)
#     for(var in var_combinations[,k]){
#       cn = colnames(P)
#       P = tidyr::expand_grid(P, unique(all_levels_no_missing[,var])) %>% magrittr::set_colnames(c(cn, var))
#     }
#     this_comb_obs_probs = obs_probs_var %>%
#       dplyr::select(state, tidyselect::all_of(var_combinations[,k]), tidyselect::all_of(paste0("prob_value_",var_combinations[,k]))) %>%
#       dplyr::distinct()
#     P = P %>% dplyr::left_join(., this_comb_obs_probs, by = c("state",var_combinations[,k]))
#     prob = apply(P %>% dplyr::select(tidyselect::starts_with("prob_value_")), 1, prod)
#     P$p0 = prob
#     model$prob_values[[this_var_combination]] = P %>% dplyr::select(state, tidyselect::all_of(var_combinations[,k]), p0)
#   }
# }




