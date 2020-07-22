
####### INITIALIZATION AND CHECKING FUNCTIONS ---------------------



#' Specifies a hidden semi-Markov model
#'
#' This function returns an object of class \code{hsmm_spec} from the parameters passed as arguments.
#' @param J an integer. The number of hidden states.
#' @param init a vector of double. The initial probabilities associated with each states. Must be a vector of length \code{J}. If values do not sum to 1 a warning message is displayed and the vector is divided by its sum.
#' @param transition a matrix of double. The transition probabilities.
#'     Must be a matrix of size \code{J x J} in which the element \code{(r, c)} provides the probability of transitioning from state \code{r} to state \code{c}.
#'     The diagonal elements must be equal to zero except for absorbing states. For absorbing states, all the other elements of that row must be zero. For other states, if there are non-zero diagonal elements, a warning message is displayed and the elements are set to zero.
#'     If the sum over the rows is different that one, a warning message is displayed and the rows are divided by their sum.
#' @param sojourn a list. The sojourn distributions. The list must have at least two elements: a \code{type} element which specifies the type of sojourn distributions and the other elements are the distribution parameters.
#' Use the function \code{available_sojourn_dist()} to see the supported sojourn distributions.
#' @param parms.emission a list. The marginal emission distributions. TODO: explain how to format them.
#' @param state_names (optional) a vector of characters. Names associated to each state. Must be of length \code{J}. If unspecified, the states are numbered from 1 to \code{J}
#' @param state_colors (optional) a vector of color-specifying characters. Colors associated to each state. Must be of length \code{J}. If unspecified, the colors are picked from the \code{viridis} palette.
#' @param augment_data_fun (optional) a function that takes as argument a matrix of observation (see xxx TODO for how observations must be formatted) and returns a matrix of augmented variables. Additionally, this function must have xxx TODO xxx explain the requirements.
#'
#' @keywords HSMM
#' @return An object of class \code{hsmm_spec} which can be used to simulate data following the model specification and which can be initialized with the function \code{initialize_hsmm()}.
#'
#' @export
#' @importFrom magrittr %>%
#' @examples
#' my_model_spec = specify_hsmm(
#'    J = 2,
#'    init = c(1,0),
#'    transition = matrix(c(0,1,1,0),2,2),
#'    sojourn = list(type = "gamma", shape = c(2,10), scale = c(10,3)),
#'    parms.emission = list(
#'        var1 = list(
#'            type = "norm",
#'            params = list(
#'                mean = c(0,1),
#'                sd = c(0.3,0.2)
#'                ),
#'            missing_prob = c(0.6,0.1)
#'            ),
#'      var2 = list(
#'          type = "binom",
#'          params = list(
#'              size = rep(1,2),
#'              prob = c(0.2,0.8)
#'              ),
#'          missing_prob = c(0,0)
#'          ),
#'      var3 = list(
#'          type = "non-par",
#'          params = list(
#'              values = c("a","b","c","d"),
#'              probs = matrix(c(0.7,0.1,0.1,0.1,
#'                             1/4,1/4,1/4,1/4), 4,2)
#'            ),
#'          missing_prob = c(0.5,0.5)
#'          )
#'      ),
#'    state_names = c("A","B"),
#'    state_colors = c("seagreen1","lightcoral")
#' )
#' class(my_model_spec)
#' Xsim = simulate_hsmm(model = my_model_spec)
#' plot_hsmm_seq(model = my_model_spec, X = Xsim)
#'
#'
specify_hsmm = function(J, init, transition, sojourn, parms.emission, state_names = NULL, state_colors = NULL, augment_data_fun = NULL){

  if((round(J) != J) | (J <= 0)) stop("J must be a strictly positive integer.")

  # initial probabilities
  if(length(init) != J) stop("initial probability vector `init` must have J elements")

  # transition probabilities
  if(NROW(transition) != J)  stop("transition matrix `transition` must have J rows")
  if(NCOL(transition) != J)  stop("transition matrix `transition` must have J columns")
  if(any(rowSums(transition) != 1)) warning("Transition matrix rows do not sum to 1. Values will be normalized such that the transition probabilities from any state sum to 1.")
  if(any(diag(transition) != 0)){
    #any absorbing states?
    j = which(diag(transition) != 0) # j are the states that have self-transitions
    k = which(rowSums(transition - diag(transition)) == 0) # k are the absorbing states
    l = setdiff(j,k) # l are the non-absorbing states with self-transitions
    if(length(l) > 0){
      warning("Some non-absorbing states had non-zero self-transition probabilities (i.e. non-zero elements on the diagonal). These self-transitions will be set to 0.")
      transition[cbind(l,l)] = 0
    }
  }
  transition = transition/rowSums(transition)

  # sojourns
  if(is.null(sojourn$type)) stop("Sojourn distribution type not specified.")
  supported_sojourns = available_sojourn_dist()$distribution_type
  if(all(sojourn$type!= supported_sojourns)) stop(paste("Invalid sojourn type specified (",sojourn$type,"). Must be one of: ", paste(supported_sojourns, collapse = ", ")))
  # TODO: check that the params of the sojourns are specified as they should.

  # emission parameters
  if(is.null(parms.emission) | (length(parms.emission) == 0)) stop("Emission parameters (parms.emission) are not specified.
                                        They should be specified as a list with one element per observed variable.
                                        Each list element should specify the name of that variable, its distribution family (e.g. 'norm', 'binom', 'non-par')
                                        and its parameters/distribution for each state.")
  if(is.null(names(parms.emission))) names(parms.emission) = 1:length(parms.emission)

  flag = FALSE
  for(var in names(parms.emission)){
    var_parem = parms.emission[[var]]
    params_names = names(var_parem$params)
    # parameters
    if(var_parem$type == "norm"){
      if(!all(c( "mean", "sd" ) %in% params_names))
        stop(paste0("Variable '",var,"' is 'norm' and must have the following params: \n",
                    "parms.emission$",var,"$params = list(\nmean = ...vector of J means...,",
                    "\nsd = ...vector of J sd...,",
                    "\nnu (optional) = ...fitting hyper-parameter for conjugate prior...,",
                    "\nalpha (optional) = ...fitting hyper-parameter for conjugate prior...,",
                    "\nbeta (optional) = ...fitting hyper-parameter for conjugate prior...)"))
      if(is.null(var_parem$params$n0)) parms.emission[[var]]$params$n0 = 200
      if(is.null(var_parem$params$alpha)) parms.emission[[var]]$params$alpha = parms.emission[[var]]$params$n0/2
      if(is.null(var_parem$params$beta)) parms.emission[[var]]$params$beta = parms.emission[[var]]$params$alpha/2 * parms.emission[[var]]$param$sd^2
      parms.emission[[var]]$params$mu_0 = parms.emission[[var]]$params$mean
      parms.emission[[var]]$breaks =
        c(-Inf,
          seq(min(parms.emission[[var]]$params$mean - 3*parms.emission[[var]]$param$sd),
              max(parms.emission[[var]]$params$mean + 3*parms.emission[[var]]$param$sd),
              len = 8),
          Inf)

    }else if(var_parem$type == "binom"){
      if(!all(c( "size", "prob" ) %in% params_names))
        stop(paste0("Variable '",var,"' is 'binom' and must have the following params: \n",
                    "parms.emission$",var,"$params = list(",
                    "\nsize = ...vector of J size...,",
                    "\nprob = ...vector of J prob...,",
                    "\nalpha (optional) = ...fitting hyper-parameter for conjugate prior...,",
                    "\nbeta (optional) = ...fitting hyper-parameter for conjugate prior...)"))

      prob_0 = parms.emission[[var]]$params$prob
      if(is.null(var_parem$params$n0)) var_parem$params$n0 = 100
      parms.emission[[var]]$params$alpha = var_parem$params$n0 * prob_0
      parms.emission[[var]]$params$beta = parms.emission[[var]]$params$alpha * (1-prob_0)/prob_0

    }else if(var_parem$type == "non-par"){
      if(!all(c( "values", "probs" ) %in% params_names))
        stop(paste0("Variable '",var,"' is 'non-par' and must have the following params: \n",
                    "parms.emission$",var,"$params = list(\nvalues = ...vector of values..., \nprobs = ...matrix (nrow = # of values, ncol = # of states)..., \nalpha (optional) = ...fitting hyper-parameter for conjugate prior...)"))
      if(is.null(var_parem$params$n0)) parms.emission[[var]]$params$n0 = 100
      parms.emission[[var]]$params$probs_0 = parms.emission[[var]]$params$probs
    }

    # missing probs
    if(is.null(var_parem$missing_prob)){
      parms.emission[[var]]$missing_prob = rep(0,J); flag = TRUE
    }else{
      lmp = length(parms.emission[[var]]$missing_prob)
      if((lmp != 1) & (lmp<J)) stop(stringr::str_c("Length of 'missing_prob' of variable ",var," is neither 1 or J."))
      if(lmp == 1) parms.emission[[var]]$missing_prob = rep(parms.emission[[var]]$missing_prob,J)
    }
  }
  if(flag) warning("Some variables are assumed to never be missing. You can change this by specifying the missing probability for each variable in each state by adding an element 'missing_prob' to the list of the variables you think may be missing sometimes.")


  if(is.null(augment_data_fun)) augment_data_fun = function(X, get_var_names_only = FALSE) if(get_var_names_only) c() else data.frame()


  if(is.null(state_names)) state_names = 1:J
  if(length(state_names) != J) stop("state_names must be of length J")

  if(is.null(state_colors)) state_colors = viridis_pal(option = "C")(J)
  if(length(state_colors) != J) stop("state_colors must be of length J")

  ans = list(J = J,
             init = init,
             transition = transition,
             sojourn = sojourn,
             parms.emission = parms.emission,
             augment_data_fun = augment_data_fun,
             state_names = state_names,
             state_colors = state_colors)

  class(ans) <- 'hsmm_spec'
  ans
}



#' Initialize a hidden semi-Markov model.
#'
#' This function initializes the joint emission probabilities given the states,
#' i.e. determines the values of \eqn{Pr(X_i, E_i | S_j)}  where \eqn{X_i} are the observations at a given time-point \eqn{i}, \eqn{E_i} are the values of the augmented variables at time-point \eqn{i} and \eqn{S_j} is the state \eqn{j}.
#'
#' @param model a \code{hsmm_spec} object. A specified hidden semi-Markov model. Use the function \code{specify_hsmm()} to specify a hidden semi-Markov model.
#' @param nmin (optional) a strictly positive integer. XXXX
#' @param seq_sim_seed (optional) an integer. XXX
#' @param verbose (optional) a logical. If TRUE, the function prints the internal steps of the function.
#'
#' @return an object of class \code{hsmm}, which can be used to decode observation sequences with the function \code{predict_states_hsmm()} or which can be fitted to observations with the function \code{fit_hsmm()}.
#'
#' @importFrom magrittr %>%
#' @export
#' @examples
#' my_model_spec = simple_model_spec
#' class(my_model_spec)
#' my_model_init = initialize_hsmm(my_model_spec)
#' class(my_model_init)
#'
initialize_hsmm = function(model,
                           nmin = NULL,
                           verbose = FALSE,
                           seq_sim_seed = NULL){

  # CHECKS
  # check if model is of class 'hsmm_spec' or 'hsmm'
  if(!(class(model) %in% c("hsmm_spec","hsmm"))) stop("model must be of class 'hsmm_spec'. Use function 'specify_hsmm' to specify your model.")
  # seeds
  if(is.null(seq_sim_seed)) seq_sim_seed = sample(1:10000,1)

  if(verbose) cat("Checks done\n")


  all_levels = .get_all_possible_levels(model = model, continuous_var_binned = TRUE)
  n_levels = apply(all_levels, 2, function(x) length(unique(x)))
  obs_probs_var = data.frame(state = 1:model$J)
  for(var in colnames(all_levels)){
    cn = colnames(obs_probs_var)
    obs_probs_var = tidyr::expand_grid(obs_probs_var, unique(all_levels[,var])) %>% magrittr::set_colnames(c(cn, var))
  }

  if(verbose) cat("Computing marginal emission distribution of model variables.\n")
  for(var in names(model$parms.emission)){
    probs_this_var = obs_probs_var %>% dplyr::select(state, dplyr::all_of(var)) %>% unique() %>%
      dplyr::rename(x = matches(var)) %>%
      dplyr::mutate(prob_missing = (!is.na(x)) + ifelse(is.na(x),1,-1)* model$parms.emission[[var]]$missing_prob[state])

    if(model$parms.emission[[var]]$type == "norm"){
      x_continuous = seq(min(model$parms.emission[[var]]$params$mean - 5 * model$parms.emission[[var]]$params$sd),
               max(model$parms.emission[[var]]$params$mean + 5 * model$parms.emission[[var]]$params$sd),
               len = 10000)
      marg_prob = tidyr::expand_grid(x_continuous = x_continuous,
                      state = 1:model$J)
      marg_prob = marg_prob %>% dplyr::mutate(x = cut(x_continuous, breaks = model$parms.emission[[var]]$breaks),
                                       d = dnorm(x = x_continuous,
                                                 mean = model$parms.emission[[var]]$params$mean[state],
                                                 sd = model$parms.emission[[var]]$params$sd[state])) %>%
        dplyr::group_by(state, x) %>% dplyr::summarize(prob = sum(d), .groups = "drop") %>%
        dplyr::group_by(state) %>% dplyr::mutate(tot = sum(prob), prob_value = prob/tot) %>% dplyr::select(-tot,-prob)
      probs_this_var = probs_this_var %>%
        dplyr::left_join(., marg_prob, by = c("state","x")) %>%
        dplyr::mutate(prob_value = prob_value %>% tidyr::replace_na(1))
    }else if(model$parms.emission[[var]]$type == "binom"){
      probs_this_var = probs_this_var %>%
        dplyr::mutate(prob_value = dbinom(x = x, size = model$parms.emission[[var]]$params$size[state], prob = model$parms.emission[[var]]$params$prob[state]))
    }else if(model$parms.emission[[var]]$type == "non-par"){
      probs_this_var = probs_this_var %>%
        dplyr::mutate(prob_value = model$parms.emission[[var]]$params$probs[cbind(match(x,model$parms.emission[[var]]$params$values),state)])

    }
    probs_this_var = probs_this_var %>%
      dplyr::mutate(prob_value = prob_value %>% tidyr::replace_na(1),
             prob = prob_missing * prob_value) %>%
      dplyr::select(state, x, prob) %>%
      magrittr::set_colnames(c("state",var,stringr::str_c("prob_",var)))

    obs_probs_var = obs_probs_var %>% dplyr::left_join(., probs_this_var, by = c("state",var))
  }


  if(length(model$augment_data_fun(get_var_names_only = TRUE)) > 0){
    if(verbose) cat("Computing marginal emission distribution of the augmented variables by simulating sequences.\n")
    nmin = 200
    Xsim = simulate_hsmm(model = model, n_state_transitions = 500, all_states = TRUE, min_tp_per_state = nmin)
    Xsima = .augment_data(model = model, X = Xsim, verbose = FALSE)

    Evar_names = model$augment_data_fun(get_var_names_only = TRUE)
    Evar_types = model$augment_data_fun(get_var_types_only = TRUE)
    for(var in Evar_names){
      probs_this_var = obs_probs_var %>% dplyr::select(state, dplyr::all_of(var)) %>% unique()
      obs_prob_this_var = Xsima %>% dplyr::select(state, dplyr::all_of(var)) %>%
        dplyr::rename(x = matches(var))
      if(Evar_types[[var]]$type != "cat") obs_prob_this_var = obs_prob_this_var %>% dplyr::mutate(x = cut(x, breaks = c(-Inf, seq(0.2,0.8,by = 0.2),Inf)))
      obs_prob_this_var = obs_prob_this_var %>%
        dplyr::group_by(state, x) %>%
        dplyr::summarize(n = n(), .groups = "drop") %>%
        dplyr::group_by(state) %>% dplyr::mutate(tot = sum(n)) %>% dplyr::ungroup() %>%
        dplyr::mutate(prob = n/tot) %>%
        dplyr::select(state, x, prob) %>%
        magrittr::set_colnames(c("state",var,"prob"))
      probs_this_var = probs_this_var %>% dplyr::left_join(., obs_prob_this_var, by = c("state",var)) %>%
        dplyr::mutate(prob = prob %>% tidyr::replace_na(0)) %>%
        magrittr::set_colnames(c("state",var,stringr::str_c("prob_",var)))
      obs_probs_var = obs_probs_var %>% dplyr::left_join(., probs_this_var, by = c("state",var))
    }
  }

  if(verbose) cat("Computing joint probabilities.\n")

  prob = apply(
    obs_probs_var %>% dplyr::select(starts_with("prob_")) %>% as.matrix(),
    1, prod)

  obs_probs_var = obs_probs_var %>%  dplyr::mutate(p = prob) %>%
    dplyr::group_by(state) %>%
    dplyr::mutate(p_max = max(p),
           p_n = p/p_max) %>% dplyr::ungroup()


  model$obs_probs = obs_probs_var %>% dplyr::select(state, dplyr::all_of(colnames(all_levels)), p, p_n)
  model$obs_probs_0 = model$obs_probs # we keep the initial one for the fitting process

  model$compute_obs_probs_fun = function(model, X, state){
    s = state
    all_levels = .get_all_possible_levels(model = model, continuous_var_binned = TRUE)
    var_names = colnames(all_levels)
    # cut continuous variables
    Xp = .cut_continuous_var(model = model, X = X)

    # join X with model$obs_probs filtered for the state
    P = model$obs_probs %>% dplyr::filter(state == s) %>%  dplyr::select(dplyr::all_of(var_names), p_n) # p_n
    Xb = dplyr::left_join(Xp, P, by = intersect(colnames(P), colnames(Xp)))
    # replace NAs by 0
    Xb = Xb %>% dplyr::mutate(p_n = p_n %>% tidyr::replace_na(0)) # p_n
    Xb$p_n # p_b
  }

  # upgrade the class of the model
  class(model) = "hsmm"
  # and return it
  model
}


.prepare_data_for_init = function(
  model,
  nmin,
  seq_sim_seed,
  verbose = FALSE){

  states = 1:model$J

  all_levels = .get_all_possible_levels(model = model, with_missing = TRUE)
  n_levels = apply(all_levels, 2, function(x) length(unique(x)))
  if(is.null(nmin)) nmin = pmin(200 * ceiling(prod(n_levels)^0.5),2000)

  if(verbose) cat("nmin = ",nmin,"\n")

  # we simulate data
  X = simulate_hsmm(model = model,
                    n_state_transitions = 500,
                    all_states = TRUE, mult_seq_allowed = TRUE, min_tp_per_state = nmin,
                    seed = seq_sim_seed)

  if(verbose) cat("data simulated \n")

  # we augment the data
  X = X %>% dplyr::arrange(seq_id, t)
  X = .augment_data(X = X, model = model, verbose = verbose)

  # weights
  X$w = 1
  X
}


.compute_observed_probabilities = function(model, X, verbose){

  all_levels = .get_all_possible_levels(model = model, with_missing = TRUE)
  var_names = colnames(all_levels)

  Xp = .cut_continuous_var(model = model, X = X)

  P = Xp %>%
    dplyr::select(dplyr::all_of(var_names), state, w) %>%
    dplyr::group_by(.dots = dplyr::all_of(c(var_names,"state"))) %>%
    dplyr::summarize(n = sum(w), .groups = "drop") %>%
    dplyr::group_by(state) %>% dplyr::mutate(Tot = sum(n)) %>% dplyr::ungroup() %>%
    dplyr::mutate(p = n / Tot) %>%
    dplyr::group_by(state) %>%  dplyr::mutate(p_max = max(p)) %>% dplyr::ungroup() %>%
    dplyr::mutate(p_n = p/p_max)

  P
}


.cut_continuous_var = function(model, X){
  Xp = X
  for(var in names(model$parms.emission)) if(model$parms.emission[[var]]$type == "norm") Xp[,var] = X[,var] %>% unlist() %>% cut(., breaks = model$parms.emission[[var]]$breaks)
  for(var in model$augment_data_fun(get_var_names_only = TRUE)) if(typeof(X[,var]) != "integer") Xp[,var] = X[,var] %>% unlist() %>% cut(., breaks = c(-Inf, 0.2,0.4,0.6,0.8,Inf))
  Xp
}


#' @export
available_sojourn_dist = function(){
  available_sojourn = rbind(
    data.frame(distribution_type = "nonparametric", parameters = "d", stringsAsFactors = FALSE),
    data.frame(distribution_type = "ksmoothed-nonparametric", parameters = "d", stringsAsFactors = FALSE),
    data.frame(distribution_type = "gamma", parameters = "shape, scale", stringsAsFactors = FALSE),
    data.frame(distribution_type = "poisson", parameters = "shift, lambda", stringsAsFactors = FALSE),
    data.frame(distribution_type = "lnorm", parameters = "meanlog, s.dlog", stringsAsFactors = FALSE),
    data.frame(distribution_type = "logarithmic", parameters = "shape", stringsAsFactors = FALSE),
    data.frame(distribution_type = "nbinom", parameters = "size, mu (or prob), shift", stringsAsFactors = FALSE)
  )
  available_sojourn
}


#' @export
.check_data = function(data, model){
  if(!("data.frame" %in% class(data))) stop("data must be a data.frame\n")
  if(!("seq_id" %in% colnames(data))){
    warning("column 'seq_id' is missing. Assuming single sequence.")
    data$seq_id = 1
  }
  if(!("t" %in% colnames(data))){
    warning("column 't' is missing. Assuming no missing time-points and ordered data.frame.")
    data = data %>% dplyr::group_by(seq_id) %>% dplyr::mutate(t = row_number())
  }
  # checking that the data has a column for each model variable
  j = which(!(names(model$parms.emission) %in% colnames(data)))
  if(length(j)>0) stop(stringr::str_c("The data must contain a column for each of the model variable. Missing variables: ",stringr::str_c(names(model$parms.emission)[j],collapse = ", ")))
  data = data %>% dplyr::arrange(seq_id, t) %>% dplyr::ungroup()

  # select only the columns we need
  var_names = names(model$parms.emission)
  if(.is_data_augmented(data = data, model = model)) var_names = c(var_names,model$augment_data_fun(X = NULL, get_var_names_only = TRUE)) # stringr::str_c(var_names,"_M")
  data = data %>% dplyr::select(seq_id, t, dplyr::all_of(var_names))

  # checking the type and values of each variable
  data = .check_types_and_values(data = data, parem = model$parms.emission)
  data
}


.check_types_and_values = function(data, parem, continuous_var_binned = FALSE){

  for(var in names(parem)){
    if(parem[[var]]$type == "non-par") data[,var] = data[,var] %>% unlist() %>% factor(., levels = parem[[var]]$params$values)

    if(parem[[var]]$type == "norm"){
      if(!continuous_var_binned){
        data[,var] = data[,var]  %>%  unlist() %>% as.character() %>% as.double()
      }else{
        data[,var] = data[,var] %>% unlist() %>% factor(., levels = levels(cut(1:10, breaks = parem[[var]]$breaks)))
      }
    }

    if(parem[[var]]$type == "binom"){
      data[,var] =  data[,var] %>%  unlist() %>% as.character() %>% as.integer()
      if(!all(is.na(data[,var]))){
        if(min(data[,var], na.rm = TRUE) < 0) stop(stringr::str_c("variable ",var," is described by a binomial distribution but the provided data has negative values for this variable."))
        if(max(data[,var], na.rm = TRUE) > max(parem[[var]]$params$size)) stop(stringr::str_c("variable ",var," is described by a binomial distribution with max size: ",max(parem[[var]]$params$size)," but the provided data has values as high as ",max(data[,var], na.rm = TRUE)))
      }
    }
  }
  data
}

#' @export
.augment_data = function(model, X, verbose = FALSE){

  if(.is_data_augmented(data = X, model = model)) return(X)

  var_names = names(model$parms.emission)
  #M = X %>% select(all_of(var_names)) %>% is.na() %>%  set_colnames(str_c(var_names, "_M")) %>%  set_rownames(rownames(X))
  #M = M*1

  if(verbose) cat("augmenting the data ...")
  E_fun = model$augment_data_fun
  E = E_fun(X)
  if(verbose) cat("... done\n")

  if(nrow(E) == nrow(X)) augmented_X = cbind(X, E) else augmented_X = X #cbind(X, M)
  augmented_X
}

.is_data_augmented = function(data, model = model){
  required_cols = c(names(model$parms.emission),  model$augment_data_fun(X = NULL, get_var_names_only = TRUE)) # str_c(names(model$parms.emission),"_M"),
  ifelse(all(required_cols %in% colnames(data)), TRUE, FALSE)
}

.check_ground_truth = function(ground_truth, model, X){

  # check the columns
  cond_col = all(c("seq_id","t","state") %in% colnames(ground_truth))
  if(!cond_col){stop("Ground truth must be a data.frame with the following columns: seq_id, t, state, weight (optional).")}
  # ground_truth$state must be numbers from 1:J or NAs
  cond_state = all(unique(ground_truth$state) %in% c(NA, 1:model$J))
  if(!cond_state){stop("The ground truth state column must be a number in 1:J where J is the number of state in the model. NAs are accepted.")}
  # there should not be any duplicated in the ground_truth
  if(any(duplicated(ground_truth %>% dplyr::select(seq_id, t)))) stop("There are duplicated time-points in the provided ground truth.")

  # we match ground truth to X
  tmp = dplyr::left_join(X %>% dplyr::select(seq_id, t) %>%  dplyr::mutate(seq_id = as.character(seq_id), t = as.numeric(t), in_X = TRUE),
                  ground_truth %>%  dplyr::mutate(seq_id = as.character(seq_id), t = as.numeric(t), in_GT = TRUE),
                  by = c("seq_id","t"))
  ground_truth = tmp; rm(tmp)
  # give a warning if none of the labels matched the X sequences
  if(all(is.na(ground_truth$state))){warning("The provided ground truth does not match any of the sequences of the observed data.")}

  # we define states as factors
  ground_truth = ground_truth %>% dplyr::mutate(state = factor(state, levels = 1:model$J))

  # we return it
  ground_truth
}



####### DECODING (predicting state sequence) ---------------------


# wrapper for "predict" generic function
predict.hsmm <- function(object, newdata, method = "viterbi", verbose = FALSE, ...) {

  if(class(object) == 'hsmm_spec')stop(stringr::str_c("\n object's class is 'hsmm_spec'.\n",
                                             "Use the function 'initialize_hsmm' to transform its class to 'hsmm'.\n",
                                             "This function will provide your model with models and functions that ",
                                             "will be used to compute the local probability of each state."))

  if(class(object) != "hsmm") stop("object must be of class 'hsmm'.")

  model = object
  X = newdata
  ans = predict_states_hsmm(model = model, X = X, method = method, verbose = verbose, ...)
  ans
}



#' Predicts the hidden state sequence from observations
#'
#' This function predict the most likely hidden states from observations.
#' Two methods are implemented:
#' the first one applies the Viterbi algorithm and predicts the most likely sequence of hidden states,
#' the second one applies the Forward-Backward algorithm and returns the probabilitiy of each state at each time-point.
#'
#' @param model a \code{hsmm} object. The model used to predict the hidden sequence of states.
#' @param X a \code{data.frame}. The observations.
#' @param method a \code{character} specifying the decoding algorithm.
#'    \code{method = "viterbi"} returns the most likely sequence of hidden states.
#'    \code{method = "smoothed"} applies the "Forward-Backward" algorithm as described by Guédon, 2003, and returns the "smoothed" (posterior) probability of each state at each time-point.
#' @param ground_truth (optional) a \code{data.frame}.
#'
#' @export
#' @examples
#' my_model_spec = simple_model_spec
#' Xsim = simulate_hsmm(my_model_spec, n_state_transitions = 20)
#' my_model_init = initialize_hsmm(my_model_spec)
#' viterbi = predict_states_hsmm(model = my_model_init, X = Xsim, method = "viterbi")
#' Xsim$state_viterbi = viterbi$state_seq$state
#' plot_hsmm_seq(X = Xsim, model = my_model_init)
#'
predict_states_hsmm <- function(model, X, method = "viterbi", ground_truth = data.frame(), trust_in_ground_truth = 0.75, verbose = FALSE, graphical = FALSE, c_code = "original", ...) {

  # CHECKS
  # model
  if(class(model) == 'hsmm_spec')stop(stringr::str_c("\n model's class is 'hsmm_spec'.\n",
                                            "Use the function 'initialize_hsmm' to transform its class to 'hsmm'.\n",
                                            "This function will provide your model with models and functions that ",
                                            "will be used to compute the local probability of each state."))

  if(class(model) != "hsmm") stop("model must be of class 'hsmm'.")
  # check model
  # TODO

  # check data
  X = .check_data(data = X, model = model)
  X = .augment_data(model = model, X = X, verbose = verbose)


  # decoding
  if(method=="viterbi") {
    ans = .predict_states_hsmm_viterbi(model = model, X = X, ground_truth = ground_truth, trust_in_ground_truth = trust_in_ground_truth, verbose = verbose, graphical = graphical)
  }else if(method == "smoothed"){
    ans = .predict_states_hsmm_forward_backward_seq_by_seq(model = model, X = X, ground_truth = ground_truth, trust_in_ground_truth = trust_in_ground_truth, verbose = verbose, graphical = graphical, c_code = c_code)
  }else{stop(paste("Unavailable prediction method",method))}

  ans
}


#' @useDynLib HiddenSemiMarkov viterbi
.predict_states_hsmm_viterbi = function(model, X = X, ground_truth = data.frame(), trust_in_ground_truth = 0.75, verbose = FALSE, graphical = FALSE){

  # check model
  # TODO
  # check data
  X = .check_data(data = X, model = model)
  if(verbose){cat("Data checked \n")}
  X = .augment_data(model = model, X = X, verbose = verbose)
  if(nrow(ground_truth) > 0) ground_truth = .check_ground_truth(ground_truth, model, X)

  # Number of states
  J = model$J
  # Sequences lengths
  N = rle(X$seq_id %>% as.character())$lengths

  # add variables N, NN, M, d, D, b + log transform
  augmented_model = .augment_model(model = model, X = X, ground_truth = ground_truth, log = TRUE, m = -1e300, trust_in_ground_truth = trust_in_ground_truth, verbose = verbose, graphical = graphical)
  NN = augmented_model$NN
  if(verbose){cat("Model augmented \n")}
  SeqIDs = unique(X$seq_id)

  # Initialization of the output variables
  state_seq = data.frame()
  loglik = data.frame()

  # decode sequence by sequence
  if(verbose){cat(stringr::str_c("Decoding sequences (x/",length(N),"):\t"))}
  for(i in 1:length(N)) {
    if(verbose){cat(i,"\t")}
    ix = (NN[i]+1):NN[i+1]
    tmp = .C("viterbi",
             a=augmented_model$trans_log %>% as.double(), # transition probabilities
             pi=augmented_model$init_log %>% as.double(), # initial probabilities
             p=augmented_model$b_log[ix,] %>% as.double(), # local state probabilities
             d=augmented_model$d_log %>% as.double(), # sojourn density
             D=augmented_model$D_log %>% as.double(), # sojourn cumulative density
             timelength=as.integer(N[i]), # sequence length
             J=as.integer(J), # number of states
             M=as.integer(rep(augmented_model$M,J)), # max sojourn
             alpha = double(N[i]*J), # initialization of the state probabilities
             statehat=integer(N[i]), # initialization of the estimated state sequence
             psi_state0=integer(N[i]*J), # initialization of the "previous-state" matrix
             psi_time0=integer(N[i]*J) # initialization of the sojourn matrix
             ,PACKAGE='HiddenSemiMarkov')

    # needed to fix the backtracking
    alpha = matrix(tmp$alpha, nrow = N[i], ncol = J) # alpha is the likelihood matrix
    psi_state = matrix(tmp$psi_state0,  nrow = N[i], ncol = J)+1 # psi_state is the "previous-state" matrix
    psi_time = matrix(tmp$psi_time0,  nrow = N[i], ncol = J) # psi_time is the sojourn matrix

    # adding results to the output variables
    s = .hsmm_viterbi_backtracking(alpha = alpha, sojourns = psi_time, previous_states = psi_state) # estimated sequence of states
    statehat = tmp$statehat+1
    state_seq = dplyr::bind_rows(state_seq,
                          data.frame(seq_id = SeqIDs[i] %>% factor(.,levels = SeqIDs),
                                     t = X$t[ix],
                                     state = s,
                                     state_deprecated = statehat,
                                     likelihood = alpha[cbind(1:N[i],s)]))
    loglik = dplyr::bind_rows(loglik,
                       data.frame(seq_id = SeqIDs[i] %>% factor(.,levels = SeqIDs),
                                  loglik = max(alpha[N[i],(1:J)])))
  }
  if(verbose) cat("\n")

  ans <- list(state_seq = state_seq,
              loglik = loglik)
  return(ans)
}


.hsmm_viterbi_backtracking = function(alpha, sojourns, previous_states){
  i = nrow(alpha)
  sequence = which.max(alpha[i,])
  while(length(sequence) < nrow(alpha)){
    this_state_sojourn = sojourns[i, dplyr::first(sequence)]
    sequence = c(rep(dplyr::first(sequence), this_state_sojourn-1), sequence)
    i = i-this_state_sojourn
    previous_state = previous_states[i,dplyr::first(sequence)]
    sequence = c(previous_state, sequence)
  }
  sequence
}



#' @useDynLib HiddenSemiMarkov backward
#' @useDynLib HiddenSemiMarkov backward_original
.predict_states_hsmm_forward_backward = function(model, X = X, ground_truth = ground_truth, trust_in_ground_truth = 0.75,  verbose = FALSE, graphical = FALSE, c_code = c("log-prob","original")){

  # check model
  # TODO

  # check data
  X = .check_data(data = X, model)
  if(verbose){cat("Data checked \n")}
  X = .augment_data(model = model, X = X, verbose = verbose)
  if((nrow(ground_truth)>0) && (!all(is.na(ground_truth$state)))) ground_truth = .check_ground_truth(ground_truth, model, X)

  # add variables N, NN, M, d, D, b
  augmented_model = .augment_model(model = model, X = X, ground_truth = ground_truth, log = FALSE, m = -1e300, trust_in_ground_truth = trust_in_ground_truth, graphical = graphical)
  if(verbose) cat("Model augmented \n")

  # shortcuts
  J = augmented_model$J
  N = augmented_model$N
  NN = augmented_model$NN
  M = augmented_model$M
  Nseq = length(augmented_model$N)
  SeqIDs = unique(X$seq_id)
  nX = sum(N)

  # local distribution on states
  local = data.frame(
    seq_id = rep(X$seq_id, J),
    t = rep(X$t, J),
    state = rep(1:J, each = nrow(X)),
    local = augmented_model$b %>% as.vector(),
    stringsAsFactors = FALSE
  )

  # run the smoothed algorithm on all sequences at once
  if(verbose) cat(stringr::str_c("Decoding sequences (",length(N),")... \n"))

  #algorithm = ifelse(c_code == "original", "backward_original", "backward")
  algorithm = "backward_original"

  fwbw_res  = .C(algorithm,
                 transition=as.double(augmented_model$transition), # transition probabilities
                 init=as.double(augmented_model$init), # initial probabilities
                 p=as.double(augmented_model$b), # joint emission probabilities
                 d=as.double(augmented_model$d), # sojourn (d)
                 D=as.double(augmented_model$D), # sojourn (D)
                 timelength=as.integer(N), # length of each sequence
                 J=as.integer(J), # number of states
                 M=as.integer(rep(M,J)), # length of longest sequence
                 L1 = double(nX*J), # ??
                 N = double(nX), # will return the likelihood
                 eta = double(M*J), # ?? M
                 F1=double(J*nX), # ??
                 si=double(J*nX), # ??
                 gamma=double(J*nX), # will return the posterior probabilities
                 nsequences=as.integer(Nseq), # number of sequences
                 totallength=as.integer(nX), # sum of lenghts of all sequences
                 G=double(J*nX), # ??
                 PACKAGE='HiddenSemiMarkov')

  # Check for errors
  if(any(is.nan(fwbw_res$gamma))) { # gamma is the smoothed probability of each state at each time-point
    stop("NaNs detected in posterior probabilities.")
  }
  if(any(fwbw_res$gamma<0)) fwbw_res$gamma = zapsmall(fwbw_res$gamma)
  if(any(fwbw_res$eta<0)) fwbw_res$eta = zapsmall(fwbw_res$eta)
  if(any(fwbw_res$N<0))  fwbw_res$N = zapsmall(fwbw_res$N)


  if(verbose) cat("Formating results \n")
  ##### Format results
  # State probabilities
  state_probs = local %>%
    dplyr::mutate(posterior = fwbw_res$gamma)
  # State sequence
  p = matrix(fwbw_res$gamma, ncol=J)
  s = apply(p, 1, which.max);
  state_seq = X %>% dplyr::select(seq_id, t) %>%
    dplyr::mutate(state = s, likelihood = fwbw_res$N)
  # Log Likelihood of each sequence
  loglik = state_seq %>%
    dplyr::mutate(ll = log(likelihood)) %>%
    dplyr::group_by(seq_id) %>%
    dplyr::summarize(loglik = sum(ll), .groups = "drop") %>%
    dplyr::ungroup()

  ##### return results
  ans = list(state_seq = state_seq,
             loglik = loglik,
             state_probs = state_probs,
             fwbw_res = fwbw_res)
  ans
}



.predict_states_hsmm_forward_backward_seq_by_seq = function(model, X = X, ground_truth = data.frame(), trust_in_ground_truth = 0.75,  verbose = FALSE, graphical = FALSE, c_code = "original"){

  # check model
  # TODO
  # check data
  X = .check_data(data = X, model = model)
  if(verbose){cat("Data checked \n")}
  X = .augment_data(model = model, X = X, verbose = verbose)
  if(verbose){cat("Data augmented \n")}
  if(nrow(ground_truth)>0) ground_truth = .check_ground_truth(ground_truth = ground_truth, model = model, X = X)

  # shortcuts
  SeqIDs = unique(X$seq_id)
  Nseq = length(SeqIDs)

  # initialize output variables
  state_seq = data.frame()
  state_probs = data.frame()
  loglik = data.frame()
  sequences_with_error = c()

  # run the smoothed algorithm sequence by sequence
  if(verbose) cat(stringr::str_c("Decoding sequences (x/",length(Nseq),"):\t"))
  for(i in 1:Nseq){
    if(verbose) cat(i,"\t")
    SeqID = SeqIDs[i]
    Xi = X %>% dplyr::filter(seq_id == SeqID) %>% dplyr::arrange(seq_id, t)
    if(nrow(ground_truth)>0) ground_truth_i = ground_truth %>% dplyr::filter(seq_id == SeqID) %>% dplyr::arrange(seq_id, t) else ground_truth_i = data.frame()

    Di = try(.predict_states_hsmm_forward_backward(model = model, X = Xi, ground_truth = ground_truth_i, verbose = FALSE, graphical = graphical, c_code = c_code))
    if(class(Di) == "try-error"){
      sequences_with_error = c(sequences_with_error, SeqID)
    }else{
      state_seq = dplyr::bind_rows(state_seq, Di$state_seq)
      state_probs = dplyr::bind_rows(state_probs, Di$state_probs)
      loglik = dplyr::bind_rows(loglik, Di$loglik)
    }
  }
  if(verbose) cat("\n done \n")

  # return results
  ans = list(state_seq = state_seq,
             loglik = loglik,
             state_probs = state_probs,
             sequences_with_error = sequences_with_error)
  ans
}



.augment_model = function(model, X, ground_truth, log = FALSE, m = -1e300, trust_in_ground_truth = 0.75, verbose = FALSE, graphical = FALSE){
  augmented_model = model

  # SEQUENCES start index
  # N
  augmented_model$N = rle(X$seq_id %>% as.character())$lengths
  # NN
  augmented_model$NN = cumsum(c(0,augmented_model$N))
  # M
  augmented_model$M = max(augmented_model$N) #
  #augmented_model$M = max(c(augmented_model$M, 3700))

  # SOJOURN
  # d
  augmented_model$d = .build_d_from_sojourn_dist(model = model, M = augmented_model$M)
  # D
  augmented_model$D = apply(augmented_model$d,2,function(x) rev(cumsum(rev(x))))

  if(verbose) cat("Sojourn matrix built\n")


  # OBSERVATION PROBABILITIES
  augmented_model$b = sapply(1:model$J, function(state) model$compute_obs_probs_fun(model = model, X = X, state = state))
  if(graphical) matplot(augmented_model$b, type = "l", col = model$state_colors, lty = 1)

  # LOCAL state probabilities
  # compute_local_state_prob = model$compute_local_state_prob_fun
  # augmented_model$b = sapply(1:model$J, function(state) compute_local_state_prob(model = model, X = X, state = state))
  # augmented_model$b = .compute_local_prob_with_KNN_imputation(
  #   X = X,
  #   ground_truth = ground_truth,
  #   model =  model,
  #   K = 20)

  if(verbose) cat("joint emission probabilities computed\n")


  # Helping the decoding with provided ground truth
  if((nrow(ground_truth)>0) && !(all(is.na(ground_truth$state))) ){
    # building the matrix GT that is a one-hot encoding of the ground-truth states
    GT = matrix(0, nrow = nrow(X), ncol = model$J)
    j = which(!is.na(ground_truth$state))
    GT[cbind(j,ground_truth$state[j])] = 1
    # modifying b (the local state probabilities) where we have a ground truth
    j = which(rowSums(GT)==1)
    augmented_model$b[j,] = (1 - trust_in_ground_truth)* augmented_model$b[j,] + trust_in_ground_truth * GT[j,]
  }

  if(log){
    # initial probabilities
    augmented_model$init_log = model$init %>% log()
    augmented_model$init_log[augmented_model$init_log == -Inf | is.nan(augmented_model$init_log)] = m

    # transition probabilities
    augmented_model$trans_log = model$transition %>% log()
    augmented_model$trans_log[augmented_model$trans_log==-Inf]=m

    # sojourns probabilities
    augmented_model$d_log = augmented_model$d %>% log()
    augmented_model$d_log[augmented_model$d_log==-Inf]=m
    augmented_model$D_log = augmented_model$D %>% log()
    augmented_model$D_log[augmented_model$D_log==-Inf]=m

    # local state probabilities
    augmented_model$b_log = augmented_model$b %>% log()
    augmented_model$b_log[augmented_model$b_log==-Inf]= -1e300

  }else{
    augmented_model$b = augmented_model$b + 0.01
  }
  augmented_model
}



.build_d_from_sojourn_dist <- function(model,M) {
  # shortcuts
  J = model$J
  sojourn.distribution=model$sojourn$type

  # Is d specified as a matrix in the sojourn specifications?
  d_exists_in_sojourn_spec = !is.null(model$sojourn$d)
  # if yes, we use it instead of the parameters.
  if(d_exists_in_sojourn_spec)
  {
    # we check that it has J columns
    if(ncol(model$sojourn$d)!=J) stop("ncol(model$d)!=J")
    # we get the longuest sojourn (runlength), i.e. the number of rows of d
    max_rl = nrow(model$sojourn$d)
    # we extend or shorten the matrix so that it matches the provided M
    if(max_rl >= M) d = head(model$sojourn$d, M)
    else d = rbind(model$sojourn$d, matrix(0, ncol = J, nrow = (M-max_rl)))
  }
  # if d is not specified in the sojourn specification, we build it from the distribution parameters
  else
  {
    # we initialize d to an empty matrix
    d = matrix(nrow=M, ncol=model$J)

    # non-par
    if(sojourn.distribution=='nonparametric' | sojourn.distribution=="ksmoothed-nonparametric") {
      stop("Sojourn distribution (model$sojourn$d) not specified.")
    }
    # poisson
    if(sojourn.distribution=="poisson") {
      if(is.null(model$sojourn$lambda)) stop('Invalid sojourn parameters supplied')
      if(is.null(model$sojourn$shift)) stop('No shift parameter provided for Poisson sojourn distribution (must be at least 1)')
      for(i in 1:J) d[,i] = .dpois.hsmm.sojourn(1:M,model$sojourn$lambda[i],model$sojourn$shift[i])
    }
    # log normal
    if(sojourn.distribution=="lnorm") {
      if(is.null(model$sojourn$meanlog) | is.null(model$sojourn$s.dlog)) stop('Invalid sojourn parameters supplied')
      for(i in 1:J) d[,i] = dlnorm(1:M,model$sojourn$meanlog[i],model$sojourn$s.dlog[i])
    }
    # gamma
    if(sojourn.distribution=="gamma") {
      if(is.null(model$sojourn$shape) | is.null(model$sojourn$scale)) stop('Invalid sojourn parameters supplied')
      for(i in 1:J) d[,i] = dgamma(1:M,shape=model$sojourn$shape[i],scale=model$sojourn$scale[i])
    }
    # log
    if(sojourn.distribution=="logarithmic") {
      if(is.null(model$sojourn$shape)) stop('Invalid sojourn parameters supplied')
      for(i in 1:J) d[,i] = .dlog(1:M,model$sojourn$shape[i])
    }
    # nbinom
    if (sojourn.distribution == "nbinom") {
      if(is.null(model$sojourn$mu) & is.null(model$sojourn$prob)) stop('Invalid sojourn parameters supplied')
      if(is.null(model$sojourn$mu))   for (i in 1:J) d[, i] = .dnbinom.hsmm.sojourn(1:M,size=model$sojourn$size[i],prob=model$sojourn$prob[i],shift=model$sojourn$shift[i])
      if(is.null(model$sojourn$prob)) for (i in 1:J) d[, i] = .dnbinom.hsmm.sojourn(1:M,size=model$sojourn$size[i],mu=model$sojourn$mu[i],shift=model$sojourn$shift[i])
    }
  }

  # finally, we need to normalize so that the sojourn distribution sums to 1 for each state
  d = apply(d,2,function(x) x/sum(x))
  d
}




####### FITTING ---------------------



#' Fit a hidden semi-Markov model to data sequences
#'
#' This function relies on a EM procedure to fit the model parameters to maximize the likelihood of the decoded hidden state sequence.
#' It returns a list whose first element is the fitted model (an object of class \code{hsmm}) and whose second elements provides information about the EM procedure (convergence, number of iteration, likelihood).
#' @param model a \code{hsmm} object. The model whose parameters will be re-estimated.
#' @param X a \code{data.frame} of observations.
#' @param ground_truth (optional) a \code{data.frame} of ground truth, _i.e._ time-points where the hidden state is known. Default is an empty \code{data.frame}, _i.e._ no ground truth.
#' @param n_iter (optional) an integer specifiying the maximal number of iterations for the EM-procedure. Default value is 10.
#' @param rel_tol (optional) a positive double specifying the tolerance at which to stop the EM. If the difference in likelihood (normalized by the total sequences length) between two iterations of the EM is smaller than \code{rel_tol}, then the EM procedure is considered to have converged to a local maximum.
#' @param lock.transition (optional) a logical. Default is \code{FALSE}. Specifies if the transition probability should be locked (kept as is) or re-estimated at the M-step of the EM.
#' @param lock.sojourn (optional) a logical. Default is \code{FALSE}. Specifies if the sojourn distributions should be locked (kept as is) or re-estimated at the M-step of the EM.
#' @param lock.emission (optional) a logical. Default is \code{FALSE}. Specifies if the emission distributions should be locked (kept as is) or re-estimated at the M-step of the EM.
#' @param trust_in_ground_truth (optional) a double between 0 and 1 specifying the trust in the ground truth. A value of 0 indicates no trust and is equivalent to not providing ground-truth. A value of 1 indicates full trust and the ground truth will not be modulated by the probability of the values of the observations.
#' @param alpha (optional) a positive number specifying the strength of the prior TODO: explain better.
#' @keywords HSMM
#' @return A list. First element of the list (\code{$model}) is a \code{hsmm} object (the fitted model) and the second element (\code{$fit_param}) provides information about the EM-procedure. The second element can be visualized by calling the function \code{plot_hsmm_fit_param()}.
#'
#' @export
#' @importFrom magrittr %>%
#' @examples
#' my_model_spec = simple_model_spec  # simple_model_spec is a model attached to the HiddenSemiMarkov package for demos
#' Xsim = simulate_hsmm(my_model_spec, n_state_transitions = 20) # a short sequence is simulated
#'
#' my_model_init = initialize_hsmm(my_model_spec) # the model is initialized
#' my_model_fit = fit_hsmm(model = my_model_init, X = Xsim) # the model is fit to the observations.
#' plot_hsmm_fit_param(model = my_model_fit)
#'
#' viterbi_init = predict_states_hsmm(model = my_model_init, X = Xsim, method = "viterbi") # predict the states with the initial model
#' viterbi_fit = predict_states_hsmm(model = my_model_fit$model, X = Xsim, method = "viterbi") # predict the states with the fit model
#' Xsim$state_viterbi_init = viterbi_init$state_seq$state
#' Xsim$state_viterbi_fit = viterbi_fit$state_seq$state
#' plot_hsmm_seq(X = Xsim, model = my_model_init)


fit_hsmm = function(model, X, ground_truth = data.frame(), n_iter = 10, rel_tol = 1/20,
                    lock.transition = FALSE, lock.sojourn = FALSE, lock.emission = FALSE,
                    trust_in_ground_truth = 0.75, alpha = 200,
                    verbose = FALSE, graphical = FALSE){

  # check/initialize the model
  # TODO
  init_model = model # keep the initial model (XXX necessary?)

  # check the data
  if(missing(X)) stop("X missing!")
  X = .check_data(data = X, model = model)
  if(verbose) cat("Data checked\n")
  X = .augment_data(model = model, X = X, verbose = verbose)
  if(verbose){cat("Data augmented \n")}

  # if any ground_truth, check it
  if(nrow(ground_truth)>0){
    ground_truth = .check_ground_truth(ground_truth, model = model, X = X)
  }else{
    ground_truth = X %>%  dplyr::select(seq_id, t) %>% dplyr::mutate(state = NA)
  }

  if((trust_in_ground_truth < 0) | (trust_in_ground_truth > 1)) stop("trust_in_ground_truth must be a number between 0 and 1.")

  # initializing the vector which keeps track of the log likelihood
  ll = c(); message = "Reached n_iter"
  for(it in 1:n_iter){
    if(verbose){cat("i: ",it,"\n")}

    # E-step
    if(verbose){cat("\t E-step \n")}
    smoothed_res = try(.predict_states_hsmm_forward_backward(model = model, X = X, ground_truth = ground_truth, trust_in_ground_truth = trust_in_ground_truth))
    if(class(smoothed_res) == "try-error"){
      message = paste0("Error in the E-step. Forward-Backward decoding threw an error at iteration ",it,".") #TODO Model from previous iteration is returned.
      break()
    }
    p = smoothed_res$state_probs %>% dplyr::select(-local) %>%
      tidyr::pivot_wider(id_cols = c("seq_id","t"), names_from = "state", values_from = "posterior") %>%
      dplyr::select(-seq_id, -t) %>% as.matrix(., ncol = model$J, nrow = nrow(X))
    p = p/rowSums(p)
    if(graphical) matplot(p, type = "l", lty = 1, col = model$state_colors)

    # M-step
    if(verbose){cat("\t M-step \n")}
    new_model = .re_estimate_parameters(model = model, X = X, p = p,
                                        fwbw_res = smoothed_res$fwbw_res,
                                        ground_truth = ground_truth,
                                        lock.transition = lock.transition,
                                        lock.sojourn = lock.sojourn,
                                        lock.emission = lock.emission,
                                        trust_in_ground_truth = trust_in_ground_truth,
                                        alpha = alpha,
                                        verbose = verbose)
    model = new_model
    # Keep iterating?
    ll[it]= sum(smoothed_res$loglik$loglik)/nrow(X) # sum(log(estep_variables$N))
    if(verbose) cat("\t ll: ",ll[it],"\n")
    if(it>=2) if( abs(ll[it]-ll[it-1]) < rel_tol*(abs(ll[it])+ abs(ll[it-1]))/2 ){message = "Converged"; break()}
  }
  if(verbose) cat("Model fitted\n")
  out = list(model = model, fit_param = list(ll = ll, message = message, n_iter = length(ll)))
  out
}

.re_estimate_parameters = function(model, X,
                                   p, fwbw_res,
                                   ground_truth,
                                   lock.transition = FALSE,
                                   lock.sojourn = FALSE,
                                   lock.emission = FALSE,
                                   trust_in_ground_truth = 0.75, alpha = 200,
                                   verbose = FALSE){
  new_model = model
  if(!lock.emission){
    if(verbose) cat("updating emission parameters\n")
    new_model = .re_estimate_emission_parameters(model = new_model, X = X, p = p, ground_truth = ground_truth, trust_in_ground_truth = trust_in_ground_truth, alpha = alpha, verbose = verbose)
  }
  if(!lock.transition){
    if(verbose) cat("updating transition matrix\n")
    new_model = .re_estimate_transition_probabilities(model = new_model, fwbw_res = fwbw_res)
  }
  if(!lock.sojourn){
    if(verbose) cat("updating sojourn distributions\n")
    new_model = .re_estimate_sojourn_distributions(model = new_model, fwbw_res = fwbw_res, graphical = verbose)
  }
  return(new_model)
}

.re_estimate_emission_parameters = function(model, X, p, ground_truth = ground_truth, trust_in_ground_truth = 0.75, alpha = 200, verbose = FALSE){

  # initialization of the weights matrix
  weights = p

  #### First, we need to combine the ground_truth with the results from the E-step
  weights_ground_truth = weights
  j = which(!is.na(ground_truth$state))
  if(length(j)>0){
    weights_ground_truth[j,] = 0
    weights_ground_truth[cbind(j,ground_truth$state[j])] = 1
  }

  weights = (1-trust_in_ground_truth) * weights + trust_in_ground_truth * weights_ground_truth

  state_seq = apply(weights,1,which.max)

  # Check for un-visited states or rarely visited states
  if(verbose){
    if(any(table(state_seq)<(nrow(X)/model$J/10))){warning("Some states are rarely visited\n")}
    if(any(table(state_seq)==0)){warning("Some states are never visited\n")}
  }

  #### Then we re-estimate the emission (and missingness) parameters (for observation simulation / variable generation)
  new_model = model
  new_model$parms.emission = .re_estimate_emission_distributions(model = model, X = X, w = weights)


  #### Finally, we re-estimate the observation probabilities
  new_model$obs_probs = .re_estimate_obs_probs(model = model, X = X, w = weights, alpha = alpha)

  ##### Return the new model with updated emission parameters
  new_model
}

.re_estimate_emission_distributions = function(model = model, X = X, w = w){

  # keeping only the observations (= variables for which we have parms.emission)
  obs_names = names(model$parms.emission)
  obs = X %>% dplyr::select(dplyr::all_of(obs_names)) %>% as.data.frame()

  new_parem = .compute_new_em_parms(obs = obs, w = w, parem = model$parms.emission)
  new_parem
}


#library(SDMTools)

.compute_new_em_parms = function(obs, w, parem){
  new.parem = parem
  for(i in 1:length(parem)){
    #cat(i, "\n")
    if(new.parem[[i]]$type == "norm"){
      # weighted mean
      x_bar = sapply(1:ncol(w),function(state) weighted.mean(obs[,i],w = w[,state], na.rm = TRUE))
      x_bar[is.na(x_bar)] = parem[[i]]$params$mean[is.na(x_bar)]
      # weighted sd
      # sd_bar = sapply(1:ncol(w),function(state) wt.sd(obs[,i],w = w[,state]))
      sd_bar = sapply(1:ncol(w), function(state) sqrt( sum( w[,state]  * (obs[,i] - x_bar[state])^2)))
      sd_bar[is.na(sd_bar)] = parem[[i]]$params$sd[is.na(sd_bar)]
      # variables and hyper-parameters
      n = apply(w,2,sum)
      mu_0 = parem[[i]]$params$mu_0
      n0 = parem[[i]]$params$n0
      new_n0 = n + n0
      new_alpha = parem[[i]]$params$alpha + n/2
      new_beta = parem[[i]]$params$beta + (n-1)/2 * sd_bar^2 + n * n0 / (n + n0) * (x_bar - mu_0)^2 / 2

      # posteriors
      new_mean = (n0 * mu_0 + n*x_bar)/(n0 + n)
      new_sd = sqrt( new_beta * (new_n0 + 1) / (new_n0 * new_alpha) )

      # updating the model
      new.parem[[i]]$params$mean = new_mean
      new.parem[[i]]$params$sd = new_sd
    }else if(new.parem[[i]]$type == "binom"){

      # observed proportions
      obs_prob = sapply(1:ncol(w),function(state) weighted.mean(obs[,i]/(parem[[i]]$param$size[state]),w = w[,state], na.rm = TRUE)) # would be better with an EM approach, but it's good enough for now
      obs_prob[is.na(obs_prob)] = parem[[i]]$params$prob[is.na(obs_prob)]
      n = apply(w,2,sum)
      obs_successes = n*obs_prob
      # hyper-parameters
      new_alpha = parem[[i]]$params$alpha + obs_successes
      new_beta = parem[[i]]$params$beta + n * parem[[i]]$param$size  - obs_successes

      # posterior
      new_prob = new_alpha / (new_alpha + new_beta)

      # update the model
      new.parem[[i]]$params$prob = new_prob

    }else if(new.parem[[i]]$type == "non-par"){
      # observed proportions
      obs_probs = sapply(
        1:ncol(w), # for each state
        function(state) {
          tt = parem[[i]]$params$probs[,state]
          j = which(w[,state]>0)
          if(length(j)>0) tt = table(sample(obs[j,i], prob = w[j,state], size = nrow(obs), replace = TRUE))
          tt = tt/sum(tt)
          return(t(tt))
        }
      )

      # variables
      n = apply(w,2,sum)
      n0 = parem[[i]]$params$n0
      probs_0 = parem[[i]]$params$probs_0
      # posterior
      new_probs = n0 * probs_0 + t(n * t(obs_probs))
      new_probs = t(t(new_probs)/colSums(new_probs))

      # updating the model
      new.parem[[i]]$params$probs = new_probs
    }else{ stop("This type of distribution has not been handled yet")}
    # we also need to impute the missingness of this variable
    # TODO
    new_missing_prob = sapply(1:ncol(w), function(state) weighted.mean(is.na(obs[,i]), w = w[,state], na.rm = FALSE))
    new_missing_prob[is.na(new_missing_prob)] = parem[[i]]$missing_prob[is.na(new_missing_prob)]
    new.parem[[i]]$missing_prob = new_missing_prob # %>% pmin(.,0.9999) %>% pmax(.,0.0001)
  }
  return(new.parem)
}


.re_estimate_obs_probs = function(model, X, w, alpha = 200){

  # first we need to format the data, i.e. have one columns specifiying the state and another giving the weight for that state
  Xx = X[rep(1:nrow(X), model$J),]
  Xx$state = rep(1:model$J, each = nrow(X))
  Xx$w = as.vector(w)

  # then check the observed probabilities
  obs_prob = .compute_observed_probabilities(model = model, X = Xx, verbose = verbose)
  #print(obs_prob %>%  dplyr::arrange(state))


  # and combine them with the prior distributions via conjugate prior formula for categorical variables
  P0 = model$obs_probs_0 %>% dplyr::select(-p_n) %>%  dplyr::rename(p0 = p)
  PX = obs_prob %>% dplyr::rename(px = p)
  new_obs_probs = dplyr::full_join(P0, PX, by = intersect(colnames(P0), colnames(PX)),
                            na_matches = "na") %>%
    dplyr::select(-p_n,-p_max) %>%
    dplyr::mutate(p0 = p0 %>% tidyr::replace_na(0),
           px = px %>% tidyr::replace_na(0),
           n = n %>% tidyr::replace_na(0),
           alpha = alpha) %>%
    dplyr::group_by(state) %>%
    dplyr::mutate(Tot = sum(n, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(p = (p0 * alpha + Tot * px)/(alpha + Tot)) %>%
    dplyr::group_by(state) %>%
    dplyr::mutate(p_max = max(p),
           p_n = p / p_max) %>%
    dplyr::ungroup()

  new_obs_probs
}




.get_all_possible_levels = function(model = model, with_missing = TRUE, continuous_var_binned = FALSE){
  lv = list()
  for(var in names(model$parms.emission)){
    if(model$parms.emission[[var]]$type == "non-par") lv[[var]] = model$parms.emission[[var]]$params$values
    if(model$parms.emission[[var]]$type == "binom") lv[[var]] = 0:max(model$parms.emission[[var]]$params$size)
    if(model$parms.emission[[var]]$type == "norm") lv[[var]] = model$parms.emission[[var]]$breaks[!is.infinite(model$parms.emission[[var]]$breaks)]
    if(continuous_var_binned & (model$parms.emission[[var]]$type == "norm")) lv[[var]] = levels(cut(lv[[var]], breaks = model$parms.emission[[var]]$breaks))
    if(with_missing) lv[[var]] = c(lv[[var]], NA)
  }
  if(length(model$augment_data_fun(get_var_names_only = TRUE))>0){
    augm_var = model$augment_data_fun(get_var_names_only = TRUE)
    augm_var_val = model$augment_data_fun(get_var_types_only = TRUE)
    for(var in augm_var){
      if(augm_var_val[[var]]$type == "cat") lv[[var]] = augm_var_val[[var]]$values
      if(augm_var_val[[var]]$type != "cat") lv[[var]] = levels(cut(seq(0,1,by = 0.2), breaks = c(-Inf, seq(0.2,0.8,by = 0.2),Inf)))
    }
  }
  max_n = max(lengths(lv))
  all_levels = data.frame(row.names = 1:max_n, stringsAsFactors = FALSE)
  for(var in names(lv)) all_levels = cbind(all_levels, var = rep(lv[[var]],max_n)[1:max_n])
  colnames(all_levels) = names(lv)
  all_levels = .check_types_and_values(data = all_levels, parem = model$parms.emission, continuous_var_binned = continuous_var_binned)
  all_levels
}



.re_estimate_transition_probabilities = function(model, fwbw_res){

  new_model = model
  new_model$init = pmax(0, fwbw_res$init)
  new_model$transition = pmax(0, fwbw_res$transition) %>% matrix(.,ncol = model$J)

  new_model
}


.re_estimate_sojourn_distributions = function(model, fwbw_res, graphical = FALSE){

  new_model = model
  # shortcuts
  sojourn_distribution = model$sojourn$type

  # the smoothed (forward/backward) algorithm returns a new "d" matrix
  # which needs to be normalized such that the density distribution of runlength for any state sums to 1.
  ds = fwbw_res$eta %>% matrix(.,ncol = model$J) %>% apply(.,2,function(x) x/sum(x))
  ds_i = ds
  M = nrow(ds)

  # We can re-use this ds as is if the sojourn distribution family is non-parametric.
  # If the sojourn distribution is defined as a parametric distribution, we can use ds to estimate the parameters of these distributions

  if(sojourn_distribution == "nonparametric"){
    new_model$sojourn$d = ds
  }else if(sojourn_distribution == "ksmoothed-nonparametric"){
    ds = ds + 1e-100
    ksmooth.thresh = 1e-20 #this is a threshold for which d(u) values to use - if we throw too many weights in the default density() seems to work quite poorly
    for(i in 1:model$J){
      u = which(ds[,i]>ksmooth.thresh)
      ds[,i] = density(u, weights = ds[u,i], from = 1, n = M) %>% approx(.,xout=1:nrow(ds)) %>% pluck("y") %>% tidyr::replace_na(0)
      ds[,i] = ds[,i]/sum(ds[,i])
    }
    new_model$sojourn$d = ds

  }else if(sojourn_distribution == "poisson"){
    new_model$sojourn$lambda = numeric(model$J)
    new_model$sojourn$shift = numeric(model$J)
    shiftthresh = 1e-20 #threshold for effective "0" when considering d(u)
    for(i in 1:model$J) {
      eta = ds[,i]
      u = which(eta>shiftthresh); maxshift =  min(u); Mtmp = max(u) ; U = maxshift:Mtmp
      shifts = sapply(1:maxshift, function(shift) .dpois.hsmm.sojourn(x = U,lambda=(U-shift)%*%eta[U],shift=shift,log=TRUE)%*%eta[U])
      shift = which.max(shifts); new_U = shift:Mtmp
      new_model$sojourn$shift[i] = shift
      new_model$sojourn$lambda[i] = (new_U-shift)%*%eta[new_U]
      ds[,i] = .dpois.hsmm.sojourn(1:M,new_model$sojourn$lambda[i],new_model$sojourn$shift[i])
    }

  }else if(sojourn_distribution == "nbinom"){
    new_model$sojourn$size = numeric(model$J)
    new_model$sojourn$shift = integer(model$J)
    new_model$sojourn$mu = numeric(model$J)
    new_model$sojourn$prob = numeric(model$J)
    for(i in 1:model$J) {
      tmp = .fitnbinom(ds[,i])
      new_model$sojourn$shift[i] = tmp['shift']
      new_model$sojourn$size[i] =  tmp['size']
      new_model$sojourn$mu[i] =  tmp['mu']
      new_model$sojourn$prob[i] =  tmp['prob']
      ds[,i] =  .dnbinom.hsmm.sojourn(1:M,tmp['size'],tmp['prob'],tmp['shift'])
    }

  }else if(sojourn_distribution == "gamma"){
    new_model$sojourn$shape = numeric(model$J)
    new_model$sojourn$scale = numeric(model$J)
    for(i in 1:model$J) {
      tmp = gammafit(1:nrow(ds),wt=ds[,i])
      new_model$sojourn$shape[i] = tmp$shape
      new_model$sojourn$scale[i] = tmp$scale
      ds[,i] = dgamma(1:M,shape=tmp$shape,scale=tmp$scale)
    }

  }else if(sojourn_distribution == "logarithmic"){
    ds = ds+1e-100
    new_model$sojourn$shape = numeric(model$J)
    for(i in 1:model$J) {
      new_model$sojourn$shape[i] = .logdistrfit(wt=ds[,i])
      ds[,i] = .dlog(1:M,new_model$sojourn$shape[i])
    }

  }else if(sojourn_distribution == "lnorm"){
    new_model$sojourn$meanlog = numeric(model$J)
    new_model$sojourn$s.dlog = numeric(model$J)
    for(i in 1:model$J) {
      new_model$sojourn$meanlog[i] = weighted.mean(log(1:M),ds[,i])
      new_model$sojourn$s.dlog[i] = sqrt(cov.wt(data.frame(log(1:M)),ds[,i])$cov)
      ds[,i] = dlnorm(1:M,new_model$sojourn$meanlog[i],new_model$sojourn$s.dlog[i])
    }
  }else{
    stop("Invalid sojourn distribution")
  }

  ds = apply(ds,2,function(x) x/sum(x))

  if(graphical){
    ds_o = .build_d_from_sojourn_dist(model, M = M)
    par(mfrow = c(5,4), mar = c(0,0,0,0)+1.8)
    for(i in 1:model$J){
      max_y = max(c(ds_o[1:100,i],ds_i[1:100,i],ds[1:100,i]), na.rm = TRUE)
      plot(ds_o[1:100,i], type = "l", col = "black", ylim = c(0,max_y), ylab = "", xlab = "", main = i)
      points(ds_i[1:100,i], type = "l", col = "blue")
      points(ds[1:100,i], type = "l", col = "red", lty = 2)
    }
    par(mfrow = c(1,1))
  }

  #new_model$d = ds
  #new_model$D = ds %>%  apply(. ,2,function(x) rev(cumsum(rev(x))))
  new_model
}



####### SUMMARY ---------------------

#' Prints the model parameters.
#'
#' @param model a \code{hsmm} or \code{hsmm_spec} object specifying the hidden semi-Markov model.
#' @export
#'
summary_hsmm <- function(model) {
  cat("\nStarting distribution = \n")
  print(model$init,2)
  cat("\nTransition matrix = \n")
  print(model$transition,2)
  cat("\nSojourn distribution parameters = \n")
  print(model$sojourn)
  cat("\nEmission distribution parameters = \n")
  print(model$parms.emission)
}


