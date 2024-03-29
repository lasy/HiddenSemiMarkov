

######### MODEL SPECIFICATION CHECKS #############

.check_J = function(J){
  if(is.null(J)) stop("J (the number of states) must be specified.\n")
  if((round(J) != J) | (J <= 0)) stop("J must be an integer, strictly positive.\n")
  J
}

.check_init = function(init, J){
  if(is.null(init)){
    warning("Initial probabilities not specified. Assuming uniform distribution across states.\n")
    init = rep(1, J)/J
  }
  if(length(init) != J) stop("Initial probability vector `init` must have J elements.\n")
  if(round(sum(init),10) != 1) warning("Initial probabilities did not sum to one. Values are normalized such that the initial probabilities sum to 1.\n")
  init = init/sum(init)
  init
}

.check_transitions = function(transition, J){
  if(is.null(transition)){
    warning("Transition probabilities not specified. Assuming uniform transition probabilities without any absorbing state.\n")
    transition = (matrix(1, J, J) - diag(1, J, J))/(J-1)
  }
  if(NROW(transition) != J)  stop("transition matrix `transition` must have J rows.\n")
  if(NCOL(transition) != J)  stop("transition matrix `transition` must have J columns.\n")
  if(any(round(rowSums(transition),10) != 1))
    warning("Transition matrix rows do not sum to 1. Values will be normalized such that the transition probabilities from any state sum to 1.")
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
  transition
}


.check_sojourn = function(sojourn, J, state_names){

  if(is.null(sojourn)){
    sojourn = list(type = "nonparametric",
                   d = matrix(1/10, 10, J))
    warning("Sojourn distributions unspecified.
        Assuming uniform distribution ('nonparametric') over 10 time-steps. \n")
  }


  # sojourn distributions can be specified in 2 different ways:
  # 1. one sojourn type for all states
  # 2. one sojourn type for each state

  if(is.null(sojourn$type) & (length(sojourn) != J)) stop("Sojourn distribution type not specified.")

  if("type"  %in% names(sojourn)){
    supported_sojourns = unique(available_sojourn$distribution_type)
    if (!(sojourn$type %in% supported_sojourns))
      stop(paste("Invalid sojourn type specified (",sojourn$type,"). Must be one of: ", paste(supported_sojourns, collapse = ", ")))

    sojourn = .check_sojourn_all_states(sojourn, J)

    sojourn_per_state =
      purrr::map(
        .x = 1:J,
        .f = function(j){
          s = list()
          s$type = sojourn$type
          for(k in 2:length(sojourn))
            if((sojourn$type %in% c("nonparametric","ksmoothed_nonparametric")) & (names(sojourn)[k] == "d"))
              s[[k]] = sojourn[[k]][,j]
          else s[[k]] = sojourn[[k]][j]
          names(s) = names(sojourn)
          s
        }
      )
  }else{
    sojourn_per_state = purrr::map(
      .x = 1:J,
      .f = function(j){
        s = .check_sojourn_single_state(sojourn[[j]])
        s
      }
    )
  }
  names(sojourn_per_state) = state_names
  sojourn = sojourn_per_state
  sojourn
}


.check_sojourn_all_states = function(sojourn, J){

  # first we check it has all the parameters
  params = .get_required_params(sojourn_type = sojourn$type)
  if(!(
    all(params$required_params %in% names(sojourn))|
    all(params$backup_required_params %in% names(sojourn))
  )) stop("Sojourn distributions are not specified as they should.
            Type ?specify_hsmm for examples and details.\n")


  if(length(params$optional_params)!=0){
    if(("shift" %in% params$optional_params) & is.null(sojourn$shift)) sojourn$shift = rep(1, J)
    if(("bw" %in% params$optional_params) & is.null(sojourn$bw)) sojourn$bw = rep("nrd0", J)
  }

  # then we check the dimensions and/or specifications
  if(sojourn$type %in% c("nonparametric","ksmoothed_nonparametric")){
    if(ncol(sojourn$d) != J)
      stop("Non-parametric sojourn distributions must be specified by a matrix d: sojourn = list(type = 'nonparametric', d = ...matrix of dim M x J...) where M is the longuest sojourn duration. \n")
    if(any(is.na(sojourn$d))) stop("NA values found in the sojourn distribution matrix d.\n")
    if(any(colSums(sojourn$d) == 0)) stop("At least one state does not have a valid sojourn distribution. Check your matrix d\n")
    if(any(round(colSums(sojourn$d),10) != 1)){
      warning("The provided sojourn distributions do not sum to 1 for at least one state. Sojourn distributions are normalized so that they sum to 1.")
      sojourn$d = t(t(sojourn$d)/colSums(sojourn$d))
    }
    if((sojourn$type == "ksmoothed_nonparametric") & (length(sojourn$bw) != J))
      if(length(sojourn$bw) == 1) sojourn$bw = rep(sojourn$bw, J) else stop("The 'bw' parameter of the non-parametric sojourn distribution must be NULL (default) or of length 1 or J.\n")
  }else{
    for(k in 2:length(sojourn)){
      if(any(is.na(sojourn[[k]])))
        stop("NA values found in the sojourn distribution parameters.\n")
      if(length(sojourn[[k]]) == 1){
        warning(paste0("Sojourn parameter '",names(sojourn)[k],"' is of length 1. Assuming the same value for all states.\n"))
        sojourn[[k]] = rep(sojourn[[k]], J)
      }
      if(length(sojourn[[k]]) != J)
        stop(paste0("Sojourn parameter '",names(sojourn)[k],"' does not have as many element as there are states in the model.\n"))
    }
  }

  sojourn
}

.get_required_params = function(sojourn_type){

  j = (available_sojourn$distribution_type == sojourn_type)
  required_params = available_sojourn$parameters[which(j & available_sojourn$required_parameter)]
  optional_params = available_sojourn$parameters[which(j & !available_sojourn$required_parameter)]

  list(required_params = required_params,
       optional_params = optional_params)
}

.check_sojourn_single_state = function(sojourn_one_state){

  # first we check it has all the parameters
  params = .get_required_params(sojourn_type = sojourn_one_state$type)
  if(!(
    all(params$required_params %in% names(sojourn_one_state))|
    all(params$backup_required_params %in% names(sojourn_one_state))
  )) stop("Sojourn distributions are not specified as they should.
            Type ?specify_hsmm for examples and details.\n")

  if(length(params$optional_params)!=0){
    if(("shift" %in% params$optional_params) & is.null(sojourn_one_state$shift)) sojourn_one_state$shift = 1
    if(("bw" %in% params$optional_params) & is.null(sojourn_one_state$bw)) sojourn_one_state$bw = "nrd0"
  }


  # then we check the dimensions and/or specifications
  if(sojourn_one_state$type %in% c("nonparametric","ksmoothed_nonparametric")){
    if(any(is.na(sojourn_one_state$d))) stop("NA values found in the sojourn distribution vector d.\n")
    if(sum(sojourn_one_state$d) == 0) stop("At least one state does not have a valid sojourn distribution. Check your matrix or vector d\n")
    if(round(sum(sojourn_one_state$d),10) != 1){
      warning("The provided sojourn distributions do not sum to 1 for at least one state. Sojourn distributions are normalized so that they sum to 1.")
      sojourn_one_state$d = sojourn_one_state$d/sum(sojourn_one_state$d)
    }
    if((sojourn_one_state$type == "ksmoothed_nonparametric") & (length(sojourn_one_state$bw) != 1)){
      warning(paste0("Sojourn parameter 'bw' is longer than one. Taking the first value only.\n"))
      sojourn_one_state$bw  = sojourn_one_state$bw[1]
    }
  }else{
    for(k in 2:length(sojourn_one_state)){
      if(length(sojourn_one_state[[k]]) != 1){
        warning(paste0("Sojourn parameter '",names(sojourn)[k],"' is longer than one. Taking the first value only.\n"))
        sojourn_one_state[[k]]  = sojourn_one_state[[k]][1]
      }
      if(is.na(sojourn_one_state[[k]]))
        stop("NA values found in the sojourn distribution parameters.\n")
    }
  }

  sojourn_one_state
}


# .check_nonparametric_sojourn_all_states = function(sojourn, J){
#   if(is.null(sojourn$d) || (ncol(sojourn$d) != J))
#     stop("Non-parametric sojourn distributions must be specified by a matrix d: sojourn = list(type = 'nonparametric', d = ...matrix of dim M x J...) where M is the longuest sojourn duration. \n")
#   if(any(is.na(sojourn$d))) stop("NA values found in the sojourn distribution matrix d.\n")
#   if(any(colSums(sojourn$d) == 0)) stop("At least one state does not have a valid sojourn distribution. Check your matrix d\n")
#   if(any(round(colSums(sojourn$d),10) != 1)){
#     warning("The provided sojourn distributions do not sum to 1 for at least one state. Sojourn distributions are normalized so that they sum to 1.")
#     sojourn$d = t(t(sojourn$d)/colSums(sojourn$d))
#   }
#   sojourn
# }




.check_marg_em_probs = function(marg_em_probs, J){
  if(is.null(marg_em_probs)){
    marg_em_probs = list(
      obs = list(
        type = "norm",
        params = list(mean = 1:J, sd = rep(1, J))
      )
    )
    warning("Marginal emission probabilities unspecified. Assuming a model for a single variable normally distributed around 'j' (the state index) and with a unitary standard deviation.\n")
  }

  if(length(marg_em_probs) == 0)
    stop("Emission probabilities (marg_em_probs) are not specified.
    They should be specified as a list with one element per observed variable.
    Each element of this list should specify the name of that variable, its distribution family (e.g. 'norm', 'binom', 'non-par')
    and its parameters/distribution for each state.
         Type ?specify_hsmm for examples.")
  if(is.null(names(marg_em_probs))) names(marg_em_probs) = 1:length(marg_em_probs)

  for(var in names(marg_em_probs)){
    var_parem = marg_em_probs[[var]]
    params_names = names(var_parem$params)

    if(var_parem$type == "norm"){
      if(!all(c( "mean", "sd" ) %in% params_names))
        stop(paste0("Variable '",var,"' is 'norm' and must have the following params: \n",
                    "marg_em_probs$",var,"$params = list(\nmean = ...vector of J means...,",
                    "\nsd = ...vector of J sd...,"))
      marg_em_probs[[var]]$breaks =
        c(-Inf,
          seq(min(marg_em_probs[[var]]$params$mean - 3*marg_em_probs[[var]]$param$sd),
              max(marg_em_probs[[var]]$params$mean + 3*marg_em_probs[[var]]$param$sd),
              len = 8),
          Inf)

      if(is.null(var_parem$viz_options)) marg_em_probs[[var]]$viz_options = list()
      if(is.null(var_parem$viz_options$color_low)) marg_em_probs[[var]]$viz_options$color_low = "steelblue2"
      if(is.null(var_parem$viz_options$color_high)) marg_em_probs[[var]]$viz_options$color_high = "indianred1"
      if(is.null(var_parem$viz_options$color_mid)) marg_em_probs[[var]]$viz_options$color_mid = "gray80"
      if(is.null(var_parem$viz_options$mid_value)) marg_em_probs[[var]]$viz_options$mid_value = 0

    }
    if(var_parem$type == "binom"){
      if(!all(c( "size", "prob" ) %in% params_names))
        stop(paste0("Variable '",var,"' is 'binom' and must have the following params: \n",
                    "marg_em_probs$",var,"$params = list(",
                    "\nsize = ...vector of J size...,",
                    "\nprob = ...vector of J prob...,"))

      if(length(var_parem$params$size) != J)
        if(length(var_parem$params$size) == 1) var_parem$params$size = rep(var_parem$params$size, J)
        else stop(paste0("Variable '",var,"' ('binom') 'size' parameter must be of length 1 or J.\n"))

        if(length(var_parem$params$prob) != J)
          stop(paste0("Variable '",var,"' ('binom') 'prob' parameter must be of length J.\n"))



        if(is.null(var_parem$viz_options)) marg_em_probs[[var]]$viz_options = list()
        if(is.null(var_parem$viz_options$color_max)) marg_em_probs[[var]]$viz_options$color_max = "indianred1"

    }
    if(var_parem$type == "non-par"){
      if(!all(c("values", "probs" ) %in% params_names))
        stop(paste0("Variable '",var,"' is 'non-par' and must have the following params: \n",
                    "marg_em_probs$",var,"$params = list(\nvalues = ...vector of values..., \nprobs = ...matrix (nrow = # of values, ncol = # of states)..., \nalpha (optional) = ...fitting hyper-parameter for conjugate prior...)"))

      if(is.null(var_parem$viz_options)) marg_em_probs[[var]]$viz_options = list()
      if(is.null(var_parem$viz_options$colors)) marg_em_probs[[var]]$viz_options$colors = rainbow(n = length(var_parem$params$values), s = 0.8, v = 0.9)
    }
    if(var_parem$type == "beta"){
      if(!all(c("shape1", "shape2" ) %in% params_names))
        stop(paste0("Variable '",var,"' is of type 'beta' and must have the following params: \n",
                    "marg_em_probs$",var,"$params = list(",
                    "\nshape1 = ...vector of size J...,",
                    "\nshape2 = ...vector of size J...)\n see ?dbeta for the interpretation of the parameters 'shape1' and 'shape2'."))

      if(is.null(marg_em_probs[[var]]$breaks)) marg_em_probs[[var]]$breaks = seq(0,1, by = 0.1) + c(-1e-300, rep(0,10))
      if(any(c(marg_em_probs[[var]]$params$shape1,marg_em_probs[[var]]$params$shape1) < 1))
        stop(paste0("Shape parameters lower than 1 for beta distributions are currently not supported (error identified for variable '",var,"' of type 'beta')"))

      if(is.null(var_parem$viz_options)) marg_em_probs[[var]]$viz_options = list()
      if(is.null(var_parem$viz_options$color_low)) marg_em_probs[[var]]$viz_options$color_low = "steelblue2"
      if(is.null(var_parem$viz_options$color_high)) marg_em_probs[[var]]$viz_options$color_high = "indianred1"
      if(is.null(var_parem$viz_options$color_mid)) marg_em_probs[[var]]$viz_options$color_mid = "gray80"
    }
    if(!(var_parem$type %in% available_marginal_emission_dist()$distribution))
      stop(paste0("This marginal emission distribution (",var_parem$type,") hasn't been implemented yet. Type available_marginal_emission_dist() to retrieve the list of available marginal emission distributions."))
  }
  marg_em_probs
}


.check_censoring_probs = function(censoring_probs, J, K){

  if(is.null(censoring_probs)){
    censoring_probs = list(
      p = rep(0, J),
      q = matrix(0.5, nrow = K, ncol = J),
      missing_prob_specified = FALSE
    )
  }else{
    censoring_probs$missing_prob_specified = TRUE

    if(is.null(censoring_probs$p)){
      censoring_probs$p = rep(0, J)
    }else{
      if(length(censoring_probs$p) == 1) censoring_probs$p = rep(censoring_probs$p, J)
      if(length(censoring_probs$p) != J) stop("censoring_probs$p should be a numeric or a vector of length J.")
      if(any((censoring_probs$p > 1)|(censoring_probs$p < 0))) stop("censoring_prob$p should have values between 0 and 1.")
    }

    if(is.null(censoring_probs$q)){
      censoring_probs$q = matrix(0, nrow = K, ncol = J)
    }else{
      if((!is.matrix(censoring_probs$q)) ||
         (nrow(censoring_probs$q) != K) ||
         (ncol(censoring_probs$q) != J)) stop("censoring_probs$q should be a matrix of dimension K x J where K is the number of variables (i.e. length(marg_em_probs))")
      if(any((censoring_probs$q > 1) | (censoring_probs$q < 0))) stop("censoring_probs$q should have values between 0 and 1.")
    }
  }
  censoring_probs = list(p = censoring_probs$p,
                         q = censoring_probs$q,
                         missing_prob_specified = censoring_probs$missing_prob_specified)
  censoring_probs
}

.check_state_names = function(state_names, J){
  if(is.null(state_names)) state_names = 1:J
  if(length(state_names) != J) stop("state_names must be of length J")
  if(!is.character(state_names)) state_names = as.character(state_names)
  state_names
}


.check_state_colors = function(state_colors, J){
  if(is.null(state_colors)) state_colors = viridis::viridis_pal(option = "C")(J)
  if(length(state_colors) != J) stop("state_colors must be of length J")
  # TODO check if it's a valid color representation
  state_colors
}


######### CHECK FOR PREDICTIONS AND FITS #############

.check_model = function(model = model){

  if(class(model) != "hsmm") stop("model must be of class 'hsmm'.")


  model$J = .check_J(model$J)
  model$state_names = .check_state_names(model$state_names, model$J)
  model$state_colors = .check_state_colors(model$state_colors, model$J)

  model$init = .check_init(model$init, model$J)
  model$transition = .check_transitions(model$transition, model$J)

  model$sojourn = .check_sojourn(model$sojourn, model$J, model$state_names)

  model$obs_probs = .check_model_obs_probs(model)
  model$b = .check_model_b(model)

  model
}

.check_model_obs_probs = function(model){
  obs_probs = model$obs_probs
  if(any(is.na(obs_probs$p))) stop("NAs found in the emission probabilities (model$obs_probs$p).\n")
  sum_per_state = ave(obs_probs$p, obs_probs$state, FUN = sum) %>% unique()
  if(any(round(sum_per_state, 10) != 1)) stop("Emission probabilities (model$obs_probs$p) do not sum to 1 for all states.\n")
  obs_probs
}

.check_model_b = function(model){
  b = model$b
  b_mat = b[,paste0("p_",1:model$J)]
  if(any(is.na(b_mat))) stop("NAs found in the emission probabilities (model$b).\n")
  sum_per_state = colSums(b_mat)
  if(any(round(sum_per_state, 10) != 1)) stop("Emission probabilities (model$obs_probs$p) do not sum to 1 for all states.\n")
  b
}


#' @export
.check_data = function(data, model){
  if(!("data.frame" %in% class(data))) stop("data must be a data.frame\n")
  if(!("seq_id" %in% colnames(data))){
    warning("column 'seq_id' is missing. Assuming a single sequence.")
    data$seq_id = 1
  }
  if(!("t" %in% colnames(data))){
    warning("column 't' is missing. Assuming no missing time-points and ordered data.frame.")
    data = data %>% dplyr::group_by(seq_id) %>% dplyr::mutate(t = row_number())
  }
  # checking that the data has a column for each model variable
  j = which(!(names(model$marg_em_probs) %in% colnames(data)))
  if(length(j)>0) stop(stringr::str_c("The data must contain a column for each of the model variable. Missing variables: ",stringr::str_c(names(model$marg_em_probs)[j],collapse = ", ")))

  # making sure the data is ordered
  data = data %>% dplyr::arrange(seq_id, t) %>% dplyr::ungroup()

  # checking that we have time-steps of duration 1 for each sequence
  if(any(round(data$t) != data$t)) stop("The column 't' of the data must be round numbers\n")
  tmp = data %>% dplyr::group_by(seq_id) %>% dplyr::mutate(time_step = t - lag(t), time_step = tidyr::replace_na(time_step, 1))
  if(any(tmp$time_step != 1)) stop("Some sequences in the provided data have non-consecutive time-points.")


  # select only the columns we need
  var_names = names(model$marg_em_probs)
  data = data %>% dplyr::select(seq_id, t, dplyr::all_of(var_names))

  # checking the type and values of each variable
  data = .check_types_and_values(data = data, parem = model$marg_em_probs)
  data
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






