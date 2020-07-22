


#' Simulate a state and observation sequence from a hidden semi-Markov model.
#'
#' This function returns a \code{data.frame} with a simulated state sequence and observations.
#' The length of the simulated sequence can be specified either as a number of state transition (with argument \code{n_state_transitions}) or as a number of time-point (with argument \code{n_timepoints})
#' Each row of the returned \code{data.frame} is a time-point and the columns specify the sequence id, the time-point, the state and the values of the observation. Each variable has its own column.
#'
#' @param model a object of class \code{hsmm_spec} or \code{hsmm}.
#' @param seq_id (optional) a character specifying the name of the sequence to simulate.
#' @param n_state_transitions (optional) XXX
#' @param n_timepoints (optional) XXX
#' @param seed (optional) XXX
#' @param all_states (optional) XXX
#' @param min_tp_per_state (optional) XXX
#' @param mult_seq_allowed (optional) XXX
#'
#' @keywords HSMM
#' @return A data.frame with the following columns: TODO xxxx
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


simulate_hsmm = function(model,
                         seq_id = NULL,
                         n_state_transitions = NULL, n_timepoints = NULL,
                         seed = NULL,
                         all_states = FALSE,  min_tp_per_state = NULL, mult_seq_allowed = TRUE,
                         state_seq = NULL ){

  ##### What this function does:
  # 0. checks the input variables
  # 1. simulates a sequence of state transitions from the initial and transition probabilities.
  # 2. samples from the sojourn distributions to generate the state sequence.
  # 3. samples from the emission parameter distributions to generate the observations.
  # 4. truncates the simulated sequence to n_timeponts if it was specified

  ##### Implementation

  # 0. CHECKING THE INPUT PARAMETERS
  if(!is.null(seed)) set.seed(seed)
  # checking the input parameters
  if(is.null(n_state_transitions) & is.null(n_timepoints)){
    warning("sequence length not specified, simulating for 100 state transitions")
    n_state_transitions = 100
  }
  if(is.null(n_state_transitions) & !is.null(n_timepoints)) n_state_transitions = n_timepoints
  # checking the model
  # TODO
  # sojourn matrix
  d = .build_d_from_sojourn_dist(model, M = ifelse(!is.null(n_timepoints),n_timepoints,get_longest_sojourn(model)))
  # seq_id
  if(is.null(seq_id) | (typeof(seq_id) != "character")) seq_id = stringr::str_c("sim_seq ",Sys.time())


  # 1 and 2. SIMULATING A SEQUENCE OF STATES and SAMPLING THE SOJOURN DISTRIBUTIONS
  s = sim.mc(init = model$init, transition = model$transition, N = n_state_transitions)
  sojourns = sapply(s, function(x) sample(1:nrow(d), size = 1, prob = d[,x], replace = TRUE))
  df = data.frame(s = s %>% factor(.,levels = 1:model$J), sojourns = sojourns, seq_id = seq_id)


  if(all_states){
    seq_nb = 0
    df$seq_id = stringr::str_c(df$seq_id,"_",seq_nb)
    n_tp_per_state = df %>% dplyr::group_by(s, .drop = FALSE) %>% dplyr::summarize(n_tp = sum(sojourns),.groups = "drop")

    if(is.null(min_tp_per_state)) min_tp_per_state = 1

    while(any(n_tp_per_state$n_tp < min_tp_per_state)){
      seq_nb = seq_nb + 1
      init = rep(0,model$J)
      if(!mult_seq_allowed) init[last(s)] = 1 else init[which.min(n_tp_per_state$n_tp)] = 1
      add_s = sim.mc(init = init, transition = model$transition, N = n_state_transitions)
      add_sojourns = sapply(add_s, function(x) sample(1:nrow(d), size = 1, prob = d[,x], replace = TRUE))
      ns = nrow(df)
      s = c(s[-ns], add_s) # we remove the last element of s in case !mult_seq_allowed
      sojourns = c(sojourns[-ns], add_sojourns)
      df = bind_rows(df[-ns,],
                     data.frame(
                       s = add_s %>% factor(.,levels = 1:model$J),
                       sojourns = add_sojourns,
                       seq_id = str_c(seq_id,"_",seq_nb)))
      n_tp_per_state = df %>% dplyr::group_by(s, .drop = FALSE) %>% dplyr::summarize(n_tp = sum(sojourns),.groups = "drop")
    }
    if(!mult_seq_allowed) df$seq_id = seq_id
  }
  n_tp_per_state = df %>% dplyr::group_by(s, .drop = FALSE) %>% dplyr::summarize(n_tp = sum(sojourns),.groups = "drop")


  # 3. GENERATING OBSERVATIONS
  sim_seq = data.frame(seq_id = df$seq_id[rep(1:nrow(df),sojourns)],state = rep(s, sojourns))
  sim_seq = sim_seq %>% dplyr::group_by(seq_id) %>% dplyr::mutate(t = dplyr::row_number())
  # generate the observations
  X = purrr::map_dfr(unique(df$s),
                     function(s){
                       n = n_tp_per_state$n_tp[s]
                       x = generate_random_obs(n = n, state = s, parem = model$parms.emission)
                       m = generate_missingness(n = n, state = s, parem = model$parms.emission)
                       x[m == 1] = NA
                       j = which(sim_seq$state == s)
                       x = x %>% dplyr::mutate(seq_id = sim_seq$seq_id[j], t = sim_seq$t[j])
                       x
                     }
  )
  sim_seq = sim_seq %>% dplyr::full_join(.,X, by = c("seq_id","t")) %>% dplyr::ungroup() %>% dplyr::arrange(seq_id, t)

  # 4. TRUNCATING SEQUENCES IF NEEDED
  if((!is.null(n_timepoints)) & (!all_states)){sim_seq = sim_seq[1:n_timepoints,]}


  sim_seq
}



get_longest_sojourn = function(model){
  # check model
  # TODO

  # shortcut
  sojourn_distribution = model$sojourn$type
  p = 0.99

  # finding M
  if(sojourn_distribution == "nonparametric" | sojourn_distribution == "ksmoothed-nonparametric" ){
    M = nrow(model$sojourn$d)
  }else if(sojourn_distribution == "poisson"){
    M = max(model$sojourn$shift + qpois(p = p, lambda = model$sojourn$lambda))
  }else if(sojourn_distribution == "lnorm"){
    M = max(qlnorm(p = p, meanlog = model$sojourn$meanlog , sdlog = model$sojourn$s.dlog))
  }else if(sojourn_distribution == "gamma"){
    M = max(qgamma(p = p, shape = model$sojourn$shape , scale = model$sojourn$scale))
  }else if(sojourn_distribution == "logarithmic"){
    M = 10
    while(all(.dlog(M,model$sojourn$shape) > 0.01)){M = M*2}
  }else if(sojourn_distribution %in% c("nbinom")){
    M = max(qnbinom(p = p, mu = model$sojourn$mu ,prob = model$sojourn$prob))
  }else{
    stop("This sojourn distribution is currently not supported.")
  }

  M
}


#' @export
#' @useDynLib HiddenSemiMarkov sim_mc
sim.mc <- function(init,transition,N) {
  if(!all.equal(rowSums(transition),rep(1,nrow(transition)))) stop("Rows of transition matrix must sum to one")
  if(!all.equal(sum(init),1)) stop("Vector of initial probabilities must sum to one")
  a0 =  t(apply(transition,1,cumsum))
  st= cumsum(init)
  state = integer(sum(N))
  .C("sim_mc",as.double(st),as.double(a0),as.integer(nrow(transition)),state=state,as.integer(N),as.integer(length(N)),PACKAGE='HiddenSemiMarkov')$state
}




# generate random observations for all variables
generate_random_obs = function(n = 10, parem, state){
  r_obs = sapply(1:length(parem), function(i) generate_random_obs_par(n, parem[[i]],state))
  r_obs = matrix(r_obs, nrow = n) # in case we generate for only 1 observation
  colnames(r_obs) = names(parem)
  r_obs = r_obs %>% as.data.frame()
  r_obs = .check_types_and_values(data = r_obs, parem = parem)
  return(r_obs)
}


# generate random observations for one variable
generate_random_obs_par = function(n, par, state){
  if(par$type == "norm"){
    x = rnorm(n, mean = par$param$mean[state], sd = par$param$sd[state])
  }else if(par$type == "binom"){
    x = rbinom(n, size = par$param$size[state], prob = par$param$prob[state])
  }else if(par$type == "non-par"){
    x = sample(rep(par$param$values, round(par$param$probs[,state]*5000)) , n, replace = TRUE)
  }else{stop("This distributions hasn't been implemented yet")}
  return(x)
}


# impute
impute_missing_data_at_random_given_state = function(model, X, state, seed = 1){
  set.seed(seed)
  n_missing_obs =  X %>% select(all_of(names(model$parms.emission))) %>% is.na() %>% colSums()
  if(any(n_missing_obs > 0)){
    for(var in names(model$parms.emission)[n_missing_obs > 0]){
      j = which(is.na(X[, var]))
      X[j, var] = generate_random_obs_par(n = n_missing_obs[var],state = state, par = model$parms.emission[[var]])
    }
  }
  set.seed(Sys.time())
  X
}

impute_missing_data_with_most_likely_given_state = function(model, X, state){
  n_missing_obs =  X %>% select(all_of(names(model$parms.emission))) %>% is.na() %>% colSums()
  if(any(n_missing_obs > 0)){
    for(var in names(model$parms.emission)[n_missing_obs > 0]){
      j = which(is.na(X[, var]))
      val = most_probable_value(par = model$parms.emission[[var]], state = state)
      X[j, var] = rep(val, n_missing_obs[var])
    }
  }
  X
}




# generate_missingness

generate_missingness = function(n = 10, parem, state){
  r_obs = sapply(1:length(parem), function(v) rbinom(n = n, size = 1, prob = parem[[v]]$missing_prob[state]))
  r_obs = matrix(r_obs, nrow = n)
  colnames(r_obs) = names(parem)
  r_obs = r_obs %>% as.data.frame()
  return(r_obs)
}


