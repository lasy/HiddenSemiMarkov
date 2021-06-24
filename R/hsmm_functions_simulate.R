


#' Simulate a state and observation sequence from a hidden semi-Markov model.
#'
#' This function returns a \code{data.frame} with a simulated state sequence and observations.
#' The length of the simulated sequence can be specified either as a number of state transition (with argument \code{n_state_transitions}) or as a number of time-point (with argument \code{n_timepoints})
#' Each row of the returned \code{data.frame} is a time-point and the columns specify the sequence id, the time-point, the state and the values of the observation. Each variable has its own column.
#'
#' @param model a object of class or \code{hsmm}.
#' @param seq_id (optional) a \code{character} specifying the name of the sequence to simulate.
#' @param n_state_transitions (optional) an \code{integer} specifying the number of state transitions that should be simulated.
#' @param n_timepoints (optional) an \code{integer} specifying the number of time-point to simulate.
#'    One can either specify \code{n_state_transitions} or \code{n_timepoints}.
#'    If both are specified, \code{n_timepoints} is ignored and a sequence with \code{n_state_transitions} is specified.
#' @param seed (optional) the seed for the random generator. Allows to reproduce a simulated sequence.
#' @param all_states (optional) a \code{logical} specifying if all states should be present in the simulated sequence.
#'    By default, this option is \code{FALSE}. When \code{TRUE}, the sequence may be longer than the specified length via \code{n_state_transitions} or \code{n_timepoints}
#' @param min_tp_per_state (optional) if \code{all_states == TRUE}, this option can be used to specify the minimum number of time-point in each state.
#' @param mult_seq_allowed (optional) if \code{all_states == TRUE}, this option can be used to specify if representation of all states in the simulated sequence can be achieved by simulated several sequences or by simulating a single sequence.
#' @param verbose (optional) whether to print additional info when running.
#' @keywords HSMM
#' @return A \code{data.frame} with the following columns: \code{seq_id}, the identifier of the sequence, \code{t} the time-stamp, \code{state} the hidden state and one additional column for each specified variable.
#'
#' @references J. O’Connell, S. Hojsgaard, Hidden semi-Markov models for multiple observation sequences: The mhsmm package for R. Journal of Statistical Software. 39, 1–22 (2011) \url{https://www.jstatsoft.org/article/view/v039i04}
#'
#' @export
#' @importFrom magrittr %>%
#' @examples
#' my_model = simple_model # simple_model is a model attached to the HiddenSemiMarkov package for demos
#' Xsim = simulate_hsmm(model = my_model, n_state_transitions = 10)
#' plot_hsmm_seq(model = my_model, X = Xsim)


simulate_hsmm = function(model,
                         seq_id = NULL,
                         n_state_transitions = NULL, n_timepoints = NULL,
                         seed = NULL,
                         all_states = FALSE,  min_tp_per_state = NULL, mult_seq_allowed = TRUE,
                         verbose = FALSE){

  ##### What this function does:
  # 0. checks the input variables
  # 1. simulates a sequence of state transitions from the initial and transition probabilities.
  # 2. samples from the sojourn distributions to generate the state sequence.
  # 3. samples from the emission parameter distributions to generate the observations.
  # 4. truncates the simulated sequence to n_timepoints if it was specified

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
  model = .check_model(model = model)
  # sojourn matrix
  d = .build_d_from_sojourn_dist(model, M = .get_longest_sojourn(model)) # FIX THIS # ifelse(!is.null(n_timepoints),n_timepoints,.get_longest_sojourn(model))
  # seq_id
  if(is.null(seq_id) | (typeof(seq_id) != "character")) seq_id = stringr::str_c("sim_seq ",Sys.time())




  # 1 and 2. SIMULATING A SEQUENCE OF STATES and SAMPLING THE SOJOURN DISTRIBUTIONS
  if(verbose) cat("simulating a sequence of state\n")
  s = .sim.mc(init = model$init, transition = model$transition, N = n_state_transitions)

  if(verbose) cat("sampling from the sojourn distributions\n")
  M = nrow(d)
  sojourns = sapply(s, function(x) sample(1:M, size = 1, prob = d[,x], replace = TRUE))
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
      add_s = .sim.mc(init = init, transition = model$transition, N = n_state_transitions)
      add_sojourns = sapply(add_s, function(x) sample(1:nrow(d), size = 1, prob = d[,x], replace = TRUE))
      ns = nrow(df)
      s = c(s[-ns], add_s) # we remove the last element of s in case !mult_seq_allowed
      sojourns = c(sojourns[-ns], add_sojourns)
      df = dplyr::bind_rows(df[-ns,],
                     data.frame(
                       s = add_s %>% factor(.,levels = 1:model$J),
                       sojourns = add_sojourns,
                       seq_id = stringr::str_c(seq_id,"_",seq_nb)))
      n_tp_per_state = df %>% dplyr::group_by(s, .drop = FALSE) %>% dplyr::summarize(n_tp = sum(sojourns),.groups = "drop")
    }
    if(!mult_seq_allowed) df$seq_id = seq_id
  }

  if(verbose) cat("generating observations\n")

  # 3. GENERATING OBSERVATIONS

  sim_seq = data.frame(seq_id = df$seq_id[rep(1:nrow(df),sojourns)],state = rep(s, sojourns))
  sim_seq = sim_seq %>% dplyr::group_by(seq_id) %>% dplyr::mutate(t = dplyr::row_number())

  # generate the observations
  X = purrr::map_dfr(unique(df$s),
                     function(s){
                       j = which(sim_seq$state == s)
                       n = length(j)
                       x = .generate_random_obs(n = n, s = s, model = model)
                       x$seq_id = sim_seq$seq_id[j]
                       x$t = sim_seq$t[j]
                       x$state = s
                       x
                     }
  )
  X = .bin_to_continuous(X, model)
  sim_seq = X[,c("seq_id", "t","state", names(model$marg_em_probs))] %>% dplyr::arrange(seq_id, t)

  # 4. TRUNCATING SEQUENCES IF NEEDED
  if((!is.null(n_timepoints)) & (!all_states)){sim_seq = sim_seq[1:n_timepoints,]}

  sim_seq
}



#' @useDynLib HiddenSemiMarkov sim_mc
.sim.mc <- function(init,transition,N) {
  if(!all.equal(as.vector(rowSums(transition)) , rep(1, nrow(transition)))) stop("Rows of transition matrix must sum to one")
  if(!all.equal(sum(init),1)) stop("Vector of initial probabilities must sum to one")
  a0 =  t(apply(transition,1,cumsum))
  st= cumsum(init)
  state = integer(sum(N))
  .C("sim_mc",as.double(st),as.double(a0),as.integer(nrow(transition)),state=state,as.integer(N),as.integer(length(N)),PACKAGE='HiddenSemiMarkov')$state
}

.generate_random_obs = function(n, s, model){
  js = which(model$obs_prob$state == s)
  this_state_obs_probs = model$obs_prob[js,]
  j = sample(1:nrow(this_state_obs_probs),
             replace = TRUE, size = n,
             prob = this_state_obs_probs$p)
  random_obs = this_state_obs_probs[j, names(model$marg_em_probs)]
  random_obs
}

