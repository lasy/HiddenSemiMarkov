

####### SPECIFICATION AND SUPPORT FUNCTIONS ---------------------

#' Specifies a hidden semi-Markov model
#'
#' This function returns an object of class \code{hsmm} from the parameters passed as arguments.
#' @param J an integer. The number of hidden states.
#' @param init a vector of double. The initial probabilities associated with each states. Must be a vector of length \code{J}. If values do not sum to 1 a warning message is displayed and the vector is divided by its sum.
#' @param transition a matrix of double. The transition probabilities.
#'     Must be a matrix of size \code{J x J} in which the element \code{(r, c)} provides the probability of transitioning from state \code{r} to state \code{c}.
#'     The diagonal elements must be equal to zero except for absorbing states. For absorbing states, all the other elements of that row must be zero. For other states, if there are non-zero diagonal elements, a warning message is displayed and the elements are set to zero.
#'     If the sum over the rows is different that one, a warning message is displayed and the rows are divided by their sum.
#' @param sojourn a list. The sojourn distributions. The list must have at least two elements: a \code{type} element which specifies the type of sojourn distributions and the other elements are the distribution parameters.
#' Use the function \code{available_sojourn_dist()} to see the supported sojourn distributions.
#' @param marg_em_probs a list. The marginal emission distributions. The list has one element per variable.
#' The name of that element must be the name of the variable. Variable names cannot contain the character '_'.
#' Each element of the list is itself a list of at least two elements.
#' The first one, named \code{type}, is used to specify the distribution family. Type \code{available_marginal_emission_dist()} to get the currently supported distribution families.
#' The second element of the list, \code{params}, is a list that provides the parameters of the model in each state.
#' The same function (\code{available_marginal_emission_dist()}) provides information on how these parameters must be specified.
#' An optional third element can be added to the variable list: \code{viz-options}.
#' This element is a list in which each element specifies a given visualization option.
#' Type \code{marginal_emission_viz_options()} to obtain the list and description of the available visualization options for each distribution type.
#' See function \code{plot_hsmm_seq()} for visualization of observation sequences.
#' @param global_censoring_prob (optional) the probabilities of observations being censored in each state. Can be specified as a vector of length J of values between 0 (never censored) and 1 (always censored) or a a single value in [0,1] if the censoring probability is assumed to be identical in each state. If unspecified, the observations are assumed to never be censored (value 0) overall (individual variables may still be censored via their individual 'missing_prob'.)
#' @param state_names (optional) a vector of characters. Names associated to each state. Must be of length \code{J}. If unspecified, the states are numbered from 1 to \code{J}
#' @param state_colors (optional) a vector of color-specifying characters. Colors associated to each state. Must be of length \code{J}. If unspecified, the colors are picked from the \code{viridis} palette.
#' @param verbose a logical (default = \code{FALSE}). Should the function print additional information?
#'
#' @keywords HSMM
#' @return An object of class \code{hsmm} which can be used to simulate time series with the \code{simulate_hsmm()} function or to decode time-series with the \code{predict_hsmm_states()} function. The returned \code{hsmm} object (model) can also be fit to specific sequences with the \code{fit_hsmm()} function.
#'
#' @seealso \code{simulate_hsmm()} to simulate a sequence of hidden states and observations following the model specifications, \code{predict_states_hsmm()} to predict the sequence of hidden states from observation time-series and \code{fit_hsmm()} to fit a model to sequences of observations.
#'
#' @export
#' @importFrom magrittr %>%
#' @examples
#' my_model = specify_hsmm(
#'    J = 2,
#'    init = c(1,0),
#'    transition = matrix(c(0,1,1,0),2,2),
#'    sojourn = list(type = "gamma", shape = c(2,10), scale = c(10,3)),
#'    marg_em_probs = list(
#'        var1 = list(
#'            type = "norm",
#'            params = list(
#'                mean = c(0,1),
#'                sd = c(0.3,0.2)
#'                )
#'            ),
#'      var2 = list(
#'          type = "binom",
#'          params = list(
#'              size = rep(1,2),
#'              prob = c(0.2,0.8)
#'              )
#'          ),
#'      var3 = list(
#'          type = "non-par",
#'          params = list(
#'              values = c("a","b","c","d"),
#'              probs = matrix(c(0.7,0.1,0.1,0.1,
#'                             1/4,1/4,1/4,1/4), 4,2)
#'            ),
#'          viz_options = list(colors = c("black","slateblue1","#E90046","#0EC290"))
#'          )
#'      ),
#'    censoring_probs = list(p = c(0.1,0.2), q = matrix(c(0.1,0.2,0.3,0.4,0.5,0.6), nrow = 3, ncol = 2)),
#'    state_names = c("A","B"),
#'    state_colors = c("seagreen1","slategray")
#' )
#' class(my_model)
#' Xsim = simulate_hsmm(model = my_model, n_state_transitions = 20)
#' plot_hsmm_seq(model = my_model, X = Xsim)
#'
#'
specify_hsmm = function(J,
                        state_names = NULL, state_colors = NULL,
                        init, transition,
                        sojourn,
                        marg_em_probs,
                        censoring_probs = NULL,
                        verbose = FALSE){

  # 1 . checks
  if(verbose) cat("Checking inputs\n")

  J = .check_J(J)
  state_names = .check_state_names(state_names, J)
  state_colors = .check_state_colors(state_colors, J)
  init = .check_init(init, J)
  transition = .check_transitions(transition, J)
  sojourn = .check_sojourn(sojourn, J, state_names)
  marg_em_probs = .check_marg_em_probs(marg_em_probs, J)
  censoring_probs = .check_censoring_probs(censoring_probs, J, length(marg_em_probs))

  model = list(J = J,
               state_names = state_names,
               state_colors = state_colors,
               init = init,
               transition = transition,
               sojourn = sojourn,
               marg_em_probs = marg_em_probs,
               censoring_probs = censoring_probs
               )

  # 2. Initialize model emission_probabilities

  if(verbose) cat("Initializing obs_probs \n")
  model$obs_probs = .initialize_obs_probs(model)

  if(verbose) cat("Format b \n")
  model$b = model$obs_probs %>%
    tidyr::pivot_wider(names_from = state, values_from = p, names_prefix = "p_")

  if(!model$censoring_probs$missing_prob_specified){
    if(verbose) cat("khdiwhdo \n")
    j = apply(model$obs_probs, 1, function(x) any(is.na(x))) %>% which()
    model$obs_probs$p[j] = 0
    model$obs_probs = model$obs_probs %>%
      dplyr::group_by(state) %>%
      dplyr::mutate(sum_prob = sum(p),
                    p = p/sum_prob) %>%
      dplyr::select(-sum_prob) %>%
      dplyr::ungroup()
  }

  model$obs_probs = model$obs_probs %>% dplyr::ungroup()
  model$obs_probs$p0 = model$obs_probs$p

  # 3. Returns specified model
  class(model) <- 'hsmm'
  model
}


.initialize_obs_probs = function(model){

  obs_prob = data.frame(state = 1:model$J)
  # 1. marginal probabilities with missing probabilities for each variable
  for(var in names(model$marg_em_probs)){
    this_var_obs_prob = .get_obs_prob(model, var)
    obs_prob = obs_prob %>%
      full_join(., this_var_obs_prob, by =  "state")
  }
  obs_prob$p = apply(obs_prob %>%
                       dplyr::select(dplyr::starts_with("probability_of_var")) %>%
                       as.data.frame(),
                     1, prod)

  # 2. probability that all variables are missing
  all_missing_prob = data.frame(state = 1:model$J,
                                all_missing_p = model$censoring_probs$p)

  obs_prob = obs_prob %>%
    dplyr::left_join(.,
              all_missing_prob, by = "state") %>%
    dplyr::mutate(p = (1-all_missing_p)*p) %>%
    dplyr::select(-all_missing_p)

  for(var in names(model$marg_em_probs)) all_missing_prob[,var] = NA
  obs_prob =  obs_prob %>%
    dplyr::left_join(.,
              all_missing_prob,
              by = c("state", names(model$marg_em_probs))) %>%
    dplyr::mutate(all_missing_p = all_missing_p %>% tidyr::replace_na(0),
                  p = p + all_missing_p) %>%
    dplyr::select(-all_missing_p)

  # 3. normalizing probabilities (sum to 1 for each state)
  obs_prob = obs_prob %>%
    dplyr::group_by(state) %>%
    dplyr::mutate(sum_prob = sum(p),
                  p = p/sum_prob) %>%
    dplyr::select(-sum_prob)

  # 4. formatting
  obs_prob = obs_prob %>%
    dplyr::select(state, dplyr::all_of(names(model$marg_em_probs)), p)

  obs_prob
}

.get_obs_prob = function(model, var){

  var_i = which(names(model$marg_em_probs) == var)

  this_var_obs_probs = .get_marginal_prob(var_name = var, model = model)
  this_var_obs_probs[,var] = this_var_obs_probs$x
  this_var_obs_probs = this_var_obs_probs %>%
    dplyr::select(-x)

  this_var_missing_prob =
    data.frame(state = 1:model$J,
               q = model$censoring_probs$q[var_i,]
    )

  this_var_obs_probs =
    dplyr::left_join(
      this_var_obs_probs,
      this_var_missing_prob %>%
        dplyr::mutate(p_m = 1-q) %>%
        dplyr::select(-q),
      by = "state") %>%
    dplyr::mutate(p = prob * p_m) %>%
    dplyr::select(-prob, -p_m)

  this_var_missing_prob[, var] = NA

  this_var_obs_probs =
    dplyr::bind_rows(
      this_var_obs_probs,
      this_var_missing_prob %>%
        dplyr::rename(p = q)
    )

  this_var_obs_probs[,paste0("probability_of_var_",var)] = this_var_obs_probs$p
  this_var_obs_probs = this_var_obs_probs %>%
    dplyr::select(-p)
  this_var_obs_probs
}


# Computes matrix b for a given sequence of observations.
.compute_obs_probs = function(model, X){
  Xp = .cut_continuous_var(model = model, X = X) # cut continuous variables so that they match values in matrix b
  Xb = dplyr::left_join(Xp, model$b, by = intersect(colnames(model$b), colnames(Xp))) # join X with model$b
  p = Xb %>% dplyr::select(dplyr::matches("p_")) %>% as.matrix() # select the prob columns and turn them into a matrix
  p
}



####### DECODING (predicting state sequence) ---------------------

#' Predicts the hidden state sequence from observations
#'
#' This function is a wrapper around the function \code{predict_states_hsmm()}.
#' It predicts the most likely hidden states from observations.
#' Two methods are implemented:
#' \code{"Viterbi"} applies the Viterbi algorithm and predicts the most likely sequence of hidden states,
#' and \code{"FwBw"} applies the Forward-Backward algorithm and returns the probability of each state at each time-point.
#'
#' @param object an \code{hsmm} model specified via the \code{specify_hsmm()} function.
#' @param newdata a \code{data.frame} with the observation sequences.
#' @param method a \code{character} specifying the method to be used, i.e. either \code{"Viterbi"} or \code{"Fwbw"}.
#'
#' @seealso see \code{predict_states_hsmm()} for the full description, options and examples.
#'
predict.hsmm <- function(object, newdata, method = "Viterbi", verbose = FALSE, ...) {
  ans = predict_states_hsmm(model = object, X = newdata, method = method, verbose = verbose, ...)
  ans
}


#' Predicts the hidden state sequence from observations
#'
#' This function predicts the most likely hidden states from observations.
#' Two methods are implemented:
#' \code{"Viterbi"} applies the Viterbi algorithm and predicts the most likely sequence of hidden states,
#' and \code{"FwBw"} applies the Forward-Backward algorithm and returns the probability of each state at each time-point.
#'
#' @param model a \code{hsmm} object. The model used to predict the hidden sequence of states.
#' @param X a \code{data.frame}. The observation sequences.
#' This \code{data.frame} must have the following columns:
#' \code{seq_id} (\code{character}) provides the sequence id, which allows the prediction of hidden states for several sequence simultaneously,
#' \code{t} (\code{numeric}) provides the timestamp. Time-points must be separated by the same time-step,
#' \code{...} one column for each variable specified for the model.
#' Additional columns will be ignored.
#' @param method a \code{character} specifying the decoding algorithm.
#'    \code{"Viterbi"} returns the most likely sequence of hidden states.
#'    \code{"FwBw"} applies the "Forward-Backward" algorithm as described by Guédon, 2003, and returns the "smoothed" (posterior) probability of each state at each time-point.
#' @param ground_truth (optional) a \code{data.frame} with three columns (\code{seq_id, t, state}) providing the ground-truth, i.e. the actual hidden state, for a given (set of) sequence(s) and time-points.
#' @param trust_in_ground_truth (optional) a double in [0,1] that indicates the reliability of the provided ground-truth. 1 means "full trust", 0 means "no trust". Default value is 0.75.
#' @param verbose logical. Should the function prints additional information?
#'
#' @export
#' @examples
#'
#'
#' my_model = simple_model
#' Xsim = simulate_hsmm(my_model, n_state_transitions = 20)
#'
#' # viterbi decoding
#' viterbi = predict_states_hsmm(model = my_model, X = Xsim, method = "Viterbi")
#' Xsim$state_viterbi = viterbi$state_seq$state
#' plot_hsmm_seq(X = Xsim, model = my_model)
#'
#' # forward backward decoding
#' smoothed = predict_states_hsmm(model = my_model, X = Xsim, method = "FwBw")
#' Xsim$state_smoothed = smoothed$state_seq$state
#' plot_hsmm_seq(X = Xsim %>% dplyr::select(-state_viterbi), model = my_model)
#'
#' ggplot2::ggplot(
#'    smoothed$probabilities,
#'    aes(x = t, y = state_prob, col = factor(state))
#'    ) +
#' ggplot2::geom_line() +
#' ggplot2::scale_color_manual(values = my_model$state_colors)
#'
predict_states_hsmm = function(model, X,
                               method = "Viterbi",
                               ground_truth = data.frame(),
                               trust_in_ground_truth = 0.75,
                               verbose = FALSE) {

  # CHECKS
  # model
  model = .check_model(model)
  if(verbose) cat("Model checked\n")
  # check data
  X = .check_data(data = X, model = model)
  if(verbose) cat("Data checked\n")
  # ground_truth
  if((nrow(ground_truth)>0) && (!all(is.na(ground_truth$state))))
    ground_truth = .check_ground_truth(ground_truth, model, X)
  if(verbose) cat("Ground truth checked\n")

  # decoding
  if(method=="Viterbi") {
    ans = .predict_states_hsmm_viterbi(model = model, X = X,
                                       ground_truth = ground_truth, trust_in_ground_truth = trust_in_ground_truth,
                                       verbose = verbose)
  }else if(method == "FwBw"){
    ans = .predict_states_hsmm_forward_backward_seq_by_seq(model = model, X = X,
                                                           ground_truth = ground_truth, trust_in_ground_truth = trust_in_ground_truth,
                                                           verbose = verbose)
  }else{stop(paste("Unavailable prediction method",method))}

  ans
}


#' @useDynLib HiddenSemiMarkov viterbi
.predict_states_hsmm_viterbi =
  function(model, X = X,
           ground_truth = data.frame(), trust_in_ground_truth = 0.75,
           verbose = FALSE){

  # Number of states
  J = model$J
  # Sequences lengths
  N = rle(X$seq_id %>% as.character())$lengths

  # add variables N, NN, M, d, D, b + log transform
  augmented_model = .augment_model(model = model, X = X, ground_truth = ground_truth, log = TRUE, m = -1e300, trust_in_ground_truth = trust_in_ground_truth, verbose = verbose)
  NN = augmented_model$NN
  if(verbose){cat("Model augmented \n")}
  SeqIDs = unique(X$seq_id)

  # Initialization of the output variables
  state_seq = data.frame()
  loglik = data.frame()
  probabilities = data.frame()

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

    probabilities =
      dplyr::bind_rows(
        probabilities,
        augmented_model$b %>% as.data.frame() %>%
          dplyr::mutate(seq_id = X$seq_id,
                        t = X$t) %>%
          tidyr::pivot_longer(cols = dplyr::starts_with("p_"),
                              names_to = "state",
                              names_prefix = "p_",
                              values_to = "obs_prob") %>%
          dplyr::mutate(state = state %>% as.numeric)
      )

    state_seq =
      dplyr::bind_rows(
        state_seq,
        data.frame(seq_id = SeqIDs[i] %>% factor(.,levels = SeqIDs),
                   t = X$t[ix],
                   state = s,
                   #state_deprecated = statehat,
                   likelihood = alpha[cbind(1:N[i],s)])
      )

    loglik = dplyr::bind_rows(
      loglik,
      data.frame(seq_id = SeqIDs[i] %>% factor(.,levels = SeqIDs),
                 loglik = max(alpha[N[i],(1:J)]))
    )
  }
  if(verbose) cat("\n")

  ans <- list(state_seq = state_seq,
              loglik = loglik,
              probabilities = probabilities)
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
.predict_states_hsmm_forward_backward =
  function(model, X = X,
           ground_truth = ground_truth, trust_in_ground_truth = 0.75,
           verbose = FALSE){

  # add variables N, NN, M, d, D, b
  augmented_model = .augment_model(model = model, X = X, ground_truth = ground_truth, log = FALSE, m = -1e300, trust_in_ground_truth = trust_in_ground_truth)
  if(verbose) cat("Model augmented \n")

  # shortcuts
  J = augmented_model$J
  N = augmented_model$N
  NN = augmented_model$NN
  M = augmented_model$M
  Nseq = length(augmented_model$N)
  SeqIDs = unique(X$seq_id)
  nX = sum(N)

  # probabilities (obs_prob)
  probabilities = data.frame(
    seq_id = rep(X$seq_id, J),
    t = rep(X$t, J),
    state = rep(1:J, each = nrow(X)),
    obs_prob = augmented_model$b %>% as.vector(),
    stringsAsFactors = FALSE
  )

  # run the smoothed algorithm on all sequences at once
  if(verbose) cat(stringr::str_c("Decoding sequences (",length(N),")... \n"))

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
    stop("NaNs detected in posterior probabilities.") #
  }
  if(any(fwbw_res$gamma<0)) fwbw_res$gamma = zapsmall(fwbw_res$gamma)
  if(any(fwbw_res$eta<0)) fwbw_res$eta = zapsmall(fwbw_res$eta)
  if(any(fwbw_res$N<0))  fwbw_res$N = zapsmall(fwbw_res$N)


  if(verbose) cat("Formating results \n")
  ##### Format results
  # State probabilities
  probabilities = probabilities %>%
    dplyr::mutate(
      state_prob = fwbw_res$gamma)
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
             probabilities = probabilities,
             fwbw_res = fwbw_res)
  ans
}



.predict_states_hsmm_forward_backward_seq_by_seq =
  function(model, X = X,
           ground_truth = data.frame(), trust_in_ground_truth = 0.75,
           verbose = FALSE){

  # shortcuts
  SeqIDs = unique(X$seq_id)
  Nseq = length(SeqIDs)

  # initialize output variables
  state_seq = data.frame()
  probabilities = data.frame()
  loglik = data.frame()
  sequences_with_error = c()

  # run the smoothed algorithm sequence by sequence
  if(verbose) cat(stringr::str_c("Decoding sequences (x/",length(Nseq),"):\t"))
  for(i in 1:Nseq){
    if(verbose) cat(i,"\t")
    SeqID = SeqIDs[i]
    Xi = X %>% dplyr::filter(seq_id == SeqID) %>% dplyr::arrange(seq_id, t)
    if(nrow(ground_truth)>0) ground_truth_i = ground_truth %>% dplyr::filter(seq_id == SeqID) %>% dplyr::arrange(seq_id, t) else ground_truth_i = data.frame()

    Di = try(.predict_states_hsmm_forward_backward(model = model, X = Xi, ground_truth = ground_truth_i, trust_in_ground_truth = trust_in_ground_truth, verbose = FALSE))
    if(class(Di) == "try-error"){
      sequences_with_error = c(sequences_with_error, SeqID)
    }else{
      state_seq = dplyr::bind_rows(state_seq, Di$state_seq)
      probabilities = dplyr::bind_rows(probabilities, Di$probabilities)
      loglik = dplyr::bind_rows(loglik, Di$loglik)
    }
  }
  if(verbose) cat("\n done \n")

  # return results
  ans = list(state_seq = state_seq,
             loglik = loglik,
             probabilities = probabilities,
             sequences_with_error = sequences_with_error)
  ans
}




.augment_model = function(model, X, ground_truth,
                          trust_in_ground_truth = 0.75,
                          log = FALSE, m = -1e300,
                          verbose = FALSE){
  augmented_model = model

  # SEQUENCES start index
  # N
  augmented_model$N = rle(X$seq_id %>% as.character())$lengths
  # NN
  augmented_model$NN = cumsum(c(0,augmented_model$N))
  # M
  augmented_model$M = max(c(augmented_model$N, .get_longest_sojourn(model = model)))  #max(augmented_model$N);

  # SOJOURN
  # d
  augmented_model$d = .build_d_from_sojourn_dist(model = model, M = augmented_model$M)
  # D
  augmented_model$D = apply(augmented_model$d,2,function(x) rev(cumsum(rev(x))))

  if(verbose) cat("Sojourn matrix built\n")


  # OBSERVATION PROBABILITIES
  augmented_model$b = .compute_obs_probs(model = model, X = X)
  if(verbose) cat("joint emission probabilities computed\n")
  if(any(rowSums(augmented_model$b) == 0)) stop("The joint emission probabilities of some observations have a value of zero for all states (which means that some observations are impossible given the specified model). Please check your model specification and make sure that the observations at all time-point have a non-zero probability in at least one state.")
  if(any(rowSums(augmented_model$b) < 1e-8)) warning("Some observations are highly unlikely in all states. This may mean that the specified model is not able to capture all of the provided observations.")
  augmented_model$b = augmented_model$b / rowSums(augmented_model$b )


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

  }

  augmented_model
}





####### FITTING ---------------------


#' Fits a hidden semi-Markov model to data sequences
#'
#' This function relies on a EM procedure to fit the model parameters to maximize the likelihood of the decoded hidden state sequence.
#' It returns a list whose first element is the fitted model (an object of class \code{hsmm}) and whose second elements provides information about the EM procedure (convergence, number of iteration, likelihood).
#' @param model a \code{hsmm} object. The model whose parameters will be re-estimated.
#' @param X a \code{data.frame} of observations.
#' @param ground_truth (optional) a \code{data.frame} of ground truth, _i.e._ time-points where the hidden state is known. Default is an empty \code{data.frame}, _i.e._ no ground truth.
#' @param n_iter (optional) an integer specifying the maximal number of iterations for the EM-procedure. Default value is 10.
#' @param rel_tol (optional) a positive double specifying the tolerance at which to stop the EM. If the difference in likelihood (normalized by the total sequences length) between two iterations of the EM is smaller than \code{rel_tol}, then the EM procedure is considered to have converged to a local maximum.
#' @param lock_emission (optional) a logical. Default is \code{FALSE}. Specifies if the emission distributions should be locked (kept as is) or re-estimated at the M-step of the EM.
#' @param lock_transition (optional) a logical. Default is \code{FALSE}. Specifies if the transition probability should be locked (kept as is) or re-estimated at the M-step of the EM.
#' @param lock_sojourn (optional) a logical. Default is \code{FALSE}. Specifies if the sojourn distributions should be locked (kept as is) or re-estimated at the M-step of the EM.
#' @param N0 (optional) a positive number specifying the strength of the prior, i.e. the number of observations (or believed number of observations) which informed the specification of the emission distributions. This number will be used to weight the specified emission distribution against the total lenght of the sequences provided for the fit.
#' @param use_sojourn_prior (optional) a logical. Default is \code{TRUE}. Specifies if the specified sojourn distributions should be used as a prior when updating the prior distributions at the M-step of the EM.
#' @param trust_in_ground_truth (optional) a double between 0 and 1 specifying the trust in the ground truth. A value of 0 indicates no trust and is equivalent to not providing ground-truth. A value of 1 indicates full trust and the ground truth will not be modulated by the probability of the values of the observations.
#' @keywords HSMM
#' @return A list. First element of the list (\code{$model}) is a \code{hsmm} object (the fitted model) and the second element (\code{$fit_param}) provides information about the EM-procedure. The second element can be visualized by calling the function \code{plot_hsmm_fit_param()}.
#'
#' @export
#' @importFrom magrittr %>%
#' @examples
#' my_model = simple_model  # simple_model is a model attached to the HiddenSemiMarkov package for demos
#' Xsim = simulate_hsmm(my_model, n_state_transitions = 20) # a short sequence is simulated
#'
#' my_model_fit = fit_hsmm(model = my_model, X = Xsim) # the model is fit to the observations.
#' plot_hsmm_fit_status(fit_output = my_model_fit)
#'
#' viterbi_init = predict_states_hsmm(model = my_model, X = Xsim, method = "Viterbi") # predict the states with the initial model
#' viterbi_fit = predict_states_hsmm(model = my_model_fit$model, X = Xsim, method = "Viterbi") # predict the states with the fit model
#' Xsim$state_viterbi_init = viterbi_init$state_seq$state
#' Xsim$state_viterbi_fit = viterbi_fit$state_seq$state
#' plot_hsmm_seq(X = Xsim, model = my_model)

fit_hsmm = function(model, X,
                    n_iter = 10, rel_tol = 1/20,
                    lock_emission = FALSE,
                    lock_transition = FALSE,
                    lock_sojourn = FALSE,
                    N0 = 0,
                    use_sojourn_prior = FALSE,
                    ground_truth = data.frame(),
                    trust_in_ground_truth = 0.75,
                    verbose = FALSE, graphical = FALSE){

  # 1. Checks
  if(missing(model)) stop("model missing")
  model = .check_model(model = model)

  if(missing(X)) stop("X missing!")
  X = .check_data(data = X, model = model)
  if(verbose) cat("Data checked\n")

  # if any ground_truth, check it
  if(nrow(ground_truth)>0){
    ground_truth = .check_ground_truth(ground_truth, model = model, X = X)
  }else{
    ground_truth = X %>%  dplyr::select(seq_id, t) %>% dplyr::mutate(state = NA)
  }

  if((trust_in_ground_truth < 0) | (trust_in_ground_truth > 1)) stop("trust_in_ground_truth must be a number between 0 and 1.")

  if(use_sojourn_prior) model$d_prior = .build_d_from_sojourn_dist(model = model, M = max(table(X$seq_id)))

  # 2. EM
  ll = c(); message = "Reached n_iter" # initializing the vector which keeps track of the log likelihood
  for(it in 1:n_iter){
    if(verbose){cat("i: ",it,"\n")}

    ###### E-step ###
    if(verbose){cat("\t E-step \n")}
    smoothed_res = try(.predict_states_hsmm_forward_backward(model = model, X = X, ground_truth = ground_truth, trust_in_ground_truth = trust_in_ground_truth),
                       silent = TRUE)

    if(class(smoothed_res) == "try-error"){
      smoothed_res = try(predict_states_hsmm(model = model, X = X, method = "FwBw", ground_truth = ground_truth, trust_in_ground_truth = trust_in_ground_truth), silent = TRUE)
      prob_seqs = smoothed_res$sequences_with_error
      ground_truth_tmp = ground_truth
      for(prob_seq in prob_seqs){
        vit_seq = predict_states_hsmm(model = model,
                                      X = X %>% dplyr::filter(seq_id == prob_seq),
                                      ground_truth = ground_truth %>% dplyr::filter(seq_id == prob_seq),
                                      trust_in_ground_truth = trust_in_ground_truth ,
                                      method = "Viterbi")
        ground_truth_tmp = ground_truth_tmp %>%
          dplyr::full_join(., vit_seq$state_seq %>% dplyr::select(seq_id, t, state) %>% dplyr::rename(state_vit = state), by = c("seq_id","t")) %>%
          dplyr::mutate(state = ifelse(!is.na(state_vit),state_vit, state)) %>% dplyr::select(-state_vit)
      }

      smoothed_res = try(.predict_states_hsmm_forward_backward(model = model, X = X, ground_truth = ground_truth_tmp, trust_in_ground_truth = trust_in_ground_truth))

      if(class(smoothed_res) == "try-error"){
        message = paste0("Error in the E-step. Forward-Backward decoding threw an error at iteration ",it,".") #TODO Model from previous iteration is returned.
        break()
      }
    }

    p = smoothed_res$probabilities %>% dplyr::select(-obs_prob) %>%
      tidyr::pivot_wider(id_cols = c("seq_id","t"), names_from = "state", values_from = "state_prob") %>%
      dplyr::select(-seq_id, -t) %>% as.matrix(., ncol = model$J, nrow = nrow(X))
    p = p/rowSums(p)
    if(graphical) matplot(p, type = "l", lty = 1, col = model$state_colors)


    ###### M-step ###
    if(verbose){cat("\t M-step \n")}
    weights = .combine_smoothed_probs_with_ground_truth(p = p, ground_truth = ground_truth, trust_in_ground_truth = trust_in_ground_truth)
    # we check for un-visited states or rarely visited states
    state_seq = apply(weights,1,which.max)
    if(verbose){
      if(any(table(state_seq)<(nrow(X)/model$J/10))) warning("Some states are rarely visited\n")
      if(any(table(state_seq)==0)) warning("Some states are never visited\n")
    }

    new_model = .re_estimate_parameters(model = model, X = X, w = weights,
                                        fwbw_res = smoothed_res$fwbw_res,
                                        lock_transition = lock_transition,
                                        lock_sojourn = lock_sojourn,
                                        lock_emission = lock_emission,
                                        use_sojourn_prior = use_sojourn_prior,
                                        N0 = N0,
                                        verbose = verbose)
    model = new_model
    # Keep iterating?
    ll[it]= sum(smoothed_res$loglik$loglik)/nrow(X) # sum(log(estep_variables$N))
    if(verbose) cat("\t ll: ",ll[it],"\n")
    if(it>=2) if( abs(ll[it]-ll[it-1]) < rel_tol*(abs(ll[it])+ abs(ll[it-1]))/2 ){message = "Converged"; break()}
  }
  if(verbose) cat("Model fitted\n")
  # TODO: update $censored_obs_probs, $marg_em_probs and $censoring_probs
  # model$censored_obs_probs = .re_estimate_censored_obs_prob(model = model, X = X, w = weights, N0 = N0)
  # model$marg_em_probs = .re_estimate_marginal_emission_probabilities(model = model, X = X, w = weights)
  # model$censoring_probs = .re_estimate_censoring_probabilities(model = model, X = X, w = weights)

  out = list(model = model, fit_param = list(ll = ll, message = message, n_iter = length(ll)))
  out
}

.combine_smoothed_probs_with_ground_truth = function(p , ground_truth , trust_in_ground_truth){
  # initialization of the weights matrix
  weights = p

  weights_ground_truth = weights
  j = which(!is.na(ground_truth$state))
  if(length(j)>0){
    weights_ground_truth[j,] = 0
    weights_ground_truth[cbind(j,ground_truth$state[j])] = 1
  }

  weights = (1-trust_in_ground_truth) * weights + trust_in_ground_truth * weights_ground_truth

  weights
}


.re_estimate_parameters = function(model, X, w,
                                   fwbw_res,
                                   lock_transition = FALSE,
                                   lock_sojourn = FALSE,
                                   lock_emission = FALSE,
                                   use_sojourn_prior = TRUE,
                                   N0 = 200,
                                   verbose = FALSE){
  new_model = model
  if(!lock_emission){
    if(verbose) cat("updating observation probabilities\n")
    new_model = .re_estimate_obs_probs(model = new_model, X = X, w = w, N0 = N0, verbose = verbose)
  }
  if(!lock_transition){
    if(verbose) cat("updating transition matrix\n")
    new_model = .re_estimate_transition_probabilities(model = new_model, fwbw_res = fwbw_res)
  }
  if(!lock_sojourn){
    if(verbose) cat("updating sojourn distributions\n")
    new_model = .re_estimate_sojourn_distributions(model = new_model, fwbw_res = fwbw_res, use_sojourn_prior = use_sojourn_prior, graphical = verbose)
  }
  return(new_model)
}



.re_estimate_obs_probs = function(model, X, w = w, N0 = 200, verbose = FALSE){

  new_model = model

  var_names = names(model$marg_em_probs)
  Xb = .cut_continuous_var(model = model, X = X)
  Xb_with_w = Xb %>%
    dplyr::select(seq_id, t, all_of(var_names)) %>% # we only keep the seq_id, the time-points and the observations
    dplyr::full_join(., # we join with the weight, but we need to transform them to long format first
                     w %>% as.data.frame() %>%
                       magrittr::set_colnames(1:model$J) %>%
                       dplyr::mutate(seq_id = X$seq_id, t = X$t) %>%
                       tidyr::pivot_longer(col = c(-seq_id, -t), names_to = "state", values_to = "p") %>%
                       dplyr::mutate(state = as.integer(state)),
                     by = c("seq_id","t"))

  observed_obs_probs = Xb_with_w %>%
    dplyr::group_by(.dots = c("state", var_names)) %>%
    dplyr::summarize(counts = sum(p), .groups = "drop")

  observed_N = observed_obs_probs %>%
    dplyr::select(state, counts) %>%
    dplyr::group_by(state) %>%
    dplyr::summarize(N = sum(counts), .groups = "drop")

  obs_probs = model$obs_probs %>%
    dplyr::left_join(., observed_N, by = "state") %>%
    dplyr::left_join(.,
                     observed_obs_probs,
                     by = c("state", var_names)) %>%
    dplyr::mutate(counts = counts %>% tidyr::replace_na(0),
                  p = (p0 * N0 + counts)/(N0 + N)) %>%
    dplyr::select(state, dplyr::all_of(var_names), p, p0)

  new_model$obs_probs = obs_probs

  new_model$b =  new_model$obs_probs %>% dplyr::select(-p0) %>%
    tidyr::pivot_wider(names_from = state, values_from = p, names_prefix = "p_")

  new_model
}

#
# # deprecated
# .re_estimate_obs_probs = function(model, X, w = w, N0 = 200, verbose = FALSE){
#
#   new_model = model
#   var_names = names(model$marg_em_dist)
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
#   observed_missing = Xb_with_w %>% dplyr::mutate(dplyr::across(tidyselect::all_of(var_names), is.na))
#
#   # we first make a list with all combination of reported variables
#   all_obs_combs = .get_all_possible_combination_of_reported_variables(model)
#   P_values = purrr::map_dfr(.x = all_obs_combs,
#                             .f = function(obs_comb){
#                               obs_names = stringr::str_split(obs_comb, "_") %>% unlist()
#                               contingency_table = Xb_with_w %>%
#                                 dplyr::select(tidyselect::all_of(obs_names), state, p) %>%
#                                 dplyr::mutate(has_NA = observed_missing %>% dplyr::select(tidyselect::all_of(obs_names)) %>% apply(.,1,any)) %>%
#                                 dplyr::filter(!has_NA) %>% dplyr::select(-has_NA) %>%
#                                 dplyr::group_by(.dots = c("state", obs_names)) %>%
#                                 dplyr::summarise(counts = sum(p), .groups = "drop") %>%
#                                 dplyr::group_by(state) %>%
#                                 dplyr::mutate(N = sum(counts)) %>%
#                                 dplyr::ungroup()
#                               contingency_table[,setdiff(var_names, obs_names)] = NA
#                               contingency_table = contingency_table %>% dplyr::select(state, dplyr::all_of(var_names), counts, N)
#                               contingency_table
#                             })
#
#
#
#   new_model$obs_probs = model$obs_probs %>%
#     dplyr::left_join(., P_values, by = c("state", var_names)) %>%
#     dplyr::mutate(counts = counts %>% replace_na(0),
#                   N = N %>% replace_na(0),
#                   p = (p0*N0 + counts)/(N0 + N),
#                   p = p %>% replace_na(1)) %>%
#     dplyr::select(-counts, - N)
#
#   new_model$b =  new_model$obs_probs %>% dplyr::select(-p0) %>%
#     tidyr::pivot_wider(names_from = state, values_from = p, names_prefix = "p_")
#
#   new_model
# }



.get_all_possible_combination_of_reported_variables = function(model){

  all_levels = .get_all_possible_levels(model = model, with_missing = TRUE, continuous_var_binned = TRUE)
  all_levels_no_missing = .get_all_possible_levels(model = model, with_missing = FALSE, continuous_var_binned = TRUE)

  var_names = colnames(all_levels)
  vars_with_missing_values = var_names[apply(all_levels, 2, function(x) any(is.na(x)))]
  vars_never_missing = var_names %>% setdiff(vars_with_missing_values)


  start_l = ifelse(length(vars_never_missing) == 0, 1, 0)
  all_var_combs = c()
  for(l in start_l:length(vars_with_missing_values)){ # the size of the set of observable variables is length(var_never_missing) + l
    var_combinations = combn(x = vars_with_missing_values, m = l)
    var_combinations = rbind(var_combinations,
                             matrix(vars_never_missing, ncol = ncol(var_combinations), nrow = length(vars_never_missing), byrow = FALSE))
    for(k in 1:ncol(var_combinations)){
      this_var_combination = stringr::str_c(var_names[var_names %in% var_combinations[,k]], collapse = "_")
      all_var_combs = c(all_var_combs, this_var_combination)
    }
  }
  all_var_combs
}






.re_estimate_transition_probabilities = function(model, fwbw_res){

  new_model = model
  new_model$init = pmax(0, fwbw_res$init)
  new_model$transition = pmax(0, fwbw_res$transition) %>% matrix(.,ncol = model$J)

  new_model
}


.re_estimate_sojourn_distributions = function(model, fwbw_res, use_sojourn_prior = TRUE, graphical = FALSE){

  new_model = model

  # the smoothed (forward/backward) algorithm returns a new "d" matrix
  ds = fwbw_res$eta %>% matrix(.,ncol = model$J)
  # if we want to use Bayesian updating, we need to multiply this new d matrix by the prior.
  if(use_sojourn_prior){
    d_prior = model$d_prior
    if(nrow(d_prior) < nrow(ds)) d_prior = rbind(d_prior, matrix(0, nrow = nrow(ds) - nrow(d_prior), ncol = model$J))
    if(nrow(d_prior) > nrow(ds)) d_prior = d_prior[1:nrow(ds),]
    ds = ds * d_prior
  }
  # This matrix now needs to be normalized such that the sojourn density distribution sums to 1 for any state.
  ds = ds %>% apply(.,2,function(x) x/sum(x))
  ds_i = ds # we keep the "initial" ds
  M = nrow(ds) # max sojourn duration

  sojourn = purrr::map(
    .x = 1:model$J,
    .f = .re_estimate_sojourn_distribution,
    model = model,
    d = ds
  )

  new_model$sojourn = sojourn

  if(graphical){
    ds_o = .build_d_from_sojourn_dist(model, M = M)
    ds_updated = .build_d_from_sojourn_dist(new_model, M = M)
    par(mfrow = c(5,4), mar = c(0,0,0,0)+1.8)
    for(i in 1:model$J){
      max_y = max(c(ds_o[1:100,i],ds_updated[1:100,i]), na.rm = TRUE)
      plot(ds_o[1:100,i], type = "l", col = "black", ylim = c(0,max_y), ylab = "", xlab = "", main = i)
      if(use_sojourn_prior) points(d_prior[1:100,i], type = "l", col = "blue")
      points(ds_updated[1:100,i], type = "l", col = "red", lty = 2)
    }
    par(mfrow = c(1,1))
  }

  new_model
}


.re_estimate_sojourn_distribution = function(j, model, d){

  this_model_sojourn = model$sojourn[[j]]
  sojourn_distribution = this_model_sojourn$type
  d = d[,j]

  M = length(d)

  if(sojourn_distribution == "nonparametric"){

    this_model_sojourn$d = d

  }else
    if(sojourn_distribution == "ksmoothed_nonparametric"){

    d = d + 1e-100
    ksmooth.thresh = 1e-20 #this is a threshold for which d(u) values to use - if we throw too many weights in the default density() seems to work quite poorly
    u = which(d>ksmooth.thresh)
    if(length(u)>1){
      d = density(u, weights = d[u], from = 1, n = M) %>%
        approx(.,xout=1:M) %>% pluck("y") %>% tidyr::replace_na(0)
      d = d/sum(d)
    }
    this_model_sojourn$d = d

  }else
    if(sojourn_distribution == "poisson"){
      #threshold for effective "0" when considering d(u)
      shiftthresh = 1e-20
      u = which(d>shiftthresh); maxshift =  min(u); Mtmp = max(u) ; U = maxshift:Mtmp
      shifts = sapply(1:maxshift, function(shift) .dpois.hsmm.sojourn(x = U,lambda=(U-shift)%*%d[U],shift=shift,log=TRUE)%*%d[U])
      shift = which.max(shifts); new_U = shift:Mtmp
      this_model_sojourn$shift = shift
      this_model_sojourn$lambda = (new_U-shift)%*%d[new_U]
  }else
    if(sojourn_distribution == "nbinom"){
      tmp = .fitnbinom(d)
      this_model_sojourn$shift = tmp['shift']
      this_model_sojourn$size =  tmp['size']
      this_model_sojourn$mu =  tmp['mu']
      this_model_sojourn$prob =  tmp['prob']
  }else
    if(sojourn_distribution == "gamma"){

      tmp = gammafit(1:length(d),wt=d)
      this_model_sojourn$shape = tmp$shape
      this_model_sojourn$scale = tmp$scale

  }else
    if(sojourn_distribution == "logarithmic"){
    d = d + 1e-100
    this_model_sojourn$shape = .logdistrfit(wt=d)

  }else
    if(sojourn_distribution == "lnorm"){
      this_model_sojourn$meanlog[i] = weighted.mean(log(1:M),d)
      this_model_sojourn$sdlog[i] = sqrt(cov.wt(data.frame(log(1:M)),d)$cov)
  }else{
    stop("Invalid sojourn distribution")
  }

  this_model_sojourn
}



.re_estimate_censored_obs_prob = function(model, X, w, N0){

  #all_levels = .get_all_possible_levels(model = model, with_missing = TRUE, continuous_var_binned = TRUE)
  var_names = names(model$marg_em_probs)
  Xb = .cut_continuous_var(model = model, X = X)
  Xb_with_w = Xb %>%
    dplyr::select(seq_id, t, all_of(var_names)) %>% # we only keep the seq_id, the time-points and the observations
    dplyr::full_join(., # we join with the weight, but we need to transform them to long format first
                     w %>% as.data.frame() %>%
                       magrittr::set_colnames(1:model$J) %>%
                       dplyr::mutate(seq_id = X$seq_id, t = X$t) %>%
                       tidyr::pivot_longer(col = c(-seq_id, -t), names_to = "state", values_to = "p") %>%
                       dplyr::mutate(state = as.integer(state)),
                     by = c("seq_id","t"))

  observed_censored_obs_probs = Xb_with_w %>%
    dplyr::group_by(.dots = c("state", var_names)) %>%
    dplyr::summarize(counts = sum(p), .groups = "drop")

  observed_N = observed_censored_obs_probs %>%
    dplyr::select(state, counts) %>%
    dplyr::group_by(state) %>%
    dplyr::summarize(N = sum(counts), .groups = "drop")

  censored_obs_probs = model$censored_obs_probs %>%
    dplyr::left_join(., observed_N, by = "state") %>%
    dplyr::left_join(.,
                     observed_censored_obs_probs,
                     by = c("state", var_names)) %>%
    dplyr::mutate(counts = counts %>% tidyr::replace_na(0),
                  p = (p0 * N0 + counts)/(N0 + N)) %>%
    dplyr::select(state, dplyr::all_of(var_names), p, p0)

  censored_obs_probs
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
  print(model$marg_em_probs)
}



