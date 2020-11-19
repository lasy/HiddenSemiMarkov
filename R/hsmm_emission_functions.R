
#' Provides the list of currently supported marginal emission distributions and the list of their parameters.
#'
#' @return A \code{data.frame} with the following columns:
#' \code{type} (continuous, discrete or categorical), \code{distribution}, \code{params}, \code{parameter_type} and \code{parameter_size}
#'
#' @export
#' @seealso \code{available_marginal_emission_viz_options()} for the list of available visualization options for each of these distributions.
#'
available_marginal_emission_dist = function(){

  available_marg_em_probs = rbind(
    data.frame(type = "continuous", distribution = "norm", params = "mean", parameter_type = "double", parameter_size = "J", stringsAsFactors = FALSE),
    data.frame(type = "continuous", distribution = "norm", params = "sd", parameter_type = "double", parameter_size = "J", stringsAsFactors = FALSE),

    data.frame(type = "continuous", distribution = "beta", params = "shape1", parameter_type = "double >=0", parameter_size = "J", stringsAsFactors = FALSE),
    data.frame(type = "continuous", distribution = "beta", params = "shape2", parameter_type = "double >=0", parameter_size = "J", stringsAsFactors = FALSE),

    data.frame(type = "discrete", distribution = "binom", params = "size", parameter_type = "int", parameter_size = "J", stringsAsFactors = FALSE),
    data.frame(type = "discrete", distribution = "binom", params = "prob", parameter_type = "double [0,1]", parameter_size = "J", stringsAsFactors = FALSE),

    data.frame(type = "categorical", distribution = "non-par", params = "values", parameter_type = "character", parameter_size = "NV", stringsAsFactors = FALSE),
    data.frame(type = "categorical", distribution = "non-par", params = "probs", parameter_type = "double [0,1]", parameter_size = "NV x J", stringsAsFactors = FALSE)
  )
  available_marg_em_probs
}


#' Provides the list of currently supported visualization options for the observations, depending on their emission distribution.
#'
#' Specifying viz options is mostly useful for the visualization of observation sequences with the function \code{plot_hsmm_seq()}.
#'
#' @return A \code{data.frame} with the following columns:
#'  \code{distribution}, \code{viz_options} and \code{option_description}.
#'
#' @export
#' @seealso \code{available_marginal_emission_dist()}, \code{plot_hsmm_seq()}
#' @examples
#' my_model_no_viz_options = specify_hsmm(
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
#' class(my_model_no_viz_options)
#' Xsim = simulate_hsmm(model = my_model_no_viz_options, n_state_transitions = 20)
#' plot_hsmm_seq(model = my_model_no_viz_options, X = Xsim)
#'
#' my_model = my_model_no_viz_options
#' my_model$marg_em_probs$var1$viz_options$color_low = "black"
#' my_model$marg_em_probs$var2$viz_options$color_max = "blue"
#' my_model$marg_em_probs$var3$viz_options$colors = c("purple","gray","blue","orange")
#' plot_hsmm_seq(model = my_model, X = Xsim)
#'
available_marginal_emission_viz_options = function(){
  marginal_em_viz_options = rbind(
    data.frame(distribution = "norm", viz_options = "color_low", option_description = "low values color", stringsAsFactors = FALSE),
    data.frame(distribution = "norm", viz_options = "color_high", option_description = "high values color", stringsAsFactors = FALSE),
    data.frame(distribution = "norm", viz_options = "color_mid", option_description = "mid (reference) value color", stringsAsFactors = FALSE),
    data.frame(distribution = "norm", viz_options = "mid_value", option_description = "value of reference. Default is 0.", stringsAsFactors = FALSE),
    data.frame(distribution = "beta", viz_options = "color_low", option_description = "low values color", stringsAsFactors = FALSE),
    data.frame(distribution = "beta", viz_options = "color_high", option_description = "high values color", stringsAsFactors = FALSE),
    data.frame(distribution = "beta", viz_options = "color_mid", option_description = "mid (reference) value color", stringsAsFactors = FALSE),
    data.frame(distribution = "binom", viz_options = "color_max", option_description = "color of the max value", stringsAsFactors = FALSE),
    data.frame(distribution = "non-par", viz_options = "colors", option_description = "colors for each variable category", stringsAsFactors = FALSE)
  )
  marginal_em_viz_options
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




.get_all_possible_levels = function(model = model, with_missing = TRUE, continuous_var_binned = FALSE){
  lv = list()
  for(var in names(model$marg_em_probs)){
    if(model$marg_em_probs[[var]]$type == "non-par") lv[[var]] = model$marg_em_probs[[var]]$params$values
    if(model$marg_em_probs[[var]]$type == "binom") lv[[var]] = 0:max(model$marg_em_probs[[var]]$params$size)
    if(model$marg_em_probs[[var]]$type == "norm") lv[[var]] = model$marg_em_probs[[var]]$breaks[!is.infinite(model$marg_em_probs[[var]]$breaks)]
    if(model$marg_em_probs[[var]]$type == "beta") lv[[var]] = model$marg_em_probs[[var]]$breaks
    if(continuous_var_binned & (model$marg_em_probs[[var]]$type %in% c("norm","beta"))) lv[[var]] = levels(cut(lv[[var]], breaks = model$marg_em_probs[[var]]$breaks))
    if(with_missing & (any(model$censoring_probs$p > 0) | any(model$censoring_probs$q > 0))) lv[[var]] = c(lv[[var]], NA) #TODO: potentially change that line because we may still want to include the NAs, even if no missing prob has been specified
  }
  max_n = max(lengths(lv))
  all_levels = data.frame(row.names = 1:max_n, stringsAsFactors = FALSE)
  for(var in names(lv)) all_levels = cbind(all_levels, var = rep(lv[[var]],max_n)[1:max_n])
  colnames(all_levels) = names(lv)
  all_levels = .check_types_and_values(data = all_levels, parem = model$marg_em_probs, continuous_var_binned = continuous_var_binned)
  all_levels
}


.get_marginal_prob = function(var_name, model){
  if(model$marg_em_probs[[var_name]]$type == "norm"){ # we need to discretize the probabilities
    # first, we define a continuous support that will cover the whole possible value range
    x_continuous = seq(min(model$marg_em_probs[[var_name]]$params$mean - 5 * model$marg_em_probs[[var_name]]$params$sd),
                       max(model$marg_em_probs[[var_name]]$params$mean + 5 * model$marg_em_probs[[var_name]]$params$sd),
                       len = 10000)
    marg_prob =
      tidyr::expand_grid(x_continuous = x_continuous,
                         state = 1:model$J) %>% # we expand for each state
      dplyr::mutate(x = cut(x_continuous, breaks = model$marg_em_probs[[var_name]]$breaks), # we discretize the support
                    d = dnorm(x = x_continuous, # we get the probability density
                              mean = model$marg_em_probs[[var_name]]$params$mean[state],
                              sd = model$marg_em_probs[[var_name]]$params$sd[state])) %>%
      dplyr::group_by(state, x) %>% dplyr::summarize(prob = sum(d), .groups = "drop") %>% # we sum the densities
      dplyr::group_by(state) %>% dplyr::mutate(tot = sum(prob), prob = prob/tot) %>% dplyr::select(-tot) # we normalize by state
  }
  if(model$marg_em_probs[[var_name]]$type == "binom"){
    marg_prob =
      tidyr::expand_grid(state = 1:model$J,
                         x = 0:max(model$marg_em_probs[[var_name]]$params$size)) %>%
      dplyr::mutate(prob = dbinom(x = x,
                                  size = model$marg_em_probs[[var_name]]$params$size[state],
                                  prob = model$marg_em_probs[[var_name]]$params$prob[state]))
  }
  if(model$marg_em_probs[[var_name]]$type == "non-par"){
    marg_prob =
      model$marg_em_probs[[var_name]]$params$probs %>%
      as.data.frame() %>%
      magrittr::set_colnames(1:model$J) %>%
      dplyr::mutate(x = model$marg_em_probs[[var_name]]$params$values) %>%
      tidyr::pivot_longer(cols = 1:model$J, names_to = "state", values_to = "prob") %>%
      dplyr::mutate(state = as.integer(state))
  }
  if(model$marg_em_probs[[var_name]]$type == "beta"){ # we need to discretize the probabilities
    # first, we define a continuous support that will cover the whole possible value range
    x_continuous = seq(0,1,len = 10000)
    marg_prob =
      tidyr::expand_grid(x_continuous = x_continuous,
                         state = 1:model$J) %>% # we expand for each state
      dplyr::mutate(x = cut(x_continuous, breaks = model$marg_em_probs[[var_name]]$breaks), # we discretize the support
                    d = dbeta(x = x_continuous, # we get the probability density
                              shape1 = model$marg_em_probs[[var_name]]$params$shape1[state],
                              shape2 = model$marg_em_probs[[var_name]]$params$shape2[state])) %>%
      dplyr::group_by(state, x) %>% dplyr::summarize(prob = sum(d), .groups = "drop") %>% # we sum the densities
      dplyr::group_by(state) %>% dplyr::mutate(tot = sum(prob), prob = prob/tot) %>% dplyr::select(-tot) # we normalize by state
  }
  marg_prob %>% dplyr::select(state, x, prob)
}


.cut_continuous_var = function(model, X){
  Xp = X
  for(var in names(model$marg_em_probs)) if(model$marg_em_probs[[var]]$type %in% c("norm","beta")) Xp[,var] = X[,var] %>% unlist() %>% cut(., breaks = model$marg_em_probs[[var]]$breaks)
  Xp
}


.bin_to_continuous = function(Xb, model){
  X = Xb
  for(var in names(model$marg_em_probs)){
    if(model$marg_em_probs[[var]]$type %in% c("norm", "beta")){
      breaks = model$marg_em_probs[[var]]$breaks
      df = data.frame(bin = cut(seq(breaks[2]-1e-10, breaks[length(breaks)-1]+1e-10, len = length(breaks)-1),
                                breaks = breaks),
                      min = breaks[1:(length(breaks)-1)],
                      max = breaks[2:length(breaks)]) %>%
        dplyr::mutate(min = ifelse(is.infinite(min),max,min),
                      max = ifelse(is.infinite(max),min,max)) %>%
        magrittr::set_colnames(c(var, "min","max"))
      Xb = Xb %>% dplyr::left_join(., df, by = var)
      j_not_na = which(!is.na(Xb[,var]))
      X[, var] = NA_real_
      X[j_not_na, var] = runif(length(j_not_na), min = Xb$min[j_not_na], max = Xb$max[j_not_na])
      Xb = Xb %>% dplyr::select(-min, -max) # for the next variable
    }
  }
  X
}
