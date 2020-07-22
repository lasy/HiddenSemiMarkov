

#' Plots a single sequence of observation.
#' @param X a \code{data.frame} specifying the sequence of observation. Should contains data for only one sequence. Each column is a variable.
#'     Any column starting with the characters \code{state} is considered to be a state sequence and will be displayed on top of the observations.
#'     If two columns start with \code{state}, then a third line will display the aggrement between these two columns. This is useful if one desires to compare a decoded sequence with the ground truth or state sequence resulting from the decoding of models with different parameters.
#' @param model a \code{hsmm} or \code{hsmm_spec} object specifying the model associated with the observation sequence.
#' @param verbose (optional) a logical specifying if the internal steps of the function should be printed.
#'
#' @return a ggplot object.
#' @import ggplot2
#' @importFrom magrittr %>%
#' @export
#' @examples
#' my_model = my_model_spec
#' Xsim = simulate_hsmm(model = my_model, n_state_transition = 20)
#' plot_hsmm_seq(X = Xsim, model = my_model)
#'
plot_hsmm_seq = function(X, model, verbose = FALSE, selection = data.frame()){

  # CHECKS

  X = .check_types_and_values(data = X, parem = model$parms.emission)

  if(nrow(selection)>0){
    if(!all(colnames(selection) %in% c("state","start","end"))) stop("The 'selection' data.frame must have the following columns: 'state','start','end'.")
  }

  if(length(unique(X$seq_id))>1) warning("More than one sequence provided.")
  state_cols = model$state_colors


  # shorcuts
  X_names = names(model$parms.emission)


  # data preparation: adding states colors
  state_columns = stringr::str_detect(colnames(X),"state")
  if(any(state_columns)){
    # for each state column, we will define a new column with the same name + "_col" which will have the colors of the state
    for(i in which(state_columns)){
      df = data.frame(state_col = state_cols[X[,i] %>% unlist()]); colnames(df) = stringr::str_c(colnames(X)[i],"_col")
      X = cbind(X, df)
    }

    # then, if there are exactly two state columns, we will create a third state column that will have the difference between the two state columns
    if(sum(state_columns) == 2){
      X = X %>%
        dplyr::mutate(
          state_diff = X[,which(state_columns)[1]] == X[,which(state_columns)[2]],
          state_diff_col = c("red","blue")[state_diff+1])
    }
  }

  if(nrow(selection)>0)
    selection = selection %>% dplyr::mutate(state_color = state_cols[state], t = start)


  g_axis = ggplot(X, aes(x = t))
  if(nrow(selection)>0){
    g_axis = g_axis +
      geom_rect(data = selection, aes(xmin = t, xmax = end, ymin = -Inf, ymax = Inf, fill = state_color), alpha = 0.3)
  }
  g_axis = g_axis +
    scale_x_continuous(limits = c(min(X$t)-1, max(X$t)+1), expand = expansion(add = 0))+
    theme_set(theme_minimal()) + theme(panel.background = element_rect(color = "transparent", fill = "transparent"))+
    guides(fill = FALSE)

  g_base = g_axis +
    theme_set(theme_minimal())+
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 10))

  plotlist = list()
  rel_heights = c()

  if(any(state_columns)){
    j = stringr::str_which(colnames(X),"state.*_col")
    for(i in j){
      df = X %>% dplyr::select(seq_id, t, matches(colnames(X)[i])) %>%
        magrittr::set_colnames(c("seq_id","t","state_col"))

      state_label = colnames(X)[i] %>%
        stringr::str_remove(.,"state") %>% stringr::str_remove(.,"_col") %>%
        stringr::str_replace(.,"_"," ") %>%
        stringr::str_remove(.,"^ ")
      if(stringr::str_length(state_label)>0) state_label = stringr::str_c("states (",state_label,")") else state_label = "states"

      g_state = g_base +
        geom_tile(data = df, aes(y = 1, width = 1, height = 1, fill = state_col))+
        scale_fill_identity()+
        ggtitle(state_label)+
        scale_y_discrete(breaks = NULL, expand = expansion(add = 0.05))

      plotlist[[state_label]] = g_state
      rel_heights = c(rel_heights,2)
    }
  }

  for(var in X_names){
    if(verbose) cat(var, "\n")
    X_var = X %>% dplyr::select(t, all_of(var)) %>% dplyr::mutate(y = X[,var] %>%  unlist())
    if(model$parms.emission[[var]]$type == "non-par"){
      if(is.numeric(model$parms.emission[[var]]$params$values)){
        X_var = X_var %>%  dplyr::mutate(y_num = y %>% as.character() %>%  as.numeric())
        y_breaks = model$parms.emission[[var]]$params$values
        y_limits = range(y_breaks)
        y_labels = model$parms.emission[[var]]$params$values
        y_ref = min(model$parms.emission[[var]]$params$values)
        add_line = TRUE
        add_point = FALSE
      }else{
        X_var = X_var %>%  dplyr::mutate(y_num = y %>% as.numeric())
        y_breaks = 1:length(model$parms.emission[[var]]$params$values)
        y_limits = c(0,max(y_breaks))
        y_labels = model$parms.emission[[var]]$params$values
        y_ref = 0
        add_line = TRUE
        add_point = FALSE
      }
      rel_heights = c(rel_heights, length(model$parms.emission[[var]]$params$values)-1)
    }else if(model$parms.emission[[var]]$type == "norm"){
      X_var = X_var %>%  dplyr::mutate(y_num = y)
      y_limits = range(X_var$y_num); y_limits = y_limits+c(-1,1)*max(abs(y_limits))*0.01
      y_breaks = seq(-10, 10, by = 1)
      y_labels = y_breaks
      y_ref = 0
      add_line = TRUE
      add_point = FALSE
      rel_heights = c(rel_heights, 5)
    }else if(model$parms.emission[[var]]$type == "binom"){
      X_var = X_var %>%  dplyr::mutate(y_num = y)
      y_breaks = 0:max(model$parms.emission[[var]]$params$size)
      y_limits = range(y_breaks)
      y_labels = y_breaks
      y_ref = 0
      add_line = FALSE
      add_point = TRUE
      rel_heights = c(rel_heights, pmax(2,length(y_breaks)))
    }

    g_var = g_base
    if(add_line) g_var = g_var + geom_segment(data = X_var, aes(xend = t, y = y_num, yend = y_ref, col = y_num))
    if(add_point) g_var = g_var + geom_point(data = X_var, aes(y = y_num, col = y_num))
    g_var = g_var +
      scale_y_continuous(breaks = y_breaks, labels = y_labels, minor_breaks = NULL, limits = y_limits)+
      scale_color_gradient2(low = "blue", high = "red", mid = "gray90", midpoint = 0)+
      scale_fill_identity()+
      guides(col = FALSE)+
      ggtitle(var)

    plotlist[[var]] = g_var
  }

  plotlist$axis = g_axis
  rel_heights = c(rel_heights, 1)

  g = cowplot::plot_grid(plotlist = plotlist, align = "v", axis = "ltb", ncol = 1, rel_heights = rel_heights)

  g


}




#' Visualize the marginal emission distribution of a \code{hsmm} or \code{hsmm_spec} model.
#'
#' @param model a \code{hsmm} or \code{hsmm_spec} object specifying the model for which the marginal emission distribution should be visualized.
#' @param verbose (optional) a logical specifying if the internal steps of the function should be printed.
#'
#' @return a ggplot object.
#' @export
#' @importFrom magrittr %>%
#' @import ggplot2
#'
plot_hsmm_em_par = function(model, verbose = FALSE, show_all_missing_probs = TRUE){

  X_names = names(model$parms.emission)
  state_names = model$state_names
  state_cols = model$state_colors

  DF_prob = data.frame()
  DF_missing_prob = data.frame()
  var_value_levels = data.frame(ix = 0, value = "init", stringsAsFactors = FALSE)
  cols = c("state","variable", "value","prob")

  for(var in X_names){
    dist_type = model$parms.emission[[var]]$type
    if(verbose) cat(var, "\n")
    if(dist_type == "non-par"){
      x = model$parms.emission[[var]]$params$values
      df_prob = model$parms.emission[[var]]$params$probs %>%
        as.data.frame() %>%
        magrittr::set_colnames(state_names) %>%
        dplyr::mutate(value = x) %>%
        tidyr::pivot_longer(cols = -value, names_to = "state", values_to = "prob")
    }else if(dist_type == "norm"){
      means = model$parms.emission[[var]]$params$mean
      sds = model$parms.emission[[var]]$params$sd
      from = min(means-2.5*sds); to = max(means+2.5*sds)
      x = axisTicks(c(from, to), log = FALSE, nint = 10)
      df_prob = data.frame(state = state_names, mean = means, sd = sds)
      df_prob = df_prob[rep(1:nrow(df_prob), each = length(x)),]
      df_prob = df_prob %>%
        dplyr::mutate(value = rep(x, model$J),
                      prob = dnorm(value, mean = mean, sd = sd))
    }else if(dist_type == "binom"){
      sizes = model$parms.emission[[var]]$params$size
      probs = model$parms.emission[[var]]$params$prob
      x = 0:max(sizes)
      df_prob = data.frame(state = state_names, size = sizes, p = probs)
      df_prob = df_prob[rep(1:nrow(df_prob), each = length(x)),]
      df_prob = df_prob %>%
        dplyr::mutate(value = rep(x, model$J),
                      prob = dbinom(value, size = size, prob = p))
    }
    # trick to keep the levels
    nix = last(var_value_levels$ix)
    ix = (nix+1):(nix+length(x))
    this_var_value_levels = data.frame(ix = ix, value = x, stringsAsFactors = FALSE)
    var_value_levels = rbind(var_value_levels, this_var_value_levels)

    df_prob = df_prob %>%
      dplyr::mutate(variable = var,
                    value = this_var_value_levels$ix[match(value, this_var_value_levels$value)]) %>%
      dplyr::select(all_of(cols))
    DF_prob = rbind(DF_prob, df_prob)

    # missing probs
    if(show_all_missing_probs | (!all(model$parms.emission[[var]]$missing_prob == 0))){
      df_missing_prob = data.frame(state = state_names %>% factor(),
                                   variable = var,
                                   missing_prob = model$parms.emission[[var]]$missing_prob,
                                   stringsAsFactors = FALSE)
    }else{
      df_missing_prob = data.frame()
    }
    DF_missing_prob = rbind(DF_missing_prob, df_missing_prob)
  }

  DF_prob = DF_prob %>%
    dplyr::mutate(state = factor(state, levels = state_names),
                  original_prob = prob)

  DF_prob = DF_prob %>%
    dplyr::group_by(variable) %>%
    dplyr::mutate(max_prob = max(original_prob),
                  prob = original_prob/max_prob)

  DF_prob = DF_prob %>% dplyr::full_join(., DF_missing_prob, by = c("state", "variable"))
  DF_prob = DF_prob %>% dplyr::mutate(observed_prob = 1-missing_prob)

  g_prob = ggplot(DF_prob, aes(x = value, y = prob, fill = state, alpha = observed_prob))
  g_prob = g_prob +
    facet_grid(state ~ variable, scale = "free_x", space = "free")+
    scale_fill_manual(values = state_cols)+
    guides(fill = FALSE)+
    ylab("")+
    scale_x_continuous(breaks = var_value_levels$ix, labels = var_value_levels$value, expand = c(0,0))+
    scale_alpha_continuous("Probability of being observed",range = c(0.2,1), limits = c(0,1))+
    geom_bar(stat = "identity", orientation = "x")+
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          panel.grid.minor = element_blank(),
          panel.spacing.x = unit(30,unit = "pt"),
          strip.text.y = element_text(angle = 0, hjust = 0),
          legend.position = "bottom")


  g = g_prob
  g
}







#' Visualizes the model graph with transition probabilities.
#'
#' The edges width is proportional to the transition probability.
#'
#' @param model a \code{hsmm} or \code{hsmm_spec} object specifying the model for which the transition probabilities should be visualized.
#' @param layout (optional) a matrix of dimension \code{model$J x 2} which provides the coordinates (x,y) of where the states should be placed on the canvas.
#' If not provided, the layout may slightly differ each time the function is ran.
#'
#' @return a ggplot object.
#' @export
#' @importFrom magrittr %>%
#' @import ggplot2
#'
plot_hsmm_transitions = function(model,layout = NULL){

  state_names = model$state_names
  state_cols = model$state_colors

  transition_mat = model$transition
  colnames(transition_mat) = state_names; rownames(transition_mat) = state_names

  g = graph_from_adjacency_matrix(transition_mat, weighted = TRUE, mode = "directed")
  if(is.null(layout)) layout = layout_nicely(g)

  #E(g)$width <- 1+E(g)$weight*5
  E(g)$width <- E(g)$weight
  plot(g,
       vertex.size=50, vertex.frame.color="white",
       vertex.label.color="white", vertex.label.family  = "sans",
       vertex.color = state_cols,
       layout = layout
  )
}




#' Visualizes the sojourn distributions of a hidden semi-Markov model.
#'
#' @param model a \code{hsmm} or \code{hsmm_spec} object specifying the model for which the sojourn distributions should be visualized.
#' @param maxt (optional) an integer specifying the upper limit for the x-axis, which is the sojourn length.
#'
#' @return a ggplot object.
#' @export
#' @importFrom magrittr %>%
#' @import ggplot2
#'
plot_hsmm_sojourn_dist = function(model, maxt = 100){

  state_names = model$state_names
  state_cols = model$state_colors

  if(is.null(model$sojourn$d)) d = .build_d_from_sojourn_dist(model = model, M = 1000) else d = model$sojourn$d
  M = nrow(d)
  df = data.frame(state = rep(1:model$J, each = M), t = rep(1:M, model$J), d = as.vector(d))
  df = df %>% filter(t <= maxt) %>%
    mutate(state_col = state_cols[state],
           state_name = state_names[state],
           state = factor(state))

  g = ggplot(df, aes(x = t, y = d, col = state, group = state))+
    geom_line()+
    scale_color_manual("States", values = state_cols, labels = state_names)
  g
}


#' Visualizes the joint emission probabilities of a hidden semi-Marvov model.
#'
#'
#' @param model a \code{hsmm} object specifying the model for which the joint emission probabilities should be visualized.
#' @param title (optional) an character specifying the title of the visualization (typically the name of the model). By default, there is no title.
#'
#' @return a ggplot object.
#' @export
#' @importFrom magrittr %>%
#' @import ggplot2
#'
plot_hsmm_joint_em_prob = function(model, title = ""){
  if(class(model) != "hsmm") stop("Needs an initialized model (class = 'hsmm'). Use the function 'initialize_hsmm()' to initialize your model.")

  P = model$obs_probs
  P = P %>% mutate(state_name = model$state_names[state])


  all_levels = .get_all_possible_levels(model)
  var_names = colnames(all_levels)
  K = length(var_names)
  if(K==1) return(plot_hsmm_em_par(model = model))

  plotlist = list()
  i = 1
  for(k in 1:(K-1)){
    for(l in (k+1):K){
      this_P = P %>% select(all_of(var_names[c(k,l)]), p_n, state_name) %>% magrittr::set_colnames(c("x","y","p_n","state_name")) %>%
        group_by(x,y,state_name) %>%
        summarize(p_n = sum(p_n), .groups = "drop")
      g = ggplot(this_P, aes(x = x, y = y, fill = p_n)) +
        geom_tile()+
        scale_fill_gradient(low = "white", high = "purple")+
        guides(fill = FALSE)+
        facet_grid(state_name ~ .) +
        xlab(var_names[k])+ylab(var_names[l])+
        theme(strip.text.y = element_text(angle = 0, hjust = 0))
      plotlist[[i]] = g
      i = i+1
    }
  }

  g = ggarrange(plotlist = plotlist , nrow = 1)

  if(title != ""){
    g_title = ggplot()+ggtitle(title)
    g = ggarrange(g_title, g, ncol = 1, heights = c(1,10))
  }

  g
}



plot_hsmm_local_prob_model = function(model){

  state_names = model$state_names
  state_cols = model$state_colors

  if(model$local_state_prob_method == "glm"){

    df = purrr::map_dfr(
      .x = 1:model$J,
      .f = function(s){
        glm_model = model$local_state_prob_models[[s]]
        coef = glm_model$coefficients
        data.frame(state = s, coef_name = names(coef), coef_value = coef)
      })

    df = df %>% mutate(state_cols = state_cols[state],
                       state_name = state_names[state] %>% factor(., levels = state_names),
                       coef_name = coef_name %>% factor(.,unique(df$coef_name)))

    g = ggplot(df, aes(x = coef_name, y = coef_value, col = state_cols))
    g = g +
      geom_hline(yintercept = 0, col = "gray")+
      geom_point()+
      geom_segment(aes(xend = coef_name, yend = 0))+
      scale_color_identity()+
      ylab("value of coefficient")+xlab("")+
      facet_grid(state_name ~ ., scale = "free")+
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            strip.text.y = element_text(angle = 0, hjust = 0))

  }else{
    stop(str_c("Visualization of the models following the model's method (",model$local_state_prob_method ,") not yet implemented."))
  }

  g
}


#' Visualizes the status of the EM-procedure.
#'
#' If a model has been fitted using the function \code{fit_hsmm()}, this function can be used to visualize the convergence of the EM.
#'
#' @param model the output of the \code{fit_hsmm()} function.
#' @param title (optional) an character specifying the title of the visualization (typically the name of the model). By default, there is no title.
#'
#' @return a ggplot object.
#' @export
#' @importFrom magrittr %>%
#' @import ggplot2
#'
plot_hsmm_fit_param = function(model, title = NULL){
  if(is.null(title)) title = "EM status"
  df = data.frame(iter = 1:model$fit_param$n_iter, log_likelihood = model$fit_param$ll)
  g = ggplot(df, aes(x = iter, y = log_likelihood))
  g = g + geom_point() + geom_line() + expand_limits(y = 0) + ggtitle(title, subtitle = model$fit_param$message) + xlab("EM Iterations") + ylab("Log Likelihood") +
    scale_x_continuous(breaks = df$iter, minor_breaks = NULL)
  g
}

### DEPRECATED ----------------------------



plot_hsmm_em_par_horizontal_states = function(model, verbose = FALSE, show_all_missing_probs = FALSE){

  X_names = names(model$parms.emission)
  state_names = model$state_names
  state_cols = model$state_colors

  DF_prob = data.frame()
  DF_missing_prob = data.frame()
  var_value_levels = data.frame(ix = 0, value = "init", stringsAsFactors = FALSE)
  cols = c("state","variable", "value","prob")

  for(var in X_names){
    dist_type = model$parms.emission[[var]]$type
    if(verbose) cat(var, "\n")
    if(dist_type == "non-par"){
      x = model$parms.emission[[var]]$params$values
      df_prob = model$parms.emission[[var]]$params$probs %>%
        as.data.frame() %>%
        magrittr::set_colnames(state_names) %>%
        dplyr::mutate(value = x) %>%
        tidyr::pivot_longer(cols = -value, names_to = "state", values_to = "prob")
    }else if(dist_type == "norm"){
      means = model$parms.emission[[var]]$params$mean
      sds = model$parms.emission[[var]]$params$sd
      from = min(means-2.5*sds); to = max(means+2.5*sds)
      x = axisTicks(c(from, to), log = FALSE, nint = 10)
      df_prob = data.frame(state = state_names, mean = means, sd = sds)
      df_prob = df_prob[rep(1:nrow(df_prob), each = length(x)),]
      df_prob = df_prob %>%
        mutate(value = rep(x, model$J),
               prob = dnorm(value, mean = mean, sd = sd))
    }else if(dist_type == "binom"){
      sizes = model$parms.emission[[var]]$params$size
      probs = model$parms.emission[[var]]$params$prob
      x = 0:max(sizes)
      df_prob = data.frame(state = state_names, size = sizes, p = probs)
      df_prob = df_prob[rep(1:nrow(df_prob), each = length(x)),]
      df_prob = df_prob %>%
        mutate(value = rep(x, model$J),
               prob = dbinom(value, size = size, prob = p))
    }
    # trick to keep the levels
    nix = last(var_value_levels$ix)
    ix = (nix+1):(nix+length(x))
    this_var_value_levels = data.frame(ix = ix, value = x, stringsAsFactors = FALSE)
    var_value_levels = rbind(var_value_levels, this_var_value_levels)

    df_prob = df_prob %>%
      mutate(variable = var,
             value = this_var_value_levels$ix[match(value, this_var_value_levels$value)]) %>%
      select(all_of(cols))
    DF_prob = rbind(DF_prob, df_prob)

    # missing probs
    if(show_all_missing_probs | (!all(model$parms.emission[[var]]$missing_prob == 0))){
      df_missing_prob = data.frame(state = state_names %>% factor(),
                                   variable = var,
                                   missing_prob = model$parms.emission[[var]]$missing_prob,
                                   stringsAsFactors = FALSE)
    }else{
      df_missing_prob = data.frame()
    }
    DF_missing_prob = rbind(DF_missing_prob, df_missing_prob)
  }

  DF_prob = DF_prob %>%
    mutate(state = factor(state, levels = state_names))


  g_prob = ggplot(DF_prob, aes(x = prob, y = value, fill = state))
  g_prob = g_prob +
    facet_grid(variable ~ state, scale = "free_y", space = "free")+
    scale_fill_manual(values = state_cols)+
    guides(fill = FALSE)+
    ylab("")+
    scale_y_continuous(breaks = var_value_levels$ix, labels = var_value_levels$value, expand = c(0,0))+
    geom_bar(stat = "identity", orientation = "y")+
    ggtitle("Emission distributions")+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.spacing.y = unit(20,unit = "pt"),
          strip.text.y = element_text(angle = 0, hjust = 0))


  if(nrow(DF_missing_prob) > 0){

    DF_missing_prob = DF_missing_prob %>%
      mutate(state = factor(state, levels = state_names))



    g_missing = ggplot(DF_missing_prob, aes(x = state, xend = state, y = 0, yend = 1 - missing_prob, col = state))
    g_missing = g_missing +
      geom_segment(size = 2)+
      scale_color_manual(values = state_cols)+
      facet_grid(variable ~ state, scale = "free_x")+
      scale_y_continuous(breaks = seq(0,1,len = 3), minor_breaks = seq(0,1,len = 5))+
      guides(col = FALSE)+
      xlab("")+ylab("")+
      ggtitle("Probabilities of being observed")+
      theme(axis.title.x  = element_blank(),
            axis.text.x = element_blank(),
            panel.grid.major.x = element_blank(),
            strip.text.y = element_text(angle = 0, hjust = 0))

    g = cowplot::plot_grid(g_prob,g_missing , align = "v", ncol = 1,
                           rel_heights = c(1+1.5*length(unique(DF_prob$variable)),1+length(unique(DF_missing_prob$variable))))
  }else{
    g = g_prob
  }

  g

}


plot_hsmm_em_par_deprecated = function(model, state_names = NULL, state_cols = NULL, verbose = FALSE){

  X_names = names(model$parms.emission)
  if(is.null(state_names)) state_names = 1:model$J
  if(is.null(state_cols)) state_cols = rainbow(n = model$J)

  plotlist = list()
  rel_heights = c()

  for(var in X_names[c(1,2,3,4,5)]){
    if(verbose) cat(var, "\n")
    if(model$parms.emission[[var]]$type == "non-par"){
      df_prob = model$parms.emission[[var]]$params$probs %>%
        as.data.frame() %>%
        magrittr::set_colnames(state_names) %>%
        dplyr::mutate(value = model$parms.emission[[var]]$params$values %>%  factor(.,levels = model$parms.emission[[var]]$params$values)) %>%
        tidyr::pivot_longer(cols = -value, names_to = "state", values_to = "prob")
      rel_heights = c(rel_heights, length(model$parms.emission[[var]]$params$values))
    }else if(model$parms.emission[[var]]$type == "norm"){
      means = model$parms.emission[[var]]$params$mean
      sds = model$parms.emission[[var]]$params$sd
      x = seq(from = min(means-3*sds),to = max(means+3*sds), len = 15)
      df_prob = data.frame(state = state_names, mean = means, sd = sds)
      df_prob = df_prob[rep(1:nrow(df_prob), each = length(x)),]
      df_prob = df_prob %>%
        mutate(value = rep(x, model$J),
               prob = dnorm(value, mean = mean, sd = sd))
      #df_prob = df_prob %>% mutate(value = value  %>%  factor(., levels = x))

      rel_heights = c(rel_heights, 5)
    }else if(model$parms.emission[[var]]$type == "binom"){
      sizes = model$parms.emission[[var]]$params$size
      probs = model$parms.emission[[var]]$params$prob
      x = 0:max(sizes)
      df_prob = data.frame(state = state_names, size = sizes, p = probs)
      df_prob = df_prob[rep(1:nrow(df_prob), each = length(x)),]
      df_prob = df_prob %>%
        mutate(value = rep(x, model$J),
               prob = dbinom(value, size = size, prob = p))
      rel_heights = c(rel_heights, length(x))
    }

    df_prob = df_prob %>% mutate(state = factor(state, levels = state_names))
    df_missing_prob = data.frame(state = state_names %>% factor(),
                                 missing_prob = model$parms.emission[[var]]$missing_prob,
                                 stringsAsFactors = FALSE)

    g_prob = ggplot(df_prob, aes(x = prob, y = value, fill = state))
    if(model$parms.emission[[var]]$type == "binom") g_prob = g_prob + scale_y_continuous(breaks = unique(df_prob$value))
    g_prob = g_prob +
      facet_grid(. ~ state)+
      scale_fill_manual(values = state_cols)+
      guides(fill = FALSE)+
      geom_bar(stat = "identity", orientation = "y")+
      ylab(var)+
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            strip.text = element_blank(),
            panel.grid.minor = element_blank())

    plotlist[[str_c(var,"_prob")]] = g_prob


    g_missing = ggplot(df_missing_prob, aes(x = 1, y = 1, col = state, size = 1-missing_prob))
    g_missing = g_missing +
      geom_point()+
      scale_size(range = c(0,2), limits = c(0,1))+
      scale_color_manual(values = state_cols)+
      facet_grid(. ~ state)+
      guides(col = FALSE, size = FALSE)+
      xlab("")+ylab("")+
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            strip.text = element_blank(),
            panel.grid = element_blank())
    rel_heights = c(rel_heights, 1)
    plotlist[[str_c(var,"_missing")]] = g_missing
  }

  g = cowplot::plot_grid(plotlist = plotlist, align = "v", ncol = 1, rel_heights = rel_heights)

  g

}






### FUNCTIONS FROM PREVIOUS PACKAGE ------------------


plot.hsmm <- function(x,...) {
  tmp = x$model$d
  plot(1:nrow(tmp),tmp[,1],type='l',...,ylab="d(u)",xlab="u",ylim=range(tmp))
  for(i in 2:x$J)
    lines(tmp[,i],type='l',col=i)
  legend("topright",legend=1:x$J,col=1:x$J,lty=1)
}


plot.hsmm.data <- function(x,...) {
  plot(ts(x$x),...)
  if(!is.null(x$s)) .add.states(x$s,ht=axTicks(2)[1],time.scale=1)
  if(length(x$N)>1) abline(v=cumsum(x$N),lty=2)
}

addStates <- function (states,x=NULL,ybot = axTicks(2)[1], ytop=ybot + (axTicks(2)[2] - axTicks(2)[1])/5,dy  = ytop - ybot,greyscale = FALSE, leg = NA,
                       J = length(unique(states)), time.scale = 1, shiftx = 0)
{

  draw.it <- function(hats, ybot, ytop, cols, greyscale){
    ##cat("ybot", ybot, "ytop", ytop, "\n")
    for (ii in 1:length(hats$state)){
      if (greyscale) {
        rect(xleft   = hats$intervals[ii],
             ybottom = ybot,
             xright  = hats$intervals[ii + 1],
             ytop    = ytop,
             col = cols[hats$state[ii]], border = 1)
      } else {
        rect(xleft   = hats$intervals[ii],
             ybottom = ybot,
             xright  = hats$intervals[ii + 1],
             ytop    = ytop,
             col = cols[hats$state[ii]], border = cols[hats$state[ii]])
      }
    }
  }


  if (is.null(states)){
    states <- x
    if (!is.list(states))
      states <- list(states)
    x <- seq_along(states[[1]])
  } else {
    if (!is.list(states))
      states <- list(states)
    if(is.null(x)) x <- seq_along(states[[1]])
  }

  ##   cat("states:\n");
  ##   print(states)
  ##   cat("x:\n");
  ##   print(x)

  x <- as.numeric(x)
  rr  <- range(x)

  J = length(unique(states))
  if (greyscale) {
    cols <- c("#FFFFFF", "#F0F0F0", "#D9D9D9", "#BDBDBD", "#969696", "#737373", "#525252", "#252525")
  } else {
    cols <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3")
  }


  st.list <- lapply(states,
                    function(st){
                      runs = rle(st)
                      cs  <- cumsum(c(0, runs$lengths))
                      hats <- list(intervals=rr[1]+ diff(rr)*cs/max(cs), states=runs$values)
                      hats
                    })



  ##cat("dy:", dy, "\n")
  for (ii in seq_along(st.list)){
    draw.it (st.list[[ii]], ybot, ytop, cols, greyscale)
    ybot <- ytop + .2*dy
    ytop <- ybot + dy
  }

  if (any(!is.na(leg)))
    legend("topleft", legend = leg, fill = cols, bg = "white")
}




.add.states <- function(states,ht=0,greyscale=FALSE,leg=NA,J=length(unique(states)),time.scale=24,shift=0) {
  J = length(unique(states))

  if(greyscale) cols=c("#FFFFFF" ,"#F0F0F0" ,"#D9D9D9", "#BDBDBD" ,"#969696", "#737373", "#525252", "#252525")
  else cols = c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#B3B3B3") #kind regards to RBrewerPal for these values

  hats = rle(states)
  hats = list(intervals=cumsum(c(0,hats$lengths))/time.scale+shift,state=hats$values)
  for(ii in 1:length(hats$state))
    if(greyscale)  rect(hats$intervals[ii],ht,hats$intervals[ii+1],ht+(axTicks(2)[2]-axTicks(2)[1])/5,col=cols[hats$state[ii]],border=1)
  else rect(hats$intervals[ii],ht,hats$intervals[ii+1],ht+(axTicks(2)[2]-axTicks(2)[1])/5,col=cols[hats$state[ii]],border=cols[hats$state[ii]])
  if(any(!is.na(leg))) legend("topleft",legend=leg,fill=cols,bg="white")
}
