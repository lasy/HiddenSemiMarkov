

#' Plots a single sequence of observation.
#'
#'
#' @param X a \code{data.frame} specifying the sequence of observation.
#'     Should contains data for only one sequence. Each column is a variable.
#'     Any column starting with the characters \code{state} is considered to be a state sequence and will be displayed on top of the observations.
#'     Any column starting with the characters \code{state_prob} is considered to be the probability associated with the corresponding state at each time-point. That column will be used to alter the transparency of the state sequence visualization.
#'     Several state columns can be provided. For example: \code{"state_ground_truth"} and \code{"state_Viterbi"}.
#'     A state probability can be specified for each of the state column by providing, for example, the columns \code{"state_prob_ground_truth"} and \code{"state_prop_Viterbi"}.
#' @param model a \code{hsmm} or \code{hsmm_spec} object specifying the model associated with the observation sequence.
#' @param title (optional) a \code{character} specifying the title of the plot.
#' @param show_state_diff (optional) a \code{logical} specifying if, in the case there are two "state" columns, a third line showing the agreement between these two columns should be displayed.
#'     Default value is \code{TRUE}. This is useful if one desires to compare a decoded sequence with the ground truth or state sequence resulting from the decoding of models with different parameters.
#' @param add_state_color_legend (optional) a \code{logical} specifying if the color legend for the model latent states should be printed. Default value is \code{FALSE}.
#' @param compact_view (optional) a \code{logical} specifying if the visualization of the observed variables should be compact,
#'     i.e. using color-coding only to display each variable on a single line.
#' @param add_color_legend_in_compact_view (optional)  a \code{logical} specifying if the color legend should be added for each variable when displaying time-series in compact view. Default is \code{TRUE}.
#' @param selection (optional) a \code{data.frame} specifying the start and end of a "selection", i.e. a part of sequence that needs to be highlighted.
#'     If not \code{NULL} (the default value), this option allows to display a transparent rectangle across all state and variable lines of a color of a given state.
#'     The \code{data.frame} must have the following columns: \code{start, end, state (integer)}.
#' @param verbose (optional) a logical specifying if the internal steps of the function should be printed.
#'
#' @return a ggplot object.
#' @import ggplot2
#' @importFrom magrittr %>%
#' @export
#' @examples
#' my_model = simple_model
#' Xsim = simulate_hsmm(model = my_model, n_state_transition = 20)
#' plot_hsmm_seq(X = Xsim, model = my_model)
#' plot_hsmm_seq(X = Xsim, model = my_model, title = "Simulated sequence", add_state_color_legend = TRUE)
#' plot_hsmm_seq(X = Xsim, model = my_model, title = "Simulated sequence (compact view)", compact_view = TRUE)


plot_hsmm_seq = function(X, model,
                         title = NULL,
                         show_state_diff = TRUE,
                         compact_view = FALSE,
                         add_color_legend_in_compact_view = TRUE,
                         add_state_color_legend = FALSE,
                         selection = data.frame(),
                         verbose = FALSE){

  # CHECKS
  # we check that the provided sequence is compatible with the model definition
  X = .check_types_and_values(data = X, parem = model$marg_em_probs)
  # the 'selection' data.frame is used for interactive labeling, to show the selected data-points.
  if(nrow(selection)>0){
    if(!all(colnames(selection) %in% c("state","start","end"))) stop("The 'selection' data.frame must have the following columns: 'state','start','end'.")
  }
  # this function can only be used to visualize one sequence at a time.
  if(length(unique(X$seq_id))>1) warning("More than one sequence provided.")

  # shortcuts
  X_names = names(model$marg_em_probs)
  state_cols = model$state_colors

  # DATA PREPARATION (in case state information is provided in the matrix X)
  # adding states colors
  state_columns = stringr::str_detect(colnames(X),"^state")
  state_prob_columns = stringr::str_detect(colnames(X),"^state_prob")
  state_columns = state_columns & !state_prob_columns
  if(any(state_columns)){
    # for each state column, we define a new column with the same name + "_col" which will have the colors of the state
    for(i in which(state_columns)){
      df = data.frame(state_col = state_cols[X[,i] %>% unlist()])
      colnames(df) = stringr::str_c(colnames(X)[i],"_col")
      X = cbind(X, df)
    }

    # then, if there are exactly two state columns, we create a third state column that has the difference between the two state columns
    if((sum(state_columns) == 2) & show_state_diff){
      X = X %>%
        dplyr::mutate(
          state_diff = (X[,which(state_columns)[1]] == X[,which(state_columns)[2]]),
          state_diff_col = c("red","blue")[state_diff+1]
          )
    }
  }

  # we also add the state color to the selection.
  if(nrow(selection)>0)
    selection = selection %>% dplyr::mutate(state_color = state_cols[state], t = start)

  # VISUALIZATIONS
  # We first define two ggplots objects: g_axis and g_base.
  # g_axis contains x-axis info + the selection + general theme
  g_axis = ggplot(X, aes(x = t))
  if(nrow(selection)>0){
    g_axis = g_axis +
      geom_rect(data = selection, aes(xmin = t, xmax = end, ymin = -Inf, ymax = Inf, fill = state_color), alpha = 0.3) +
      scale_fill_identity() +
      ggnewscale::new_scale_fill()
  }
  g_axis = g_axis +
    scale_x_continuous(limits = c(min(X$t)-1, max(X$t)+1), expand = expansion(add = 0)) +
    theme_set(theme_minimal()) + theme(panel.background = element_rect(color = "transparent", fill = "transparent"))+
    guides(alpha = FALSE) # fill = FALSE,

  # g_base is based on g_axis but all elements of x axes are removed.
  g_base = g_axis +
    theme_set(theme_minimal())+
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 10))


  # PLOTLIST
  # plotlist holds the time-series viz starting from the top to the bottom:
  # title + states + variables + x-axis
  plotlist = list()
  rel_heights = c() # the relative height of each plot.

  # title
  if(!is.null(title)){
    if(!is.character(title)) stop("title must be of type character \n")
    g_title = ggplot() + ggtitle(title)
    plotlist[["title"]] = g_title
    rel_heights = c(rel_heights, 1)
  }

  # states
  if(any(state_columns)){
    j = stringr::str_which(colnames(X),"state.*_col")
    for(i in j){

      # state name
      state_label = colnames(X)[i] %>%
        stringr::str_remove(.,"state") %>% stringr::str_remove(.,"_col") %>%
        stringr::str_replace(.,"_"," ") %>%
        stringr::str_remove(.,"^ ")

      cols_to_select = c("seq_id","t",colnames(X)[i]); cols_names = c("seq_id","t","state_col")

      # checking if the state has a probability attached to it.
      state_prob_column = ifelse(state_label == "",
                                 "state_prob",
                                 stringr::str_c("state_prob",state_label, sep = "_"))
      k = which(colnames(X) == state_prob_column)
      if(length(k)>0){ cols_to_select = c(cols_to_select, state_prob_column); cols_names = c(cols_names, "state_alpha") }
      # selecting and formatting the columns of interest
      df = X %>% dplyr::select(dplyr::all_of(cols_to_select)) %>% magrittr::set_colnames(cols_names)
      if(!("state_alpha" %in% cols_names)) df = df %>% dplyr::mutate(state_alpha = 1)

      # plot title
      if(stringr::str_length(state_label)>0) state_label = stringr::str_c("states (",state_label,")") else state_label = "states"

      g_state = g_base +
        geom_tile(data = df,
                  aes(y = 1, fill = state_col, alpha = state_alpha),
                  na.rm = TRUE)+ #  width = 1, height = 1,
        scale_alpha(limits = c(0,1), range = c(0, 1)) +
        scale_fill_identity("",
                            guide = "legend",
                            labels = model$state_names,
                            breaks = model$state_colors) +
        theme(legend.position = "right",
              legend.direction = "horizontal",
              legend.title = element_blank())

      if(compact_view) g_state = g_state + scale_y_continuous(breaks = 1, labels = state_label)
      if(!compact_view) g_state = g_state + scale_y_continuous(breaks = NULL) + ggtitle(state_label)

      if(!(compact_view & add_color_legend_in_compact_view)) g_state = g_state + guides(fill = FALSE)
      if(i > j[1]) g_state = g_state + guides(fill = FALSE)

      plotlist[[state_label]] = g_state
      rel_heights = c(rel_heights,ifelse(compact_view,1,2))
    }
  }

  # variables
  for(var in X_names){
    if(verbose) cat(var, "\n")
    X_var = X %>%
      dplyr::select(t, all_of(var)) %>%
      dplyr::mutate(Y = X[,var] %>%  unlist(),
                    y = Y)
    if(compact_view) X_var = X_var %>% dplyr::mutate(y = 1)

    g_var = g_base + guides(col = FALSE) # , fill = FALSE

    if(model$marg_em_probs[[var]]$type == "non-par"){
      X_var = X_var %>%
        dplyr::mutate(color = model$marg_em_probs[[var]]$viz_options$colors[Y %>% as.numeric()])

      g_var = g_var +
        geom_tile(data = X_var, aes(y = y, fill = color),
                  na.rm = TRUE) +
        scale_fill_identity("",
                            guide = "legend",
                            labels = model$marg_em_probs[[var]]$params$values,
                            breaks = model$marg_em_probs[[var]]$viz_options$colors)

      this_var_rel_height = 1 + sqrt(length(model$marg_em_probs[[var]]$params$values))
    }

    if(model$marg_em_probs[[var]]$type == "binom"){
      color_max = model$marg_em_probs[[var]]$viz_options$color_max
      max_y = max(model$marg_em_probs[[var]]$params$size)
      g_var = g_var +
        geom_tile(data = X_var,
                  aes(y = y, fill = Y),
                  na.rm = TRUE) +
        scale_fill_gradient("", low = "gray", high = color_max, na.value = "transparent",
                            breaks = 0:max_y,
                            guide = guide_colourbar(nbin = max_y+1, raster = FALSE, frame.colour = "white"))

      if(!compact_view) g_var = g_var +
        scale_y_continuous(limits = c(-0.51,max_y+0.51), breaks = 0:max_y, minor_breaks = NULL)

      this_var_rel_height = 1 + sqrt(max_y)
    }

    if(model$marg_em_probs[[var]]$type == "norm"){
      color_high = model$marg_em_probs[[var]]$viz_options$color_high
      color_low = model$marg_em_probs[[var]]$viz_options$color_low
      color_mid = model$marg_em_probs[[var]]$viz_options$color_mid
      mid_value = model$marg_em_probs[[var]]$viz_options$mid_value

      if(compact_view){
        g_var = g_var +
          geom_tile(data = X_var,
                    aes(y = 1, fill = Y),
                    na.rm = TRUE) +
          scale_fill_gradient2(low = color_low, high = color_high, mid = color_mid,
                               midpoint = mid_value, na.value = "transparent")
        this_var_rel_height = 1
      }else{
        g_var = g_var +
          geom_line(data = X_var,
                    aes(y = y, col = Y),
                    na.rm = TRUE) +
          geom_point(data = X_var,
                     aes(y = y, col = Y),
                     size = 0.5,
                     na.rm = TRUE) +
          scale_color_gradient2(low = color_low, high = color_high, mid = color_mid,
                                midpoint = mid_value, na.value = "transparent")

        this_var_rel_height = 4
      }
    }

    if(model$marg_em_probs[[var]]$type == "beta"){
      color_high = model$marg_em_probs[[var]]$viz_options$color_high
      color_low = model$marg_em_probs[[var]]$viz_options$color_low
      color_mid = model$marg_em_probs[[var]]$viz_options$color_mid
      mid_value = 0.5

      if(compact_view){
        g_var = g_var +
          geom_tile(data = X_var,
                    aes(y = 1, fill = Y),
                    na.rm = TRUE) +
          scale_fill_gradient2(low = color_low, high = color_high, mid = color_mid,
                               midpoint = mid_value, na.value = "transparent")
        this_var_rel_height = 1
      }else{
        g_var = g_var +
          geom_line(data = X_var,
                    aes(y = y, col = Y),
                    na.rm = TRUE) +
          geom_point(data = X_var,
                     aes(y = y, col = Y),
                     size = 0.5,
                     na.rm = TRUE) +
          scale_y_continuous(limits = c(0,1))+
          scale_color_gradient2(low = color_low, high = color_high, mid = color_mid,
                                midpoint = mid_value, na.value = "transparent")

        this_var_rel_height = 4
      }
    }

    if(compact_view) g_var = g_var + scale_y_continuous(breaks = 1, labels = var)
    if(!compact_view) g_var = g_var + ggtitle(var)
    if(compact_view) this_var_rel_height = 1
    if(compact_view) g_var = g_var +
      theme(legend.position = "right",
            legend.direction = "horizontal",
            legend.title = element_blank())
    if(!(compact_view & add_color_legend_in_compact_view)) g_var = g_var + guides(fill = "none")




    plotlist[[var]] = g_var
    rel_heights = c(rel_heights, this_var_rel_height)
  }


  plotlist$axis = g_axis
  rel_heights = c(rel_heights, 1)

  g = cowplot::plot_grid(plotlist = plotlist, align = "v", axis = "lrtb", ncol = 1, rel_heights = rel_heights)
  if(add_state_color_legend & (!(compact_view & add_color_legend_in_compact_view)))
    g = cowplot::plot_grid(g, plot_hsmm_state_colors(model = model), nrow = 1, rel_widths = c(10,1))

  g
}



#' Visualization of the model's state colors.
#'
#' @param model a \code{hsmm} model.
#'
#' @export
#' @importFrom magrittr %>%
#' @import ggplot2
plot_hsmm_state_colors = function(model){
  df = data.frame(state_num = 1:model$J,
                  state_names = model$state_names %>% factor(., levels = model$state_names),
                  state_cols = model$state_colors)
  ggplot(df, aes(x = 1, y = 1, fill = state_cols))+
    geom_tile()+
    scale_fill_identity()+
    xlab("")+ylab("")+
    facet_wrap(state_names ~ ., dir = "v", nrow = 6)+
    theme_set(theme_minimal())+
    theme(panel.grid = element_blank(),
          axis.text = element_blank(),
          strip.text = element_text(face = 2))

}




#' Visualizes the marginal emission distribution of a \code{hsmm} model.
#'
#' @param model a \code{hsmm} object specifying the model for which the marginal emission distributions should be visualized.
#' @param show_missing_probs (optional) a \code{logical} specifying if transparency should be used to reflect how likely variables are going to be missing in each state.
#' @param verbose (optional) a \code{logical} specifying if the internal steps of the function should be printed.
#'
#' @return a ggplot object.
#' @export
#' @importFrom magrittr %>%
#' @import ggplot2
#' @examples
#' my_model = simple_model
#' plot_hsmm_marg_dist(model = my_model)
plot_hsmm_marg_dist = function(model, show_missing_probs = TRUE, verbose = FALSE){

  X_names = names(model$marg_em_probs)
  state_names = model$state_names
  state_cols = model$state_colors

  obs_probs = model$obs_probs

  df = purrr::map_dfr(
    .x = X_names,
    .f = function(var){
      this_var_marg_dist =
        obs_probs %>%
        dplyr::select(state, dplyr::all_of(var), p) %>%
        dplyr::group_by(.dots = c("state", dplyr::all_of(var))) %>%
        dplyr::summarise(p = sum(p), .groups = "drop")
      this_var_marg_dist$var_name = var
      this_var_marg_dist$var_value = this_var_marg_dist[,var] %>% unlist()
      this_var_marg_dist = this_var_marg_dist %>%
        dplyr::arrange(state, var_value) %>%
        dplyr::mutate(var_value = var_value %>% as.character()) %>%
        dplyr::select(state, var_name, var_value, p)

      this_var_df = this_var_marg_dist %>%
        dplyr::filter(!is.na(var_value)) %>%
        dplyr::left_join(.,
                         this_var_marg_dist %>%
                           dplyr::filter(is.na(var_value)) %>%
                           dplyr::rename(missing_prob = p) %>%
                           dplyr::select(-var_value),
                         by = c("state", "var_name")) %>%
        dplyr::group_by(state) %>%
        dplyr::mutate(sum_p = sum(p),
                      p = p/sum_p) %>%
        dplyr::select(-sum_p) %>%
        dplyr::ungroup()

      this_var_df
    }
  )

  df = df %>%
    dplyr::mutate(var_value_fct = paste0(var_name,'_',var_value),
                  var_value_fct = var_value_fct %>%
                    factor(., levels = unique(var_value_fct)),
                  state_name = state_names[state] %>% factor(., levels = state_names))

  x_scale = df %>%
    dplyr::select(var_value, var_value_fct) %>%
    dplyr::distinct()

  if(!show_missing_probs)
    df$missing_prob = 0

  g = ggplot(df, aes(x = var_value_fct, y = p, alpha = 1-missing_prob, fill = state_name)) +
    geom_bar(stat = "identity") +
    scale_x_discrete(breaks = x_scale$var_value_fct, labels = x_scale$var_value) +
    xlab("") +
    ylab("probability") +
    scale_fill_manual(values = state_cols, guide = "none") +
    scale_alpha("probability of being observed", range = c(0,1), limits = c(0,1),
                guide = ifelse(show_missing_probs,"legend","none")) +
    facet_grid(state_name ~ var_name, scales = "free", space = "free") +
    theme_set(theme_minimal()) +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 90, hjust = 1),
          strip.background = element_rect(color = NA, fill = "gray90"),
          strip.text.y = element_text(angle = 0,hjust = 0))
  g
}






# plot_hsmm_marg_dist = function(model, show_missing_probs = TRUE, verbose = FALSE){
#
#   X_names = names(model$marg_em_probs)
#   state_names = model$state_names
#   state_cols = model$state_colors
#
#   DF_prob = data.frame()
#   DF_missing_prob = data.frame()
#   var_value_levels = data.frame(ix = 0, value = "init", stringsAsFactors = FALSE)
#   cols = c("state","variable", "value","prob")
#
#   for(var in X_names){
#     dist_type = model$marg_em_probs[[var]]$type
#     if(verbose) cat(var, "\n")
#     if(dist_type == "non-par"){
#       x = model$marg_em_probs[[var]]$params$values
#       df_prob = model$marg_em_probs[[var]]$params$probs %>%
#         as.data.frame() %>%
#         magrittr::set_colnames(state_names) %>%
#         dplyr::mutate(value = x) %>%
#         tidyr::pivot_longer(cols = -value, names_to = "state", values_to = "prob")
#     }else if(dist_type == "norm"){
#       means = model$marg_em_probs[[var]]$params$mean
#       sds = model$marg_em_probs[[var]]$params$sd
#       from = min(means-2.5*sds); to = max(means+2.5*sds)
#       x = axisTicks(c(from, to), log = FALSE, nint = 10)
#       df_prob = data.frame(state = state_names, mean = means, sd = sds)
#       df_prob = df_prob[rep(1:nrow(df_prob), each = length(x)),]
#       df_prob = df_prob %>%
#         dplyr::mutate(value = rep(x, model$J),
#                       prob = dnorm(value, mean = mean, sd = sd))
#     }else if(dist_type == "beta"){
#       shape1s = model$marg_em_probs[[var]]$params$shape1
#       shape2s = model$marg_em_probs[[var]]$params$shape2
#       from = 0; to = 1
#       x = axisTicks(c(from, to), log = FALSE, nint = 10)
#       df_prob = data.frame(state = state_names, shape1 = shape1s, shape2 = shape2s)
#       df_prob = df_prob[rep(1:nrow(df_prob), each = length(x)),]
#       df_prob = df_prob %>%
#         dplyr::mutate(value = rep(x, model$J),
#                       prob = dbeta(value, shape1 = shape1, shape2 = shape2))
#     }else if(dist_type == "binom"){
#       sizes = model$marg_em_probs[[var]]$params$size
#       probs = model$marg_em_probs[[var]]$params$prob
#       x = 0:max(sizes)
#       df_prob = data.frame(state = state_names, size = sizes, p = probs)
#       df_prob = df_prob[rep(1:nrow(df_prob), each = length(x)),]
#       df_prob = df_prob %>%
#         dplyr::mutate(value = rep(x, model$J),
#                       prob = dbinom(value, size = size, prob = p))
#     }
#     # trick to keep the levels
#     nix = last(var_value_levels$ix)
#     ix = (nix+1):(nix+length(x))
#     this_var_value_levels = data.frame(ix = ix, value = x, stringsAsFactors = FALSE)
#     var_value_levels = rbind(var_value_levels, this_var_value_levels)
#
#     df_prob = df_prob %>%
#       dplyr::mutate(variable = var,
#                     value = this_var_value_levels$ix[match(value, this_var_value_levels$value)]) %>%
#       dplyr::select(all_of(cols))
#     DF_prob = rbind(DF_prob, df_prob)
#
#     # missing probs
#     if(show_missing_probs){
#       df_missing_prob = data.frame(state = state_names %>% factor(),
#                                    variable = var,
#                                    missing_prob = model$censoring_probs$p + (1 - model$censoring_probs$p) *  model$censoring_probs$q[which(X_names == var), ],
#                                    stringsAsFactors = FALSE)
#     }else{
#       df_missing_prob = data.frame()
#     }
#     DF_missing_prob = rbind(DF_missing_prob, df_missing_prob)
#   }
#
#   DF_prob = DF_prob %>%
#     dplyr::mutate(state = factor(state, levels = state_names),
#                   original_prob = prob)
#
#   DF_prob = DF_prob %>%
#     dplyr::group_by(variable) %>%
#     dplyr::mutate(max_prob = max(original_prob),
#                   prob = original_prob/max_prob)
#
#   if(show_missing_probs){
#     DF_prob = DF_prob %>% dplyr::full_join(., DF_missing_prob, by = c("state", "variable"))
#     DF_prob = DF_prob %>% dplyr::mutate(observed_prob = 1-missing_prob)
#   }else{
#     DF_prob = DF_prob %>% dplyr::mutate(observed_prob = 1)
#   }
#
#
#
#   g_prob = ggplot(DF_prob, aes(x = value, y = prob, fill = state, alpha = observed_prob))
#   g_prob = g_prob +
#     facet_grid(state ~ variable, scales = "free_x", space = "free")+
#     scale_fill_manual(values = state_cols)+
#     guides(fill = FALSE)+
#     ylab("")+
#     scale_x_continuous(breaks = var_value_levels$ix, labels = var_value_levels$value, expand = c(0,0))+
#     scale_alpha_continuous("Probability of being observed", range = c(0.2,1), limits = c(0,1))+
#     geom_bar(stat = "identity", orientation = "x") +
#     theme_set(theme_minimal()) +
#     theme(axis.title.y = element_blank(),
#           axis.text.y = element_blank(),
#           axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
#           panel.grid.minor = element_blank(),
#           panel.spacing.x = unit(30,unit = "pt"),
#           strip.background = element_rect(fill = "gray90", color = "transparent"),
#           strip.text.y = element_text(angle = 0, hjust = 0),
#           legend.position = "bottom")
#
#
#   g = g_prob
#   g
# }



#' Visualizes the model graph with transition probabilities.
#'
#' The edges width is proportional to the transition probability.
#'
#' @param model a \code{hsmm} or \code{hsmm_spec} object specifying the model for which the transition probabilities should be visualized.
#' @param size (optional) the size of the nodes (hidden states).
#' @param label_size (optional) the size of the labels, i.e. hidden state names.
#' @param label_color (optional) the color of the labels. Default is white. Any valid color specification can be used here.
#' @param arrow_gap (optional) a \code{double} value in [0,1] specifying the gap between the arrow tips and the nodes.
#'    If \code{NULL} (default value), the function adapts the value of this parameter to the specified size.
#'
#' @return a ggplot object.
#' @export
#' @importFrom magrittr %>%
#' @import ggplot2 dplyr geomnet
#' @examples
#' my_model = simple_model
#' plot_hsmm_transitions(model = my_model)
plot_hsmm_transitions = function(model,
                                 size = 20,
                                 label_size = NULL,
                                 label_color = "white",
                                 arrow_gap = NULL){

  state_names = model$state_names
  state_cols = model$state_colors

  if(is.null(label_size)) label_size = 0.3*size
  if(is.null(arrow_gap)) arrow_gap = 0.025+0.0015*size


  transition_mat = model$transition
  colnames(transition_mat) = state_names; rownames(transition_mat) = state_names


  t_long = transition_mat %>%
    as.data.frame() %>%
    dplyr::mutate(from_id = rownames(transition_mat), from_nb = 1:nrow(transition_mat)) %>%
    tidyr::pivot_longer(cols = tidyselect::all_of(state_names), names_to = "to_id", values_to = "transition_prob") %>%
    dplyr::filter(transition_prob > 0) %>%
    dplyr::mutate(state_col = state_cols[from_nb])

  g = ggplot(t_long, aes(from_id = from_id, to_id = to_id)) +
    geom_net(
      aes(linewidth = transition_prob,
          color = state_col),
      fontsize = label_size,
      labelon = TRUE, labelcolour = label_color, vjust = 0.5,
      directed = TRUE,
      size = size,
      arrowgap = arrow_gap,
      layout.alg = "fruchtermanreingold"
    ) +
    theme_net() +
    theme(text = element_text(vjust = 0.5, face = "bold"))+
    scale_color_identity()

  g

}




#' Visualizes the sojourn distributions of a hidden semi-Markov model.
#'
#' @param model a \code{hsmm} object specifying the model for which the sojourn distributions should be visualized.
#' @param maxt (optional) an \code{integer} specifying the upper limit for the x-axis, which is the sojourn length.
#' @param one_panel_per_state (optional) a \code{logical} specifying if the sojourn distribution of each state should be displayed in separate vertically stacked panels. Default value is \code{FALSE}.
#'
#' @return a ggplot object.
#' @export
#' @importFrom magrittr %>%
#' @import ggplot2
#' @examples
#' my_model = simple_model
#' plot_hsmm_sojourn_dist(model = my_model)
plot_hsmm_sojourn_dist = function(model, maxt = 100, one_panel_per_state = FALSE){

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

  if(one_panel_per_state) g = g + facet_grid(state ~ ., scales = "free")

  g
}


#' Visualizes the joint emission probabilities (2x2) of a hidden semi-Markov model.
#'
#' @param model a \code{hsmm} object specifying the model for which the joint emission probabilities should be visualized.
#' @param title (optional) an character specifying the title of the visualization (typically the name of the model). By default, there is no title.
#'
#' @return a ggplot object.
#' @export
#' @importFrom magrittr %>%
#' @import ggplot2
#' @examples
#' my_model = simple_model
#' plot_hsmm_joint_em_prob(model = my_model)
plot_hsmm_joint_em_prob = function(model, title = ""){
  if(class(model) != "hsmm") stop("Needs an hsmm model (class = 'hsmm').")

  P = model$b %>%
    tidyr::pivot_longer(cols = tidyselect::starts_with("p_"), names_to = "state", values_to = "p") %>%
    dplyr::mutate(state = stringr::str_remove(state,"p_") %>% as.integer())
  P = P %>% mutate(state_name = model$state_names[state])

  all_levels = .get_all_possible_levels(model)
  var_names = colnames(all_levels)
  K = length(var_names)
  if(K==1) return(plot_hsmm_em_par(model = model))

  plotlist = list()
  i = 1
  for(k in 1:(K-1)){
    for(l in (k+1):K){
      this_P = P %>% select(all_of(var_names[c(k,l)]), p, state_name) %>% magrittr::set_colnames(c("x","y","p","state_name")) %>%
        group_by(x,y,state_name) %>%
        summarize(p = sum(p), .groups = "drop")
      g = ggplot(this_P, aes(x = x, y = y, fill = p)) +
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

  g = ggpubr::ggarrange(plotlist = plotlist , nrow = 1)

  if(title != ""){
    g_title = ggplot()+ggtitle(title)
    g = ggpubr::ggarrange(g_title, g, ncol = 1, heights = c(1,10))
  }

  g
}


#' Visualizes the status of the EM-procedure.
#'
#' If a model has been fitted using the function \code{fit_hsmm()}, this function can be used to visualize the convergence of the EM.
#'
#' @param fit_output the output of the \code{fit_hsmm()} function.
#' @param title (optional) an character specifying the title of the visualization (typically the name of the model). By default, there is no title.
#' @param y_axis_limits (optional) a 2-element vector specifying the limits of the y-axis.
#' @return a ggplot object.
#' @export
#' @importFrom magrittr %>%
#' @import ggplot2
#' @examples
#' my_model = simple_model
#' X_sim = simulate_hsmm(model = my_model, n_state_transitions = 20)
#' fit_results = fit_hsmm(model = my_model, X = X_sim)
#' plot_hsmm_fit_status(fit_output = fit_results)
plot_hsmm_fit_status = function(fit_output, title = NULL, y_axis_limits = NULL){
  if(is.null(title)) title = "EM status"
  df = data.frame(iter = 1:fit_output$fit_param$n_iter, log_likelihood = fit_output$fit_param$ll)

  color = "black"
  if(stringr::str_detect(fit_output$fit_param$message,"n_iter")) color = "orange"
  if(stringr::str_detect(fit_output$fit_param$message,"error")) color = "red"

  g = ggplot(df, aes(x = iter, y = log_likelihood))
  g = g + geom_point(col = color) + geom_line(col = color) +
    ggtitle(title, subtitle = fit_output$fit_param$message) +
    xlab("EM Iterations") + ylab("Log Likelihood") +
    scale_x_continuous(breaks = df$iter, minor_breaks = NULL)

  if(!is.null(y_axis_limits)) g = g + ylim(y_axis_limits)

  g
}

