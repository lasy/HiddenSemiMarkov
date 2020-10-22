

#' Interactive app for the labeling of sequences with mode states.
#'
#' The function \code{label_sequences()} launches an interactive app in your default web-browser that allows you to
#' (1) visualize observations and existing ground truth or other state sequences,
#' (2) label time-points with states of the provided \code{hsmm} or \code{hsmm_spec} model,
#' (3) delete any existing ground truth, and
#' (4) validate state sequences provided with the observation so that they are included in the ground truth.
#'
#' The app lets you select the sequence if the argument \code{X} contains several sequences and lets you zoom on particular sections of the sequences.
#' Another slider allows you to select specific time-points of the sequence for labeling them.
#' Closing the app window will stop the app and returns the updated ground truth.
#'
#' @param model a \code{hsmm} or \code{hsmm_spec} object: the hidden semi-Markov model defined to decode the observations
#' @param X a \code{data.frame} of observations.
#'    In addition to the sequence ID (\code{seq_id}), time-point (\code{t}) and variable columns, state sequences that the user would like to validate or compare to the ground truth can be specified as additional columns.
#'    These columns names should be \code{'state_xyz'} where \code{'xyz'} is a string identifying the decoding. \code{'xyz'} cannot be \code{'ground_truth'}, any other string is accepted.
#' @param ground_truth (optional) a \code{data.frame} providing the initial set of labeled time-points. If not provided, the ground truth is initialized to an empty \code{data.frame}.
#'
#' @return a \code{"ground truth" data.frame} which combines the initial ground truth (if provided) and the newly labeled time-points.
#' @importFrom magrittr %>%
#' @import shiny
#' @export
#' @examples
#' my_model = simple_model_spec
#' Xsim = simulate_hsmm(model = my_model, n_state_transition = 20, seq_id = "test_seq")
#' Xsim = dplyr::rename(Xsim, state_simulated_sequence = state)
#' ground_truth = label_sequences(model = my_model, X = Xsim)
#' if(nrow(ground_truth)>0) table(model$state_names[ground_truth$state])

label_sequences = function(model, X, ground_truth = data.frame(), verbose = FALSE){

  # Copying original X and ground truth
  provided_X = X
  provided_ground_truth = ground_truth

  # CHECKS
  model = .check_model(model = model)

  X = .check_data(data = X, model = model)
  if(nrow(ground_truth)>0)
    ground_truth = .check_ground_truth(ground_truth = ground_truth, model = model, X = X)

  # CHECKING FOR DECODINGS TO BE VALIDATED
  state_columns = colnames(provided_X)[grep("state", colnames(provided_X))]
  if(length(state_columns)>0){
    decodings = stringr::str_remove(state_columns, "state_?")
    if(any(stringr::str_length(decodings)== 0) | any(decodings == "ground_truth")) stop("If state sequences are provided as columns of X, their column name must start with 'state_name' where 'name' is specifying the type of state sequence. 'name' cannot be 'ground_truth'.")
  }else{decodings = c()}
  X = X %>% dplyr::left_join(., provided_X %>% dplyr::select(seq_id, t, dplyr::all_of(state_columns)), by = c("seq_id","t"))

  # adding the ground truth to X
  if(nrow(ground_truth)>0){
    X = X %>% dplyr::left_join(., ground_truth %>% dplyr::rename(state_ground_truth = state), by = c("seq_id","t"))
  }else{
    X$state_ground_truth = NA_integer_
  }

  # INITIALIZATION
  all_seq = unique(X$seq_id)
  selected_seq = all_seq[1]
  this_seq_X = X %>% dplyr::filter(seq_id == selected_seq)
  max_t = max(this_seq_X$t)
  if(length(decodings) > 0) decoding = decodings[1] else decoding = c()

  # INTERACTIVE GUI (SHINY APP)
  app <- shinyApp(
    ui =
      fluidPage(

        selectInput(width = 380, inputId = "Sequence", label = "Select sequence",choices = all_seq, selected = selected_seq),
        hr(),

        h4("Manual labelling"),
        fluidRow(title = "Manual labelling",
                 column(width = 3, selectInput(inputId = "State",label = "States", choices = model$state_names)),
                 column(width = 2, style = "margin-top: 25px;", actionButton(inputId = "Save", label = "Save this label")),
                 column(width = 4, style = "margin-top: 25px;", actionButton(inputId = "Erase", label = "Erase manual labels at the selected positions"))
        ),
        hr(),


        if(length(decodings) > 0) h4("Decoding validation"),
        if(length(decodings) > 0)
          fluidRow(title = "Decoding validation",
                   column(width = 3, selectInput(inputId = "Decoding_type", label = "Select decoding",
                                                 choices = decodings, selected = decoding)),
                   column(width = 9, style = "margin-top: 25px;", actionButton(inputId = "Validate",
                                                                               label = "Validate the decoding at the selected positions"))
          ),
        if(length(decodings) > 0) hr(),

        plotOutput(outputId = "Sequence_viz"),

        sliderInput(inputId = "Selection_range", label = "Selection range",min = 1, max = min(max_t/2, 365), value = c(5,20), step = 1, round = TRUE, width = '100%'),

        sliderInput(inputId = "Zoom_range", label = "Zoom Range",min = 1, max = max_t, value = c(1,min(max_t/2, 365)), step = 1, round = TRUE, width = '100%')
      ),

    server =
      function(input, output, session){

        # Updating sequence
        observeEvent(input$Sequence,{
          selected_seq <<- input$Sequence
        })

        # Update the zoom sliders when we change sequence
        observe({
          max_t = max(X$t[which(X$seq_id == input$Sequence)], na.rm = TRUE)
          maxi_slider = min(c(365,max_t,max_t/5))
          updateSliderInput(session, "Zoom_range", max = max_t, value = c(1, maxi_slider))
        })

        # Update the selection_range slider in function of the zoom
        observe({
          min_zoom = input$Zoom_range[1]
          max_zoom = input$Zoom_range[2]
          updateSliderInput(session, "Selection_range", min = min_zoom, max = max_zoom)
        })

        # Adding a manual label (Save selection)
        observeEvent(input$Save,{
          isolate({
            X <<- X %>%
              dplyr::mutate(state_ground_truth =
                              ifelse(
                                (seq_id == input$Sequence) & (t %in% input$Selection_range[1]:input$Selection_range[2]),
                                match(input$State, model$state_names),
                                state_ground_truth));
          })
        })

        # Erasing manual labels (Erasing the current selection)
        observeEvent(input$Erase,{
          isolate({
            X <<- X %>%
              dplyr::mutate(state_ground_truth =
                              ifelse(
                                (seq_id == input$Sequence) & (t %in% input$Selection_range[1]:input$Selection_range[2]),
                                NA,
                                state_ground_truth));
          })
        })


        # Validating decoded labels
        observeEvent(input$Validate,{
          isolate({
            j = which((X$seq_id == input$Sequence) & (X$t %in% input$Selection_range[1]:input$Selection_range[2]))
            X[j,"state_ground_truth"] <<- as.integer(unlist(X[j,paste0("state_",input$Decoding_type)]))
          })

        })

        # VIZ
        output$Sequence_viz = renderPlot({
          # the viz is changed any time:
          input$Save # a manual label is saved
          input$Erase # a manual label is erased
          input$Validate # a decoding is validated
          input$Sequence # a different sequence is selected
          input$Decoding_type # a different decoding is selected
          input$Zoom_range # the zoom range is changed
          input$Selection_range # the user selects a different region

          # we select the
          t1 = input$Zoom_range[1]; t2 = input$Zoom_range[2]
          cols_to_select = c("seq_id","t",names(model$marg_em_probs), "state_ground_truth")
          if(length(decodings)>0) cols_to_select = c(cols_to_select, paste0("state_",input$Decoding_type))

          this_seq_X = X %>%
            dplyr::filter(seq_id == selected_seq, t %in% t1:t2) %>%
            dplyr::select(all_of(cols_to_select))

          current_selection = data.frame(state = match(input$State, model$state_names),
                                         start = input$Selection_range[1]-0.5,
                                         end = input$Selection_range[2]+0.5)
          g = suppressWarnings(plot_hsmm_seq(X = this_seq_X, model = model, selection = current_selection)) #
          print(g)
        })

        session$onSessionEnded(function() {
          stopApp()
        })

      }
  )

  runApp(app, display.mode = "normal", launch.browser = TRUE)


  ground_truth = X %>% dplyr::select(seq_id, t, state_ground_truth) %>% dplyr::rename(state = state_ground_truth) %>% dplyr::filter(!is.na(state))
  ground_truth
}
