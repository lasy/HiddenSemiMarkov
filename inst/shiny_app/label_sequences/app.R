
ui <- fluidPage(

  selectInput(width = 380, inputId = "Sequence", label = "Select sequence",choices = all_seq, selected = selected_seq),
  hr(),

  h4("Manual labelling"),
  fluidRow(title = "Manual labelling",
           column(width = 3, selectInput(inputId = "State",label = "States", choices = model$state_names)),
           column(width = 2, style = "margin-top: 25px;", actionButton(inputId = "Save", label = "Save this label")),
           column(width = 4, style = "margin-top: 25px;", actionButton(inputId = "Erase", label = "Erase manual labels at the selected positions"))
  ),
  hr(),

  h4("Decoding validation"),
  fluidRow(title = "Decoding validation",
           column(width = 3, selectInput(inputId = "Decoding_type", label = "Select decoding",
                                         choices = decodings, selected = decoding)),
           column(width = 9, style = "margin-top: 25px;", actionButton(inputId = "Validate",
                                                                       label = "Validate the decoding at the selected positions"))
  ),
  hr(),

  plotOutput(outputId = "Sequence_viz"),

  fluidRow(
    column(
      sliderInput(inputId = "Selection_range", label = "Selection range",min = 1, max = min(max_t/2, 365), value = c(5,20), step = 1, round = TRUE, width = '100%'),
      offset = 1, width = 11
    )
  ),

  sliderInput(inputId = "Zoom_range", label = "Zoom Range",min = 1, max = max_t, value = c(1,min(max_t/2, 365)), step = 1, round = TRUE, width = '100%')
)


server <- function(input, output, session){

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
    this_label = isolate({
      data.frame(seq_id = input$Sequence,
                 t = input$Selection_range[1]:input$Selection_range[2],
                 state = match(input$State, model$state_names),
                 stringsAsFactors = FALSE);
    })
    # remove existing labels at the same position
    j = which((ground_truth$seq_id == this_label$seq_id) &
                (ground_truth$t %in% this_label$t))
    if(length(j)>0){ ground_truth <<- ground_truth[-j,] }
    # add the new labels
    ground_truth <<- rbind(ground_truth, this_label)
  })

  # Erasing manual labels (Erasing the current selection)
  observeEvent(input$Erase,{
    j = isolate({
      which((ground_truth$seq_id == input$Sequence)&
              ground_truth$t %in% (input$Selection_range[1]:input$Selection_range[2]))
    })
    if(length(j)>0){ ground_truth <<- ground_truth[-j,] }
  })


  # Validating decoded labels
  observeEvent(input$Validate,{
    this_label = isolate({
      X %>%
        dplyr::select(seq_id, t, all_of(state_columns[match(input$Decoding_type, decodings)])) %>%
        dplyr::filter(
          seq_id == input$Sequence,
          t %in% input$Selection_range[1]:input$Selection_range[2]
        ) %>%
        magrittr::set_colnames(c("seq_id","t","state"))
    })
    # remove existing labels at the same position
    j = which((ground_truth$seq_id == this_label$seq_id) &
                (ground_truth$t %in% this_label$t))
    if(length(j)>0){ ground_truth <<- ground_truth[-j,] }
    # add the new labels
    ground_truth <<- rbind(ground_truth, this_label)
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
    this_seq_X = X %>%
      dplyr::filter(seq_id == selected_seq, t %in% t1:t2) %>%
      dplyr::select(seq_id, t, all_of(names(model$parms.emission)), state_ground_truth, all_of(paste0("state_",input$Decoding_type)) )

    current_selection = data.frame(seq_id = input$Sequence,
                                   state = input$State,
                                   start = input$Selection_range[1]-0.5, end = input$Selection_range[2]+0.5)
    g = plot_hsmm_seq(X = this_seq_X, model = model) # , selection = current_selection
    print(g)
  })

}


shinyApp(ui = ui, server = server)
