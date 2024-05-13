library(chron)
library(data.table)
library(repifriends)
library(ggmap)
library(ggplot2)
library(gridExtra)
library(magick)
library(shiny)
library(shinyjs)
library(RANN)


ui <- fluidPage(
  useShinyjs(),
  titlePanel("EpiFRIenDs"),
  tabsetPanel(
    tabPanel("Distribution Analysis",
             #############################################
             #########   DISTRIBUTION ANALYSIS    ########
             #############################################
             ##### SIDEBAR LAYOUT
             sidebarLayout(
               sidebarPanel(
                 numericInput(inputId = "min_neighbours",
                              label = "Min. number of neighbors:",
                              value = 2, min = 0, max = 10),
                 numericInput(inputId = "link_d",
                             label = "Select linking distance:",
                             value = 0,min = 0, step = 0.1),
                 selectInput(inputId = "prevalence",
                             label = "Select column on data refering to prevalence:",
                             choices = c("None", "prevalence"),
                             selected = "None"),
                 selectInput(inputId = "keep_null_tests",
                             label = "How to deal with missing data",
                             choices = c("Remove", "Keep", "Impute"),
                             selected = "Centroid"),
                 conditionalPanel(
                   condition = "input.keep_null_tests == 'Impute'",
                   numericInput(inputId = "imp_value",
                                label = "Imputation value:",
                                value = 0)
                 ),
                 selectInput(inputId = "in_latlon",
                             label = "Treatment of coordinates:",
                             choices = c("Longitude/Latitude", "Cartesian"),
                             selected = "Longitude/Latitude"),
                 conditionalPanel(
                   condition = "input.in_latlon == 'Longitude/Latitude'",
                   textInput(inputId = "to_epsg",
                             label = "EPSG number for the projection to use:",
                             value = "32736")),
                 selectInput(inputId = "data_type",
                             label = "Type of data:",
                             choices = c("Mock up", "Real"),
                             selected = "Mock up"),
                 # Collapsible section
                 tags$details(
                   tags$summary(
                     tags$span(
                       id = "advancedToggle",
                       style = "cursor: pointer;",
                       HTML("&nbsp;&nbsp;&#43;")
                     ),
                     "Advanced Options"
                   ),
                   tags$div(
                     class = "advanced-options",
                     id = "advancedOptions",
                     selectInput(inputId = "optimize_link_d",
                                 label = "Optimize linking distance: ",
                                 choices = c("Yes", "No"),
                                 selected = "No"),
                     selectInput(inputId = "false_detect",
                                 label = "Consider False Detections (plus exec time): ",
                                 choices = c("Yes", "No"),
                                 selected = "No"),
                     conditionalPanel(
                       condition = "input.false_detect == 'Yes'",
                       numericInput(inputId = "n_simulations",
                                    label = "Number of simulations to perform",
                                    value = 10000, step = 50)
                     ),
                     selectInput(inputId = "link_neighbors",
                                 label = "Use linking neighbors instead of linking distance: ",
                                 choices = c("Yes", "No"),
                                 selected = "No"),
                     selectInput(inputId = "method",
                                 label = "Select method to account fot local prevalence:",
                                 choices = c("base", "kmeans", "centroid"),
                                 selected = "base"),
                     conditionalPanel(
                       condition = "input.method == 'centroid'",
                       numericInput(inputId = "max_epi_cont",
                                    label = "Percentage of epifriends contribution to the cluster:",
                                    value = 0.5, min = 0, max = 1, step = 0.1),
                       numericInput(inputId = "max_thr_data",
                                    label = "Percentage of all data for local prevalence calculus:",
                                    value = 0.1, min = 0, max = 1, step = 0.1)
                     ),
                     numericInput(inputId = "max_p",
                                  label = "Maximum p-value:",
                                  value = 1, min = 0, max = 1, step = 0.1),
                     numericInput(inputId = "min_pos",
                                  label = "Minimum number of positives in a cluster:",
                                  value = 2, step = 1),
                     numericInput(inputId = "min_total",
                                  label = "Minimum number of entries in a cluster:",
                                  value = 2, step = 1),
                     numericInput(inputId = "min_pr",
                                  label = "Minimum positivity rate:",
                                  value = 0, step = 0.1)
                   )
                 ),
                 fileInput(inputId = "data",
                           label = "Load data (CSV format):"),
                 actionButton(inputId = "load_data_button", label = "Load Data"),
                 actionButton(inputId = "epifriends_run", label = "Run Spatial Analysis")
               ),
               ##### DIFFERENT PANELS INSIDE DISTRIBUTION ANALYSIS
               mainPanel(
                 tabsetPanel(
                   tabPanel("Data",
                            tableOutput(outputId = "load_data")
                   ),
                   tabPanel("Distribution",
                            plotOutput(outputId = "distribution")
                   ),
                   tabPanel("Summary-EpiFRIenDs",
                            h3("Summary table with detected clusters and significance"),
                            br(),
                            tableOutput(outputId = "epifriends_summary"),
                            br(),
                            h3("Summary plots"),
                            br(),
                            plotOutput(outputId = "epifriends_v1"),
                            plotOutput(outputId = "epifriends_v3"),
                            plotOutput(outputId = "epifriends_v2")
                   ),
                   tabPanel("Clusters-EpiFRIenDs",
                            h3("Coordinates & clusters detected by Epifriends"),
                            br(),
                            tableOutput(outputId = "epifriends_results")
                   )
                 )
               )
             )
    ),
    
    tabPanel("Temporal Analysis", 
             #############################################
             #########     TEMPORAL ANALYSIS      ########
             #############################################
             sidebarLayout(
               sidebarPanel(
                 numericInput(inputId = "min_neighbours_temp",
                              label = "Min. number of neighbors:",
                              value = 2, min = 0, max = 10),
                 selectInput(inputId = "prevalence_temp",
                             label = "Select column on data refering to prevalence:",
                             choices = c("None", "prevalence"),
                             selected = "None"),
                 numericInput(inputId = "link_d_temp",
                              label = "Select linking distance:",
                              value = 0,min = 0, step = 0.1),
                 numericInput(inputId = "time_width_temp",
                              label = "Width of the time window (in number of days)",
                              value = 180),
                 numericInput(inputId = "time_steps_temp",
                              label = "Number of days moved forward in each time step",
                              value = 90),
                 numericInput(inputId = "linking_time_temp",
                              label = "Time steps considered to link temporal clusters",
                              value = 3),
                 numericInput(inputId = "linking_dist_temp",
                              label = "Spatial distance (to link clusters from different timeframes)",
                              value = 0.15),
                 selectInput(inputId = "keep_null_tests_temp",
                             label = "How to deal with missing data",
                             choices = c("Remove", "Keep", "Impute"),
                             selected = "Centroid"),
                 conditionalPanel(
                   condition = "input.keep_null_tests_temp == 'Impute'",
                   numericInput(inputId = "imp_value_temp",
                                label = "Imputation value:",
                                value = 0)
                 ),
                 selectInput(inputId = "in_latlon_temp",
                             label = "Treatment of coordinates:",
                             choices = c("Longitude/Latitude", "Cartesian"),
                             selected = "Longitude/Latitude"),
                 conditionalPanel(
                   condition = "input.in_latlon_temp == 'Longitude/Latitude'",
                   textInput(inputId = "to_epsg_temp",
                             label = "EPSG number for the projection to use:",
                             value = "32736")),
                 selectInput(inputId = "data_type_temp",
                             label = "Type of data:",
                             choices = c("Mock up", "Real"),
                             selected = "Mock up"),
                 selectInput(inputId = "get_plots_temp",
                             label = "Store individual timeframe plots (increase execution time):",
                             choices = c("Yes", "No"),
                             selected = "Yes"),
                 # Collapsible section
                 tags$details(
                   tags$summary(
                     tags$span(
                       id = "advancedToggle",
                       style = "cursor: pointer;",
                       HTML("&nbsp;&nbsp;&#43;")
                     ),
                     "Advanced Options"
                   ),
                   tags$div(
                     class = "advanced-options",
                     id = "advancedOptions_temp",
                     selectInput(inputId = "optimize_link_d_temp",
                                 label = "Optimize linking distance: ",
                                 choices = c("Yes", "No"),
                                 selected = "No"),
                     selectInput(inputId = "false_detect_temp",
                                 label = "Consider False Detections (plus exec time): ",
                                 choices = c("Yes", "No"),
                                 selected = "No"),
                     conditionalPanel(
                       condition = "input.false_detect == 'Yes'",
                       numericInput(inputId = "n_simulations_temp",
                                    label = "Number of simulations to perform",
                                    value = 10000, step = 50)
                     ),
                     selectInput(inputId = "method_temp",
                                 label = "Select method to account fot local prevalence:",
                                 choices = c("base", "kmeans", "centroid"),
                                 selected = "base"),
                     conditionalPanel(
                       condition = "input.method_temp == 'centroid'",
                       numericInput(inputId = "max_epi_cont_temp",
                                    label = "Percentage of epifriends contribution to the cluster:",
                                    value = 0.5, min = 0, max = 1, step = 0.1),
                       numericInput(inputId = "max_thr_data_temp",
                                    label = "Percentage of all data for local prevalence calculus:",
                                    value = 0.1, min = 0, max = 1, step = 0.1)
                     ),
                     numericInput(inputId = "max_p_temp",
                                  label = "Maximum p-value:",
                                  value = 1, min = 0, max = 1, step = 0.1),
                     numericInput(inputId = "min_pos_temp",
                                  label = "Minimum number of positives in a cluster:",
                                  value = 2, step = 1),
                     numericInput(inputId = "min_total_temp",
                                  label = "Minimum number of entries in a cluster:",
                                  value = 2, step = 1),
                     numericInput(inputId = "min_pr_temp",
                                  label = "Minimum positivity rate:",
                                  value = 0, step = 0.1)
                   )
                 ),
                 fileInput(inputId = "data_temp",
                           label = "Load data (CSV format):"),
                 actionButton(inputId = "load_data_button_temp", label = "Load Data"),
                 actionButton(inputId = "epifriends_run_temp", label = "Run Temporal Analysis")
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Data",
                            tableOutput(outputId = "load_data_temp")
                   ),
                   tabPanel("Summary-EpiFRIenDs",
                            h3("Summary plots"),
                            br(),
                            plotOutput(outputId = "epifriends_v1_temp"),
                            plotOutput(outputId = "epifriends_v2_temp"),
                            plotOutput(outputId = "epifriends_v3_temp"),
                            
                   ),
                   tabPanel("Animation-EpiFRIenDs",
                            h3("Store animated video of the detected hotspots"),
                            br(),
                            textInput(inputId = "gif_input",
                                      label = "Path to store the Animated GIF file",
                                      value = getwd()),
                            textInput(inputId = "filename",
                                      label = "Name of the file",
                                      value = "temporal_analysis.gif"),
                            actionButton("download_gif", "Download GIF")
                   ),
                 )
               )
             )
    )
  )
)

server <- function(input, output) {
  #output$result1 <- renderText({
  #  paste("You selected", input$option1, "on Page 1")
  #})
  
  algorithm_run <- reactiveVal(FALSE)
  epifriends <- reactiveVal(NULL)
  
  algorithm_run_temp <- reactiveVal(FALSE)
  epifriends_temp <- reactiveVal(NULL)
  
  #########################
  ### LOAD DATA
  ## LOAD SPATIAL DATA
  #########################
  data <- eventReactive(input$load_data_button, {
    req(input$data)
    read.csv(input$data$datapath)
  })
  
  # Check if the column names match the required column names
  observeEvent(data(), {
    required_names <- c("x", "y", "test")
    if (!all(required_names %in% colnames(data()))) {
      showModal(modalDialog(
        title = "Error",
        "The uploaded file does not have the required column names.",
        easyClose = TRUE
      ))
    }
  })
  
  output$load_data <- renderTable({
    data()
  })
  
  #########################
  ## LOAD TEMPORAL DATA
  # Load & show data
  #########################
  data_temp <- eventReactive(input$load_data_button_temp, {
    req(input$data_temp)
    read.csv(input$data_temp$datapath)
  })
  
  # Check if the column names match the required column names
  observeEvent(data_temp(), {
    required_names <- c("x", "y", "test", "date")
    if (!all(required_names %in% colnames(data_temp()))) {
      showModal(modalDialog(
        title = "Error",
        "The uploaded file does not have the required column names.",
        easyClose = TRUE
      ))
    }
  })
  
  output$load_data_temp <- renderTable({
    data_temp()
  })
  
  #########################
  ### RUN ALGO
  ## SPATIAL ALGORITHM
  #########################
  observeEvent(input$epifriends_run, {
    
    # DEFINE PARAMETERS
    df <- as.data.table(data())
    link_d <- input$link_d
    prevalence <- input$prevalence
    if(prevalence == "None"){
      prevalence = NULL
    }else{
      prevalence <- df[[prevalence]]
    }
    cluster_id <- NULL
    min_neighbours <- input$min_neighbours
    link_neighbors <- input$link_neighbors
    if (link_neighbors == 'Yes'){
      use_link_d <- FALSE
    }else{
      use_link_d <- TRUE
    }
    max_p <- input$max_p
    min_pos <- input$min_pos
    min_total <- input$min_total
    min_pr <- input$min_pr
    method <- input$method
    keep_null_tests <- input$keep_null_tests
    if(keep_null_tests == 'Remove'){
      keep_null_tests <- FALSE
    }else if (keep_null_tests == 'Keep'){
      keep_null_tests <- TRUE
    }else{
      keep_null_tests <- input$imp_value
    }
    in_latlon <- input$in_latlon
    if(in_latlon == "Longitude/Latitude"){
      in_latlon <- TRUE
    }else{
      in_latlon <- FALSE
    }
    to_epsg <- as.numeric(input$to_epsg)
    max_epi_cont <- input$max_epi_cont
    max_thr_data <- input$max_thr_data
    consider_fd <- input$false_detect
    if(consider_fd == 'Yes'){
      consider_fd <- TRUE
    }else{
      consider_fd <- FALSE
    }
    n_simulations <- as.numeric(input$n_simulations)
    optimize_link_d <- input$optimize_link_d
    if(optimize_link_d == 'Yes'){
      optimize_link_d <- TRUE
    }else{
      optimize_link_d <- FALSE
    }
    verbose <- FALSE
    
    # RUN CATALOGUE
    epifriends(catalogue(
      x = df$x,
      y = df$y,
      test_result = df$test,
      link_d = link_d,
      prevalence = prevalence,
      cluster_id = cluster_id,
      min_neighbours = min_neighbours,
      use_link_d=use_link_d,
      max_p = max_p,
      min_pos = min_pos,
      min_total = min_total,
      min_pr = min_pr,
      method = method,
      keep_null_tests = keep_null_tests,
      in_latlon = in_latlon,
      to_epsg = to_epsg,
      max_epi_cont = max_epi_cont,
      max_thr_data = max_thr_data,
      consider_fd = consider_fd,
      n_simulations = n_simulations,
      optimize_link_d = optimize_link_d,
      verbose = verbose))
    algorithm_run(TRUE)
  })
  
  observeEvent(input$epifriends_run, {
    algorithm_run(TRUE)
  })
  
  #########################
  ## RUN ALGO
  ## TEMPORAL ALGORITHM
  #########################
  observeEvent(input$epifriends_run_temp, {
    
    # DEFINE PARAMETERS
    df <- as.data.table(data_temp())
    link_d <- input$link_d_temp
    prevalence <- input$prevalence_temp
    if(prevalence == "None"){
      prevalence = NULL
    }else{
      prevalence <- df[[prevalence]]
    }
    min_neighbours <- input$min_neighbours_temp
    time_width <- input$time_width_temp
    link_neighbors <- input$link_neighbors
    if (link_neighbors == 'Yes'){
      use_link_d <- FALSE
    }else{
      use_link_d <- TRUE
    }
    max_p <- input$max_p_temp
    min_pos <- input$min_pos_temp
    min_total <- input$min_total_temp
    min_pr <- input$min_pr_temp
    linking_time = input$linking_time_temp
    linking_dist = input$linking_dist_temp
    time_steps <- input$time_steps_temp
    method <- input$method_temp
    keep_null_tests <- input$keep_null_tests_temp
    if(keep_null_tests == 'Remove'){
      keep_null_tests <- FALSE
    }else if (keep_null_tests == 'Keep'){
      keep_null_tests <- TRUE
    }else{
      keep_null_tests <- input$imp_value_temp
    }
    in_latlon <- input$in_latlon_temp
    if(in_latlon == "Longitude/Latitude"){
      in_latlon <- TRUE
    }else{
      in_latlon <- FALSE
    }
    to_epsg <- as.numeric(input$to_epsg_temp)
    max_epi_cont <- input$max_epi_cont_temp
    max_thr_data <- input$max_thr_data_temp
    consider_fd <- input$false_detect_temp
    if(consider_fd == 'Yes'){
      consider_fd <- TRUE
    }else{
      consider_fd <- FALSE
    }
    n_simulations <- as.numeric(input$n_simulations_temp)
    optimize_link_d <- input$optimize_link_d_temp
    if(optimize_link_d == 'Yes'){
      optimize_link_d <- TRUE
    }else{
      optimize_link_d <- FALSE
    }
    verbose <- FALSE
    get_plot <- input$get_plots_temp
    if(get_plot == 'Yes'){
      get_plot = TRUE
    }else{
      get_plot <- FALSE
    }
    
    # Add hour, minute & second to date column so that it can be converted to chron
    data_type <- input$data_type_temp
    if(data_type == 'Real'){
      use_geom_map <- TRUE
    }else{
      use_geom_map <- FALSE
    }
    
    # Process dates
    parsed_dates <- strptime(df$date, format = c("%y:%m:%d %H:%M:%S", "%y:%m:%d"))
    is_full_datetime <- !is.na(parsed_dates)
    if (all(is_full_datetime)) {
      format_type <- "YY:MM:DD HH:MM:SS"
    } else if (any(is_full_datetime)) {
      format_type <- "Mixed formats"
    } else {
      format_type <- "YY:MM:DD alone"
    }
    
    if(format_type == 'YY:MM:DD alone'){
      df$date <- paste0(df$date, " 00:00:00")
    }
    
    #Get temporal IDs
    epifriends_temp(temporal_catalogue(
      x = df$x,
      y = df$y, 
      test_result = df$test,
      dates = df$date,
      link_d = link_d,
      prevalence = prevalence,
      use_link_d=use_link_d,
      min_neighbours =  min_neighbours,
      time_width = time_width,
      min_date = min(df$date),
      max_date = max(df$date),
      time_steps = time_steps,
      max_p = max_p,
      min_pos = min_pos,
      min_total = min_total,
      min_pr = min_pr,
      add_temporal_id = TRUE,
      linking_time = linking_time,
      linking_dist = linking_dist,
      get_timelife = TRUE,
      optimize_link_d = optimize_link_d,
      method = method,
      keep_null_tests = keep_null_tests, 
      in_latlon = in_latlon,
      to_epsg = to_epsg, 
      consider_fd = consider_fd,
      n_simulations = n_simulations,
      verbose = FALSE,
      store_gif = get_plot,
      use_geom_map = use_geom_map))
    algorithm_run_temp(TRUE)
  })
  
  observeEvent(input$epifriends_run_temp, {
    algorithm_run_temp(TRUE)
  })
  
  
  #########################
  ### TABLES RENDERING
  ## SPATIAL ALGO TABLES
  #########################
  output$epifriends_summary <- renderTable({
    if(algorithm_run()){
      df <- as.data.table(data())
      epi <- epifriends()
      general <- data.table("cluster_id" = 1:length(epi$epifriends_catalogue$p))
      general$pvalue <- epi$epifriends_catalogue$p
      general <- cbind(general,rbindlist(epi$epifriends_catalogue$mean_position_all))
      general <- general[,.(cluster_id = cluster_id, p_value = pvalue, mean_x = x, mean_y = y)]
      general
    }      
  }, 
  table.attr = "class='centered'"
  )
  
  output$epifriends_results <- renderTable({
    if(algorithm_run()){
      df <- as.data.table(data())
      epi <- epifriends()
      # Plot significant clusters
      coords <- data.table::copy(df)
      coords[, cluster := 0]
      coords[, index := 1:nrow(coords)]
      for(clusters in which(epi$epifriends_catalogue$p <= 0.05)){
        coords[index %in% epi$epifriends_catalogue$indeces[[clusters]], cluster := clusters]
      }
      coords <- coords[cluster != 0]
      coords <- coords[order(cluster)]
      coords[, ':='(X = NULL, id = NULL, index = NULL)]
      
    }      
  }, 
  table.attr = "class='centered'"
  )
  
  #########################
  ### PLOTS RENDERING
  ## SPATIAL ALGO PLOTS
  #########################
  output$distribution <- renderPlot({
    df <- as.data.table(data())
    keep_null_tests <- input$keep_null_tests
    if(keep_null_tests == 'Remove'){
      keep_null_tests <- FALSE
    }else if (keep_null_tests == 'Keep'){
      keep_null_tests <- TRUE
    }else{
      keep_null_tests <- input$imp_value_temp
    }
    use_geom_map <- input$data_type
    if(use_geom_map == 'Mock up'){
      use_geom_map <- FALSE
    }else{
      use_geom_map <- TRUE
    }
    
    pos = clean_unknown_data(df, cols_impute = c("test"), keep_null_tests,FALSE)
    positions = data.table(x = pos$x, y = pos$y, test = pos$test)
    
    if(use_geom_map){
      my_location <- c(min(positions$x),min(positions$y), max(positions$x), max(positions$y))
      my_map <- get_map(location = my_location, source = "stamen", maptype = "terrain")
      map_plot <- ggmap::ggmap(my_map)
      graph <- map_plot +
        geom_point(data = positions, aes(x = x, y = y, color = as.factor(test))) +
        geom_point(size = 2.5) +
        labs(title = "Distribution of Positive and Negative Cases") +
        coord_equal()
    }else{
      graph <- ggplot(data(), aes(x = x, y = y, color = as.factor(test))) +
        geom_point(size = 2.5) +
        labs(title = "Distribution of Positive and Negative Cases") +
        coord_equal()
      
    }
    
    graph
    
  })
  
  output$epifriends_v1 <- renderPlot({
    if(algorithm_run()){
      df <- as.data.table(data())
      epi <- epifriends()
      use_geom_map <- input$data_type
      if(use_geom_map == 'Mock up'){
        use_geom_map <- FALSE
      }else{
        use_geom_map <- TRUE
      }
      max_epi_cont <- input$max_epi_cont
      max_thr_data <- input$max_thr_data
      
      keep_null_tests <- input$keep_null_tests
      if(keep_null_tests == 'Remove'){
        keep_null_tests <- FALSE
      }else if (keep_null_tests == 'Keep'){
        keep_null_tests <- TRUE
      }else{
        keep_null_tests <- input$imp_value_temp
      }
      
      if (input$method == 'kmeans'){
        # Plot of KMeans-Identified Clusters
        kmeans_prev <- compute_kmeans(
          clean_data(df[,.(x, y, test)])[,.(x,y)], 
          clean_data(df[,.(x, y, test)])$test)
        prevalence <- kmeans_prev$prevalence
      }else if (input$method == 'centroid'){
        copy_position <- data.table::copy(df[,.(x,y)])	
        copy_position[, id := 1:nrow(copy_position)]	
        total_friends_indeces <- epi$epifriends_catalogue$indeces	
        copy_position[, prevalence := 0]	
        radial <- c()	
        centroid_df <- list()	
        for(i in 1:length(total_friends_indeces)){
          result <- compute_centroid(	
            positions = copy_position,
            test_result = data.table("test" = df$test),	
            total_friends_indeces = total_friends_indeces[[i]],	
            thr_data = max_thr_data,
            max_epi_cont = max_epi_cont)	
          
          # Assign computed prevalence	
          rows <- dim(copy_position[id %in% result$local_id])[1]	
          copy_position[id %in% result$local_id,
                        prevalence := rep(result$prevalence, rows)]	
        }
        prevalence <- copy_position$prevalence
      }else{
        prevalence <- rep(sum(df$test) / length(df$test), length(df$test))
      }
      
      pos = clean_unknown_data(df, cols_impute = c("test"), keep_null_tests,FALSE)
      positions = data.table(x = pos$x, y = pos$y, test = pos$test)
      
      plot1 <- scatter_pval(
        coordinates = positions[,.(x,y)], 
        id_data = epi$cluster_id, 
        positive = (positions$test == 1), 
        prevalence = prevalence, 
        epi_catalogue = epi$epifriends_catalogue,
        use_geom_map = use_geom_map,
        xlims = NULL,
        ylims = NULL,
        title = NULL)
      grid.arrange(plot1, ncol=1)
    }      
  })
  
  output$epifriends_v2 <- renderPlot({
    if(algorithm_run()){
      df <- as.data.table(data())
      epi <- epifriends()
      plot2 <- size_histogram(epi)
      grid.arrange(plot2, ncol=1)
    }      
  })
  
  output$epifriends_v3 <- renderPlot({
    if(algorithm_run()){
      df <- as.data.table(data())
      epi <- epifriends()
      use_geom_map <- input$data_type
      if(use_geom_map == 'Mock up'){
        use_geom_map <- FALSE
      }else{
        use_geom_map <- TRUE
      }
      keep_null_tests <- input$keep_null_tests
      if(keep_null_tests == 'Remove'){
        keep_null_tests <- FALSE
      }else if (keep_null_tests == 'Keep'){
        keep_null_tests <- TRUE
      }else{
        keep_null_tests <- input$imp_value_temp
      }
      # Plot significant clusters
      coords = clean_unknown_data(df, cols_impute = c("test"), keep_null_tests,FALSE)
      coords = data.table(x = coords$x, y = coords$y, test = coords$test)
      coords[, cluster := 0]
      coords[, index := 1:nrow(coords)]
      for(clusters in which(epi$epifriends_catalogue$p <= 0.05)){
        coords[index %in% epi$epifriends_catalogue$indeces[[clusters]], cluster := clusters]
      }
      
      if(use_geom_map){
        my_location <- c(min(coords$x),min(coords$y), max(coords$x), max(coords$y))
        my_map <- get_map(location = my_location, source = "stamen", maptype = "terrain")
        map_plot <- ggmap::ggmap(my_map)
        plot3 <- map_plot +
          geom_point(data = coords[cluster != 0], aes(x = x, y = y, color = as.factor(cluster))) +
          geom_point(size = 2.5)  + geom_point(data = coords[cluster == 0], aes(x=x, y=y, color = "#FFFFFF"), 
                                               shape=21, stroke = 1) +
          labs(title = "Distribution of Significant Clusters") + coord_equal()
      }else{
        plot3 <- ggplot(coords[cluster != 0], aes(x = x, y = y, color = as.factor(cluster))) +
          geom_point(size = 2.5)  + geom_point(data = coords[cluster == 0], aes(x=x, y=y, color = "#FFFFFF"), shape=21, stroke = 1) +
          labs(title = "Distribution of Significant Clusters") + coord_equal()
        
      }
      
      grid.arrange(plot3, ncol=1)
    }      
  })
  
  #########################
  ### PLOTS RENDERING
  ## TEMPORAL ALGO PLOTS
  #########################
  output$epifriends_v1_temp <- renderPlot({
    if(algorithm_run_temp()){
      epi_catalogue_list <- epifriends_temp()
      num_clusters <- c()
      for(i in epi_catalogue_list$temporal_catalogues){
        num_clusters <- append(num_clusters,length(i$id))
      }
      
      df_all = data.table("mean_date" = epi_catalogue_list$mean_date,
                          "num_clusters" = num_clusters)
      graph <- ggplot(df_all, aes(x = as.Date(mean_date), y = num_clusters)) +
        geom_line() + labs(x = "Dates", y = "Number of cluster",
                           title = "Number of clusters across dates")
      grid.arrange(graph, ncol=1)
    }
  })
  
  output$epifriends_v2_temp <- renderPlot({
    if(algorithm_run_temp()){
      epi_catalogue_list <- epifriends_temp()
      graph <- hist_timelifes(epi_catalogue_list$temporal_catalogues)
      grid.arrange(graph, ncol=1)
    }
  })
  
  output$epifriends_v3_temp <- renderPlot({
    if(algorithm_run_temp()){
      epi_catalogue_list <- epifriends_temp()
      graph <- lifetime_timeline(epi_catalogue_list$temporal_catalogues, epi_catalogue_list$mean_date, input$time_steps_temp)
      grid.arrange(graph, ncol=1)
    }
  })
  
  file_path <- reactive({
    # Your process to determine the file path of the file to be downloaded
    # For example:
    file_path <- paste0(getwd(),"/www/animated_gif.gif")
    return(file_path)
  })
  
  observeEvent(input$download_gif, {
    folder_path <- paste0(input$gif_input, "/")
    
    # Use the file.copy() function to copy the file to the user's desired folder
    file.copy(from = file_path(), to = folder_path)
    file.rename(
      from = paste0(folder_path, basename(file_path())),
      to = paste0(folder_path, basename(input$filename)))
    
    showNotification("File has been successfully stored!", type = "message")
    
  })
}

shinyApp(ui = ui, server = server)
