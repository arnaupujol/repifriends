library(shiny)
library(epifriends)
library(data.table)
library(gridExtra)
library(RANN)
library(ggplot2)
library(chron)
library(shinyjs)

ui <- fluidPage(
  useShinyjs(),
  titlePanel("EpiFRIenDs"),
  tabsetPanel(
    tabPanel("Distribution Analysis",
             sidebarLayout(
               sidebarPanel(
                 
               numericInput(inputId = "min_neighbours",
                            label = "Min. number of neighbors:",
                            value = 2, min = 0, max = 10),
               sliderInput(inputId = "link_d",
                           label = "Select linking distance:",
                           min = 0.05, max = 0.2, value = 0.1, step = 0.05),
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
               fileInput(inputId = "data",
                         label = "Load data (CSV format):"),
               actionButton(inputId = "load_data_button", label = "Load Data"),
               actionButton(inputId = "epifriends_run", label = "Run Spatial Analysis")
               ),
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
             sidebarLayout(
               sidebarPanel(
                 
                 numericInput(inputId = "min_neighbours_temp",
                              label = "Min. number of neighbors:",
                              value = 2, min = 0, max = 10),
                 sliderInput(inputId = "link_d_temp",
                             label = "Select linking distance:",
                             min = 0.05, max = 0.2, value = 0.1, step = 0.05),
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
                   numericInput(inputId = "imp_value",
                                label = "Imputation value:",
                                value = 0)
                 ),
                 selectInput(inputId = "in_latlon_temp",
                             label = "Treatment of coordinates:",
                             choices = c("Longitude/Latitude", "Cartesian"),
                             selected = "Longitude/Latitude"),
                 conditionalPanel(
                   condition = "input.in_latlon_temp == 'Longitude/Latitude'",
                   textInput(inputId = "to_epsg",
                             label = "EPSG number for the projection to use:",
                             value = "32736")),
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
                   tabPanel("Distribution",
                            plotOutput(outputId = "distribution_temp")
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
  
  ### LOAD DATA
  ## LOAD SPATIAL DATA
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
  
  ## LOAD TEMPORAL DATA
  # Load & show data
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
  
  ### RUN ALGO
  ## SPATIAL ALGORITHM
  observeEvent(input$epifriends_run, {
    df <- as.data.table(data())
    keep_null_tests <- input$keep_null_tests
    in_latlon <- input$in_latlon
    if(keep_null_tests == 'Remove'){
      keep_null_tests <- FALSE
    }else if (keep_null_tests == 'Keep'){
      keep_null_tests <- TRUE
    }else{
      keep_null_tests <- input$imp_value
    }
    
    if(in_latlon == "Longitude/Latitude"){
      in_latlon <- FALSE
    }
    epifriends(catalogue(
      positions = df[,.(x,y)], 
      test_result = df[,.(test)], 
      link_d = input$link_d,
      min_neighbours = input$min_neighbours, 
      max_p = 1, 
      min_pos = 2, 
      min_total = 2,
      min_pr = 0, 
      keep_null_tests = keep_null_tests, 
      in_latlon = in_latlon,
      to_epsg = as.numeric(input$to_epsg), 
      verbose = FALSE))
    algorithm_run(TRUE)
  })
  
  observeEvent(input$epifriends_run, {
    algorithm_run(TRUE)
  })
  
  ## TEMPORAL ALGORITHM
  observeEvent(input$epifriends_run_temp, {
    df <- as.data.table(data_temp())
    keep_null_tests <- input$keep_null_tests_temp
    in_latlon <- input$in_latlon_temp
    if(keep_null_tests == 'Remove'){
      keep_null_tests <- FALSE
    }else if (keep_null_tests == 'Keep'){
      keep_null_tests <- TRUE
    }else{
      keep_null_tests <- input$imp_value_temp
    }
    
    if(in_latlon == "Longitude/Latitude"){
      in_latlon <- FALSE
    }
    
    dtparts = t(as.data.frame(strsplit(df$date," ")))
    row.names(dtparts) = NULL
    datesform <- chron(dates=dtparts[,1],times=dtparts[,2],format=c('y-m-d','h:m:s'))
    dates <- as.numeric(datesform)
    
    #Get temporal IDs
    epifriends_temp(temporal_catalogue(
      positions = df[,.(x,y)], 
      test_result = df[,.(test)],
      dates = df$date,
      link_d = input$link_d_temp,
      min_neighbours = input$min_neighbours_temp,
      time_width = input$time_width_temp,
      min_date = min(datesform),
      max_date = max(datesform),
      time_steps = input$time_steps_temp,
      add_temporal_id = TRUE,
      linking_time = input$linking_time_temp,
      linking_dist = input$linking_dist_temp,
      get_timelife = TRUE,
      cluster_id = NULL,
      keep_null_tests = keep_null_tests, 
      in_latlon = in_latlon,
      to_epsg = as.numeric(input$to_epsg), 
      verbose = FALSE,
      store_gif = TRUE))
    algorithm_run_temp(TRUE)
  })
  
  observeEvent(input$epifriends_run_temp, {
    algorithm_run_temp(TRUE)
  })
  
  
  ### TABLES RENDERING
  ## SPATIAL ALGO TABLES
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
  
  ### PLOTS RENDERING
  ## SPATIAL ALGO PLOTS
  output$distribution <- renderPlot({
    ggplot(data(), aes(x = x, y = y, color = as.factor(test))) +
      geom_point(size = 2.5) +
      labs(title = "Distribution of Positive and Negative Cases") +
      coord_equal()
  })
  
  output$epifriends_v1 <- renderPlot({
    if(algorithm_run()){
      df <- as.data.table(data())
      epi <- epifriends()
      plot1 <- scatter_pval(df[,.(x,y)], epi$cluster_id, (df$test == 1), epi$epifriends_catalogue)
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
      # Plot significant clusters
      coords <- data.table::copy(df)
      coords[, cluster := 0]
      coords[, index := 1:nrow(coords)]
      for(clusters in which(epi$epifriends_catalogue$p <= 0.05)){
        coords[index %in% epi$epifriends_catalogue$indeces[[clusters]], cluster := clusters]
      }
      
      plot3 <- ggplot(coords[cluster != 0], aes(x = x, y = y, color = as.factor(cluster))) +
        geom_point(size = 2.5)  + geom_point(data = coords[cluster == 0], aes(x=x, y=y, color = "#FFFFFF"), shape=21, stroke = 1) +
        labs(title = "Distribution of Significant Clusters") + coord_equal()
      grid.arrange(plot3, ncol=1)
    }      
  })
  
  ## TEMPORAL ALGO PLOTS
  output$distribution_temp <- renderPlot({
    df <- data_temp()
    dtparts = t(as.data.frame(strsplit(df$date," ")))
    row.names(dtparts) = NULL
    datesform <- chron(dates=dtparts[,1],times=dtparts[,2],format=c('y-m-d','h:m:s'))
    dates <- as.numeric(datesform)
    graph <- ggplot(df,aes(x=x, y=y))+
        geom_point(aes(fill = dates, colour = as.factor(test)), shape = 21, size = 2.5)+
        scale_fill_gradientn(colors = c("midnight blue", "turquoise", "yellow"))+
        scale_color_manual(values = c("grey", "red"), name = "Test", labels = c("Negative","Positive")) + coord_equal()
    grid.arrange(graph, ncol=1)
  })
  
  output$epifriends_v1_temp <- renderPlot({
    if(algorithm_run_temp()){
      epi_catalogue_list <- epifriends_temp()
      num_clusters <- c()
      for(i in epi_catalogue_list$temporal_catalogues){
        num_clusters <- append(num_clusters,length(i$id))
      }
      
      df_all = data.table("mean_date" = epi_catalogue_list$mean_date,
                          "num_clusters" = num_clusters)
      graph <- ggplot(df_all, aes(x = mean_date, y = num_clusters)) +
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

