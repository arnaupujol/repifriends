library(shiny)
library(epifriends)
library(data.table)
library(gridExtra)

scatter_pval <- function(coordinates, id_data, positive, epi_catalogue){
  # This method shows a scatter plot of the data distribution, showing in 
  # colour the positive cases that belong to foci (with the colour 
  # representing their p-value) and in grey the rest of the data. 
  # 
  # Parameters:
  # -----------
  # coordinates: data.frame
  #     data frame with the values of the coordinates
  # id_data: data.frame
  #     data frame with the cluster ID associated to each positive case, o for 
  #     no associated cluster
  # positive: data.frame
  #     Boolean vector that indicates if the case is infected. 
  # epi_catalogue: list
  #     List of the EpiFRIenDs catalogue
  # 
  # Returns:
  # --------
  # Scatter plot of the data distribution, showing in 
  # colour the positive cases that belong to foci (with the colour 
  # representing their p-value) and in grey the rest of the data.
  pos <- data.frame(coordinates[positive,],id_data)
  p_vals <- c()
  for(i in id_data[id_data > 0]){
    p_vals <- append(p_vals, epi_catalogue$p[epi_catalogue$id == i])
  }
  # plot(coordinates$x, coordinates$y, pch = 19, col = "grey")
  # points(pos$x[id_data > 0], pos$y[id_data>0], pch = 19, col = rainbow(100)[factor(p_vals)])
  graph <- ggplot(coordinates,aes(x=x, y=y))+
    geom_point(color = "grey", size = 2.5)+
    geom_point(data = pos[id_data >0,], aes(colour = p_vals), size = 2.5)+
    scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07"))+
    ggtitle("P-value of hotspots")
  return(graph)
}

size_histogram <- function(catalogue){
  # This method shows a histogram of the size (total number of cases) in foci 
  # for foci with p>0.05 and with p<0.05. 
  #   
  # Parameters:
  # -----------
  # catalogue: list
  #   List of the EpiFRIenDs catalogue
  #   
  # Returns:
  # --------
  # Histogram of number of foci per total number of cases with p<0.05 in red 
  # and p>0.05 in blue.
  cases <- data.frame(catalogue$epifriends_catalogue$total, catalogue$epifriends_catalogue$p)
  colnames(cases) <- c("total","p")
  #auxiliar vector to avoid converting doubles to strings
  aux <- cases$p
  for(i in 1:length(cases$p)){
    if(aux[i] > 0.05){
      cases$p[i] <- "p > 0.05"
    }else{
      cases$p[i] <- "p < 0.05"
    }
  }
  histo <- ggplot(cases,aes(x=total, fill=p))+
    geom_histogram()+
    scale_fill_manual(values = c("red", "blue"))+
    xlab("Number of cases in focci")+
    ylab("Number of focci")
  return(histo)
}

ui <- fluidPage(
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
               selectInput(inputId = "method",
                         label = "Methodology for local prevalence calculus:",
                         choices = c("KMeans", "Radial", "Centroid"),
                         selected = "Centroid"),
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
               actionButton(inputId = "epifriends_run", label = "Run Algorithm")
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Data",
                            tableOutput(outputId = "load_data")
                   ),
                   tabPanel("Distribution",
                            plotOutput(outputId = "distribution")
                   ),
                   tabPanel("Summary EpiFRIenDs",
                            plotOutput(outputId = "epifriends"))
                 )
               )
             )
    ),
             
    tabPanel("Temporal Analysis", 
             selectInput(inputId = "option2",
                         label = "Select an option for page 2:",
                         choices = c("Option X", "Option Y", "Option Z"),
                         selected = "Option X"),
             textOutput(outputId = "result2"))
  )
)

server <- function(input, output) {
  #output$result1 <- renderText({
  #  paste("You selected", input$option1, "on Page 1")
  #})
  
  algorithm_run <- reactiveVal(FALSE)
  epifriends <- reactiveVal(NULL)
  
  # Load & show data
  data <- eventReactive(input$load_data_button, {
    req(input$data)
    read.csv(input$data$datapath)
  })
  
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
      to_epsg = input$to_epsg, 
      verbose = FALSE))
    algorithm_run(TRUE)
  })
  
  observeEvent(input$epifriends_run, {
    algorithm_run(TRUE)
  })
  
  output$load_data <- renderTable({
    data()
  })
  
  output$distribution <- renderPlot({
    ggplot(data(), aes(x = x, y = y, color = as.factor(test))) +
      geom_point(size = 2.5) +
      labs(title = "Distribution of Positive and Negative Cases")
  })
  
  output$epifriends <- renderPlot({
    if(algorithm_run()){
      df <- as.data.table(data())
      epi <- epifriends()
      plot1 <- scatter_pval(df[,.(x,y)], epi$cluster_id, (df$test == 1), epi$epifriends_catalogue)
      
      plot2 <- histo_rand <- size_histogram(epi)
      
      # Arrange the plots in a 2x1 grid layout
      grid.arrange(plot1, plot2, ncol = 1)
    }      
  })

}

shinyApp(ui = ui, server = server)
