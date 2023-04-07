
library(data.table)
library(readr)
library(factoextra)
library(NbClust)
library(doParallel)
library(reshape2)
library(pdftools)
library(clValid)
library(epifriends)
library(ggplot2)
library(ggforce)
library(useful)
library(RANN)
library(patchwork)
library(knitr)
library(kableExtra)
library(grid)

setwd("/Users/ericmatamoros/Desktop/TFM/repifriends/")
source("./utils.R")
# Load data
# Prevalence defined as number of positives versus total number, where total number
# can either be a cluster, a radious or all the population size.
files <- c("example_multiple.csv", "example_one_isolated_one_populated_v2.csv", "example_one_isolated_one_populated.csv")

for(file in files){
  data_v1 <- as.data.table(read_csv(paste0("./data/", file)))
  data_v1[, id := 1:nrow(data_v1)]
  
  position_rand <- data_v1[,3:4]
  positive_rand = (data_v1$test == 1)
  test_rand <- data_v1$test
  
  # Plot data
  graph_base <- ggplot(data_v1, aes(x = x, y = y, color = as.factor(test))) +
    geom_point(size = 2.5) +
    labs(
      title = "Distribution of Positive and Negative Cases",
      subtitle = file
      )
  
  print(graph_base)
  
}

link_d_to_analyse = c(0.1)
min_neighbours = 2
methods_list <- c("kmeans", "radial", "centroid", "base")
store_PDF = TRUE

for(file in files[1]){
  
  catalogue_methods <- c()
  for(method in methods_list){
    if (store_PDF){
      pdf_name <- paste0("pdfs/", strsplit(file, ".csv")[[1]],"_",method,".pdf")
      pdf(pdf_name, width = 11, height = 8.5)
    }
    
    for(link_d in link_d_to_analyse){
      
      #########################
      ####### READ DATA #######
      #########################
      data_v1 <- as.data.table(read_csv(paste0("data/", file)))
      data_v1[, id := 1:nrow(data_v1)]
      
      position_rand <- data_v1[,3:4]
      positive_rand = (data_v1$test == 1)
      test_rand <- data_v1$test
      
      ####################################################################
      ####### GENERATE PLOTS OF DISTRIBUTION AND DETECTED CLUSTERS #######
      ####################################################################
      # Plot data
      graph_base <- ggplot(data_v1, aes(x = x, y = y, color = as.factor(test))) +
        geom_point(size = 2.5) +
        labs(title = "Distribution of Positive and Negative Cases")
      
      #print(graph_base)
      
      # Compute Epifriends
      categories <- catalogue(
        position_rand, test_rand, link_d, cluster_id = NULL, min_neighbours = min_neighbours,
        method = method)
      catalogue_methods[[method]] <-  categories
      
      # Plot significant clusters
      coords <- data.table::copy(position_rand)
      coords[, cluster := 0]
      coords[, index := 1:nrow(coords)]
      for(clusters in which(categories$epifriends_catalogue$p <= 0.05)){
        coords[index %in% categories$epifriends_catalogue$indeces[[clusters]], cluster := clusters]
      }
      
      graph_clusters <- ggplot(coords[cluster != 0], aes(x = x, y = y, color = as.factor(cluster))) +
        geom_point(size = 2.5)  + geom_point(data = coords[cluster == 0], aes(x=x, y=y, color = "#FFFFFF"), shape=21, stroke = 1) +
        labs(title = "Distribution of Significant Clusters")
      
      if(method != 'radial'){
        graph_clusters <- graph_clusters + geom_encircle(
          data = get_region_prev(data_v1, categories$epifriends_catalogue), 
          aes(x=x, y=y, group = epifriends, color = as.factor(epifriends)), expand = 0.1 , size = 2)
      }
      
      ###################################
      ####### TABLE WITH INSIGHTS #######
      ###################################
      # Table general overview
      general <- data.table("cluster_id" = 1:length(categories$epifriends_catalogue$p))
      general$pvalue <- categories$epifriends_catalogue$p
      general$mean_pr <- categories$epifriends_catalogue$mean_pr
      general$mean_pr <- categories$epifriends_catalogue$mean_local_prev
      general$n_positives <- categories$epifriends_catalogue$positives
      general$n_negatives <- categories$epifriends_catalogue$negatives
      general$n_total <- categories$epifriends_catalogue$total
      
      general_two = data.table("cluster_id" = 1:length(categories$epifriends_catalogue$p))
      general_two <- cbind(general_two,rbindlist(categories$epifriends_catalogue$mean_position_all))
      
      #grid::grid.newpage()
      grid.arrange(tableGrob(general), tableGrob(general_two), ncol = 1)
      print(graph_base + graph_clusters)

      if(method == 'kmeans'){
        kmeans_prev <- compute_kmeans(
          clean_data(data_v1[,.(x, y, test)])[,.(x,y)], 
          clean_data(data_v1[,.(x, y, test)])$test)
        prevalence <- kmeans_prev$prevalence
        
        k_results <- ggplot(kmeans_prev, aes(x = x, y = y, color = as.factor(clusters))) + 
          geom_point(size = 2.5) + labs(title = "Detected Clusters using KMeans")

        print(graph_clusters + k_results)
        
      }else if(method == "radial"){
        index_friends <- unique(unlist(categories$epifriends_catalogue$indeces))
        radial_prev <- prevalence_radial(data_v1, position_rand)
        prevalence <- radial_prev$prevalence$prevalence
      }else if (method == "base"){
        prevalence <- rep(sum(test_rand) / length(test_rand), length(test_rand))
      }else if (method == "centroid"){
        copy_position <- data.table::copy(position_rand)
        copy_position[, id := 1:nrow(copy_position)]
        total_friends_indeces <- categories$epifriends_catalogue$indeces
        copy_position[, prevalence := 0]
        copy_position[, radious := 0]
        radial <- c()
        centroid_df <- list()
        for(i in 1:length(total_friends_indeces)){
          result <- get_radious(
            positions = copy_position,
            test_result = data.table("test" = test_rand),
            total_friends_indeces = total_friends_indeces[[i]],
            thr_data = 0.1, 
            max_epi_cont = 0.5)
          
          # Assign computed prevalence
          rows <- dim(copy_position[id %in% result$id])[1]
          copy_position[id %in% result$id,
                        prevalence := rep(result$prevalence, rows)]
          
          prevalence <- copy_position$prevalence
          centroid_df[[i]] <- result$centroid
          radial[i] <- result$radious
        }
        
        centroid_df <- rbindlist(centroid_df)
        centroid_df[, radious := radial]
            
      }else{
        prevalence <- NULL
      }
      
      graphs_prev <- scatter_pval(
          data.table::copy(position_rand), 
          categories$cluster_id, 
          data.table::copy(positive_rand), 
          data.table::copy(prevalence),
          categories$epifriends_catalogue,
          paste0(method,"-Link_d: ",link_d)) 
      
      graphs_no_prev <- scatter_pval(
        position_rand, 
        categories$cluster_id, 
        positive_rand, 
        NULL,
        categories$epifriends_catalogue,
        paste0(method,"-Link_d: ",link_d)) 
      histo <- size_histogram(categories)
      
      print(graph_clusters+graphs_prev)
      #print(graph_clusters + graphs_no_prev)
      print(histo)
      
    }
    
    if (store_PDF){
      dev.off()
    }
  }
  
  pdf_name <- paste0("pdfs/", strsplit(file, ".csv")[[1]],"_GENERAL_INSIGHTS.pdf")
  pdf(pdf_name, width = 11, height = 8.5)
  
  general <- get_all_results(catalogue_methods, methods_list)
  print(general$chart)
  if (store_PDF){
    dev.off()
  }
}
