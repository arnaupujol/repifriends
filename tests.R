
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

# Load data
# Prevalence defined as number of positives versus total number, where total number
# can either be a cluster, a radious or all the population size.

files <- c("example_multiple.csv", "example_one_isolated_one_populated_v2.csv", "example_one_isolated_one_populated.csv")

for(file in files){
  data_v1 <- as.data.table(read_csv(paste0("data/", file)))
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

link_d_to_analyse = c(0.05)
min_neighbours = 2
methods_list <- c("kmeans", "radial", "centroid", "base")
store_PDF = TRUE

for(file in files[2]){
  
  catalogue_methods <- c()
  for(method in methods_list){
    if (store_PDF){
      pdf_name <- paste0("pdfs/", strsplit(file, ".csv")[[1]],"_",method,".pdf")
      pdf(pdf_name, width = 11, height = 8.5)
    }
    
    for(link_d in link_d_to_analyse){
      data_v1 <- as.data.table(read_csv(paste0("data/", file)))
      data_v1[, id := 1:nrow(data_v1)]
      
      position_rand <- data_v1[,3:4]
      positive_rand = (data_v1$test == 1)
      test_rand <- data_v1$test
      
      # Plot data
      graph_base <- ggplot(data_v1, aes(x = x, y = y, color = as.factor(test))) +
        geom_point(size = 2.5) +
        labs(title = "Distribution of Positive and Negative Cases")
      
      print(graph_base)
      
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
      
      # Table general overview
      general <- data.table("cluster_id" = 1:length(categories$epifriends_catalogue$p))
      general$pvalue <- categories$epifriends_catalogue$p
      general <- cbind(general,rbindlist(categories$epifriends_catalogue$mean_position_all))
      
      grid::grid.newpage()
      grid.table(general)
      print(graph_base + graph_clusters)
      
      
      
      if(method == 'kmeans'){
        kmeans_prev <- compute_kmeans(position_rand, test_rand)
        prevalence <- kmeans_prev$prevalence
        
        k_results <- kmeans(position_rand[,.(x,y)], centers = max(kmeans_prev$clusters), nstart = 500)
        
        print(graph_clusters + plot(k_results, data = position_rand[,.(x,y)]))
        
      }else if(method == "radial"){
        index_friends <- unique(unlist(categories$epifriends_catalogue$indeces))
        radial_prev <- prevalence_radial(data_v1, position_rand)
        prevalence <- radial_prev$prevalence$prevalence
        
        radial_prev$distances
        centroid_df <- merge(data_v1[id %in% index_friends, .(id, x, y)],
              radial_prev$distances,
              by.x = "id",
              by.y = "coord1", how = "left")
        names(centroid_df) <- c("id", "x", "y", "radious")
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
      
      if(method != "centroid"){
        graphs_prev <- scatter_pval(
          data.table::copy(position_rand), 
          categories$cluster_id, 
          data.table::copy(positive_rand), 
          data.table::copy(prevalence),
          categories$epifriends_catalogue,
          NULL,
          paste0(method,"-Link_d: ",link_d)) 
      }else{
        graphs_prev <- scatter_pval(
          data.table::copy(position_rand), 
          categories$cluster_id, 
          data.table::copy(positive_rand), 
          prevalence,
          categories$epifriends_catalogue,
          data.table::copy(centroid_df),
          paste0(method,"-Link_d: ",link_d)) 
        
      }
      
      graphs_no_prev <- scatter_pval(
        position_rand, 
        categories$cluster_id, 
        positive_rand, 
        NULL,
        categories$epifriends_catalogue,
        paste0(method,"-Link_d: ",link_d)) 
      histo <- size_histogram(categories)
      
      print(graph_clusters+graphs_prev)
      print(graph_clusters + graphs_no_prev)
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

##############################
####  COMPUTE PREVALENCES ####
##############################
## RADIAL APPROACH
# Create a table with coord-1, coord-2 and distance
prevalence_radial <- function(df, positions){
  df <- data.table::copy(df[,.(x,y, id, test)])
  n <- nrow(positions)
  radial_dist <- data.frame(coord1 = rep(1:n, each = n),
                      coord2 = rep(1:n, n),
                      distance = as.vector(as.matrix(proxy::dist(positions[,.(x, y)], pairwise = TRUE))))
  thr_dist <- link_d_to_analyse * 2
  radial_dist <- merge(radial_dist, df[,.(id, test)], by.x = 'coord2', by.y = 'id', how = "left")
  radial_dist <- as.data.table(radial_dist)
  distances <- radial_dist[distance <= thr_dist, .(max_distance = max(distance)), by = 'coord1']
  radial_dist <- radial_dist[distance <= thr_dist, .(prevalence = sum(test) / .N), by = 'coord1']
  
  radial_prev <- data.table::copy(df)
  radial_prev <- merge(radial_prev, radial_dist, by.x = 'id', by.y = 'coord1', how = "left")
  return(list("prevalence" = radial_prev, "distances" = distances))
}

## CENTROID APPROACH
# Compute euclidean distance

get_radious <- function(
    positions, 
    test_result, 
    total_friends_indeces,
    thr_data = 0.1, 
    max_epi_cont = 0.5){
  # Find centroid of X & Y
  centroid_x <- mean(positions[total_friends_indeces]$x)
  centroid_y <- mean(positions[total_friends_indeces]$y)
  centroid_df <- data.table("x" = centroid_x, "y" = centroid_y)
  
  combs <- calc_ind_prev(centroid_df,positions, test_result)

  # Filter by 10% of overall data closest to epifriends
  combs_eval <- combs[1: (thr_data * nrow(combs))]
  combs_eval[, is_epifriends := ifelse(id %in% total_friends_indeces, 1, 0)]
  
  # Evaluate if there is more than 50% of obs belonging to epifriends.
  # In that case, increase sample size to satisfy condition
  if( (nrow(combs_eval[is_epifriends == 1]) / nrow(combs_eval)) > max_epi_cont){
    n_total_eval = nrow(combs_eval[is_epifriends == 1]) / max_epi_cont
    combs_eval <- combs[1:n_total_eval]
    combs_eval <- na.omit(combs_eval)
  }
  
  return(list(
    "centroid" = centroid_df,
    "radious" = max(combs_eval$distance), 
    "prevalence" = sum(combs_eval$test) / nrow(combs_eval),
    "id" = combs_eval$id))
}

get_all_results <- function(catalogue, methods_list){
  results <- list()
  counter = 1
  for(method in methods_list){
    mean_pos <- rbindlist(catalogue_methods[[method]]$epifriends_catalogue$mean_position_all)
    mean_pos[, pvalue := catalogue_methods[[method]]$epifriends_catalogue$p]
    mean_pos[, method := method]
    mean_pos[, epifriends := 1:nrow(mean_pos)]
    results[[counter]] <- mean_pos
    counter <- counter +1
  }
  results <- rbindlist(results)
  
  results[, method_epifriends := paste0(epifriends, "_", method)]
  
  graph_all <- ggplot(results, aes(x = method_epifriends, y = pvalue,  fill = epifriends)) +
    geom_bar(stat = "identity") +
    ylim(0, 0.2) + 
    labs(title = "P-Value for each Epifriends-Methodology Combination", 
         x = "Epifriends-Methodology", y = "P-Value") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(list("table" = results, "chart" = graph_all))
}

# Plot significant clusters
coords <- data.table::copy(position_rand)
coords[, cluster := 0]
coords[, index := 1:nrow(coords)]
for(clusters in which(categories$epifriends_catalogue$p <= 0.05)){
  coords[index %in% categories$epifriends_catalogue$indeces[[clusters]], cluster := clusters]
}

