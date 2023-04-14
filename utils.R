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
  
  # Compute euclidean distance
  combs <- calc_ind_prev(centroid_df,positions, test_result)
  
  # Filter by 10% of overall data closest to epifriends
  combs_eval <- combs[1: (thr_data * nrow(combs))]
  combs_eval[, is_epifriends := ifelse(id %in% total_friends_indeces, 1, 0)]
  
  rows_epi <-  nrow(combs_eval[is_epifriends == 1])
  if(rows_epi < length(total_friends_indeces)){
    n_total_eval <- length(total_friends_indeces)/max_epi_cont
    combs_eval <- combs[1:n_total_eval]
    combs_eval <- na.omit(combs_eval)
    
  }else{
    
    # Evaluate if there is more than 50% of obs belonging to epifriends.
    # In that case, increase sample size to satisfy condition
    if( (rows_epi / nrow(combs_eval)) > max_epi_cont){
      n_total_eval = nrow(combs_eval[is_epifriends == 1]) / max_epi_cont
      combs_eval <- combs[1:n_total_eval]
      combs_eval <- na.omit(combs_eval)
    }
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
    labs(title = "P-Value for each Epifriends-Methodology Combination", 
         x = "Epifriends-Methodology", y = "P-Value") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  graph_all <- graph_all + geom_hline(yintercept = 0.05, linetype = "dashed")
  
  return(list("table" = results, "chart" = graph_all))
}

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

get_region_prev <- function(data, epifriends_catalogue, significant = FALSE){
  epi_df <- data.table('id' = as.character(), 'x' = as.numeric(), 'y' = as.numeric(), 'epifriends' = as.numeric())
  
  if(significant){
    levels <- which(epifriends_catalogue$p <= 0.05)
  }else{
    levels <- 1:length(epifriends_catalogue$indeces_local_prev)
  }
  for(idx in levels){
    df <- data.table(
      'id' = epifriends_catalogue$indeces_local_prev[[idx]],
      'epifriends' = rep(idx, length(epifriends_catalogue$indeces_local_prev[[idx]]))
    )
    df = merge(df, data[,.(x,y,id)], on = 'id', how= "left")
    epi_df <- rbind(epi_df, df)
  }
  
  return(epi_df)
}

clean_data <- function(positions, keep_null_tests = FALSE){
  to_impute <- colnames(positions)[!(colnames(positions) %in% c("x", "y"))]
  positions = clean_unknown_data(positions,to_impute,keep_null_tests,FALSE)
  return(positions)
}

get_min_distances <- function(position, positive, min_neighbors){
  
  n <- dim(position[positive, .(x,y)])[1]
  radial_dist <- data.table(coord1 = rep(1:n, each = n),
                            coord2 = rep(1:n, n),
                            distance = as.vector(as.matrix(proxy::dist(position[positive, .(x,y)], pairwise = TRUE))))
  # Remove identities
  radial_dist <- radial_dist[coord1 != coord2]
  
  # Determine minimal distance based on min_neighbors
  min_distances <- c()
  counter <- 1
  for(coords in unique(radial_dist$coord1)){
    rad_test <- radial_dist[coord1 == coords]
    rad_test <- rad_test[order(distance)]
    rad_test[, index := 1:dim(rad_test)[1]]
    min_distances[counter] <- min(rad_test[index >= min_neighbors]$distance)[1]
    counter <- counter +1 
  }
  
  return(min_distances)
}

best_link_d <- function(position, positive, test, min_neighbours){
  
  # Get minimum distance & quantiles of it
  min_distances <- get_min_distances(position[,.(x,y)], positive, min_neighbours)
  quantiles_dist <- quantile(min_distances,probs = c(0.25,0.5,0.75),na.rm = TRUE)
  
  # Determine the level of significance of p-values based on the distribution of distances
  sign_level <- c()
  for(quant in 1:length(quantiles_dist)){
    categories <- catalogue(
      data.table::copy(position), test, quantiles_dist[quant], cluster_id = NULL, min_neighbours = min_neighbours,
      method = method)
    
    sign_level[quant] <- length(which(categories$pval_cluster <= 0.05))
  }
  
  if(all(sign_level == 0)){
    return(quantiles_dist[2])
  }else{
    return(quantiles_dist[min(which(sign_level == max(sign_level)))])
  }
}


get_link <- function(file){
  data_v1 <- as.data.table(read_csv(paste0("data/", file)))
  data_v1[, id := 1:nrow(data_v1)]
  
  position_rand <- data_v1[,3:4]
  positive_rand = (data_v1$test == 1)
  test_rand <- data_v1$test
  
  return(as.numeric(best_link_d(position_rand, positive_rand, test_rand, min_neighbours)))
}

store_indiv <- function(pdf_path, files){
  for( output in as.character(sapply(files, function(x) strsplit(x, ".csv")[[1]][1] ))){
    
    path_specific <- paste0(pdf_path, output)
    # Create directory if not exists
    if(!(dir.exists(paste0(pdf_path, output)))){
      dir.create(path_specific)
    }
    
    # Remove files if not empty
    if (length(dir(path_specific)) != 0) {
      unlink(path_specific, recursive = TRUE)
    }
    
    files_specific <- list.files(pdf_path, pattern = output)
    files_specific <- files_specific[!(files_specific %in% output)]
    
    for(file_specific in files_specific){
      file.copy(
        paste0(pdf_path, file_specific), 
        paste0(pdf_path, output, "/", file_specific))
      file.remove(paste0(pdf_path, file_specific))
    }
  }
}

