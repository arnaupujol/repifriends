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
    ylim(0, 0.2) + 
    labs(title = "P-Value for each Epifriends-Methodology Combination", 
         x = "Epifriends-Methodology", y = "P-Value") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
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

get_region_prev <- function(data, epifriends_catalogue){
  epi_df <- data.table('id' = as.character(), 'x' = as.numeric(), 'y' = as.numeric(), 'epifriends' = as.numeric())
  for(idx in 1:length(epifriends_catalogue$indeces_local_prev)){
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
