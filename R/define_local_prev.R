# This source code is licensed under the GNU GENERAL PUBLIC LICENSE license

#' This function performs a KMeans clustering on the coordinates data.
#'
#' @param positions data.frame with the positions of parameters we want to query 
#' with shape (n,2) where n is the number of positions.
#' @param test_result data.frame with the test results (0 or 1).
#' @param seed Seed for reproducibility.
#' @param coord_cols Vector with the name of the columns corresponding to the 
#' x,y coordinates.
#'
#' @return  Data.table with new columns corresponding to the clusters & the 
#' individual prevalence of each coordinate.
#'   
#' @importFrom clValid clValid optimalScores
#' @importFrom stats kmeans
#' 
#' @export
#' 
#' @author Eric Matamoros.
#'
#' @examples
#' # Creation of x vector of longitude coordinates, y vector of latitude coordinates 
#' # and finally merge them on a position data frame.
#' x <- c(1,2,3,4,7.5,8,8.5,9,10,13,13.1,13.2,13.3,14,15,30)
#' y <- c(1,2,3,4,7.5,8,8.5,9,10,13,13.1,13.2,13.3,14,15,30)
#' pos <- data.table("x" = x,"y" = y)
#'
#' # Creation of test data frame with 0 for negative cases and 1 for positive 
#' # classes for each position.
#' test <- c(0,1,1,0,1,0,1,1,0,0,0,0,1,0,1,1)
#'
#' # Creation of catalogue for this positions, linking distance 2 and default values.
#' pos_rate <- compute_kmeans(pos, test)
#' 
compute_kmeans <- function(positions, test, seed = 123, coord_cols = c("x", "y")){
  
  df = data.table(positions)
  df[, test := test]
  df[, id := 1:nrow(df)]

  # Leverage clValid to automatic detect optimal number of clusters
  set.seed(seed)
  
  suppressWarnings({
    # Optimal cluster based on connectivity
    intern <- clValid(as.data.frame(df[,..coord_cols]),
                      nClust = 2:8,clMethods = c("kmeans"),
                      validation = "internal",
                      maxitems = nrow(df))
  })
  
  rows <- rownames(optimalScores(intern))
  optimal <- as.data.table(optimalScores(intern))
  optimal[, approach := rows]
  n_clusters <- as.numeric(optimal[approach == "Dunn"]$Clusters[1])
  
  k_results <- kmeans(df[,..coord_cols], centers = n_clusters, nstart = 500)
  
  df[, clusters := k_results$cluster]
  
  df[, prevalence := sum(test) / .N, by = "clusters"]
  
  return(df)
}

#' This function determines the prevalence based on the closest groups to a centroid.
#'
#' @param positions data.frame with the positions.
#' @param total_friends_indeces indices of the events detected by the epifriends.
#' @param test_result data.frame with the test results (0 or 1).
#' @param max_epi_cont Maximum percentage contribution of the number of epifriends 
#' indices with respect to all the ones.
#' @param thr_data Percentage of overall data considered for the local prevalence calculus.
#'
#' @return Local prevalence based on the positive cases for the closest datapoints 
#' to the centroid defined by epifriends indices
#' 
#' @export
#' 
#' @author Eric Matamoros.
#'
#' @examples
#' # Creation of x vector of longitude coordinates, y vector of latitude coordinates 
#' # and finally merge them on a position data frame.
#' x <- c(1,2,3,4,7.5,8,8.5,9,10,13,13.1,13.2,13.3,14,15,30)
#' y <- c(1,2,3,4,7.5,8,8.5,9,10,13,13.1,13.2,13.3,14,15,30)
#' pos <- data.table("x" = x,"y" = y)
#' pos$id <- 1:length(x)
#' friends_ind <- c(1,2,3)
#'
#' # Creation of test data frame with 0 for negative cases and 1 for positive clases 
#' # for each position.
#' test <- c(0,1,1,0,1,0,1,1,0,0,0,0,1,0,1,1)
#' test_result <- data.table(test = test)
#'
#' # Calculus of prevalence based on the centroid
#' pos_rate <- compute_centroid(pos, friends_ind, test_result)
#'
compute_centroid <- function(
    positions, 
    total_friends_indeces,
    test_result, 
    max_epi_cont = 0.5, 
    thr_data = 0.1){
  
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
    'prevalence' = sum(combs_eval$test) / nrow(combs_eval),
    'local_id' = sort(combs_eval$id))
  )
}


#' Recursively merge KMean clusters based on proximity if percentage of epifriend 
#' data is above specific threshold.
#'
#' @param pos_clusters data.frame with the positions.
#' @param total_friends_indeces indices of the events detected by the epifriends.
#' @param thr_perc_epi percentage allowed contribution of epifriends to the overall cluster.
#'
#' @return Local prevalence based on the positive cases for the closest datapoints to the centroid 
#' defined by epifriends indices
#' 
#' @importFrom proxy dist
#'
#' @export
#' 
#' @author Eric Matamoros.
#'
#' @examples
#' # Creation of x vector of longitude coordinates, y vector of latitude coordinates 
#' # and finally merge them on a position data frame.
#' x <- c(1,2,3,4,7.5,8,8.5,9,10,13,13.1,13.2,13.3,14,15,30)
#' y <- c(1,2,3,4,7.5,8,8.5,9,10,13,13.1,13.2,13.3,14,15,30)
#' test <- c(1,1,0,1,0,0,0,1,1,0,1,0,0,1,1,0)
#' test_result <- data.table(test = test)
#' pos <- data.table("x" = x,"y" = y, "test" = test)
#' pos$id <- 1:length(x)
#' pos$cluters <- 1
#' friends_ind <- c(1,2,3)
#'
#' # Calculus of prevalence based on the centroid
#' # Compute KMeans
#' pos_clusters <- compute_kmeans(pos, test_result)
#' merged_clusters <- merge_kmeans_clusters(pos_clusters,friends_ind, 0.5)
#'
merge_kmeans_clusters<- function(pos_clusters,total_friends_indeces, thr_perc_epi = 0.5){
  
  epi_clusters <- unique(pos_clusters[total_friends_indeces]$clusters)
  perc_epi <- dim(pos_clusters[total_friends_indeces])[1] /
    dim(pos_clusters[clusters %in% epi_clusters])[1]
  
  if(perc_epi >= thr_perc_epi){
    # Merge clusters sequentially and re-compute prevalence based on new added clusters
    counter <- 0
    while( (perc_epi >= thr_perc_epi) & (counter < 3)){
      mean_epis <- pos_clusters[clusters %in% epi_clusters, .(x = mean(x), y = mean(y))]
      remain_epis <- pos_clusters[!(clusters %in% epi_clusters),
                                  .(x = mean(x), y = mean(y)), by = c("clusters")]
      # If no clusters rema
      if(dim(remain_epis)[1] == 0){
        break
      }
      combs <- t(proxy::dist(mean_epis,remain_epis[,.(x, y)], pairwise = FALSE))
      remain_epis$distance <- as.vector(combs)
      
      # Merge new closest cluster to all
      epi_clusters <- c(epi_clusters, remain_epis[distance == min(distance)]$clusters[1])
      
      pos_clusters[clusters %in% epi_clusters, prevalence := sum(test) / .N]
      counter <- counter + 1
      perc_epi <- dim(pos_clusters[total_friends_indeces])[1] /
        dim(pos_clusters[clusters %in% epi_clusters])[1]
    }
  }
  
  return(list("pos_clusters" = pos_clusters, "epi_clusters" = epi_clusters))
  
}


#' This function is a helper of calc_distance and helps compute the distance of 
#' one ID versus the rest
#'
#' @param df Dataframe with asked coordinates.
#' @param positions data.frame with the positions of parameters we want to query 
#' with shape (n,2) where n is the number of positions.
#' @param test_result data.frame with the test results (0 or 1).
#'
#' @return  Data.table computed distances resulting from df mapped versus all 
#' coordinates in positions.
#'   
#' @importFrom proxy dist
#' @import data.table
#'
#' @export
#' 
#' @author Eric Matamoros.
#'
#' @examples
#' # Creation of x vector of longitude coordinates, y vector of latitude coordinates 
#' # and finaly merge them on a position data frame.
#' x <- c(1,2,3,4,7.5,8,8.5,9,10,13,13.1,13.2,13.3,14,15,30)
#' y <- c(1,2,3,4,7.5,8,8.5,9,10,13,13.1,13.2,13.3,14,15,30)
#' pos <- data.table("x" = x,"y" = y)
#'
#' # Creation of test data frame with 0 for negative cases and 1 for positive clases 
#' # for each position.
#' test <- data.frame(c(0,1,1,0,1,0,1,1,0,0,0,0,1,0,1,1))
#' test_result <- data.table(test = test)
#' 
#' df <- data.table("x" = 1, "y" = 4)
#'
#' # Creation of catalogue for this positions, linking distance 2 and default values.
#' ind_dist <- calc_ind_prev(df, pos, test)
#'
calc_ind_prev <- function(df, positions, test_result){
  combs <- t(proxy::dist(df[,.(x,y)], positions[,.(x, y)], pairwise = FALSE))
  combs <- data.table(
    "id" = positions$id, "distance" = combs[,1], "test" = test_result$test)
  setorderv(combs, cols = "distance", order=1L)
  return(combs)
}
