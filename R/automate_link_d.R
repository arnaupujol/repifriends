# This source code is licensed under the GNU GENERAL PUBLIC LICENSE license found in the
# LICENSE file in the root directory of this source tree.

#' This method optimizes the linking distance based on the distribution of the minimum distances that satisfy a given
#' minimum neighbor restriction.
#'
#' @param df data.table with coordinates and test results. 
#' @param cluster_id Numeric vector with the cluster IDs of each position, with 0 for those without a cluster. Give NULL if cluster_id must be calculated.
#' @param min_neighbours Minium number of neighbours in the radius < link_d needed to link cases as friends.
#' @param dist_prop Distance threshold to merge two quantiles due to proximity.
#' @param thr_impr Percentage of improvement of the optimization metric to keep iterating.
#' @param diff_quantile Step increase/decrease of best quantile to avoid local minima.
#' @param quantile_est If 'fixed', quantiles used will be from 0 to 1 with a step of 0.125. If 'random', a random estimation of the step will be done.
#' @param keep_null_tests Whether to remove or not missings. If numeric, provide value to impute.
#' @param in_latlon:  If True, x and y coordinates are treated as longitude and latitude respectively, otherwise they are treated as cartesian coordinates.
#' @param to_epsg: If in_latlon is True, x and y are reprojected to this EPSG.
#' @param verbose: If TRUE, print information of the process; else, do not print.
#'
#' @return Best linking distance based on optimization metrics
#'   
#' @export
#' 
#' @author Eric Matamoros.
#'
#' @examples
#' # Required packages
#' if(!require("RANN")) install.packages("RANN")
#' if(!require("doParallel")) install.packages("doParallel")
#' library("RANN")
#' library("doParallel")
#'
#' # Creation of x vector of longitude coordinates, y vector of latitude coordinates and finaly merge them on a position data frame.
#' x <- c(1,2,3,4,7.5,8,8.5,9,10,13,13.1,13.2,13.3,14,15,30)
#' y <- c(1,2,3,4,7.5,8,8.5,9,10,13,13.1,13.2,13.3,14,15,30)
#' test <-c(0,1,1,0,1,0,1,1,0,0,0,0,1,0,1,1)
#' 
#' df <- data.table(x = x, y = y, test = test)
#'
#' # Creation of catalogue for this positions, linking distance 2 and default values.
#' link_d <- opt_link_d(df, 2)
opt_link_d <- function(df, min_neighbors, cluster_id=NULL, dist_prop = 0.1,
                       thr_impr = 0.05, diff_quantile = 0.05, quantile_est =  "fixed", 
                       keep_null_tests = FALSE, in_latlon = FALSE, to_epsg = NULL, 
                       verbose = FALSE){
  
  # Obtain position coordinates & test results
  position <- df[,.(x,y)]
  positive = (df$test == 1)
  test_result <- df[,.(test)]
  
  # Obtain minimum distances that satisfy the min_neighbors constrain
  if(verbose){print("Obtain minimum distances satisfying the min_neighbor restriction")}
  min_distances <- get_min_distances(position[,.(x,y)], positive, min_neighbors)

  # Get quantiles from whose distance will be evaluated
  if(verbose){print("Estimate quantiles based on distances")}
  quantiles <- quantile_estimation(quantile_est, verbose)
  
  if(verbose){print("Simplify quantiles based on distance distributions")}
  # Normalize and compute quantiles dist
  min_distances_norm <- (min_distances - min(min_distances)) / (max(min_distances) - min(min_distances))
  quantiles_dist <- quantile(min_distances_norm,probs = c(quantiles),na.rm = TRUE)
  dt = data.table('quantiles' = quantiles, 'distances'= as.numeric(quantiles_dist))
  # Simplify estimated quantiles based on distance proxy
  dt <- simplify_distributions(dt[!(quantiles %in% c(0, 1))], dist_prop)
  dt[, distances := quantile(min_distances,probs = c(dt$quantiles),na.rm = TRUE)]
  
  if(verbose){print("Optimize linking distance based on cases of real positives")}
  opt <- optimize_positives(position, test_result, min_neighbors,cluster_id,min_distances, dt$quantiles, 
                            thr_impr, diff_quantile, keep_null_tests, in_latlon,to_epsg)
  
  
  return(opt$optimal_link_d)
}

#' Obtain desired quantiles that will be used to evaluate distances.
#'
#' @param quantile_est If 'fixed', quantiles used will be from 0 to 1 with a step of 0.125. If 'random', a random estimation of the step will be done.
#' @param verbose: If TRUE, print information of the process; else, do not print.
#'
#' @return Vector with quantiles
#'   
#' @export
#' 
#' @author Eric Matamoros.
#'
#' @examples
#' quantiles <- quantile_estimation('fixed', FALSE)
quantile_estimation <- function(quantile_est, verbose){
  if(quantile_est == 'random'){
    incr <- 0.05
    quantiles <- c(0, 0.25,0.5,0.75, 1)
    # Generate 2 random values between each of the standard quantiles
    quantiles <- sort(as.vector(sapply(1:(length(quantiles) - 1), function(x) 
      runif(2, min = quantiles[x], max = quantiles[x+1])
    )))
    if(verbose){print("Random selection method of quantiles")}
  }else if(quantile_est == 'fixed'){
    quantiles <- seq(0, 1, by = 0.125)
    if(verbose){print("Fixed selection method of quantiles by intervals of 0.125")}
  }else{
    stop("Any valid selected method for quantile estimation. Please check documentation")
  }
  
  return(quantiles)
}

#' This method optimizes the linking distance based on the distribution of the minimum distances that satisfy a given
#' minimum neighbor restriction.
#'
#' @param positions data.frame with the positions of parameters we want to query with shape (n,2) where n is the number of positions.
#' @param positive Array with location of positive tests.
#' @param min_neighbours Minium number of neighbours in the radius < link_d needed to link cases as friends.
#'
#' @return Minimum distance for each positive test where a given contrain of min_neighbors is satisfied.
#'   
#' @export
#' 
#' @author Eric Matamoros.
#'
#' @examples
#' # Creation of x vector of longitude coordinates, y vector of latitude coordinates and finaly merge them on a position data frame.
#' x <- c(1,2,3,4,7.5,8,8.5,9,10,13,13.1,13.2,13.3,14,15,30)
#' y <- c(1,2,3,4,7.5,8,8.5,9,10,13,13.1,13.2,13.3,14,15,30)
#' test <-c(0,1,1,0,1,0,1,1,0,0,0,0,1,0,1,1)
#' 
#' df <- datatable(x = x, y = y, test = test)
#' position <- df[,.(x,y)]
#' positive = (df$test == 1)
#' min_neighbors <- 2
#'
#' # Get minimum distances based on neighbors
#' min_distances <- get_min_distances(position, positive, min_neighbors)
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
    
    if( nrow(rad_test[index >= min_neighbors]) == 0){
      min_distances[counter] <- 0
    }else{
      min_distances[counter] <- min(rad_test[index >= min_neighbors]$distance)[1]
    }
    counter <- counter +1 
  }
  
  return(min_distances)
}

#' Merge quantiles based on proximity of distances given its distribution.
#'
#' @param dt data.table with quantiles and minimum distances.
#' @param dist_prop Distance threshold to merge two quantiles due to proximity.
#'
#' @return data.table with simplified quantiles based on proximity.
#'   
#' @export
#' 
#' @author Eric Matamoros.
#'
#' @examples
#'
#' # Creation of x vector of longitude coordinates, y vector of latitude coordinates and finaly merge them on a position data frame.
#' quantiles <- c(0, 0.125, 0.250, 0.5)
#' y <- c(1,2,3,4,7.5,8,8.5,9,10,13,13.1,13.2,13.3,14,15,30)
#' test <-c(0,1,1,0,1,0,1,1,0,0,0,0,1,0,1,1)
#' 
#' dt <- datatable( quantiles = c(0, 0.125, 0.250, 0.5), distances = c(0.01, 0.02, 0.03, 0.04)
#' 
#' # Simplify distributions based on proximity
#' dt <- simplify_distributions(dt, 0.1)
simplify_distributions<- function(dt, dist_prop){
  dt_prev <- dt[1:2,]
  dt_post <- dt[1,]
  counter <- 1
  while( (dim(dt_prev)[1] != dim(dt_post)[1]) & (counter < 4)){
    # Loop over consecutive pairs of rows
    dt_prev <- data.table::copy(dt)
    for (i in 1:(nrow(dt)-1)) {
      
      if(i >= nrow(dt)){
        break
      }
      # Check if distances are lower than threshold (dist_thr)
      if (abs(dt$distances[i] - dt$distances[i+1])  < dist_prop) {
        # Compute mean of distances and update first row
        dt[i, distances := mean(c(dt$distances[i], dt$distances[i+1]))]
        dt[i, quantiles := mean(c(dt$quantiles[i], dt$quantiles[i+1]))]
        # Delete second row
        dt <- dt[-(i+1), ]
        # Update row index
        i <- i - 1
      }
    }
    dt_post <- data.table::copy(dt)
    counter <- counter +1 
  }
  
  return(dt)
}

#' This method optimizes the linking distance based on the distribution of the minimum distances that satisfy a given
#' minimum neighbor restriction.
#'
#' @param df data.table with coordinates and test results. 
#' @param cluster_id Numeric vector with the cluster IDs of each position, with 0 for those without a cluster. Give NULL if cluster_id must be calculated.
#' @param min_neighbours Minium number of neighbours in the radius < link_d needed to link cases as friends.
#' @param dist_prop Distance threshold to merge two quantiles due to proximity.
#' @param thr_impr Percentage of improvement of the optimization metric to keep iterating.
#' @param diff_quantile Step increase/decrease of best quantile to avoid local minima.
#' @param quantile_est If 'fixed', quantiles used will be from 0 to 1 with a step of 0.125. If 'random', a random estimation of the step will be done.
#' @param keep_null_tests Whether to remove or not missings. If numeric, provide value to impute.
#' @param in_latlon:  If True, x and y coordinates are treated as longitude and latitude respectively, otherwise they are treated as cartesian coordinates.
#' @param to_epsg: If in_latlon is True, x and y are reprojected to this EPSG.
#' @param verbose: If TRUE, print information of the process; else, do not print.
#'
#' @return List with linking_distance & trace of iterations
#'   
#' @export
#' 
#' @author Eric Matamoros.
#'
#' @examples
#' # Required packages
#' if(!require("RANN")) install.packages("RANN")
#' if(!require("doParallel")) install.packages("doParallel")
#' library("RANN")
#' library("doParallel")
#'
#' # Creation of x vector of longitude coordinates, y vector of latitude coordinates and finaly merge them on a position data frame.
#' x <- c(1,2,3,4,7.5,8,8.5,9,10,13,13.1,13.2,13.3,14,15,30)
#' y <- c(1,2,3,4,7.5,8,8.5,9,10,13,13.1,13.2,13.3,14,15,30)
#' pos <- data.frame(x,y)
#'
#' # Creation of test data frame with 0 for negative cases and 1 for positive clases for each position.
#' test <- data.frame(c(0,1,1,0,1,0,1,1,0,0,0,0,1,0,1,1))
#' 
#' min_distances <- c(0.01, 0.02, 0.05, 0.02, 0.07, 0.08, 0.01, 0.02)
#'
#' # Creation of catalogue for this positions, linking distance 2 and default values.
#' link_d_list <- optimize_positives(pos, test, 2, NULL, min_distances)
#' 
optimize_positives<- function(position, test_result, min_neighbors,cluster_id,
                              min_distances, quantiles, thr_impr = 0.05, 
                              diff_quantile = 0.05, keep_null_tests = FALSE, 
                              in_latlon = FALSE, to_epsg = NULL){
  
  opt_metric <- 0
  counter <- 1
  opt_quantile <- c()
  
  traces <- list()
  while( (thr_impr >= 0.05) & (counter < 4)){
    
    quantiles_dist <- quantile(min_distances,unique(quantiles),na.rm = TRUE)
    
    cl <- makeCluster(detectCores() - 2)
    registerDoParallel(cl)
    
    ## Compute scoring metric using parallel or for loop depending on size of the evaluation inputs
    if(length(quantiles_dist) > 2){
      
      # Use foreach to parallelize the for loop
      metric <- foreach(quant = 1:length(quantiles_dist), .combine = "c") %dopar% {
        
        library(epifriends)
        library(RANN)
        # Code to execute in parallel
        categories <- catalogue(
          data.table::copy(position), test_result, quantiles_dist[quant], cluster_id = NULL, 
          min_neighbours = min_neighbors, keep_null_tests = keep_null_tests, in_latlon=in_latlon,
          to_epsg=to_epsg)
        
        kpi <- 0
        for(j in 1:length(categories$cluster_id)){
          if(categories$cluster_id[j]!=0){
            kpi <- kpi + (1-categories$pval_cluster[j])
          }
        }
        
        return(kpi)
      }
      # Stop the parallel backend
      stopCluster(cl)
    }else{
      metric <- c()
      for(quant in 1:length(quantiles_dist)){
        categories <- catalogue(
          data.table::copy(position), test_result, quantiles_dist[quant], cluster_id = NULL, 
          min_neighbours = min_neighbors, keep_null_tests = keep_null_tests, in_latlon=in_latlon,
          to_epsg=to_epsg)
        
        kpi <- 0
        for(j in 1:length(categories$cluster_id)){
          if(categories$cluster_id[j]!=0){
            kpi <- kpi + (1-categories$pval_cluster[j])
          }
        }
        
        metric[quant] <- kpi
      }
      
    }
    
    # Obtain best quantile and best metric
    best_quantile <- quantiles[which(metric == max(metric))]
    best_metric <- max(metric)
    
    # Store for tracing and plots
    traces[[counter]] <- data.table('quantiles' = unique(quantiles), 'metrics' = metric, 'iteration' = counter)
    
    # Check if metric has been improved by more than a 5%
    perc_impr <- (best_metric - (opt_metric + 0.00000001)) / (opt_metric + 0.00000001)
    if(perc_impr > thr_impr){
      opt_metric <- best_metric
      opt_quantile <- best_quantile
      counter <- counter + 1
    }else{
      break
    }
    
    # Create list of vectors and cap
    quantiles <- c(opt_quantile - diff_quantile, opt_quantile, opt_quantile + diff_quantile)
    quantiles <- pmin(quantiles, 1)
    quantiles <- pmax(quantiles, 0)
  }
  
  opt_link_d <- as.numeric(quantile(min_distances,opt_quantile,na.rm = TRUE))
  
  # Obtain traces
  traces <- unique(rbindlist(traces))
  setorderv(traces, 'quantiles')
  traces[, link_d := quantile(min_distances,quantiles,na.rm = TRUE)]
  
  return(list('optimal_link_d' = opt_link_d, 'traces' = traces))
  
}