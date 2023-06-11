# This source code is licensed under the GNU GENERAL PUBLIC LICENSE license found in the
# LICENSE file in the root directory of this source tree.

#' Remove false detections from the clusters identified by the catalogue
#'
#' @param positions data.frame with the positions of parameters we want to query with shape (n,2) where n is the number of positions.
#' @param test_result data.frame with the test results (0 or 1).
#' @param link_d The linking distance to connect cases. Should be in the same scale as the positions.
#' @param use_link_d: If True, use linking distance to determine the closest neighbours. If False, use default linking neighbours based on proximity.
#' @param min_neighbours: Minium number of neighbours in the radius < link_d needed to link cases as friends.
#' @param link_neighbours: Number of surrounding neighbors to link. 
#' @param n_simulations Number of random simulations to be performed.
#' @param cluster_id Numeric vector with the cluster IDs of each position, with 0 for those without a cluster. Give NULL if cluster_id must be calculated.
#' @param min_neighbours Minium number of neighbours in the radius < link_d needed to link cases as friends.
#' @param max_p Maximum value of the p-value to consider the cluster detection.
#' @param min_pos Threshold of minimum number of positive cases in clusters applied.
#' @param min_total Threshold of minimum number of cases in clusters applied.
#' @param min_pr Threshold of minimum positivity rate in clusters applied.
#' @param keep_null_tests Whether to remove or not missings. If numeric, provide value to impute.
#'
#' @return  data.frame with the number of isolated positives detected per iteration and its frequency.
#'   
#' @export
#' 
#' @author Eric Matamoros
#'
#' @examples
#' # Required packages
#' if(!require("RANN")) install.packages("RANN")
#' library("RANN")
#'
#' # Creation of x vector of longitude coordinates, y vector of latitude coordinates and finaly merge them on a position data frame.
#' x <- c(1,2,3,4,7.5,8,8.5,9,10,13,13.1,13.2,13.3,14,15,30)
#' y <- c(1,2,3,4,7.5,8,8.5,9,10,13,13.1,13.2,13.3,14,15,30)
#' pos <- data.frame(x,y)
#'
#' # Creation of test data frame with 0 for negative cases and 1 for positive clases for each position.
#' test <- c(0,1,1,0,1,0,1,1,0,0,0,0,1,0,1,1)
#'
#' # Creation of catalogue for this positions, linking distance 2 and default values.
#' cat <- rm_false_detections(pos, test, 2, 100, filter_fd = TRUE)
#' catalogues <- cat$catalogue
#' 
get_false_detection <- function(positions, test_result, link_d,use_link_d = TRUE, 
                                link_neighbours = 10, n_simulations = 500,
                                min_neighbours = 2, max_p = 1, min_pos = 2, min_total = 2,
                                min_pr = 0,  keep_null_tests = FALSE){
  
  ## Generate simulations
  scenarios <- generate_simulations(positions = positions, test_result = test_result, 
                                    link_d = link_d, use_link_d = use_link_d, 
                                    link_neighbours = link_neighbours, n_simulations = n_simulations,
                                    min_neighbours = min_neighbours, max_p=max_p,
                                    min_pos=min_pos, min_total=min_total, min_pr=min_pr, keep_null_tests=keep_null_tests)
  
  # Obtain mean epifriends per number of positives
  random_dist <- scenarios[,.(mean_epi = mean(epifriends)), by = "num_pos"]
  random_dist[, num_pos := as.numeric(num_pos)]
  setorderv(random_dist, "num_pos")
  
  ## Compute catalogue with real infections
  cat <- catalogue(x = positions$x, y = positions$y,test_result = test_result, link_d=link_d, 
                   use_link_d = use_link_d, link_neighbours = link_neighbours,
                   cluster_id = NULL, min_neighbours = min_neighbours, max_p=max_p,
                   min_pos=min_pos, min_total=min_total, min_pr=min_pr, keep_null_tests=keep_null_tests,
                   in_latlon=FALSE, to_epsg = NULL, verbose = FALSE)
  cat_real <- data.table(
    'id' = cat$epifriends_catalogue$id,
    'positive' = cat$epifriends_catalogue$positives,
    'negative' = cat$epifriends_catalogue$negatives
  )
  cat_real <- cat_real[, .N, by = "positive"]
  names(cat_real) <- c("num_pos", "epifriends")
  cat_real[, num_pos := as.numeric(num_pos)]
  
  # Merge information and compute probability of false detection
  fp <- merge(cat_real, random_dist, by = c("num_pos"), how = "left")
  fp[, prob_fd := mean_epi / epifriends]
  
  # Pending to think about how 
  return(fp)
}

#' Generate random distributions based on the global prevalence of infections and compute its Epifriends to
#' undermine false positive detections.
#'
#' @param positions data.frame with the positions of parameters we want to query with shape (n,2) where n is the number of positions.
#' @param test_result data.frame with the test results (0 or 1).
#' @param link_d The linking distance to connect cases. Should be in the same scale as the positions.
#' @param use_link_d: If True, use linking distance to determine the closest neighbours. If False, use default linking neighbours based on proximity.
#' @param min_neighbours: Minium number of neighbours in the radius < link_d needed to link cases as friends.
#' @param link_neighbours: Number of surrounding neighbors to link. 
#' @param n_simulations Number of random simulations to be performed.
#' @param cluster_id Numeric vector with the cluster IDs of each position, with 0 for those without a cluster. Give NULL if cluster_id must be calculated.
#' @param min_neighbours Minium number of neighbours in the radius < link_d needed to link cases as friends.
#' @param max_p Maximum value of the p-value to consider the cluster detection.
#' @param min_pos Threshold of minimum number of positive cases in clusters applied.
#' @param min_total Threshold of minimum number of cases in clusters applied.
#' @param min_pr Threshold of minimum positivity rate in clusters applied.
#' @param keep_null_tests Whether to remove or not missings. If numeric, provide value to impute.
#' @param in_latlon:  If True, x and y coordinates are treated as longitude and latitude respectively, otherwise they are treated as cartesian coordinates.
#' @param to_epsg: If in_latlon is True, x and y are reprojected to this EPSG.
#' @param verbose: If TRUE, print information of the process; else, do not print.
#'
#' @return  data.frame with the number of isolated positives detected per iteration and its frequency.
#'   
#' @export
#' 
#' @author Eric Matamoros
#'
#' @examples
#' # Required packages
#' if(!require("RANN")) install.packages("RANN")
#' library("RANN")
#'
#' # Creation of x vector of longitude coordinates, y vector of latitude coordinates and finaly merge them on a position data frame.
#' x <- c(1,2,3,4,7.5,8,8.5,9,10,13,13.1,13.2,13.3,14,15,30)
#' y <- c(1,2,3,4,7.5,8,8.5,9,10,13,13.1,13.2,13.3,14,15,30)
#' pos <- data.frame(x,y)
#'
#' # Creation of test data frame with 0 for negative cases and 1 for positive clases for each position.
#' test <- c(0,1,1,0,1,0,1,1,0,0,0,0,1,0,1,1)
#'
#' # Creation of catalogue for this positions, linking distance 2 and default values.
#' cat <- generate_simulations(pos, test, 2, 100)
#' 
generate_simulations <- function(positions, test_result, link_d, use_link_d = TRUE, 
                                 link_neighbours = 10, n_simulations = 500,
                                 min_neighbours = 2, max_p = 1, min_pos = 2, min_total = 2,
                                 min_pr = 0, keep_null_tests = FALSE){
  
  # Generate data.table
  df <- data.table(x = positions$x, y = positions$y, test = test_result)
  
  prevalence <- sum(df$test) / nrow(df)
  
  cl <- makeCluster(detectCores() - 2)
  registerDoParallel(cl)
  
  # Use foreach to parallelize the for loop
  scenarios <- foreach(i = 1:n_simulations) %dopar% {
    
    library(epifriends)
    library(RANN)
    df[, test := stats::rbinom(nrow(df), 1, prevalence)]
    cat <- catalogue(x = df$x, y = df$y, test_result = df$test, link_d = link_d, 
                     link_d = link_d, use_link_d = use_link_d, 
                     link_neighbours = link_neighbours,
                     cluster_id = NULL, min_neighbours = min_neighbours,
                     max_p = max_p, min_pos = min_pos, min_total = min_total,
                     min_pr = min_pr, keep_null_tests = keep_null_tests, in_latlon = FALSE, 
                     to_epsg = NULL, verbose = FALSE)

    cat_count <- data.table(
      'positive' = cat$epifriends_catalogue$positives,
      'negative' = cat$epifriends_catalogue$negatives
    )
    
    cat_count <- cat_count[, .N, by = "positive"]
    names(cat_count) <- c("num_pos", "epifriends")
    cat_count[, iteration := i]
    
    return(cat_count)
  }
  # Stop the parallel backend
  stopCluster(cl)
  
  cat_df <- rbindlist(scenarios)
  
  return(cat_df)
}