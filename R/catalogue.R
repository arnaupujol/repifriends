# This source code is licensed under the GNU GENERAL PUBLIC LICENSE license found in the
# LICENSE file in the root directory of this source tree.

#' This method runs the DBSCAN algorithm (if cluster_id is NULL) and obtains the mean positivity rate (PR) of each cluster extended with the non-infected cases closer than the link_d.
#'
#' @param x Vector of x positions.
#' @param y Vector of y positions.
#' @param test_result vector of test results (0 or 1).
#' @param link_d The linking distance to connect cases. Should be in the same scale as the positions.
#' @param prevalence Probability of having an infected case for each individual.
#' @param cluster_id Numeric vector with the cluster IDs of each position, with 0 for those without a cluster. Give NULL if cluster_id must be calculated.
#' @param min_neighbours Minium number of neighbours in the radius < link_d needed to link cases as friends.
#' @param max_p Maximum value of the p-value to consider the cluster detection.
#' @param min_pos Threshold of minimum number of positive cases in clusters applied.
#' @param min_total Threshold of minimum number of cases in clusters applied.
#' @param min_pr Threshold of minimum positivity rate in clusters applied.
#' @param method Method that wants to be used to compute the local prevalence - either 'kmeans', 'centroid', or 'base'. Defaults to 'base'.
#' @param keep_null_tests Whether to remove or not missings. If numeric, provide value to impute.
#' @param in_latlon:  If True, x and y coordinates are treated as longitude and latitude respectively, otherwise they are treated as cartesian coordinates.
#' @param to_epsg: If in_latlon is True, x and y are reprojected to this EPSG.
#' @param max_epi_cont: Maximum contribution of the detected Epifriends with respect to the total local data selected. Only applies for "centroid" method.
#' @param max_thr_data: Percentage of data used to compute the local prevalence. Only applies for "centroid" method.
#' @param consider_fd: If True, consider false detections and adjust p-value of that.
#' @param n_simulations: Numeric value with the number of desired iterations to compute the false-detected clusters.
#' @param sampling_sim: Numeric value with the number of random simulations to be performed to determine the real p-value.
#' @param optimize_link_d: If True, optimize the linking distance based on minimum distribution of distances between neighbors.
#' @param verbose: If TRUE, print information of the process; else, do not print.
#'
#' @details The epifriends package uses the RANN package which can gives the exact nearest neighbours using the friends of friends algorithm. For more information on the RANN library please visit https://cran.r-project.org/web/packages/RANN/RANN.pdf
#'
#' @return  List with the next objects:
#' cluster_id: numeric vector
#'   Vector of the cluster IDs of each position, with 0 for those
#'   without a cluster.
#' mean_pr_cluster: numeric vector
#'   Mean PR corresponding to cluster_id.
#' pval_cluster: numeric vector
#'   P-value corresponding to cluster_id.
#' epifriends_catalogue: List
#'   Catalogue of the epifriends clusters and their main characteristics.
#'   
#' @export
#' 
#' @author Mikel Majewski Etxeberria based on earlier python code by Arnau Pujol.
#'
#' @examples
#' # Required packages
#' if(!require("RANN")) install.packages("RANN")
#' library("RANN")
#'
#' # Creation of x vector of longitude coordinates, y vector of latitude coordinates and finaly merge them on a position data frame.
#' x <- c(1,2,3,4,7.5,8,8.5,9,10,13,13.1,13.2,13.3,14,15,30)
#' y <- c(1,2,3,4,7.5,8,8.5,9,10,13,13.1,13.2,13.3,14,15,30)
#'
#' # Creation of test data frame with 0 for negative cases and 1 for positive clases for each position.
#' test <- c(0,1,1,0,1,0,1,1,0,0,0,0,1,0,1,1)
#'
#' # Creation of catalogue for this positions, linking distance 2 and default values.
#' cat <- catalogue(x, y, test, 2)
#' 
catalogue <- function(x, y, test_result, link_d,  prevalence = NULL,  cluster_id = NULL,
                      min_neighbours = 2, max_p = 1, min_pos = 2, min_total = 2,
                      min_pr = 0, method = "base",keep_null_tests = FALSE,in_latlon = FALSE,
                      to_epsg = NULL, max_epi_cont = 0.5, 
                      max_thr_data = 0.1, consider_fd = FALSE, n_simulations= 500,
                      sampling_sim = 1000, optimize_link_d = FALSE, verbose = FALSE){

  # Create data.frame
  suppressWarnings({
    positions <- data.table(
      "x" = x,
      "y" = y,
      "test" = test_result,
      "prevalence" = prevalence)
  })
  to_impute <- colnames(positions)[!(colnames(positions) %in% c("x", "y"))]
  positions = clean_unknown_data(positions,to_impute,keep_null_tests,verbose)
  test_result = data.table("test" = positions$test)
  if(!is.null(prevalence)){ prevalence = positions$prevalence}
  
  #Defining 2d-positions
  positions = get_2dpositions(
    x = positions$x, 
    y = positions$y, 
    in_latlon = in_latlon, 
    to_epsg = to_epsg,
    verbose = verbose)
  positions[, id := 1:nrow(positions)]

  # Optimize linking distance if not provided or specified to do so
  if(optimize_link_d | is.null(link_d)){
    if(verbose){print("Automating the calculus of the linking distance")}
    link_d <- opt_link_d(df=data.table(x = positions$x, y = positions$y, test = pos$test), 
                         min_neighbors=min_neighbors, cluster_id=cluster_id, 
                         keep_null_tests = keep_null_tests, in_latlon = in_latlon, 
                         to_epsg = to_epsg, verbose = verbose)
  }

  #Define positions of positive cases
  positive_positions <- positions[which(test_result == 1), .(x,y)]
  #Computing cluster_id if needed
  if(is.null(cluster_id)){
    cluster_id = dbscan(x = positive_positions$x, y = positive_positions$y, link_d, min_neighbours = min_neighbours)
  }

  #Define total number of positive cases
  total_positives = sum(test_result)
  #Define total number of cases
  total_n = dim(test_result)[1]
  #Initialising mean PR and p-value for the positive cases in clusters
  mean_pr_cluster <- numeric(length(cluster_id))
  pval_cluster <- rep(1, length(cluster_id))
  #EpiFRIenDs cluster catalogue
  epifriends_catalogue <- vector(mode='list', length=9)
  names(epifriends_catalogue) <- c("id",                  #EpiFRIenDs id
                                   "mean_position_pos",   #Mean position of positive cases
                                   "mean_position_all",   #Mean position of all cases
                                   "mean_pr",             #Positivity rate
                                   "positives",           #Number of positive cases
                                   "negatives",           #Number of negative cases
                                   "total",               #Total number of positions
                                   "indeces",             #Indeces of all positions
                                   "p"                    #p-value of detection
  )
  next_id <- 1
  sort_unici <- sort(unique(cluster_id[cluster_id>0]))
  #Case without any clusters
  if(length(sort_unici) == 0){
    return(epifriends_catalogue)
  }
  
  # Generate simulations for false positive detection
  if(consider_fd){
    false_det <- get_false_detection(positions=positions, test_result=test_result$test,
                                     link_d=link_d, n_simulations=n_simulations,
                                     min_neighbours = min_neighbours, max_p = max_p, 
                                     min_pos = min_pos, min_total = min_total,
                                     min_pr = min_pr, keep_null_tests = keep_null_tests, 
                                     in_latlon = in_latlon,to_epsg = to_epsg, verbose = verbose)
  }else{
    false_det <- NULL
  }
  
  
  for(i in 1:length(sort_unici)){
    
    #get all indeces with this cluster id
    cluster_id_indeces <- which(cluster_id == sort_unici[i])
    #for all these indeces, get list of friends from all positions
    all_friends_indeces <- find_indeces(positive_positions[cluster_id_indeces,.(x,y)], link_d, positions[,.(x,y)])
    #get unique values of such indeces
    total_friends_indeces <- sort(unique(unlist(all_friends_indeces)))
    #get positivity rate from all the unique indeces
    mean_pr <- test_result[total_friends_indeces, .(mean(test))][[1]]
    npos <- sum(test_result[total_friends_indeces,])
    ntotal <- length(total_friends_indeces)

    if(!is.null(prevalence)){
      if(verbose){print("Using user-given prevalence.")}
      ind_pos_rate = prevalence[positive_positions]
      trials <- simulate_trial(sampling_sim, 1, ind_pos_rate)
      pos_rate <- length(which(trials >= (npos / ntotal))) / sampling_sim
      
    }else{
      # Approaches to estimate the p-value
      if(method == "kmeans"){
        if(verbose){print("Using KMeans method to account for local prevalence.")}
        
        # Compute KMeans
        pos_clusters <- compute_kmeans(positions, test_result)
        ind_pos_rate <- pos_clusters[total_friends_indeces]
        
        # Merge recursively clusters if epifriend percentage over local population is above thr
        merged_clusters <- merge_kmeans_clusters(pos_clusters,total_friends_indeces, max_epi_cont)
        epi_clusters <- merged_clusters$epi_clusters
        pos_clusters <- merged_clusters$pos_clusters
        indices_local <- pos_clusters[clusters %in% epi_clusters]$id
        
        # Simulate trials
        trials <- simulate_trial(sampling_sim, 1, ind_pos_rate$prevalence)
        pval <- length(which(trials >= (npos / ntotal))) / sampling_sim
        mean_prev <- mean(ind_pos_rate$prevalence)
        
      }else if(method == "centroid"){
        if(verbose){print("Using Centroid method to account for local prevalence.")}
        prev_indices <- compute_centroid(positions, total_friends_indeces,
                                     test_result,max_epi_cont, max_thr_data)
        
        # Obtain prevalence and local indices used for that
        ind_pos_rate <- prev_indices$prevalence
        indices_local <- prev_indices$local_id
        
        pval <- 1 - pbinom(npos - 1, ntotal, ind_pos_rate)
        mean_prev <- mean(ind_pos_rate)
      }else if(method == "base"){
        if(verbose){print("Accounting for global prevalence.")}
        total_positives = sum(test_result)
        ntotal <- length(total_friends_indeces)
        pos_rate <- total_positives/total_n
        pval <- 1 - pbinom(npos - 1, ntotal, pos_rate)
        
        # For base approach global indices == local indices
        indices_local <- positions$id
        mean_prev <- pos_rate
      }else{
        stop("None of the methods specified is valid. Please check the documentation.")
      }
    }
      
    # Adjust p-value based on: adj-pval = (1-p-val) * (1-p-fd)
    if(!is.null(false_det)){
      if(nrow(fp[num_pos ==npos]) != 0){
        prob_fd <- max(1-fp[num_pos ==npos]$prob_fd, 0)
        adj_pval <- (1-pval)*(prob_fd)
        pval <- 1 - adj_pval
      }
    }
    
    #setting EpiFRIenDs catalogue
    if(pval <= max_p && npos >= min_pos && ntotal >= min_total && mean_pr >= min_pr){
      epifriends_catalogue[['id']] <- append(epifriends_catalogue[['id']], next_id)
      cluster_id[cluster_id_indeces] <- next_id
      next_id = next_id + 1
      
      mean_pr_cluster[cluster_id_indeces] <- mean_pr
      pval_cluster[cluster_id_indeces] <- pval
      mean_pos <- as.data.table(t(colMeans(positive_positions[cluster_id_indeces,])))
      epifriends_catalogue[["mean_position_pos"]] <- append(epifriends_catalogue[["mean_position_pos"]], list(mean_pos))
      mean_pos_ext <-  as.data.table(t(colMeans(positions[total_friends_indeces,])))
      epifriends_catalogue[["mean_position_all"]] <- append(epifriends_catalogue[["mean_position_all"]], list(mean_pos_ext))
      epifriends_catalogue[["mean_pr"]] <- append(epifriends_catalogue[["mean_pr"]], mean_pr)
      epifriends_catalogue[["mean_local_prev"]] <- append(epifriends_catalogue[["mean_local_prev"]], mean_prev)
      epifriends_catalogue[["positives"]] <- append(epifriends_catalogue[["positives"]], as.integer(npos))
      epifriends_catalogue[["negatives"]] <- append(epifriends_catalogue[["negatives"]], as.integer(ntotal - npos))
      epifriends_catalogue[["total"]] <- append(epifriends_catalogue[["total"]], as.integer(ntotal))
      epifriends_catalogue[["indeces"]] <- append(epifriends_catalogue[["indeces"]], list(total_friends_indeces))
      epifriends_catalogue[["indeces_local_prev"]] <- append(epifriends_catalogue[["indeces_local_prev"]], list(indices_local))
      epifriends_catalogue[["p"]] <- append(epifriends_catalogue[["p"]], pval)
    }else{
      cluster_id[cluster_id_indeces] <- 0
    }
  }
  #Creation of list with return variables
  return <- list(cluster_id, mean_pr_cluster, pval_cluster, epifriends_catalogue)
  names(return) <- c("cluster_id", "mean_pr_cluster", "pval_cluster", "epifriends_catalogue")
  return(return)
}

#' Perform binomial simulations with different independent probabilities of success.
#'
#' @param n Number of observations.
#' @param size Number of trials (zero or more).
#' @param probs Vector of probabilities for each simulation to be performed.
#'
#' @return Ratio between the sum of each n and the total number of trials.
#'   
#' @export
#' 
#' @author Eric Matamoros.
#'
#' @examples
#' # Required packages
#' if(!require("RANN")) install.packages("RANN")
#' library("RANN")
#'
#' # Creation of x vector of longitude coordinates, y vector of latitude coordinates and finaly merge them on a position data frame.
#' n <- 10000
#' size <- 1
#' probs <- c(0.3, 0.4)
#'
#' # Do simulations
#' trials <- simulate_trial(n, size, probs)
#'
simulate_trial <- function(n, size, probs){

  # Create data.table with individual simulations
  df <- data.table("id" = 1:n)
  counter <- 1
  for(prob in probs){
    df[, (as.character(counter)) := rbinom(n = n, size = size, prob = prob)]
    counter <- counter +1
    }
  df[, id := NULL]
  
  # Compute ratio
  return(rowSums(df) / ncol(df))
}


