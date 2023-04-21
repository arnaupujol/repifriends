# This source code is licensed under the GNU GENERAL PUBLIC LICENSE license found in the
# LICENSE file in the root directory of this source tree.

#' This method runs the DBSCAN algorithm (if cluster_id is NULL) and obtains the mean positivity rate (PR) of each cluster extended with the non-infected cases closer than the link_d.
#'
#' @param positions data.frame with the positions of parameters we want to query with shape (n,2) where n is the number of positions.
#' @param test_result data.frame with the test results (0 or 1).
#' @param link_d The linking distance to connect cases. Should be in the same scale as the positions.
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
#' pos <- data.frame(x,y)
#'
#' # Creation of test data frame with 0 for negative cases and 1 for positive clases for each position.
#' test <- data.frame(c(0,1,1,0,1,0,1,1,0,0,0,0,1,0,1,1))
#'
#' # Creation of catalogue for this positions, linking distance 2 and default values.
#' cat <- catalogue(pos, test, 2)
#' 
catalogue <- function(positions, test_result, link_d=NULL, cluster_id = NULL,
                      min_neighbours = 2, max_p = 1, min_pos = 2, min_total = 2,
                      min_pr = 0, keep_null_tests = FALSE, in_latlon = FALSE, 
                      to_epsg = NULL, optimize_link_d = FALSE, verbose = FALSE){
  
  # Remove or impute missings
  pos = clean_unknown_data(positions,test_result[[1]],keep_null_tests,verbose)
  positions = pos$position
  test_result = data.frame("test_result" = pos$test)
  
  #Defining 2d-positions
  positions = get_2dpositions(
    x = positions$x, 
    y = positions$y, 
    in_latlon = in_latlon, 
    to_epsg = to_epsg,
    verbose = verbose)
  
  if(optimize_link_d | is.null(link_d)){
    if(verbose){print("Automating the calculus of the linking distance")}
    link_d <- opt_link_d(df=data.table(x = positions$x, y = positions$y, test = pos$test), 
                         min_neighbors=min_neighbors, cluster_id=cluster_id, 
                         keep_null_tests = keep_null_tests, in_latlon = in_latlon, 
                         to_epsg = to_epsg, verbose = verbose)
  }

  #Define positions of positive cases
  positive_positions <- positions[which(test_result == 1),]
  #Computing cluster_id if needed
  if(is.null(cluster_id)){
    cluster_id = dbscan(positive_positions, link_d, min_neighbours = min_neighbours)
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
  for(i in 1:length(sort_unici)){
    #get all indeces with this cluster id
    cluster_id_indeces <- which(cluster_id == sort_unici[i])
    #for all these indeces, get list of friends from all positions
    all_friends_indeces <- find_indeces(positive_positions[cluster_id_indeces,], link_d, positions)
    #get unique values of such indeces
    total_friends_indeces <- sort(unique(unlist(all_friends_indeces)))
    #get positivity rate from all the unique indeces
    mean_pr <- mean(test_result[total_friends_indeces,])
    npos <- sum(test_result[total_friends_indeces,])
    ntotal <- length(total_friends_indeces)
    pval <- 1 - pbinom(npos - 1, ntotal, total_positives/total_n)
    #setting EpiFRIenDs catalogue
    if(pval <= max_p && npos >= min_pos && ntotal >= min_total && mean_pr >= min_pr){
      epifriends_catalogue[['id']] <- append(epifriends_catalogue[['id']], next_id)
      cluster_id[cluster_id_indeces] <- next_id
      next_id = next_id + 1

      mean_pr_cluster[cluster_id_indeces] <- mean_pr
      pval_cluster[cluster_id_indeces] <- pval
      mean_pos <- colMeans(positive_positions[cluster_id_indeces,])
      epifriends_catalogue[["mean_position_pos"]] <- append(epifriends_catalogue[["mean_position_pos"]], list(mean_pos))
      mean_pos_ext <-  colMeans(positions[total_friends_indeces,])
      epifriends_catalogue[["mean_position_all"]] <- append(epifriends_catalogue[["mean_position_all"]], list(mean_pos_ext))
      epifriends_catalogue[["mean_pr"]] <- append(epifriends_catalogue[["mean_pr"]], mean_pr)
      epifriends_catalogue[["positives"]] <- append(epifriends_catalogue[["positives"]], as.integer(npos))
      epifriends_catalogue[["negatives"]] <- append(epifriends_catalogue[["negatives"]], as.integer(ntotal - npos))
      epifriends_catalogue[["total"]] <- append(epifriends_catalogue[["total"]], as.integer(ntotal))
      epifriends_catalogue[["indeces"]] <- append(epifriends_catalogue[["indeces"]], list(total_friends_indeces))
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
