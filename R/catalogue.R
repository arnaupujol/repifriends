# Mikel Majewski Etxeberria(Ver 28-11-2022)

##################################################################################################
# CREATION OF CATALOGUE WITH TEST VALUES.
#################################################################################################
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
#' @param verbose If TRUE, print information of the process; else, do not print.
#'
#' @details The epifriends package uses the RANN package which can gives the exact nearest neighbours using the friends of friends algorithm. For more information on the RANN library please visit https://cran.r-project.org/web/packages/RANN/RANN.pdf
#'#'
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
#' @export
#' #'
#'
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

catalogue <- function(positions, test_result, link_d, cluster_id = NULL,
                      min_neighbours = 2, max_p = 1, min_pos = 2, min_total = 2,
                      min_pr = 0, keep_null_tests = FALSE, verbose = FALSE){
  # This method runs the DBSCAN algorithm (if cluster_id is NULL) and obtains the mean
  # positivity rate (PR) of each cluster extended with the non-infected cases
  # closer than the link_d.
  #
  # Parameters:
  #   -----------
  # positions: List of class data.frame
  #   A list with the position parameters we want to query with shape (n,2),
  #   where n is the number of positions.
  # test_result: List of class data.frame
  #   A list with the test results (0 or 1).
  # link_d: double
  #   The linking distance to connect cases.
  # cluster_id: numeric vector
  #   Vector of the cluster IDs of each position, with 0 for those
  #   without a cluster. Give NULL if cluster_id must be calculated.
  # min_neighbours: integer
  #   Minium number of neighbours in the radius < link_d needed to link cases
  #   as friends.
  # max_p: double
  #   Maximum value of the p-value to consider the cluster detection.
  # min_pos: integer
  #   Threshold of minimum number of positive cases in clusters applied.
  # min_total: integer
  #   Threshold of minimum number of cases in clusters applied.
  # min_pr: double
  #   Threshold of minimum positivity rate in clusters applied.
  # keep_null_tests: numeric of logical
  #   Whether to remove or not missings. If numeric, provide value to impute.
  # verbose: logical
  #   If TRUE, print information of the process; else, do not print.
  #
  # Returns:
  #   --------
  # return: list
  #   List which contains all the next objects. If there are not clusters an empty catalogue is returned with NULL values.
  # cluster_id: numeric vector
  #   Vector of the cluster IDs of each position, with 0 for those
  #   without a cluster.
  # mean_pr_cluster: numeric vector
  #   Mean PR corresponding to cluster_id.
  # pval_cluster: numeric vector
  #   P-value corresponding to cluster_id.
  # epifriends_catalogue: List
  #   Catalogue of the epifriends clusters and their main characteristics.
  
  # Remove or impute missings
  pos = clean_unknown_data(positions,test_result[[1]],keep_null_tests,verbose)
  positions = pos$position
  test_result = data.frame("test_result" = pos$test)

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
