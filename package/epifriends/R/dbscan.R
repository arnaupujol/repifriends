# Mikel Majewski Etxeberria (Ver 28-11-2022)

##################################################################################################
# CLUSTERING POINTS USING FRIENDS OF FRIENDS ALGORITHM
#################################################################################################
#' This method finds the DBSCAN clusters from a set of positions and returns their cluster IDs.
#'
#' @param positions data.frame with the positions of parameters we want to query with shape (n,2) where n is the number of positions.
#' @param link_d: The linking distance to connect cases. Should be in the same scale as the positions.
#' @param min_neighbours: Minium number of neighbours in the radius < link_d needed to link cases as friends.
#'
#' @details The epifriends package uses the RANN package which can gives the exact nearest neighbours using the friends of friends algorithm. For more information on the RANN library please visit https://cran.r-project.org/web/packages/RANN/RANN.pdf
#'
#' @return cluster_id: Vector of the cluster IDs of each position, with 0 for those without a cluster. Returns empty numeric vector if positions vector is empty.
#' @export
#' #'
#'
#' @author Mikel Majewski Etxeberria based on earlier python code by Arnau Pujol.
#'
#' @examples
#' Required packages
#' if(!require("RANN")) install.packages("RANN")
#' library("RANN")
#'
#' # Creation of x vector of longitude coordinates, y vector of latitude coordinates and finaly merge them on a position data frame.
#' x <- c(1,2,3,4,7.5,8,8.5,9,10,13,13.1,13.2,13.3,14,15,30)
#' y <- c(1,2,3,4,7.5,8,8.5,9,10,13,13.1,13.2,13.3,14,15,30)
#' pos <- data.frame(x,y)
#'
#' # Computation of clusters of hotspots for positions with dbscan algorithm using linking distance 2 and minimum 3 neighbours.
#' db <- dbscan(pos, 2 ,3)

dbscan <- function(positions, link_d, min_neighbours = 2){
  # This method finds the DBSCAN clusters from a set of positions and
  # returns their cluster IDs.
  #
  # Parameters:
  #   -----------
  # positions: List of class data.frame
  #   A list with the position parameters we want to query with shape (n,2),
  #   where n is the number of positions.
  # link_d: double
  #   The linking distance to connect cases.
  # min_neighbours: integer
  #   Minium number of neighbours in the radius < link_d needed to link cases
  #   as friends.
  #
  # Returns:
  #   --------
  # cluster_id: numeric vector
  #   Vector of the cluster IDs of each position, with 0 for those
  #   without a cluster. Returns empty numeric vector if positions vector is empty.

  #Create cluster id
  cluster_id <- integer(dim(positions)[1])
  #Query KDTree
  indeces <- find_indeces(positions, link_d, positions)
  #inicialize ID variable
  last_cluster_id = 0
  #Have in count the posible exception of not having any positions
  if(dim(positions)[1] == 0){
    return(cluster_id)
  }
  for(i in 1:dim(positions)[1]){
    #check if ith position has the minimum neighbours
    if(length(indeces[[i]]) >= min_neighbours){
      indeces_friends <- indeces[[i]]
      #cluster_ids of these friends
      cluster_id_friends <- cluster_id[indeces_friends]
      #Unique values of cluster_ids
      unique_cluster_ids <- unique(cluster_id_friends)
      #check values of cluster_id in these neighbours
      if(length(unique_cluster_ids) == 1){
        if(unique_cluster_ids[1] == 0){
          #assign to ith and friends last_cluster_id
          cluster_id[indeces_friends] = last_cluster_id + 1
          last_cluster_id = last_cluster_id + 1
        }else{
          #if one cluster_id different than 0, assign it to ith and friends
          cluster_id[indeces_friends] = unique_cluster_ids[1]
        }
      }else{
        #Define the cluster_id to assign for merging several clusters
        min_cluster_id = min(unique_cluster_ids[unique_cluster_ids != 0])
        #Assign this cluster_id to ith and its friends
        cluster_id[indeces_friends] = min_cluster_id
        for(j in unique_cluster_ids[unique_cluster_ids != 0]){
          cluster_id[cluster_id == j] = min_cluster_id
        }
      }
    }
  }
  #Rename cluster_id to continuous integers
  sort_unici <- sort(unique(cluster_id[cluster_id>0]))
  for(i in 1:length(sort_unici)){
    cluster_id[cluster_id == sort_unici[i]] = i
  }
  return(cluster_id)
}

find_indeces <- function(positions, link_d, tree){
  # This method returns the indeces of all the friends
  # of each position from positions given a KDTree.
  #
  # Parameters:
  # -----------
  # positions: List of class data.frame
  #   A list with the position parameters we want to query with shape (n,2),
  #   where n is the number of positions.
  # link_d: double
  #   The linking distance to label friends.
  # tree: List of class data.frame
  #   A list build with the positions of the target data.
  #
  # Returns:
  # --------
  # indeces: list
  #   List with an array of the indeces of the friends of each
  #   position.
  #Creation of empty list where the clusters of points will be saved
  indeces <- list()
  #Have in count the posible exception of not having any positions
  if(dim(positions)[1] == 0){
    return(indeces)
  }
  #Case if there is only one position
  if(dim(positions)[1] == 1){
    indeces[[1]] <- 1
    return(indeces)
  }
  #loop for each position
  for(i in 1:dim(positions)[1]){
    #creatios of list where the linked positions of the given position will be saved
    indecesaux <- c()
    #inicialitation of distance and k number of nearest positions
    dist = 0
    kth = 0
    #loop which stops when the maximum linking distance is overcomed
    while((dist <= link_d) && (kth < dim(tree)[1])){
      #Aplication of KDTree method for the k nearest neighbors while the linking distance is not overcomed
      kth = kth + 1
      query <- nn2(tree, positions[i,], kth)
      index <- query$nn.idx[length(query$nn.idx)]
      dist <- query$nn.dists[length(query$nn.dists)]
      #Addition of the last point to the indeces list if it is in the wanted range
      if((dist[length(dist)] <= link_d) && (kth < dim(tree)[1])){
        indecesaux <- append(indecesaux,index)
      }else{
        #Final list for the position i
        indeces[[i]] <- indecesaux
        break
      }
    }
  }
  return(indeces)
}
