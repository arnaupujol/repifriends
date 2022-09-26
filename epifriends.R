#This module contains methods to identify EpiFRIenDs clusters.

install.packages("RANN")  #KDTree method package for R. Documentation: https://cran.r-project.org/web/packages/RANN/RANN.pdf

library("RANN")
library("readxl")

datarand <- read.csv("mock_data_rand.csv")
datasin <- read.csv("mock_data_sin.csv")
dataclus <- read.csv("mock_data_clustered.csv")

position_rand <- datarand[3:4]
test <- datarand[5]

typeof(position_rand)

#Prouebas
# dim(datarand)[2]
# indecesp <- list()
# indecesp[0] <- c()
# position_rand[2:4,] <- position_rand[1,]
# query <- nn2(position_rand, position_rand[1,], 4)
# query$nn.idx[length(query$nn.idx)]

#Functions

find_indeces <- function(positions, link_d, tree){
    # This method returns the indeces of all the friends
    # of each position from positions given a KDTree.
    # 
    # Parameters:
    # -----------
    # positions: list
    #   A list with the position parameters we want to query with shape (n,2),
    #   where n is the number of positions
    # link_d: float
    #   The linking distance to label friends
    # tree: list
    #   A list build with the positions of the target data
    # 
    # Returns:
    # --------
    # indeces: list
    #   List with an array of the indeces of the friends of each
    #   position
  #Creation of empty list where the clusters of points will be saved
  indeces <- list()
  #loop for each position
  for(i in 1:dim(positions)[1]){
    #creatios of list where the linked positions of the given position will be saved
    indecesaux <- c()
    #inicialitation of distance and k number of nearest positions
    dist = 0
    kth = 0
    #loop which stops when the maximum linking distance is overcomed
    while(dist <= link_d){
      #Aplication of KDTree method for the k nearest neighbors while the linking distance is not overcomed
      kth = kth + 1
      query <- nn2(tree, positions[i,], kth)
      index <- query$nn.idx[length(query$nn.idx)]
      dist <- query$nn.dists[length(query$nn.dists)]
      #Addition of the last point to the indeces list if its in the wanted range
      if(dist[length(dist)] <= link_d){
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

ind <- find_indeces(position_rand, 0.01, position_rand)

dbscan <- function(positions, link_d, min_neighbours = 2){
  # This method finds the DBSCAN clusters from a set of positions and
  # returns their cluster IDs.
  # 
  # Parameters:
  #   -----------
  # positions: list
  #   A list with the position parameters we want to query with shape (n,2),
  #   where n is the number of positions
  # link_d: float
  #   The linking distance of the DBSCAN algorithm
  # min_neighbours: int
  #   Minium number of neighbours in the radius < link_d needed to link cases
  #   as friends
  # 
  # Returns:
  #   --------
  # cluster_id: list
  #   List of the cluster IDs of each position, with 0 for those
  #   without a cluster.
  
  #Create cluster id
  cluster_id <- integer(dim(positions)[1])
  #Query KDTree
  indeces <- find_indeces(positions, link_d, positions)
  #inicialize ID variable
  last_cluster_id = 0
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

#Pruebas
# position_rand
# last_cluster_id = 0
# indeces_friends <- ind[[1]]
# cluster_id = integer(dim(position_rand)[1])
# cluster_id[ind[[1]]]
# cluster_id_friends <- cluster_id[indeces_friends]
# unique_cluster_ids <- unique(cluster_id_friends)
# length(unique_cluster_ids)
# cluster_id[indeces_friends] = last_cluster_id + 1
# last_cluster_id = last_cluster_id + 1

db <- dbscan(position_rand, 0.05 ,2)

#Catalogue

catalogue <- function(positions, test_result, link_d, cluster_id = NULL,
                      min_neighbours = 2, max_p = 1, min_pos = 2, min_total = 2,
                      min_pr = 0){
  # This method runs the DBSCAN algorithm (if cluster_id is None) and obtains the mean
  # positivity rate (PR) of each cluster extended with the non-infected cases
  # closer than the link_d.
  # 
  # Parameters:
  #   -----------
  # positions: list
  #   A list with the position parameters we want to query with shape (n,2),
  #   where n is the number of positions
  # test_result: list
  #   A list with the test results (0 or 1)
  # link_d: float
  #   The linking distance to connect cases
  # cluster_id: list
  #   A list with the cluster ids of the positive cases
  # min_neighbours: int
  #   Minium number of neighbours in the radius < link_d needed to link cases
  #   as friends
  # max_p: float
  #   Maximum value of the p-value to consider the cluster detection
  # min_pos: int
  #   Threshold of minimum number of positive cases in clusters applied
  # min_total: int
  #   Threshold of minimum number of cases in clusters applied
  # min_pr: float
  #   Threshold of minimum positivity rate in clusters applied
  # 
  # Returns:
  #   --------
  # cluster_id: list
  #   List of the cluster IDs of each position, with 0 for those
  #   without a cluster.
  # mean_pr_fof: list
  #   Mean PR corresponding to cluster_id
  # pval_fof: list
  #   P-value corresponding to cluster_id
  # epifriends_catalogue: geopandas.DataFrame
  #   Catalogue of the epifriends clusters and their main characteristics
  
  #Define positions of positive cases
  positive_positions <- positions[test_result == 1,]
  #Computing cluster_id if needed
  if(is.null(cluster_id)){
    cluster_id = dbscan(positive_positions, link_d, min_neighbours = min_neighbours)
  }
  #Define total number of positive cases
  total_positives = sum(test_result)
  #Define total number of cases
  total_n = dim(test_result)[1]
  #Initialising mean PR and p-value for the positive cases in clusters
  mean_pr_cluster <- numeric(length(cluster_id))                       #!!!!!!!!!!!!!
  pval_cluster <- rep(1, length(cluster_id))                           #!!!!!!!!!!!!!
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
    if(pval < max_p && npos >= min_pos && ntotal >= min_total && mean_pr >= min_pr){
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
  #Creation of list of return variables
  return <- list(cluster_id, mean_pr_cluster, pval_cluster, epifriends_catalogue)
  names(return) <- c("cluster_id", "mean_pr_cluster", "pval_cluster", "epifriends_catalogue")
  return(return)
}

cat <- catalogue(position_rand, test, 0.05)

#Pruebas
# positive_positions <- position_rand[test == 1,]
# cluster_id = dbscan(positive_positions, 0.01, min_neighbours = 2)
# total_positives = sum(test)
# total_n = dim(test)[1]
# epifriends_catalogue <- vector(mode='list', length=9)
# names(epifriends_catalogue) <- c("id", "mean_position_pos", "mean_position_all", "mean_pr", "positives", "negatives", "total", "indeces", "p")
# cluster_id_indeces = which(cluster_id == sort_unici[1])
# all_friends_indeces = find_indeces(positive_positions[cluster_id_indeces,], 0.05, position_rand)
# total_friends_indeces <- unique(unlist(all_friends_indeces))
# mean_pr <- mean(test[total_friends_indeces,])
# npos <- sum(test[total_friends_indeces,])
# total <- length(total_friends_indeces)
append(c(c(1,2)), c(c(2,3)))
