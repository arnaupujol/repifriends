#This module contains methods to identify EpiFRIenDs clusters.

install.packages("RANN")  #KDTree method package for R. Documentation: https://cran.r-project.org/web/packages/RANN/RANN.pdf

library("RANN")
library("readxl")

datarand <- read.csv("mock_data_rand.csv")
datasin <- read.csv("mock_data_sin.csv")
dataclus <- read.csv("mock_data_clustered.csv")

position_rand <- datarand[3:4]

typeof(position_rand)

#Prouebas
dim(datarand)[2]
indecesp <- list()
indecesp[0] <- c()
position_rand[2:4,] <- position_rand[1,]
query <- nn2(position_rand, position_rand[1,], 4)
query$nn.idx[length(query$nn.idx)]

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
position_rand
last_cluster_id = 0
indeces_friends <- ind[[1]]
cluster_id = integer(dim(position_rand)[1])
cluster_id[ind[[1]]]
cluster_id_friends <- cluster_id[indeces_friends]
unique_cluster_ids <- unique(cluster_id_friends)
length(unique_cluster_ids)
cluster_id[indeces_friends] = last_cluster_id + 1
last_cluster_id = last_cluster_id + 1

db <- dbscan(position_rand, 0.05 ,2)