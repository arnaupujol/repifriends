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
