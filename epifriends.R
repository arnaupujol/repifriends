#This module contains methods to identify EpiFRIenDs clusters.

install.packages("RANN")  #KDTree method package for R. Documentation: https://cran.r-project.org/web/packages/RANN/RANN.pdf

library("RANN")
library("readxl")

datarand <- read.csv("mock_data_rand.csv")
datasin <- read.csv("mock_data_sin.csv")
dataclus <- read.csv("mock_data_clustered.csv")

position_rand <- datarand[3:4]

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
    # positions: np.ndarray
    #     An array with the position parameters with shape (n,2),
    #     where n is the number of positions
    # link_d: float
    #     The linking distance to label friends
    # tree: scipy.spatial.KDTree
    #     A KDTree build from the positions of the target data
    # 
    # Returns:
    # --------
    # indeces: list
    #     List with an array of the indeces of the friends of each
    #     position
  indeces <- list()
  for(i in 1:dim(positions)[1]){
    indeces[i] <- c()
    dist = 0
    kth = 0
    while(dist <= link_d){
      kth = kth + 1
      query <- nn2(tree, positions[i,], kth)
      index <- query$nn.idx[length(query$nn.idx)]
      dist <- query$nn.dists[length(query$nn.dists)]
      if(dist[length(dist)] <= link_d){
        append(indeces[i],index)
      }else{
        break
      }
    }
    return(indeces)
  }
  
}