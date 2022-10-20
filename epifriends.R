#This module contains methods to identify EpiFRIenDs clusters.

install.packages("RANN")  #KDTree method package for R. Documentation: https://cran.r-project.org/web/packages/RANN/RANN.pdf
install.packages("chron")

library("RANN")
library("readxl")
library("chron")

datarand <- read.csv("data\\mock_data_rand.csv")
datasin <- read.csv("data\\mock_data_sin.csv")
dataclus <- read.csv("data\\mock_data_clustered.csv")

position_rand <- datarand[3:4]
test <- datarand[5]

typeof(position_rand)

#Pruebas
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
    # positions: DataFrame
    #   A list with the position parameters we want to query with shape (n,2),
    #   where n is the number of positions
    # link_d: float
    #   The linking distance to label friends
    # tree: DataFrame
    #   A list build with the positions of the target data
    # 
    # Returns:
    # --------
    # indeces: list
    #   List with an array of the indeces of the friends of each
    #   position
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

ind <- find_indeces(position_rand, 0.05, position_rand)

dbscan <- function(positions, link_d, min_neighbours = 2){
  # This method finds the DBSCAN clusters from a set of positions and
  # returns their cluster IDs.
  # 
  # Parameters:
  #   -----------
  # positions: DataFrame
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
  # positions: DataFrame
  #   A list with the position parameters we want to query with shape (n,2),
  #   where n is the number of positions
  # test_result: DataFrame
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
  # return: list 
  #   List which contains all the values to be returned
  # cluster_id: list
  #   List of the cluster IDs of each position with positive test, with 0 for those
  #   without a cluster.
  # mean_pr_fof: list
  #   Mean PR corresponding to cluster_id
  # pval_fof: list
  #   P-value corresponding to cluster_id
  # epifriends_catalogue: List of lists
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
  #Case where there are not positive positions
  if(dim(positive_positions)[1] == 0){
    return(epifriends_catalogue)
  }
  
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

datarandtemp <- read.csv("data\\mock_data_rand_temp.csv")

position_randtemp <- datarandtemp[3:4]
testtemp <- datarandtemp[5]
datetemp <- datarandtemp[7]

distance <- function(pos_a, pos_b){
    # This method calculates the Euclidean distance between two positions.
    # 
    # Parameters:
    # -----------
    # pos_a: vector
    #     First position
    # pos_b: vector
    #     Second position
    # 
    # Returns:
    # --------
    # dist: float
    #     Distance between positions
  dist = sqrt(sum((pos_a - pos_b)**2))
  return(dist)
}

a <- c(1,2)
b <- c(4,1)

distance(a,b)

temporal_catalogue <- function(positions, test_result, dates, link_d, min_neighbours = 2,
                   time_width, min_date = NULL, max_date = NULL, time_steps = 1,
                   max_p = 1, min_pos = 2, min_total = 2, min_pr = 0){
  # This method generates a list of EpiFRIenDs catalogues representing different time frames
  # by including only cases within a time window that moves within each time step.
  # 
  # Parameters:
  #   -----------
  # positions: DataFrame
  #   A list with the position parameters we want to query with shape (n,2),
  #   where n is the number of positions
  # test_result: list
  #   A list with the test results (0 or 1)
  # dates: list
  #   List of the date times of the corresponding data
  # link_d: float
  #   The linking distance to connect cases
  # min_neighbours: int
  #   Minium number of neighbours in the radius < link_d needed to link cases
  #   as friends
  # time_width: int
  #   Number of days of the time window used to select cases in each time step
  # min_date: pd.DateTimeIndex
  #   Initial date used in the first time step and time window selection
  # max_date: pd.DateTimeIndex
  #   Final date to analyse, defining the last time window as the one fully overlapping
  #   the data
  # time_steps: int
  #   Number of days that the time window is shifted in each time step
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
  # returns: list
  #   List that contais the two return variables
  # temporal_catalogues: list of catalogues
  #   List of EpiFRIenDs catalogues, where each element contains the catalogue in each
  #   time step
  # mean_date: list
  #   List of dates corresponding to the median time in each time window
  #   Save dates in a temporal format
  #Case of empty list of dates
  if(length(dates) == 0){
    print("There are no dates")
    return()
  }
  dtparts = t(as.data.frame(strsplit(dates," ")))
  row.names(dtparts) = NULL
  dates <- chron(dates=dtparts[,1],times=dtparts[,2],format=c('y-m-d','h:m:s'))
  #Defining temporal range
  if(is.null(min_date)){
    min_date <- min(dates)
  } 
  if(is.null(max_date)){
    max_date <- max(dates)
  }
  #temporal loop until the last time frame that fully overlaps the data
  temporal_catalogues <- list()
  #Mean dates defines as the median time in each time window
  mean_date <- list()
  step_num <- 0
  while(min_date + time_steps*step_num + time_width <= max_date){
    #select data in time window
    selected_data <- (dates >= min_date + time_steps*step_num)&
                     (dates <= min_date + time_steps*step_num + time_width)
    selected_positions <- positions[selected_data,]
    selected_test_results <- as.data.frame(test_result[selected_data,])
    
    #get catalogue
    Newcatalogue <- catalogue(selected_positions, selected_test_results,
                                     link_d, min_neighbours = min_neighbours,
                                     max_p = max_p, min_pos = min_pos,
                                     min_total = min_total, min_pr = min_pr)
    #we only add the mean date if the catalogue is not empty
    if(FALSE %in% lapply(Newcatalogue$epifriends_catalogue, is.null)){
      mean_date <- append(mean_date, min_date + time_steps*step_num + 0.5*time_width)
      Newcatalogue$epifriends_catalogue["Date"] <- toString(min_date + time_steps*step_num + 0.5*time_width)
      temporal_catalogues[[step_num+1]] <- Newcatalogue$epifriends_catalogue
    }
    step_num = step_num + 1
  }
  returns <- list(temporal_catalogues, mean_date)
  names(returns) <- c("temporal_catalogues", "mean_date")
  return(returns)
}

tcat <- temporal_catalogue(position_randtemp, testtemp, datarandtemp$date, 0.05, time_width = 180, time_steps = 90)
Newcatalogue <- catalogue(selected_positions, selected_test_results,0.05)

tcat$temporal_catalogues[[2]]

#Proves
# dtimes <- datarandtemp$date
# dtparts = t(as.data.frame(strsplit(dtimes," ")))
# row.names(dtparts) = NULL
# dates = chron(dates=dtparts[,1],times=dtparts[,2],format=c('y-m-d','h:m:s'))
# 
# min_date <- min(dates)
# max_date <- max(dates)
# 
# temporal_catalogues <- list()
# mean_date <- list()
# 
# min(dates) + 10*70 + 31 <= max(dates)
# max(dates) - min(dates)
# selected_data <- (dates >= min_date)&(dates <= min_date + 180)
# selected_positions <- position_randtemp[selected_data,]
# selected_test_results <- as.data.frame(testtemp[selected_data,])
# Newcatalogue <- catalogue(selected_positions, selected_test_results,
#                           0.05, min_neighbours = 2,
#                           max_p = 1, min_pos = 2,
#                           min_total = 2, min_pr = 0)
# mean_date <- append(mean_date, min_date + 0.5*180)
# Newcatalogue$epifriends_catalogue["Date"] <- min_date + 0.5*180
# chron(18892.42,format=c('y-m-d','h:m:s'))
#
# temporal_catalogues <- append(temporal_catalogues, Newcatalogue$epifriends_catalogue)

get_label_list <- function(df_list, label = "tempID"){
  # This method gives the unique values of a column in a list
  # of data frames.
  # 
  # Parameters:
  #   -----------
  # df_list: list
  #   List of dataframes
  # label: str
  #   Name of column to select
  # 
  # Returns:
  #   --------
  # label_list: list
  #   List of unique values of the column over all dataframes from the list
  for(i in 1:length(df_list)){
    mask = df_list[[i]][label][[1]][df_list[[i]][label][[1]] != 0]
    if(i == 1){
      label_list <- unique(mask)
    }else{
      label_list <- unique(c(label_list, unique(mask)))
    }
  }
  return(label_list)
}

lablis <- get_label_list(tcatid)

a <- unique(tcatid[[1]]["tempID"][[1]][tcatid[[1]]["tempID"][[1]] != 0])
b <- unique(tcatid[[2]]["tempID"][[1]][tcatid[[2]]["tempID"][[1]] != 0])
unique(c(a,b))

get_lifetimes <- function(catalogue_list){
  # This method obtains the first and last time frames for each
  # temporal ID from a list of EpiFRIenDs catalogues and the corresponding
  # timelife.
  # 
  # Parameters:
  # -----------
  # catalogue_list: list
  #     List of EpiFRIenDs catalogues, each element of the list
  #     corresponding to the EpiFRIenDs catalogue of each timestep
  # 
  # Returns:
  # --------
  # catalogue_list: list
  #     List of hotspot catalogues with the added fields 'first_timestep',
  #     'last_timestep' and 'lifetime'
  
  #getting list of temporal IDs appearing in catalogue_list
  tempid_list <- get_label_list(catalogue_list, label = "tempID")
  #Creating empty columns for first timestep, last timestep and lifteime
  for(t in 1:length(catalogue_list)){
    catalogue_list[[t]]["first_timestep"] <- c()
    catalogue_list[[t]]["last_timestep"] <- c()
    catalogue_list[[t]]["lifetime"] <- c()
  }
  for(tempid_num in tempid_list){
    appearances = c()
    for(i in 1:length(catalogue_list)){
      if(tempid_num %in% unique(catalogue_list[[i]]["tempID"][[1]])){
        appearances <- append(appearances,i)
      }
    }
    min_appearance = min(appearances)
    max_appearance = max(appearances)
    lifetime = max_appearance - min_appearance
    for(i in min_appearance:max_appearance){
      catalogue_list[[i]]["first_timestep"][[1]][catalogue_list[[i]]["tempID"][[1]] == tempid_num] <- min_appearance
      catalogue_list[[i]]["last_timestep"][[1]][catalogue_list[[i]]["tempID"][[1]] == tempid_num] <- max_appearance
      catalogue_list[[i]]["lifetime"][[1]][catalogue_list[[i]]["tempID"][[1]] == tempid_num] <- lifetime
    }
  }
  return(catalogue_list)
}

catlif <- get_lifetimes(tcatid)
catlif

# # c <- c()
# # c[1] <- 1
# # c[3] <- 2
# 
# catlif[[1]]["id"][[1]][catlif[[1]]["tempID"][[1]] == 1]
# 
# lablis[1] %in% unique(tcatid[[4]]["tempID"][[1]])

add_temporal_id <- function(catalogue_list, linking_time, linking_dist, get_timelife = TRUE){
  #Case of empty catalogue list
  if(length(catalogue_list) == 0){
    print("There are no catalogues")
    return()
  }
  #setting empty values of temp_id
  for(t in 1:length(catalogue_list)){
    print(catalogue_list[[t]]$id)
    aux <- data.frame(matrix(0,length(catalogue_list[[t]]$id)))
    colnames(aux) <- "tempID"
    catalogue_list[[t]] <- append(catalogue_list[[t]],aux)
    #catalogue_list[[t]]["tempID"] = vector(mode="list", length=length(catalogue_list[[t]]$id))
  }
  #Initialising tempID value to assign
  next_temp_id = 0
  #Loop over all timesteps
  for(t in 1:(length(catalogue_list)-1)){
    #Loop over all timesteps within linking_time
    for (t2 in (t + 1):min(t + linking_time, length(catalogue_list))){
      #Loop over all points of catalogue number 1
      for(f in 1:length(catalogue_list[[t]]$id)){
        #Loop over all points of catalogue number 2
        for(f2 in 1:length(catalogue_list[[t2]]$id))
          #Calculating distance between clusters
          dist <- distance(catalogue_list[[t]]["mean_position_pos"][[1]][[f]], catalogue_list[[t2]]["mean_position_pos"][[1]][[f2]])  #To improve. Better not to use [[1]]
        if(dist <= linking_dist){
          temp_id1 <- catalogue_list[[t]]["tempID"][[1]][[f]] 
          temp_id2 <- catalogue_list[[t2]]["tempID"][[1]][[f2]] 
          #Assign tempIDs to linked clusters
          if((temp_id1 == 0) && (temp_id2 == 0)){
            catalogue_list[[t]]["tempID"][[1]][[f]]  <- next_temp_id + 1
            catalogue_list[[t2]]["tempID"][[1]][[f2]] <- next_temp_id + 1
            next_temp_id = next_temp_id + 1
          }else if((temp_id1 == 0)){
            catalogue_list[[t]]["tempID"][[1]][[f]]  <- temp_id2
          }else if((temp_id2 == 0)){
            catalogue_list[[t2]]["tempID"][[1]][[f2]] <- temp_id1 
          }else if(temp_id1 != temp_id2){
            for(t3 in 1:length(catalogue_list)){
              catalogue_list[[t3]]["tempID"][catalogue_list[[t3]]["tempID"] == temp_id2] <- temp_id1
            }
          }
        }
      }
    }
  }
  if(get_timelife){
    catalogue_list <- get_lifetimes(catalogue_list)
  }
  return(catalogue_list)
}

tcatid <- add_temporal_id(tcat$temporal_catalogues, 3, 0.15, get_timelife = TRUE)
tcatid

# tcat$temporal_catalogues[[1]]["tempID"][[1]][1]
# t(tcat$temporal_catalogues[1])[[1]]$mean_position_pos[[1]]
# length(tcat$temporal_catalogues[[1]]$id)
# length(tcat$temporal_catalogues)

#Pruebas replica add_temporal_id
# for(t in 1:length(tcat$temporal_catalogues)){
#   tcat$temporal_catalogues[1]["tempID"] = vector(mode="list", length=length(tcat$temporal_catalogues[[1]]["id"]))
# }
# 
# for(t in 1:length(tcat$temporal_catalogues)){
#   tcat$temporal_catalogues[1]["tempID"] = vector(mode="list", length=length(tcat$temporal_catalogues[[1]]["id"]))
# }

# for(t in 1:length(tcat$temporal_catalogues)){
#   aux <- data.frame(matrix(NA,length(tcat$temporal_catalogues[[t]]$id)))
#   colnames(aux) <- "tempID"
#   tcat$temporal_catalogues[[t]] <- append(tcat$temporal_catalogues[[t]],aux)
#   #catalogue_list[[t]]["tempID"] = vector(mode="list", length=length(catalogue_list[[t]]$id))
# }

#Initialising tempID value to assign

#tcat$temporal_catalogues[[1]]["mean_position_pos"][[1]][[8]]

# c <- c()
# c[1] <- 1
# c[3] <- 2

catlif[[1]]["id"][[1]][catlif[[1]]["tempID"][[1]] == 1]

lablis[1] %in% unique(tcatid[[4]]["tempID"][[1]])


#VALIDATIONS

x <- c(1,2,3,4,7.5,8,8.5,9,10,13,13.1,13.2,13.3,14,15,16)
y <- c(1,2,3,4,7.5,8,8.5,9,10,13,13.1,13.2,13.3,14,15,16)
length(x)
pos <- data.frame(x,y)

ind <- find_indeces(pos, 2 ,pos)
ind

db <- dbscan(pos, 2 ,2)
db

db <- dbscan(pos, 10 ,2)
db

db <- dbscan(pos, 2 ,5)
db

db <- dbscan(pos, 2 ,7)
db

test <- data.frame(c(0,1,1,0,1,0,1,1,0,0,0,0,1,0,0,1))
length(test)
cat <- catalogue(pos, test, 2)
cat

dates <- c("2020-11-03 05:33:07","2021-05-19 10:29:59","2021-02-09 14:53:20","2021-11-21 02:35:38","2020-11-19 05:57:24",
           "2021-06-09 07:50:30","2021-09-18 05:53:26","2020-03-19 17:16:56","2021-06-08 12:40:46","2020-06-26 05:01:31",
           "2020-10-15 04:40:27","2021-05-28 15:23:23","2020-05-01 02:56:54","2020-08-19 22:45:35","2021-10-23 18:56:35",
           "2020-10-19 00:01:25")

typeof(dates)
typeof(datarandtemp$date)
tcat <- temporal_catalogue(pos, test, dates ,2, time_width = 180, time_steps = 90)
tcat$temporal_catalogues
length(tcat$temporal_catalogues)

tcatid <- add_temporal_id(tcat$temporal_catalogues, 3, 0.15, get_timelife = FALSE)
tcatid

#VALIDATION with no positions
x <- c()
y <- c()
pos <- data.frame(x,y)

ind <- find_indeces(pos, 2 ,pos)
ind
#ind[1,]
sort(unique(unlist(ind)))

db <- dbscan(pos, 2 ,2)
db

test <- data.frame(c())
cat <- catalogue(pos, test, 2)
cat

dates <- c()
tcat <- temporal_catalogue(pos, test, dates ,2, time_width = 180, time_steps = 90)

tcatid <- add_temporal_id(tcat$temporal_catalogues, 3, 0.15, get_timelife = TRUE)
tcatid

#VALIDATION with only one position
x <- c(1)
y <- c(1)
pos <- data.frame(x,y)

ind <- find_indeces(pos, 2 ,pos)
ind

db <- dbscan(pos, 2 ,2)
db

test <- data.frame(c(1))
cat <- catalogue(pos, test, 2)
cat

dates <- c("2020-11-03 05:33:07")
tcat <- temporal_catalogue(pos, test, dates ,2, time_width = 180, time_steps = 90)
tcat

tcatid <- add_temporal_id(tcat$temporal_catalogues, 3, 0.15, get_timelife = TRUE)
tcatid

#VALIDATION with two positions
x <- c(1,3)
y <- c(1,3)
pos <- data.frame(x,y)

ind <- find_indeces(pos, 4 ,pos)
ind

db <- dbscan(pos, 4 ,2)
db

test <- data.frame(c(1,0))
cat <- catalogue(pos, test, 4)
cat

dates <- c("2020-11-03 05:33:07", "2021-05-19 10:29:59")
tcat <- temporal_catalogue(pos, test, dates ,2, time_width = 180, time_steps = 90)
tcat$temporal_catalogues

tcatid <- add_temporal_id(tcat$temporal_catalogues, 3, 0.15, get_timelife = TRUE)
tcatid
#VALIDATION with repetition of positions

x <- c(1,1,1,3)
y <- c(1,1,1,3)
pos <- data.frame(x,y)

ind <- find_indeces(pos, 4 ,pos)
ind

ind <- find_indeces(pos, 2 ,pos)
ind

db <- dbscan(pos, 4 ,2)
db

db <- dbscan(pos, 2 ,2)
db

test <- data.frame(c(1,1,1,0))
cat <- catalogue(pos, test, 4)
cat

dates <- c("2020-11-03 05:33:07","2020-11-03 05:33:07","2020-11-03 05:33:07","2021-05-19 10:29:59")
tcat <- temporal_catalogue(pos, test, dates ,2, time_width = 180, time_steps = 90)
tcat$temporal_catalogues  

tcatid <- add_temporal_id(tcat$temporal_catalogues, 3, 0.15, get_timelife = TRUE)
tcatid
