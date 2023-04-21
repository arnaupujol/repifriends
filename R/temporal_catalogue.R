# This source code is licensed under the GNU GENERAL PUBLIC LICENSE license found in the
# LICENSE file in the root directory of this source tree.

#' This method generates a list of EpiFRIenDs catalogues representing different time frames by including only cases within a time window that moves within each time step.
#'
#' @param positions data.frame with the positions of parameters we want to query with shape (n,2) where n is the number of positions.
#' @param test_result data.frame with the test results (0 or 1).
#' @param dates Vector of the date times of the corresponding data in format 'y-m-d h:m:s'.
#' @param use_link_d: If True, use linking distance to determine the closest neighbours. If False, use default linking neighbours based on proximity.
#' @param min_neighbours: Minium number of neighbours in the radius < link_d needed to link cases as friends.
#' @param link_neighbours: Number of surrounding neighbors to link. 
#' @param time_width Number of days of the time window used to select cases in each time step.
#' @param min_date Initial date used in the first time step and time window selection. In format 'y-m-d h:m:s'.
#' @param max_date Final date to analyse, defining the last time window as the one fully overlapping the data. In format 'y-m-d h:m:s'.
#' @param time_steps Number of days that the time window is shifted in each time step.
#' @param max_p Maximum value of the p-value to consider the cluster detection.
#' @param min_pos Threshold of minimum number of positive cases in clusters applied.
#' @param min_total Threshold of minimum number of cases in clusters applied.
#' @param min_pr Threshold of minimum positivity rate in clusters applied.
#' @param add_temporal_id Boolean that indicates if we want to add a temporal id to the hotspots indicating which hotspots are near in time and space. Only posible if there 2 or more temporal windows.
#' @param linking_time Maximum number of timesteps of distance to link hotspots with the same temporal ID. Only necesary if add_temporal_id is TRUE.
#' @param linking_dist Linking distance used to link the clusters from the different time frames. Only necesary if add_temporal_id is TRUE.
#' @param get_timelife It specifies if the time periods and timelife of clusters are obtained. Only necesary if add_temporal_id is TRUE.
#'
#' @details The epifriends package uses the RANN package which can gives the exact nearest neighbours using the friends of friends algorithm. For more information on the RANN library please visit https://cran.r-project.org/web/packages/RANN/RANN.pdf
#'
#' @return  List with the next objects:
#' temporal_catalogues: list of catalogues
#'   List of EpiFRIenDs catalogues, where each element contains the catalogue in each
#'   time step.
#' mean_date: vector of chorn date times.
#'  List of dates corresponding to the median time in each time window
#'  Save dates in a temporal format.
#'
#' @export
#' 
#' @author Mikel Majewski Etxeberria based on earlier python code by Arnau Pujol.
#'
#' @examples
#' # Required packages
#' if(!require("RANN")) install.packages("RANN")
#' if(!require("chron")) install.packages("chron")
#' library("RANN")
#' library("chron")
#'
#' # Creation of x vector of longitude coordinates, y vector of latitude coordinates and finaly merge them on a position data frame.
#' x <- c(1,2,3,4,7.5,8,8.5,9,10,13,13.1,13.2,13.3,14,15,30)
#' y <- c(1,2,3,4,7.5,8,8.5,9,10,13,13.1,13.2,13.3,14,15,30)
#' pos <- data.frame(x,y)
#'
#' # Creation of test data frame with 0 for negative cases and 1 for positive clases for each position.
#' test <- data.frame(c(0,1,1,0,1,0,1,1,0,0,0,0,1,0,1,1))
#'
#' #Creation of chron dates vector in format 'y-m-d h:m:s'.
#' dates <- c("2020-11-03 05:33:07","2021-05-19 10:29:59","2021-02-09 14:53:20","2021-11-21 02:35:38","2020-11-19 05:57:24",
#' "2021-06-09 07:50:30","2021-09-18 05:53:26","2020-03-19 17:16:56","2021-06-08 12:40:46","2020-06-26 05:01:31",
#' "2020-10-15 04:40:27","2021-05-28 15:23:23","2020-05-01 02:56:54","2020-08-19 22:45:35","2021-10-23 18:56:35",
#' "2020-10-19 00:01:25")
#'
#' # Creation of temporal catalogue for this data.
#' tcat <- tcat <- temporal_catalogue(pos, test, dates ,link_d = 2, time_width = 305, time_steps = 305, linking_time = 3, linking_dist = 2)
#' 
temporal_catalogue <- function(positions, test_result, dates, link_d=NULL, min_neighbours = 2,
                               use_link_d = TRUE, link_neighbours = 10,time_width, 
                               min_date = NULL, max_date = NULL, time_steps = 1,
                               max_p = 1, min_pos = 2, min_total = 2, min_pr = 0,
                               add_temporal_id = TRUE, linking_time, linking_dist, get_timelife = TRUE){

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
  }else{
    #If min_date is passed as a character
    if(is.character(min_date)){
      #Convert string format to date.
      dtparts <- strsplit(min_date," ")
      min_date <- chron(dates=dtparts[[1]][1],times=dtparts[[1]][2],format=c('y-m-d','h:m:s'))
    }
  }
  if(is.null(max_date)){
    max_date <- max(dates)
  }else{
    #If max_date is passed as a character
    if(is.character(max_date)){
      #Convert string format to date.
      dtparts <- strsplit(max_date," ")
      max_date <- chron(dates=dtparts[[1]][1],times=dtparts[[1]][2],format=c('y-m-d','h:m:s'))
    }
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
    Newcatalogue <- catalogue(positions=selected_positions,test_result= selected_test_results,
                              link_d=link_d, min_neighbours = min_neighbours,
                              use_link_d=use_link_d, link_neighbours=link_neighbours,
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
  if(add_temporal_id == TRUE && length(temporal_catalogues)>1){
    temporal_catalogues <- add_temporal_id(temporal_catalogues,linking_time,linking_dist,get_timelife)
  }
  returns <- list(temporal_catalogues, mean_date)
  names(returns) <- c("temporal_catalogues", "mean_date")
  return(returns)
}

#' This method calculates the Euclidean distance between two positions.
#'
#' @param pos_a First position.
#' @param pos_b: Second position.
#'
#' @return Distance between positions.
#' 
#' @export
#' 
#' @author Eric Matamoros Morales based on earlier python code by Arnau Pujol.
#'
#' @examples
#' # Creation of x vector of longitude coordinates, y vector of latitude coordinates and finaly merge them on a position data frame.
#' pos_a <- c(1,2)
#' pos_b <- c(4,7)
#'
#' # Compute distances between two positions
#' x <- distance(pos_a, pos_b)
#' 
distance <- function(pos_a, pos_b){
  dist = sqrt(sum((pos_a - pos_b)**2))
  return(dist)
}

#' This method gives the unique values of a variable in a list of lists(temporal catalogues with temporal_id in our case).
#'
#' @param df_list List of lists(temporal catalogues with temporal_id in our case).
#' @param label Name of variable to select.
#'
#' @return List of unique values of the selected variable over all lists (temporal catalogues with temporal_id in our case) from the list.
#' 
#' @export
#' 
#' @author Eric Matamoros Morales based on earlier python code by Arnau Pujol.
#'
get_label_list <- function(df_list, label = "tempID"){
  for(i in 1:length(df_list)){
    mask = df_list[[i]][label][[1]][!is.na(df_list[[i]][label][[1]])]
    if(i == 1){
      label_list <- unique(mask)
    }else{
      label_list <- unique(c(label_list, unique(mask)))
    }
  }
  return(label_list)
}

#' This method obtains the first and last time frames for each temporal ID from a list of EpiFRIenDs catalogues and the corresponding timelife.
#'
#' @param catalogue_list List of lists (temporal catalogues with temporal_id in our case).
#'
#' @return List of hotspot catalogues with the added fields 'first_timestep', 'last_timestep' and 'lifetime'.
#' 
#' @export
#' 
#' @author Eric Matamoros Morales based on earlier python code by Arnau Pujol.
#'
get_lifetimes <- function(catalogue_list){
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

#' This method generates the temporal ID of EpiFRIenDs clusters by linking clusters from different time frames, assigning the same temporal ID to them when they are close enough in time and space.
#'
#' @param catalogue_list List of EpiFRIenDs catalogues, each element of the list corresponding to the catalogue of each timestep.
#' @param linking_time Maximum number of timesteps of distance to link hotspots with the same temporal ID.
#' @param linking_dist Linking distance used to link the clusters from the different time frames.
#' @param get_timelife It specifies if the time periods and timelife of clusters are obtained.
#' 
#' @return List of EpiFRIenDs catalogues with the added variable 'tempID' (and optionally the variables 'first_timestep', 'last_timestep' and 'lifetime').
#' 
#' @export
#' 
#' @author Eric Matamoros Morales based on earlier python code by Arnau Pujol.
#'
add_temporal_id <- function(
    catalogue_list, 
    linking_time, 
    linking_dist, 
    get_timelife = TRUE){
  
  #Case of empty catalogue list
  if(length(catalogue_list) == 0){
    return()
  }
  #setting empty values of temp_id
  for(t in 1:length(catalogue_list)){
    aux <- data.frame(matrix(NA,length(catalogue_list[[t]]$id)))
    colnames(aux) <- "tempID"
    catalogue_list[[t]] <- append(catalogue_list[[t]],aux)
    #catalogue_list[[t]]["tempID"] = vector(mode="list", length=length(catalogue_list[[t]]$id))
  }
  #Initialising tempID value to assign
  next_temp_id <- 0
  #Loop over all timesteps
  for(t in 1:(length(catalogue_list)-1)){
    #Loop over all timesteps within linking_time
    for (f in 1:length(catalogue_list[[t]]$id)){
      #Loop over all points of catalogue number 1
      for(t2 in (t + 1):min(t + linking_time, length(catalogue_list))){
        #Loop over all points of catalogue number 2
        for(f2 in 1:length(catalogue_list[[t2]]$id)){
          dist <- distance(catalogue_list[[t]]["mean_position_pos"][[1]][[f]], catalogue_list[[t2]]["mean_position_pos"][[1]][[f2]])
          if(dist <= linking_dist){
            temp_id1 <- catalogue_list[[t]]["tempID"][[1]][[f]]
            temp_id2 <- catalogue_list[[t2]]["tempID"][[1]][[f2]]
            #Assign tempIDs to linked clusters
            if((is.na(temp_id1)) && is.na((temp_id2))){
              catalogue_list[[t]]["tempID"][[1]][[f]]  <- next_temp_id
              catalogue_list[[t2]]["tempID"][[1]][[f2]] <- next_temp_id
              next_temp_id = next_temp_id + 1
            }else if(is.na(temp_id1)){
              catalogue_list[[t]]["tempID"][[1]][[f]]  <- temp_id2
            }else if(is.na(temp_id2)){
              catalogue_list[[t2]]["tempID"][[1]][[f2]] <- temp_id1
            }else if(temp_id1 != temp_id2){
              for(t3 in 1:length(catalogue_list)){
                catalogue_list[[t3]]$tempID <- replace(catalogue_list[[t3]]$tempID, catalogue_list[[t3]]$tempID == temp_id2, temp_id1)
              }
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
