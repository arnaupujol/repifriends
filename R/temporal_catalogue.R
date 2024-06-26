# This source code is licensed under the GNU GENERAL PUBLIC LICENSE license 

#' This method generates a list of EpiFRIenDs catalogues representing different 
#' time frames by including only cases within a time window that moves within 
#' each time step.
#'
#' @param x Vector of x positions.
#' @param x Vector of y positions.
#' @param test_result vector of test results (0 or 1).
#' @param dates Date times of the corresponding data in format 'y-m-d h:m:s'.
#' @param link_d the linking distance of the DBSCAN algorithm. Should be in the
#'same scale as the positions.
#' @param use_link_d: If True, use linking distance to determine the closest 
#' neighbors. If False, use default linking neighbours based on proximity.
#' @param prevalence Probability of having an infected case for each individual.
#' @param link_neighbours: Number of surrounding neighbors to link. 
#' @param min_neighbours Minium number of neighbours in the radius < link_d 
#' needed to link cases as friends.
#' @param time_width Number of days of the time window used to select cases in 
#' each time step.
#' @param min_date Initial date used in the first time step and time window 
#' selection. In format 'y-m-d h:m:s'.
#' @param max_date Final date to analyse, defining the last time window as the 
#' one fully overlapping the data. In format 'y-m-d h:m:s'.
#' @param time_steps Number of days that the time window shifts in each time step.
#' @param max_p Maximum value of the p-value to consider the cluster detection.
#' @param min_pos Threshold of minimum number of positive cases in clusters applied.
#' @param min_total Threshold of minimum number of cases in clusters applied.
#' @param min_pr Threshold of minimum positivity rate in clusters applied.
#' @param add_temporal_id Boolean that indicates if we want to add a temporal id
#' to the hotspots indicating which hotspots are near in time and space. 
#' Only posible if there 2 or more temporal windows.
#' @param linking_time Maximum number of timesteps of distance to link hotspots 
#' with the same temporal ID. Only necesary if add_temporal_id is TRUE.
#' @param linking_dist Linking distance used to link the clusters from the different 
#' time frames. Only necesary if add_temporal_id is TRUE.
#' @param get_timelife It specifies if the time periods and timelife of clusters 
#' are obtained. Only necesary if add_temporal_id is TRUE.
#' @param optimize_link_d: If True, optimize the linking distance based on minimum 
#' distribution of distances between neighbors.
#' @param method Method that wants to be used to compute the local prevalence - 
#' either 'kmeans', 'centroid', or 'base'. Defaults to 'base'.
#' @param keep_null_tests Whether to remove or not missings. If numeric, 
#' provide value to impute.
#' @param in_latlon:  If True, x and y coordinates are treated as longitude and 
#' latitude respectively, otherwise they are treated as cartesian coordinates.
#' @param to_epsg: If in_latlon is True, x and y are reprojected to this EPSG.
#' @param consider_fd: If True, consider false detections and adjust p-value of that.
#' @param n_simulations: Numeric value with the number of desired iterations to 
#' compute the false-detected clusters.
#' @param verbose: If TRUE, print information of the process; else, do not print.
#' @param store_gif: If TRUE, store the different time-frame images in an animated 
#' GIF format.
#' @param use_geom_map: If TRUE, the generated plot will have a geo-map.
#' @param out_gif_path: Output directory of the GIF animated video. Only useful 
#' if store_gif parameter is set to TRUE. By default a new folder called /gif will 
#' be created in the working directory.
#'
#' @return  List with the next objects:
#' temporal_catalogues: list of catalogues
#'   List of EpiFRIenDs catalogues, where each element contains the catalogue in each
#'   time step.
#' mean_date: vector of chorn date times.
#'  List of dates corresponding to the median time in each time window
#'  Save dates in a temporal format.
#'  
#' @importFrom magick image_read image_join image_animate image_write
#' @importFrom chron chron
#'
#' @export
#' 
#' @author Mikel Majewski Etxeberria based on earlier python code by Arnau Pujol.
#'
#' @examples
#' # Creation of x vector of longitude coordinates, y vector of latitude 
#' # coordinates and finaly merge them on a position data frame.
#' x <- c(1,2,3,4,7.5,8,8.5,9,10,13,13.1,13.2,13.3,14,15,30)
#' y <- c(1,2,3,4,7.5,8,8.5,9,10,13,13.1,13.2,13.3,14,15,30)
#' pos <- data.frame(x,y)
#'
#' # Creation of test data frame with 0 for negative cases and 1 for positive 
#' # clases for each position.
#' test <-c(0,1,1,0,1,0,1,1,0,0,0,0,1,0,1,1)
#'
#' #Creation of chron dates vector in format 'y-m-d h:m:s'.
#' dates <- c("2020-11-03 05:33:07","2021-05-19 10:29:59","2021-02-09 14:53:20",
#' "2021-11-21 02:35:38","2020-11-19 05:57:24", "2021-06-09 07:50:30","2021-09-18 05:53:26",
#' "2020-03-19 17:16:56","2021-06-08 12:40:46","2020-06-26 05:01:31", "2020-10-15 04:40:27",
#' "2021-05-28 15:23:23","2020-05-01 02:56:54","2020-08-19 22:45:35","2021-10-23 18:56:35",
#' "2020-10-19 00:01:25")
#'
#' # Creation of temporal catalogue for this data.
#' tcat <- temporal_catalogue(x = x, y = y, test = test, dates = dates ,link_d = 2, 
#' time_width = 305, time_steps = 305, linking_time = 3, linking_dist = 2, get_timelife=FALSE)
#' 
temporal_catalogue <- function(x, y, test_result, dates, link_d, prevalence = NULL, 
                               use_link_d = TRUE,  link_neighbours = NULL,min_neighbours = 2, 
                               time_width, min_date = NULL, max_date = NULL, time_steps = 1,
                               max_p = 1, min_pos = 2, min_total = 2, min_pr = 0,
                               add_temporal_id = TRUE, linking_time, linking_dist, 
                               get_timelife = TRUE, optimize_link_d = FALSE, method = "base",
                               keep_null_tests = FALSE, in_latlon = FALSE, to_epsg = NULL,
                               consider_fd = FALSE, n_simulations= 500, verbose = FALSE, 
                               store_gif = FALSE, use_geom_map = FALSE, 
                               out_gif_path = paste0(getwd(),"/www/")
                               ){
  
  #Create data.table with all coordinates & test
  positions = data.table("x" = x, "y" = y, "test" = test_result, "prevalence" = prevalence)
  
  #Case of empty list of dates
  if(length(dates) == 0){
    print("There are no dates")
    return()
  }
  dtparts = t(as.data.frame(strsplit(dates," ")))
  row.names(dtparts) = NULL
  date_array <- chron(dates=dtparts[,1],times=dtparts[,2],format=c('y-m-d','h:m:s'))
  #Defining temporal range
  if(is.null(min_date)){
    min_date <- min(date_array)
  }else{
    #If min_date is passed as a character
    if(is.character(min_date)){
      #Convert string format to date.
      dtparts <- strsplit(min_date," ")
      min_date <- chron(dates=dtparts[[1]][1],times=dtparts[[1]][2],format=c('y-m-d','h:m:s'))
    }
  }
  if(is.null(max_date)){
    max_date <- max(date_array)
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
    selected_data <- (date_array >= min_date + time_steps*step_num)&
      (date_array <= min_date + time_steps*step_num + time_width)
    selected_positions <- positions[selected_data,]
    selected_test_results <- as.data.frame(test_result[selected_data])
    
    # Remove or impute missings
    to_impute <- colnames(selected_positions)[!(colnames(selected_positions) %in% c("x", "y"))]
    selected_positions = clean_unknown_data(selected_positions,to_impute,keep_null_tests,verbose)
    selected_test_results = data.frame("test" = selected_positions$test)

    #get catalogue
    Newcatalogue <- catalogue(x = selected_positions$x, y = selected_positions$y, 
                              test_result = selected_test_results$test,
                              link_d = link_d, use_link_d = use_link_d,
                              prevalence=selected_positions$prevalence, cluster_id = NULL,
                              method=method, min_neighbours = min_neighbours,
                              max_p = max_p, min_pos = min_pos, min_total = min_total, 
                              min_pr = min_pr,  optimize_link_d=optimize_link_d, 
                              keep_null_tests = keep_null_tests, in_latlon = in_latlon, 
                              to_epsg = to_epsg,verbose = verbose)

    #we only add the mean date if the catalogue is not empty
    if(FALSE %in% lapply(Newcatalogue$epifriends_catalogue, is.null)){
      mean_date <- append(mean_date, min_date + time_steps*step_num + 0.5*time_width)
      
      # Check if GIF wants to be created
      if(store_gif){
        minimum_date <- (min_date + time_steps*step_num + 0.5*time_width)
        plots <- scatter_pval(
          coordinates = selected_positions[,.(x,y)], 
          id_data = Newcatalogue$cluster_id, 
          positive = (selected_test_results$test == 1), 
          prevalence = prevalence,
          epi_catalogue = Newcatalogue$epifriends_catalogue,
          use_geom_map = use_geom_map,
          c(min(na.omit(positions)$x), max(na.omit(positions)$x)),
          c(min(na.omit(positions)$y), max(na.omit(positions)$y)),
          paste0("Date: ", as.character(minimum_date))
        )
        
        out_temp_path = paste0(getwd(),"/tmp/")
        if(!dir.exists(out_temp_path)){ # create temp file for storing png if doesn't exist
          dir.create(out_temp_path)
        }
        ggsave(paste0(out_temp_path, "/",minimum_date, ".png"), plot = plots, dpi = 300)
      }
      Newcatalogue$epifriends_catalogue$Date <- (min_date + time_steps*step_num + 0.5*time_width)
      temporal_catalogues[[step_num+1]] <- Newcatalogue$epifriends_catalogue
    }
    step_num = step_num + 1
  }
  if(add_temporal_id == TRUE && length(temporal_catalogues)>1){
    temporal_catalogues <- add_temporal_id(
      temporal_catalogues,linking_time,linking_dist,get_timelife
    )
  }
  returns <- list(temporal_catalogues, mean_date)
  names(returns) <- c("temporal_catalogues", "mean_date")
  
  # Convert PNG to GIF output
  if(store_gif){
    if(!dir.exists(out_gif_path)){ # create temp file for storing png if doesn't exist
      dir.create(out_gif_path)
    }
    gif_output_name <- "animated_gif.gif"
    gif_output_name <- paste0(out_gif_path, "/", gif_output_name)
    list.files(path=out_temp_path, pattern = '*.png', full.names = TRUE) %>% 
      image_read() %>% # reads each path file
      image_join() %>% # joins image
      image_animate(fps=50) %>% # animates, can opt for number of loops
      image_write(gif_output_name) # write to current dir
  }
  
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
#' # Creation of x vector of longitude coordinates, y vector of latitude coordinates and 
#' # finally merge them on a position data frame.
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

#' This method gives the unique values of a variable in a list of 
#' lists(temporal catalogues with temporal_id in our case).
#'
#' @param df_list List of lists(temporal catalogues with temporal_id in our case).
#' @param label Name of variable to select.
#'
#' @return List of unique values of the selected variable over all lists (temporal 
#' catalogues with temporal_id in our case) from the list.
#' 
#' @export
#' 
#' @author Eric Matamoros Morales based on earlier python code by Arnau Pujol.
#' 
#' 
#' # Creation of x vector of longitude coordinates, y vector of latitude 
#' # coordinates and finaly merge them on a position data frame.
#' x <- c(1,2,3,4,7.5,8,8.5,9,10,13,13.1,13.2,13.3,14,15,30)
#' y <- c(1,2,3,4,7.5,8,8.5,9,10,13,13.1,13.2,13.3,14,15,30)
#' pos <- data.frame(x,y)
#'
#' # Creation of test data frame with 0 for negative cases and 1 for positive 
#' # clases for each position.
#' test <-c(0,1,1,0,1,0,1,1,0,0,0,0,1,0,1,1)
#'
#' #Creation of chron dates vector in format 'y-m-d h:m:s'.
#' dates <- c("2020-11-03 05:33:07","2021-05-19 10:29:59","2021-02-09 14:53:20",
#' "2021-11-21 02:35:38","2020-11-19 05:57:24", "2021-06-09 07:50:30","2021-09-18 05:53:26",
#' "2020-03-19 17:16:56","2021-06-08 12:40:46","2020-06-26 05:01:31", "2020-10-15 04:40:27",
#' "2021-05-28 15:23:23","2020-05-01 02:56:54","2020-08-19 22:45:35","2021-10-23 18:56:35",
#' "2020-10-19 00:01:25")
#'
#' # Creation of temporal catalogue for this data.
#' tcat <- temporal_catalogue(x = x, y = y, test = test, dates = dates ,link_d = 2, 
#' time_width = 305, time_steps = 305, linking_time = 3, linking_dist = 2, get_timelife=FALSE)
#'
#' tlist <- get_label_list(tcat, label = "tempID")
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

#' This method obtains the first and last time frames for each temporal ID from 
#' a list of EpiFRIenDs catalogues and the corresponding timelife.
#'
#' @param cat_list List of lists (temporal catalogues with temporal_id in our case).
#'
#' @return List of hotspot catalogues with the added fields 'first_timestep', 
#' 'last_timestep' and 'lifetime'.
#' 
#' @export
#' 
#' @author Eric Matamoros Morales based on earlier python code by Arnau Pujol.
#'
get_lifetimes <- function(cat_list){
  #getting list of temporal IDs appearing in cat_list
  tempid_list <- get_label_list(cat_list, label = "tempID")
  #Creating empty columns for first timestep, last timestep and lifteime
  for(t in 1:length(cat_list)){
    cat_list[[t]]["first_timestep"] <- c()
    cat_list[[t]]["last_timestep"] <- c()
    cat_list[[t]]["lifetime"] <- c()
  }
  for(tempid_num in tempid_list){
    appearances = c()
    for(i in 1:length(cat_list)){
      if(tempid_num %in% unique(cat_list[[i]]["tempID"][[1]])){
        appearances <- append(appearances,i)
      }
    }
    min_appearance = min(appearances)
    max_appearance = max(appearances)
    lifetime = max_appearance - min_appearance
    for(i in min_appearance:max_appearance){
      cat_list[[i]]["first_timestep"][[1]][cat_list[[i]]["tempID"][[1]] == tempid_num] <- 
        min_appearance
      cat_list[[i]]["last_timestep"][[1]][cat_list[[i]]["tempID"][[1]] == tempid_num] <-
        max_appearance
      cat_list[[i]]["lifetime"][[1]][cat_list[[i]]["tempID"][[1]] == tempid_num] <-
        lifetime
    }
  }
  return(cat_list)
}

#' This method generates the temporal ID of EpiFRIenDs clusters by linking clusters
#' from different time frames, assigning the same temporal ID to them when they 
#' are close enough in time and space.
#'
#' @param catalogue_list List of EpiFRIenDs catalogues, each element of the list 
#' corresponding to the catalogue of each timestep.
#' @param linking_time Maximum number of timesteps of distance to link hotspots 
#' with the same temporal ID.
#' @param linking_dist Linking distance used to link the clusters from the different 
#' time frames.
#' @param get_timelife It specifies if the time periods and timelife of clusters are 
#' obtained.
#' 
#' @return List of EpiFRIenDs catalogues with the added variable 'tempID' 
#' (and optionally the variables 'first_timestep', 'last_timestep' and 'lifetime').
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
    min_length <- max(1, length(catalogue_list[[t]]$id))
    aux <- data.frame(matrix(NA,min_length))
    colnames(aux) <- "tempID"
    catalogue_list[[t]] <- append(catalogue_list[[t]],aux)
  }
  
  #Initialising tempID value to assign
  next_temp_id <- 0
  #Loop over all timesteps
  for(t in 1:(length(catalogue_list)-1)){
    #Loop over all timesteps within linking_time
    for (f in 1:length(catalogue_list[[t]]$id)){
      if(is.null(catalogue_list[[t]]$id)){
        next
      }
      #Loop over all points of catalogue number 1
      for(t2 in (t + 1):min(t + linking_time, length(catalogue_list))){
        if(is.null(catalogue_list[[t2]]$id)){
          next
        }
        #Loop over all points of catalogue number 2
        for(f2 in 1:length(catalogue_list[[t2]]$id)){
          dist <- distance(catalogue_list[[t]]["mean_position_pos"][[1]][[f]], 
                           catalogue_list[[t2]]["mean_position_pos"][[1]][[f2]])
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
                catalogue_list[[t3]]$tempID <- replace(
                  catalogue_list[[t3]]$tempID, 
                  catalogue_list[[t3]]$tempID == temp_id2, 
                  temp_id1)
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
