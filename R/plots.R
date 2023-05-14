# This source code is licensed under the GNU GENERAL PUBLIC LICENSE license found in the
# LICENSE file in the root directory of this source tree.


#' This method shows a scatter plot of the data distribution, showing in colour the positive cases that belong to foci (with the colour 
#' representing their p-value) and in grey the rest of the data. 
#'
#' @param coordinates data frame with the values of the coordinates
#' @param id_data data frame with the cluster ID associated to each positive case, o for no associated cluster
#' @param positive Boolean vector that indicates if the case is infected. 
#' @param epi_catalogue List of the EpiFRIenDs catalogue
#' @param xlims: Vector with the limits of the X axis.
#' @param ylims: Vector with the limits of the Y axis
#' @param title: Title of the plot. If NULL, it defaults to 'P-value of hotspots'
#' 
#' @return  Scatter plot of the data distribution, showing in colour the positive cases that belong to foci (with the colour 
#' representing their p-value) and in grey the rest of the data.
#'
#' @export
#' 
#' @author Eric Matamoros Morales based on earlier python code by Arnau Pujol.
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
#' test <- c(0,1,1,0,1,0,1,1,0,0,0,0,1,0,1,1)
#'
#' # Creation of catalogue for this positions, linking distance 2 and default values.
#' cat <- catalogue(pos, test, 2)
#' 
#' scatter_pval(pos, cat$cluster_id, (test == 1), cat$epifriends_catalogue)
#' 
scatter_pval <- function(
    coordinates, 
    id_data, 
    positive, 
    prevalence,
    epi_catalogue,
    xlims = NULL,
    ylims = NULL,
    title = NULL){
  pos <- data.frame(coordinates[positive,],id_data)
  p_vals <- c()
  for(i in id_data[id_data > 0]){
    p_vals <- append(p_vals, epi_catalogue$p[epi_catalogue$id == i])
  }
 
  title <- "P-value on positive cases of hotspots"
  
  if(!is.null(prevalence)){
    
    pos <- as.data.table(pos)
    coordinates[, prevalence := prevalence]
    coordinates[, id := paste0(x, "_", y)]
    
    pos[, id := paste0(x, "_", y)]
    pos_filt <- pos[id_data >0,]
    pos_filt$p_vals <- p_vals
    pos_filt <- merge(pos_filt, coordinates[,.(id, prevalence)], by = "id", how = "left")
    graph <- ggplot(coordinates,aes(x=x, y=y, size = prevalence))+
      geom_point(color = "#F38B8B", shape=21, stroke = 1) +
      geom_point(
        data =pos_filt, 
        aes(color = "#F38B8B", fill = p_vals, size = prevalence), shape=21, stroke = 1)+
      scale_fill_continuous(limits = c(0, 0.2), low = "grey", high = "blue") +
      scale_size_continuous(range = c(1, 4), limits = c(0,1)) +
      ggtitle(title) + coord_equal()
    
  }else{
    
    pos_filt <- pos[id_data >0,]
    pos_filt$p_vals <- p_vals
    graph <- ggplot(coordinates,aes(x=x, y=y))+
      geom_point(color = "grey", size = 2.5)+
      geom_point(data = pos_filt, aes(colour = p_vals), size = 2.5)+
      scale_color_continuous(limits = c(0, 0.2), low = "grey", high = "red") +
      #scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07"))+
      ggtitle(title) + coord_equal()
    
  }

  return(graph)
}

#' This method shows a histogram of the size (total number of cases) in foci for foci with p>0.05 and with p<0.05. 
#'
#' @param catalogue Catalogue of the epifriends clusters and their main characteristics as output of the catalogue() function.
#' 
#' @return  Histogram of number of foci per total number of cases with p<0.05 in red and p>0.05 in blue.
#'
#' @export
#' 
#' @author Eric Matamoros Morales based on earlier python code by Arnau Pujol.
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
#' test <- c(0,1,1,0,1,0,1,1,0,0,0,0,1,0,1,1)
#'
#' # Creation of catalogue for this positions, linking distance 2 and default values.
#' cat <- catalogue(pos, test, 2)
#' 
#' size_histogram(cat)
#' 
size_histogram <- function(catalogue){
  cases <- data.frame(catalogue$epifriends_catalogue$total, catalogue$epifriends_catalogue$p)
  colnames(cases) <- c("total","p")
  #auxiliar vector to avoid converting doubles to strings
  aux <- cases$p
  for(i in 1:length(cases$p)){
    if(aux[i] > 0.05){
      cases$p[i] <- "p > 0.05"
    }else{
      cases$p[i] <- "p < 0.05"
    }
  }
  histo <- ggplot(cases,aes(x=total, fill=p))+
    geom_histogram()+
    scale_fill_manual(values = c("red", "blue"))+
    xlab("Number of cases in focci")+
    ylab("Number of focci")
  return(histo)
}

#' This method plots a histogram of the lifetime of the detected clusters. 
#'
#' @param list_catalogues List of EpiFRIenDs catalogues, where each element contains the catalogue in each time step.
#' 
#' @return  Histogram of the lifetime of the detected clusters
#'
#' @export
#' 
#' @author Eric Matamoros Morales based on earlier python code by Arnau Pujol.
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
#' test <- c(0,1,1,0,1,0,1,1,0,0,0,0,1,0,1,1)
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
#' hist_timelifes(list_catalogues$temporal_catalogues)
#' 
hist_timelifes <- function(list_catalogues) {
  all_tempids <- epifriends::get_label_list(list_catalogues, label = 'tempID')
  if (length(all_tempids) > 0) {
    lifetimes <- c()
    for (tempid in all_tempids) {
      for (t in list_catalogues) {
        if (tempid %in% unique(t[['tempID']])) {
          mask <- t[['tempID']] == tempid
          lifetimes <- c(lifetimes, mean(t[['lifetime']][mask]))
          break
        }
      }
    }
    graph <- ggplot(data = data.frame(lifetimes), aes(x = lifetimes)) +
      geom_histogram(bins = 50) +
      labs(x = "Time duration of clusters (number of time steps)",
           y = "Number of clusters",
           title = "Histogram of the lifetime of the detected clusters")
    
    return(graph)
  }
}

# This method visualises the size evolution of the clusters in time and colour-codes them with their total lifetime.
#'
#' @param list_catalogues List of EpiFRIenDs catalogues, where each element contains the catalogue in each time step.
#' @param mean_dates List of dates corresponding to the mean time in each time window
#' @param time_steps Number of days that the time window is shifted in each time step
#' 
#' @return  ggplot object with the size evolution and lifetimes of clusters
#'
#' @export
#' 
#' @author Eric Matamoros Morales based on earlier python code by Arnau Pujol.
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
#' time_steps = 305
#'
#' # Creation of test data frame with 0 for negative cases and 1 for positive clases for each position.
#' test <- c(0,1,1,0,1,0,1,1,0,0,0,0,1,0,1,1)
#'
#' #Creation of chron dates vector in format 'y-m-d h:m:s'.
#' dates <- c("2020-11-03 05:33:07","2021-05-19 10:29:59","2021-02-09 14:53:20","2021-11-21 02:35:38","2020-11-19 05:57:24",
#' "2021-06-09 07:50:30","2021-09-18 05:53:26","2020-03-19 17:16:56","2021-06-08 12:40:46","2020-06-26 05:01:31",
#' "2020-10-15 04:40:27","2021-05-28 15:23:23","2020-05-01 02:56:54","2020-08-19 22:45:35","2021-10-23 18:56:35",
#' "2020-10-19 00:01:25")
#'
#' # Creation of temporal catalogue for this data.
#' tcat <- tcat <- temporal_catalogue(pos, test, dates ,link_d = 2, time_width = 305, time_steps = time_steps, linking_time = 3, linking_dist = 2)
#' 
#' lifetime_timeline(list_catalogues$temporal_catalogues, list_catalogues$mean_dates, time_steps)
#' 
lifetime_timeline <- function(list_catalogues, mean_dates, time_steps) {
  tempids <- epifriends::get_label_list(list_catalogues, label = 'tempID')
  lifetimes <- epifriends::get_label_list(list_catalogues, label = 'lifetime')
  num_cases_df <- data.frame()
  
  num_cases_df <- data.table("mean_dates" = NULL, "tempID" = NULL, "num_cases" = NULL)
  for (i in seq_along(tempids)) {
    num_cases <- numeric(length(mean_dates))
    dates <- c()
    min_date <- list_catalogues[[1]]$Date
    tempid = tempids[i]
    
    for (f in seq_along(list_catalogues)) {
      dates[f] <- min_date + f*time_steps
      if (tempid %in% unique(list_catalogues[[f]]$tempID)) {
        mask <- list_catalogues[[f]]$tempID == tempid
        mask[is.na(mask)] <- FALSE
        num_cases[f] <- sum(list_catalogues[[f]]$total[mask])
        lifetime <- list_catalogues[[f]]$lifetime[mask][1]
      } else {
        num_cases[f] <- 0
      }
    }
    dates <- as.chron(dates)
    
    num_cases_df <- rbind(
      num_cases_df, 
      data.table("mean_dates" = as.Date(mean_dates), 
                 "tempID" = tempid, "num_cases" = num_cases, "lifetime" = lifetime)
    )
  }
  
  graph <- ggplot(num_cases_df, aes(x = mean_dates, y = num_cases, color = as.factor(lifetime))) +
    geom_line(size = 2) +
    labs(x = "Date", y = "Number of cases in cluster",
         title = "Size evolution of the clusters colored by their total lifetime")
  
  return(graph)
}


# Plot the coordinates with different colors for each of the clusters identified by the EpiFRIenDs.
#'
#' @param coordinates data frame with the values of the coordinates. 
#' @param catalogue List of EpiFRIenDs catalogues.
#' @param only_significant If True, color only the clusters identified as significant by the algorithm.
#' 
#' @return  ggplot object with the clusters identified by the EpiFRIenDs algorithm.
#'
#' @export
#' 
#' @author Eric Matamoros Morales based on earlier python code by Arnau Pujol.
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
#' test <- c(0,1,1,0,1,0,1,1,0,0,0,0,1,0,1,1)
#'
#' # Creation of catalogue for this positions, linking distance 2 and default values.
#' cat <- catalogue(pos, test, 2)
#' 
#' plot_clusters(pos, cat)
plot_clusters <- function(coordinates, catalogue, only_significant = FALSE){
  coords <- data.table::copy(coordinates)
  coords[, cluster := 0]
  coords[, index := 1:nrow(coords)]
  
  if(only_significant){
    indexes <- which(categories$epifriends_catalogue$p <= 0.05)
  }else{
    indexes <- 1:length(catalogue$epifriends_catalogue$p)
  }
  for(clusters in indexes){
    coords[index %in% catalogue$epifriends_catalogue$indeces[[clusters]], cluster := clusters]
  }
  
  graph_clusters <- ggplot(coords[cluster != 0], aes(x = x, y = y, color = as.factor(cluster))) +
    geom_point(size = 2.5)  + geom_point(data = coords[cluster == 0], aes(x=x, y=y, color = "#FFFFFF"), shape=21, stroke = 1) +
    labs(title = "Distribution of Clusters")
  
  return(graph_clusters)
}
