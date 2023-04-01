#' @export
scatter_pval <- function(coordinates, id_data, positive, epi_catalogue){
  # This method shows a scatter plot of the data distribution, showing in 
  # colour the positive cases that belong to foci (with the colour 
  # representing their p-value) and in grey the rest of the data. 
  # 
  # Parameters:
  # -----------
  # coordinates: data.frame
  #     data frame with the values of the coordinates
  # id_data: data.frame
  #     data frame with the cluster ID associated to each positive case, o for 
  #     no associated cluster
  # positive: data.frame
  #     Boolean vector that indicates if the case is infected. 
  # epi_catalogue: list
  #     List of the EpiFRIenDs catalogue
  # 
  # Returns:
  # --------
  # Scatter plot of the data distribution, showing in 
  # colour the positive cases that belong to foci (with the colour 
  # representing their p-value) and in grey the rest of the data.
  pos <- data.frame(coordinates[positive,],id_data)
  p_vals <- c()
  for(i in id_data[id_data > 0]){
    p_vals <- append(p_vals, epi_catalogue$p[epi_catalogue$id == i])
  }
  # plot(coordinates$x, coordinates$y, pch = 19, col = "grey")
  # points(pos$x[id_data > 0], pos$y[id_data>0], pch = 19, col = rainbow(100)[factor(p_vals)])
  graph <- ggplot(coordinates,aes(x=x, y=y))+
    geom_point(color = "grey", size = 2.5)+
    geom_point(data = pos[id_data >0,], aes(colour = p_vals), size = 2.5)+
    scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07"))+
    ggtitle("P-value of hotspots") + coord_equal()
  return(graph)
}

#' @export
size_histogram <- function(catalogue){
  # This method shows a histogram of the size (total number of cases) in foci 
  # for foci with p>0.05 and with p<0.05. 
  #   
  # Parameters:
  # -----------
  # catalogue: list
  #   List of the EpiFRIenDs catalogue
  #   
  # Returns:
  # --------
  # Histogram of number of foci per total number of cases with p<0.05 in red 
  # and p>0.05 in blue.
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

#' @export
hist_timelifes <- function(list_catalogues) {
  # This method plots a histogram of the lifetime of the detected clusters.
  #
  # Parameters:
  # -----------
  # list_catalogues:  list of data.frame
  #     List of EpiFRIenDs catalogues, each element of the list
  #     corresponding to the EpiFRIenDs catalogue of each timestep
  #
  # Returns:
  # --------
  # ggplot object with histogram
  
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

lifetime_timeline <- function(list_catalogues, mean_dates, time_steps) {
  # This method visualised the size evolution of the clusters in time
  # and colour-codes them with their total lifetime.
  #
  # Parameters:
  # -----------
  # list_catalogues:  list of data.frame
  #     List of EpiFRIenDs catalogues, each element of the list
  #     corresponding to the EpiFRIenDs catalogue of each timestep
  # mean_dates: list
  #     List of dates corresponding to the median time in each time window
  # time_steps: int
  #     Number of days that the time window is shifted in each time step
  #
  # Returns:
  # --------
  # ggplot object with the size evolution and lifetimes of clusters
  
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