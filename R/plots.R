scatter_pval <- function(coordinates, id_data, positive, epi_catalogue, method = NULL){
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
  
  if(!is.null(method)){
    title <- paste0("P-value of hotspots - Method: ", method)
  }else{
    title <- "P-value of hotspots"
  }
  graph <- ggplot(coordinates,aes(x=x, y=y))+
    geom_point(color = "grey", size = 2.5)+
    geom_point(data = pos[id_data >0,], aes(colour = p_vals), size = 2.5)+
    scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07"))+
    ggtitle(title)
  return(graph)
}


scatter_pr <- function(coordinates, id_data, positive, epi_catalogue, method = NULL){
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
  pr <- c()
  for(i in id_data[id_data > 0]){
    pr <- append(pr, epi_catalogue$mean_pr[epi_catalogue$id == i])
  }
  # plot(coordinates$x, coordinates$y, pch = 19, col = "grey")
  # points(pos$x[id_data > 0], pos$y[id_data>0], pch = 19, col = rainbow(100)[factor(p_vals)])
  if(!is.null(method)){
    title <- paste0("P-value of hotspots - Method: ", method)
  }else{
    title <- "P-value of hotspots"
  }
  graph <- ggplot(coordinates,aes(x=x, y=y))+
    geom_point(color = "grey", size = 2.5)+
    geom_point(data = pos[id_data >0,], aes(colour = pr), size = 2.5)+
    scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07"))+
    ggtitle(title)
  return(graph)
}


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