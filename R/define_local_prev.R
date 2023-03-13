define_loc_prev <- function(positions, test, method = 'kmeans', seed = 123, thr_radial = 0.2, 
                            link_d = NULL,coord_cols = c("x", "y"), plot_kmeans = FALSE){
  
  df = data.table(positions)
  df[, test := test]
  df[, id := 1:nrow(df)]

  if(method == 'kmeans'){
    ########################################
    ###     APPROACH 2: KMEANS METHOD    ###
    ########################################
    
    # Compute k-means with k = 4
    set.seed(seed)
    num_clusters <- ceiling(max(df$x) / link_d)
    km <- kmeans(df[,..coord_cols], centers=num_clusters, nstart=100)
    
    clusters <- dbscan::dbscan(df[,..coord_cols], eps = link_d, minPts = 2)
    
    clusters <- dbscan(df[,..coord_cols], link_d = link_d, min_neighbours = 2)
    
    #km_detect <- NbClust(data = df[,..coord_cols], distance = "euclidean", 
    #                     method = "kmeans", min.nc = 2, max.nc = 8)
    
    #print(paste0("Best split of clusters is: ", max(km_detect$Best.partition)))
    #df[, clusters := km_detect$Best.partition]
    #df[, clusters := km$cluster]
    df[, clusters := clusters]
    df[, prevalence := sum(test) / .N, by = "clusters"]
    
    if(plot_kmeans == TRUE){
      graph <- ggplot(df,aes(x=x, y=y, group = test))+
        geom_point(aes(col = clusters, size = test, shape = as.factor(test))) + 
        scale_size_continuous(range = c(1.5,2))
      graph = graph+scale_size(guide=FALSE)
      graph
    }
    
  }else if(method == 'radial'){
    
    ########################################
    ###     APPROACH 2: RADIAL METHOD    ###
    ########################################

    combs <- dist(df[,.(x,y)])
    combs <- as.data.table(melt(as.matrix(combs), 
                                varnames = c("element_1", "element_2")))
    combs[, element_2 := as.numeric(element_2)]
    setnames(combs, "value", "distance")
    
    tests <- data.table(
      "element_2" = as.numeric(rownames(df_radial)),
      "test" = df_radial$test
    )
    combs_dist <- merge(combs, tests, by = 'element_2', how = "left")
    combs_dist <- combs_dist[distance <= thr, .(prevalence = sum(test) / .N), by = "element_1"]
    combs_dist[, element_1 := as.numeric(element_1)]
    setnames(combs_dist, "element_1", "id")
    
    df = merge(df, combs_dist, by = "id", how = "left")
    
  }else{
    stop("Method specified does not exist, please check documentation")
  }
  
  return(df)
}
