library(data.table)
library(readr)
library(factoextra)
library(NbClust)
library(doParallel)
library(reshape2)
library(pdftools)
library(clValid)
library(epifriends)
library(ggplot2)
library(ggforce)
library(useful)
library(RANN)
library(patchwork)
library(knitr)
library(kableExtra)
library(grid)

####################################################################
####### AUTOMATIC DETECTION OF LINK_D BASED ON MIN NEIGHBORS #######
####################################################################

setwd("/Users/ericmatamoros/Desktop/TFM/repifriends/")
min_neighbors <- 2 ## User specified minimum number of neighbors

df <- as.data.table(read.csv("./data/mock_data_clustered.csv"))
df[, id := 1:dim(df)[1]]

link_d <- opt_link_d(df, min_neighbors, cluster_id=NULL, dist_prop = 0.1,
                       thr_impr = 0.05, diff_quantile = 0.05, quantile_est =  "fixed", 
                       keep_null_tests = FALSE, in_latlon = FALSE, to_epsg = NULL, 
                       verbose = TRUE)

position <- df[,.(x,y)]
positive = (df$test == 1)
test_result <- df[,.(test)]

# Plot data
graph_base <- ggplot(data, aes(x = x, y = y, color = as.factor(test))) +
  geom_point(size = 3) +
  labs(title = "Distribution of Positive and Negative Cases")
print(graph_base)

position = position_rand; positive = positive_rand; test = test_rand$test; min_neighbours = min_neighbours

# Get minimum distance & quantiles of it
min_distances <- get_min_distances(position[,.(x,y)], positive, min_neighbors)
hist(min_distances)

##################################################################################
### DETECTION OF TOP QUANTILES TO EVALUATE BASED ON DISTRIBUTION OF DISTANCES ####
##################################################################################
selection <- "fixed"

if(selection == 'random'){
  incr <- 0.05
  quantiles <- c(0, 0.25,0.5,0.75, 1)
  # Generate 2 random values between each of the standard quantiles
  quantiles <- sort(as.vector(sapply(1:(length(quantiles) - 1), function(x) 
    runif(2, min = quantiles[x], max = quantiles[x+1])
    )))
  print("Selection method random")
  print(quantiles)
}else if(selection == 'fixed'){
  quantiles <- seq(0, 1, by = 0.125)
  print("Selection method fixed")
  print(quantiles)
}else{
  stop("Any valid selected method for quantile estimation. Please check documentation")
}


min_distances_norm <- (min_distances - min(min_distances)) / (max(min_distances) - min(min_distances))
quantiles_dist <- quantile(min_distances_norm,probs = c(quantiles),na.rm = TRUE)
dist_thr <- 0.05 # think of proportions

dt = data.table('quantiles' = quantiles, 'distances'= as.numeric(quantiles_dist))
dt = dt[!(quantiles %in% c(0, 1))]
print(dt)

dt_prev <- dt[1:2,]
dt_post <- dt[1,]
counter <- 1
while( (dim(dt_prev)[1] != dim(dt_post)[1]) & (counter < 4)){
  # Loop over consecutive pairs of rows
  dt_prev <- data.table::copy(dt)
  for (i in 1:(nrow(dt)-1)) {
    
    if(i >= nrow(dt)){
      break
    }
    # Check if distances are lower than threshold (dist_thr)
    if (abs(dt$distances[i] - dt$distances[i+1])  < dist_thr) {
      # Compute mean of distances and update first row
      dt[i, distances := mean(c(dt$distances[i], dt$distances[i+1]))]
      dt[i, quantiles := mean(c(dt$quantiles[i], dt$quantiles[i+1]))]
      # Delete second row
      dt <- dt[-(i+1), ]
      # Update row index
      i <- i - 1
    }
  }
  dt_post <- data.table::copy(dt)
  counter <- counter +1 
}
print(dt)


# Determine the level of significance of p-values based on the distribution of distances
opt_metric <- 0
opt_quantile <- c()
thr_impr <- 0.05 # Difference in 
diff_quantile <- 0.05 # Difference in quantiles to search for new locals
quantiles <- dt$quantiles

counter <- 1

traces <- list()
while( (thr_impr >= 0.05) & (counter < 4)){
  
  quantiles_dist <- quantile(min_distances,unique(quantiles),na.rm = TRUE)
  
  cl <- makeCluster(detectCores() - 2)
  registerDoParallel(cl)
  
  ## Compute scoring metric using parallel or for loop depending on size of the evaluation inputs
  if(length(quantiles_dist) > 2){
    
    # Use foreach to parallelize the for loop
    metric <- foreach(quant = 1:length(quantiles_dist), .combine = "c") %dopar% {
      
      library(epifriends)
      library(RANN)
      # Code to execute in parallel
      categories <- catalogue(
        data.table::copy(position), test, quantiles_dist[quant], cluster_id = NULL, min_neighbours = min_neighbors,
        method = 'base')
      
      kpi <- 0
      for(j in 1:length(categories$cluster_id)){
        if(categories$cluster_id[j]!=0){
          kpi <- kpi + (1-categories$pval_cluster[j])
        }
      }
      
      return(kpi)
    }
    # Stop the parallel backend
    stopCluster(cl)
  }else{
    metric <- c()
    for(quant in 1:length(quantiles_dist)){
      categories <- catalogue(
        data.table::copy(position), test_result, quantiles_dist[quant], cluster_id = NULL, min_neighbours = min_neighbors,
        method = 'base')
      
      kpi <- 0
      for(j in 1:length(categories$cluster_id)){
        if(categories$cluster_id[j]!=0){
          kpi <- kpi + (1-categories$pval_cluster[j])
        }
      }
      
      metric[quant] <- kpi
    }
    
  }
  
  # Obtain best quantile and best metric
  best_quantile <- quantiles[which(metric == max(metric))]
  best_metric <- max(metric)
  
  # Store for tracing and plots
  traces[[counter]] <- data.table('quantiles' = unique(quantiles), 'metrics' = metric)
  
  # Check if metric has been improved by more than a 5%
  perc_impr <- (best_metric - (opt_metric + 0.00000001)) / (opt_metric + 0.00000001)
  if(perc_impr > thr_impr){
    opt_metric <- best_metric
    opt_quantile <- best_quantile
    counter <- counter + 1
  }else{
    break
  }
  
  # Create list of vectors and cap
  quantiles <- c(opt_quantile - diff_quantile, opt_quantile, opt_quantile + diff_quantile)
  quantiles <- pmin(quantiles, 1)
  quantiles <- pmax(quantiles, 0)
}

opt_link_d <- as.numeric(quantile(min_distances,opt_quantile,na.rm = TRUE))

print(paste0("Optimal quantile: ", opt_quantile))
print(paste0("Best maximization metric obtained: ", best_metric))
print(paste0("Linking Distance: ", opt_link_d))

traces <- unique(rbindlist(traces))
setorderv(traces, 'quantiles')

traces[, link_d := quantile(min_distances,quantiles,na.rm = TRUE)]
plot1 <- ggplot(traces, aes(x = link_d, y = metrics)) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  labs(title = "Best metric based on link_d")
print(plot1)
      