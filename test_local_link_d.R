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
#######         AUTOMATIC LOCAL DETECTION OF LINK_D          #######
####################################################################

setwd("/Users/ericmatamoros/Desktop/TFM/repifriends/")
source("./utils.R")
min_neighbors <- 2 ## User specified minimum number of neighbors

data_v1 <- as.data.table(read.csv("./data/example_5.csv"))
data_v1[, id := 1:dim(data)[1]]

position_rand <- data_v1[,.(x,y)]
positive_rand = (data$test == 1)
test_rand <- data_v1[,5]

kmeans_prev <- compute_kmeans(
  clean_data(data_v1[,.(x, y, test)])[,.(x,y)], 
  clean_data(data_v1[,.(x, y, test)])$test)
prevalence <- kmeans_prev$prevalence

# Plot data
graph_base <- ggplot(data_v1, aes(x = x, y = y, color = as.factor(test))) +
  geom_point(size = 3) +
  labs(title = "Distribution of Positive and Negative Cases")

k_results <- ggplot(kmeans_prev, aes(x = x, y = y, color = as.factor(clusters))) + 
  geom_point(size = 2.5) + labs(title = "Detected Clusters using KMeans")

print(graph_base + k_results)

##################################################################################
# Approach 1. Determine percentiles based on distribution. For each positive item in the
# cluster evaluate the individual prevalence and finally make the mean of the individual one.
# The percentile selected will be the one whose mean individual prevalence is closer to the
# local prevalence of the cluster chosen. 
##################################################################################

cluster_n <- 2
min_distances <- get_min_distances(
  kmeans_prev[clusters == cluster_n, .(x,y)], 
  kmeans_prev[clusters == cluster_n]$test == 1,
  min_neighbors)
min_distances

quantiles <- seq(0, 1, by = 0.125)

quantiles_dist <- quantile(min_distances,probs = c(quantiles),na.rm = TRUE)
dist_thr <- 0.01

dt = data.table('quantiles' = quantiles, 'distances'= as.numeric(quantiles_dist))
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
    if ( abs(dt$distances[i] - dt$distances[i+1]) < dist_thr) {
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

prev_vect <- c()
counter <- 1
for(distance in dt$distances){
  indeces <- find_indeces(kmeans_prev[clusters == cluster_n, .(x,y)], distance, kmeans_prev[clusters == cluster_n, .(x,y)])
  
  kmeans_filter <- kmeans_prev[clusters == cluster_n,]
  kmeans_filter[, id := 1:dim(kmeans_filter)[1]]
  pos_indeces <- which(kmeans_prev[clusters == cluster_n,]$test == 1)
  
  prevalences <- sapply(pos_indeces, function(value, min_neighbors){
      summed <- sum(kmeans_filter[id %in% indeces[[value]],]$test) / nrow(kmeans_filter[id %in% indeces[[value]],])
    return(summed)
  })
  #if(nrow(kmeans_filter[id %in% indeces[[value]],]) == min_neighbors)
  prev_vect[counter] <- median(prevalences)
  counter <- counter +1 
}
dt$mean_prev <- prev_vect
dt$cluster_prev <- sum(kmeans_prev[clusters == cluster_n]$test) / nrow(kmeans_prev[clusters == cluster_n])

dt[, diff := abs(cluster_prev - mean_prev)] 

distance_selected <- dt[diff == min(diff)]$distances
print(paste0("Distance selected is: ", distance_selected))
print(past<e0("The distance belong to the percentile: ", dt[diff == min(diff)]$quantiles))
print(paste0("The prevalence of the cluster is: ", unique(dt$cluster_prev)))
print(paste0("The mean individual prevalence of the positive cases given the link_d provided is: ",dt[diff == min(diff)]$mean_prev))

# PROBLEM: LIKELY THAT PERCENTILE 1 IS ALWAYS SELECTED


##################################################################################
# Approach 2. Generate random distributions of local prevalence for each cluster
# and compute the mean minimum linking distance that matches the min_neighbors
##################################################################################

## Generate simulations and compute the mean of the minimum distances in each of them 

df_insights <- data.table('cluster_n' = c(), 'mean_init' = c(), 'mean_random' = c())

random_init <- 500
for(cluster_n in unique(kmeans_prev$clusters)){
  print(cluster_n)
  mean_min_vect <- c()
  kmeans_filt <-  kmeans_prev[clusters == cluster_n,]
  min_distances <- get_min_distances(
    kmeans_filt[, .(x,y)], 
    kmeans_filt$test == 1,
    min_neighbors)


  mean_vect <- c()
  for(i in 1:random_init){
    kmeans_sample <- data.table::copy(kmeans_filt)
    kmeans_sample[, test := stats::rbinom(nrow(kmeans_sample), 1, prev_cluster)]
    min_distances <- get_min_distances(
      kmeans_sample[, .(x,y)], 
      kmeans_sample$test == 1,
      min_neighbors)
    
    if(length(min_distances) == 0){
      mean_vect[i] <- 0
    }else{
      mean_vect[i] <- mean(min_distances)
    }
  }

  df_sub <- data.table(
    'cluster_n' = cluster_n, 
    'mean_init' = mean(min_distances), 
    'mean_random' = mean(mean_vect))
  df_insights <- rbind(df_insights, df_sub)
}
df_insights


