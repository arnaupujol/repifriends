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
# Approach 2. Generate random distributions of local prevalence for each cluster
# and compute the mean minimum linking distance that matches the min_neighbors
##################################################################################

## Generate simulations and compute the mean of the minimum distances in each of them 

df_insights <- data.table('cluster_n' = c(), 'mean_init' = c(), 'mean_random' = c())

random_init <- 100
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

kmeans_filt[,.(sum(test)), .(clusters)]

