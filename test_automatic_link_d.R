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
source("./utils.R")
min_neighbors <- 2

data <- as.data.table(read.csv("./data/example_1.csv"))
data[, id := 1:dim(data)[1]]

position_rand <- data[,.(x,y)]
positive_rand = (data$test == 1)
test_rand <- data[,5]

# Plot data
graph_base <- ggplot(data, aes(x = x, y = y, color = as.factor(test))) +
  geom_point(size = 2.5) +
  labs(title = "Distribution of Positive and Negative Cases")
print(graph_base)

position = position_rand; positive = positive_rand; test = test_rand$test; min_neighbours = min_neighbours

# Get minimum distance & quantiles of it
min_distances <- get_min_distances(position[,.(x,y)], positive, min_neighbors)
hist(min_distances)


quantiles_dist <- quantile(min_distances,probs = c(0.25,0.5,0.75),na.rm = TRUE)
print(quantiles_dist)

# Determine the level of significance of p-values based on the distribution of distances
sign_level <- c()
for(quant in 1:length(quantiles_dist)){
  categories <- catalogue(
    data.table::copy(position), test, quantiles_dist[quant], cluster_id = NULL, min_neighbours = min_neighbors,
    method = 'base')
  
  sign_level[quant] <- length(which(categories$pval_cluster <= 0.05))
}
print("Number of focci with significance level lower than 0.05 by using the 25th, 50th & 75th percentile:")
print(sign_level)

if(all(sign_level == 0)){
  print(quantiles_dist[2])
  result <- quantiles_dist[2]
}else{
  print(quantiles_dist[min(which(sign_level == max(sign_level)))])
  result <- quantiles_dist[min(which(sign_level == max(sign_level)))]
}
