install.packages("patchwork")

library("RANN")
library("readxl")
library("chron")
library("ggplot2")
library("patchwork")
library("ggpubr")

datarand <- read.csv("data\\mock_data_rand.csv")
datasin <- read.csv("data\\mock_data_sin.csv")
dataclus <- read.csv("data\\mock_data_clustered.csv")


plot(datarand$x, datarand$y, pch = 19, col = as.factor(datarand$test))
plot(datasin$x, datasin$y, pch = 19, col = as.factor(datasin$test))
plot(dataclus$x, dataclus$y, pch = 19, col = as.factor(dataclus$test))

#EpiFRIenDs parameters
link_d = 0.01
min_neighbours = 2

position_rand <- datarand[3:4]
positive_rand = (datarand$test == 1)

position_sin <- datasin[3:4]
positive_sin = (datasin$test == 1)

position_clus <- dataclus[3:4]
positive_clus = (dataclus$test == 1)

#Running DBSCAN for the mock data with random distribution
cluster_id = dbscan(position_rand[positive_rand,], link_d, min_neighbours)

neg <- data.frame(position_rand[!positive_rand,])
pos <- data.frame(position_rand[positive_rand,],cluster_id)

neggraph <- ggplot(neg,aes(x=x, y=y))+
            geom_point(color = "grey")
posgraph <- ggplot(pos,aes(x=x,y=y))+
            geom_point(aes(colour = cluster_id))+
            scale_colour_gradientn(colors=rainbow(7))

neggraph + posgraph

ggarrange(neggraph,posgraph)

plot(datarand$x[!positive_rand], datarand$y[!positive_rand], pch = 19, col = "grey")
plot(datarand$x[positive_rand], datarand$y[positive_rand], pch = 19, col = as.factor(cluster_id))

#Running DBSCAN for the mock data with sinusoidal distribution
cluster_id = dbscan(position_sin[positive_sin,], link_d, min_neighbours)

neg <- data.frame(position_sin[!positive_sin,])
pos <- data.frame(position_sin[positive_sin,],cluster_id)

neggraph <- ggplot(neg,aes(x=x, y=y))+
  geom_point(color = "grey")
posgraph <- ggplot(pos,aes(x=x,y=y))+
  geom_point(aes(colour = cluster_id))+
  scale_colour_gradientn(colors=rainbow(40)) 

ggarrange(neggraph,posgraph)

#Running DBSCAN for the mock data with cluster distribution
cluster_id = dbscan(position_clus[positive_clus,], link_d, min_neighbours)

neg <- data.frame(position_clus[!positive_clus,])
pos <- data.frame(position_clus[positive_clus,],cluster_id)

neggraph <- ggplot(neg,aes(x=x, y=y))+
  geom_point(color = "grey")
posgraph <- ggplot(pos,aes(x=x,y=y))+
  geom_point(aes(colour = cluster_id))+
  scale_colour_gradientn(colors=rainbow(40)) 

ggarrange(neggraph,posgraph)

