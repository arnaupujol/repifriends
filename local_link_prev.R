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
#######        INTEGRATION OF LOCAL DISTRIBUTION             #######
####################################################################
setwd("/Users/ericmatamoros/Desktop/TFM/repifriends/")
df <- as.data.table(read.csv("./data/example_5.csv"))
df[, id := 1:dim(df)[1]]

link_d = 0.05
min_neighbours = 2
link_neighbours = 5

storePDF = TRUE
pdf_path = "./pdfs/"
if(storePDF){
  pdf_name <- paste0(pdf_path, "link_neighbors.pdf")
  pdf(pdf_name, width = 11, height = 8.5)
}

# Base plot
graph_base <- ggplot(df, aes(x = x, y = y, color = as.factor(test))) +
  geom_point(size = 3) +
  labs(title = "Distribution of Positive and Negative Cases")
print(graph_base)

# RUSE LINKING DISTANCE METHOD
cat_link <- catalogue(df[,.(x,y)], df[,.(test)], link_d = link_d,
                      min_neighbours=min_neighbours, use_link_d = TRUE)

graph_link <- scatter_pval(
  df[,.(x,y)], 
  cat_link$cluster_id, 
  (df$test == 1), 
  cat_link$epifriends_catalogue,
  c(0,1),
  c(0,1),
  paste0('P-value each cluster - link_dist: ', link_d))

print(graph_base + graph_link)

## USE MINIMUM NEIGHBORS METHOD
cat_neighbours <- catalogue(df[,.(x,y)], df[,.(test)], link_d = link_d,
                      min_neighbours=min_neighbours,link_neighbours= link_neighbours,
                      use_link_d = FALSE)

graph_neigh <- scatter_pval(
  df[,.(x,y)], 
  cat_neighbours$cluster_id, 
  (df$test == 1), 
  cat_neighbours$epifriends_catalogue,
  c(0,1),
  c(0,1),
  paste0('P-value each cluster - link_neigh: ', link_neighbours))

print(graph_base + graph_neigh)

# Plot both detected clusters
print(graph_link + graph_neigh)


# Plot histogram distribution
#hist_link <- size_histogram(cat_link)
#hist_neigh <- size_histogram(cat_neighbours)


mean_link <- rbindlist(cat_link$epifriends_catalogue$mean_position_all)
mean_link[, pvalue := cat_link$epifriends_catalogue$p]
mean_link[, epifriends := 1:nrow(mean_link)]
mean_link[, method := paste0("link_d_", epifriends)]

mean_neigh <- rbindlist(cat_neighbours$epifriends_catalogue$mean_position_all)
mean_neigh[, pvalue := cat_neighbours$epifriends_catalogue$p]
mean_neigh[, epifriends := 1:nrow(mean_neigh)]
mean_neigh[, method := paste0("link_neigh_", epifriends)]

graph_link <- ggplot(mean_link, aes(x = method, y = pvalue,  fill = pvalue)) +
  geom_bar(stat = "identity") +
  labs(title = "P-Value for Linking Distance method", 
       x = "Epifriends-Methodology", y = "P-Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

graph_link <- graph_link + geom_hline(yintercept = 0.05, linetype = "dashed")

graph_neigh <- ggplot(mean_neigh, aes(x = method, y = pvalue,  fill = pvalue)) +
  geom_bar(stat = "identity") +
  labs(title = "P-Value for Linking Neighbors method", 
       x = "Epifriends-Methodology", y = "P-Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

mean_neigh <- graph_neigh + geom_hline(yintercept = 0.05, linetype = "dashed")

print(graph_link + mean_neigh)

if(storePDF){
 dev.off()
}
