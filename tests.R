
library(data.table)
library(readr)
library(factoextra)
library(NbClust)
library(doParallel)
library(reshape2)

link_d_to_analyse = 0.02
min_neighbours = 2

# Load data
# Prevalence defined as number of positives versus total number, where total number
# can either be a cluster, a radious or all the population size.
datarand <- as.data.table(read_csv("data/mock_data_sin.csv"))

position_rand <- datarand[,3:4]
positive_rand = (datarand$test == 1)
test_rand <- datarand[,5]


#df <- define_loc_prev(position_rand, test_rand, method = 'kmeans')
cat_rand <- catalogue(
  position_rand, test_rand, link_d, cluster_id = NULL, min_neighbours = min_neighbours)