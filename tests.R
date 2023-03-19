
library(data.table)
library(readr)
library(factoextra)
library(NbClust)
library(doParallel)
library(reshape2)
library(pdftools)

store_PDF = TRUE

# Load data
# Prevalence defined as number of positives versus total number, where total number
# can either be a cluster, a radious or all the population size.

files <- c("mock_dataframe_clustered_v1.csv", "mock_dataframe_clustered_v2.csv", "mock_dataframe_clustered_v3.csv")

link_d_to_analyse = 0.05
min_neighbours = 2
methods_list <- c("kmeans", "radial", "centroid", "base")
store_PDF = TRUE

for(file in files[1]){
  if (store_PDF){
    pdf_name <- paste0("pdfs/", strsplit(file, ".csv")[[1]], ".pdf")
    pdf(pdf_name, width = 11, height = 8.5)
  }
  for(method in methods_list){
    data_v1 <- as.data.table(read_csv(paste0("data/", file)))
    
    position_rand <- data_v1[,3:4]
    positive_rand = (data_v1$test == 1)
    test_rand <- data_v1[,5]
    
    # Plot data
    graph_base <- ggplot(data_v1, aes(x = x, y = y, color = as.factor(test))) +
      geom_point(size = 2.5) +
      labs(title = "Distribution of Positive and Negative Cases")
    
    # Compute Epifriends
    categories <- catalogue(
      position_rand, test_rand, link_d, cluster_id = NULL, min_neighbours = min_neighbours,
      method = method)
    graphs <- scatter_pval(
      position_rand, 
      categories$cluster_id, 
      positive_rand, 
      categories$epifriends_catalogue,
      method) 
    histo <- size_histogram(categories)
    print(graph_base+graphs+histo)
  }
  
  if (store_PDF){
    dev.off()
  }
}
