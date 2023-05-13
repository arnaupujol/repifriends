
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

setwd("/Users/ericmatamoros/Desktop/TFM/repifriends/")
source("./utils.R")
# Load data
# Prevalence defined as number of positives versus total number, where total number
# can either be a cluster, a radious or all the population size.

##############################################################################
##############         OVERVIEW AT DISTRIBUTIONS OF FILES.      ##############
##############################################################################
files <- list.files("./data/", pattern = "example")

counter <- 1
for(file in files){
  data_v1 <- as.data.table(read_csv(paste0("./data/", file)))
  data_v1[, id := 1:nrow(data_v1)]
  
  position_rand <- data_v1[,3:4]
  positive_rand = (data_v1$test == 1)
  test_rand <- data_v1$test
  
  # Plot data
  graph_base <- ggplot(data_v1, aes(x = x, y = y, color = as.factor(test))) +
    geom_point(size = 3) +
    labs(
      title = "Distribution of Positive and Negative Cases",
      subtitle = paste0("example_", counter, ".csv")
      )
  
  print(graph_base)
  
  file.rename(
    paste0("./data/", file), 
    paste0("./data/example_", counter, ".csv")
    )
  counter <- counter + 1
  
}

##############################################################################
##############                     RUN ANALYSIS                 ##############
##############################################################################
files <- list.files("./data/", pattern = "example")
link_d_to_analyse = c(0.01, 0.02, 0.05, 0.1)
min_neighbours = 2
methods_list <- c("kmeans", "radial", "centroid", "base")
store_PDF = TRUE
store_individually = TRUE
pdf_path = "./pdfs/"
automatic_link_d = FALSE

for(file in files){
  
  catalogue_methods <- list()
  for(method in methods_list){
    
    # Get linking distance
    if(automatic_link_d){
      link_d_to_analyse <- get_link(file)
    }
    
    for(link_d in link_d_to_analyse){
      if (store_PDF){
        pdf_name <- paste0(pdf_path, strsplit(file, ".csv")[[1]],"_",method,"_link_d_", link_d, ".pdf")
        pdf(pdf_name, width = 11, height = 8.5)
      }
      
      #########################
      ####### READ DATA #######
      #########################
      data_v1 <- as.data.table(read_csv(paste0("data/", file)))
      data_v1[, id := 1:nrow(data_v1)]
      
      position_rand <- data_v1[,3:4]
      positive_rand = (data_v1$test == 1)
      test_rand <- data_v1$test
      
      ####################################################################
      ####### GENERATE PLOTS OF DISTRIBUTION AND DETECTED CLUSTERS #######
      ####################################################################
      # Plot data
      graph_base <- ggplot(data_v1, aes(x = x, y = y, color = as.factor(test))) +
        geom_point(size = 2.5) +
        labs(title = "Distribution of Positive and Negative Cases")
      
      #print(graph_base)
      
      # Compute Epifriends
      categories <- catalogue(
        position_rand, test_rand, link_d, cluster_id = NULL, min_neighbours = min_neighbours,
        method = method)
      catalogue_methods[[paste0(method,"_", link_d)]] <-  categories
      
        
        graphs_prev <- scatter_pval(
            data.table::copy(position_rand), 
            categories$cluster_id, 
            data.table::copy(positive_rand), 
            data.table::copy(prevalence),
            categories$epifriends_catalogue,
            paste0(method,"-Link_d: ",link_d)) 
        
        graphs_no_prev <- scatter_pval(
          position_rand, 
          categories$cluster_id, 
          positive_rand, 
          NULL,
          categories$epifriends_catalogue,
          paste0(method,"-Link_d: ",link_d)) 
        histo <- size_histogram(categories)
        
        print(graph_base + graph_clusters)
        print(graph_clusters + graph_sign_clusters)
        print(graph_clusters+graphs_prev)
        #print(graph_clusters + graphs_no_prev)
        print(histo)
        
        if (store_PDF){
          dev.off()
        }
      }else{
        file.remove(pdf_name)
      }
    }
    
  }
  
  for(link_d in link_d_to_analyse){
    if (store_PDF){
      pdf_name <- paste0(pdf_path, strsplit(file, ".csv")[[1]],"_link_d_", link_d, "_GENERAL_INSIGHTS.pdf")
      pdf(pdf_name, width = 11, height = 8.5)
    }
    all_combs <- names(catalogue_methods)
    links_combs <- as.numeric(sapply(all_combs, function(x){ strsplit(x, "_")[[1]][2]}))
    general <- get_all_results(catalogue_methods, names(catalogue_methods)[which(links_combs == link_d)])
    
    if(is.null(general$table$pvalue)){
      if (store_PDF){
        dev.off()
        file.remove(pdf_name)
      }
    }else{
      print(general$chart)
      if (store_PDF){
        dev.off()
      }
    }
  }
}

if(store_individually & store_PDF){
  store_indiv(pdf_path, files)
}

