# Epidemiological Foci Relating Infections by Distance (EpiFRIenDs)

> Author: **Eric Matamoros, Mikel Majewski, Arnau Pujol**  
> Year: **2023**  
> Version: **2.0**  

This repository contains the R package for Epidemiological Foci Relating Infections by Distance (EpiFRIenDs), a software to detect and analyse foci (clusters, outbreaks or hotspots) of infections from a given disease.

This software is fully open source and all are welcome to use or modify it for any purpose. We would kindly request that any scientific publications making use of this software cite **Pujol A., Brokhattingen N., Matambisso G., et al (in prep.)**.

Software requirements:
----------------------
All the packages that are required so that all the codes can run are listed below and can be found in the
DESCRIPTION file:
- R (>= 3.4.0)
- RANN (>= 2.6.1)
- data.table (>= 1.14.8)
- ggplot2 (>= 3.4.1)

The following are packages that need to be imported:
- cluster (>= 2.1.4)
- clValid (>= 0.7)
- chron (>= 2.3-61)
- doParallel (>= 1.0.17)
- foreach (>= 1.5.2)
- ggmap (>= 3.0.2)
- magick (>= 2.7.4)
- sf (>= 1.0-12)
- stats (>= 4.2.2)


Installation:
----------------------
To install the repository, open the RStudio and move to the working directory where the project is located. Actions needed are:
1. Build -> Build Source Package. This will install the binaries required for the installation. It will ask for dependencies installation if missing.
2. Build -> Install Package. This will install the EpiFRIenDs package. 

A successful message will automatically appear after the package is installed resulting from the command "library(repifriends)".

A detailed guide on the full setup & local deployment of the RShiny application can be found in app/run_shiny_app.pdf file


Structure of the repository:
----------------------------
The repository contains 4 main directories:
- *data*: Folder where data is stored. Contains several mock-up examples used across the examples.
- *vignettes*: Folder with some examples of the software implementation are shown in Rmarkdowns. 
- *R*: Contains the core-code of EpiFRIenDs, with functions segmented based on their utilities.
- *app*: Folder with the front-end R-shiny application.

Structure of the different functions:
----------------------------
The different functions used across EpiFRIenDs are distributed between different files with the extension .R & located in the R folder. These files are segmented based to the topic they cover . Find below a description of each file:

- *dbscan*: from some input position, linking distance and minimum number of neighbours, this function finds DBSCAN clusters and assigns a cluster ID for each position, with 0 meaning that the position does not belong to any cluster.
- *catalogue*: from some input positions, test results, linking distance and minimum number of neighbours, this function detects the EpiFRIenDs foci and outputs a catalogue of them with its associated data.
- *temporal_catalogue*: from a EpiFRIenDs catalogue and a dates vector this method creates a list of EpiFRIenDs catatalogues (each element of the   list representing a time frame) and the posiblity of adding some other  variables like a temporal id using some linking times and distances,  and in this way assigning a temporal ID to the foci, assigning the same temporal ID to foci form different time frames when they are close in time and space. Some other possible variables are those relate with the lifetimes of those foci.
- *automate_link_d*: from some input positions & test_results, it automatically computes the best linking distance based on an optimization of the maximum number of real positive tests found under different statistical representative linking distances.
- *convert_coords*: converts longitude & latitude positions into real coordinates based on a defined projection.
- *define_local_prev*: different methodologies used to account for the local prevalence of a specific region. If method is kmeans it will use the prevalence detected by the different k-clusters. If method is centroid, it will expand the centroid detected by the dbscan significant clusters until a certain criteria is met. If method is base, global prevalence will be used.
- *plots*: different plotting functionalities used accross the code and specially in the vignettes & shiny application.
- *false_detections*: code used to generate random distributions of a given prevalence with the goal to detect random distributions and help in correcting the p-value through accounting for this randomness.
- *clean_data*: code used to impute or remove missing values in the input data.

#### Examples:

Examples can be found in the directory `vignettes`, in the following Rmarkdowns:

- *epifriends_on_different_distributions.Rmd*: this Rmarkdown uses three sets of artificial data generated in the python version and shows how to run EpiFRIenDs to detect foci on them.

- *temporal_analysis.Rmd*: this Rmarkdown uses one artificial data set generates in the python version and shows how to work the temporal catalogues on them.

- *detect_false_positives.Rmd*: this Rmarkdown intends to guide the user on how the approach to detect False Positives work and how it can be leveraged to obtain to adjust the p-value of the different clusters.

- *local_prevalence.Rmd*: this Rmarkdown uses 2 different methodologies to account for the local prevalence of a certain region and compares it with the base methodology, which only considers the global prevalence.

- *optimize_link_d.Rmd*: this Rmarkdown shows how to optimize the linking_distance (distance used to detect closest neighbors) through a statistical optimization of the number of real positives detected.

- *link_neighbours.Rmd*: this Rmarkdown shows how to adapt to the local density of population by using individual linking distances based on the proximity of the neighbors instead of a local linking distance across all coordinates. 
