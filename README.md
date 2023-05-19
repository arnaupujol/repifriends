# Epidemiological Foci Relating Infections by Distance (EpiFRIenDs)

> Author: **Mikel Majewski Etxeberria, Arnau Pujol**  
> Year: **2022**  
> Version: **1.0**  

This repository contains the R package for Epidemiological Foci Relating Infections by Distance (EpiFRIenDs), a software 
to detect and analyse foci (clusters, outbreaks or hotspots) of infections from a given disease.

This software is fully open source and all are welcome to use or modify it for any purpose.
We would kindly request that any scientific publications making use of this software cite **Pujol A., Brokhattingen N., Matambisso G., et al (in prep.)**.

The package is still work in progress, consider this a Beta version. The first official release is comming very soon.
You can find [here](https://github.com/arnaupujol/epifriends) a version of EpiFRIenDs in Python. 


Software requirements:
----------------------
All the packages that are required so that all the codes can run are:
- RANN
- chron

Installation:
----------------------
To install the repository, you open your RStudio. Then in the packages window 
in the bottom right corner you click on install and you select the option 
"Package Archive File (.zip; .tar.gz)".

Then you select the "epifriends_0.1.0.tar.gz" or "epifriends_0.1.0.zip" on 
the package folder and epifriends will be instaled.

This version have been proven to work for R 4.1.1 version.

Structure of the repository:
----------------------------
The repository contains 4 main directories:
- data: where data is stored
- vignettes:where some examples of the software implementation are shown in
Rmarkdowns. 
- R: core code of EpiFRIenDs, with functions segmented by their functionalities.
- app: front-end R-shiny application

How to use it:
----------------------------
The files dbscan.R, catalogue.R, temporal_catalogue.R from the directory 
packages/epifriends/R contains all the methods (functions) that can be called within EpiFRIenDs. 
The main methods of the algorithm are:
- dbscan: from some input position, linking distance and minimum number of
neighbours, this function finds DBSCAN clusters and assigns a cluster ID for
each position, with 0 meaning that the position does not belong to any
cluster.
- catalogue: from some input positions, test results, linking distance and
minimum number of neighbours, this function detects the EpiFRIenDs foci and
outputs a catalogue of them with its associated data.
- temporal_catalogue: from a EpiFRIenDs catalogue and a dates vector
this method creates a list of EpiFRIenDs catatalogues (each element of the  
list representing a time frame) and the posiblity of adding some other 
variables like a temporal id using some linking times and distances, 
and in this way assigning a temporal ID to the foci, assigning the same temporal ID
to foci form different time frames when they are close in time and space. 
Some other posible variables are those relate with the lifetimes of those foci.
- automate_link_d: from some input positions & test_results, it automatically computes
the best linking distance based on an optimization of the maximum number of real positive
tests found under different statistical representative linking distances.
- convert_coords: converts longitude & latitude positions into real coordinates based on
a defined projection.
- define_local_prev: different methodologies used to account for the local prevalence of
a specific region. If method is kmeans it will use the prevalence detected by the different
k-clusters. If method is centroid, it will expand the centroid detected by the dbscan
significant clusters until a certain criteria is met. If method is base, global prevalence will
be used.
- plots: different plotting functionalities used accross the code and specially in the vignettes &
shiny application.
- false_detections: code used to generate random distributions of a given prevalence with the goal
to detect random distributions and help in correcting the p-value through accounting for this
randomness.
- clean_data: code used to impute or remove missing values in the input data.

#### Examples:

Examples can be found in the directory `vignettes`, in the following Rmarkdowns:

- epifriends_on_different_distributions.Rmd: this Rmarkdown uses three
sets of artificial data generated in the python version 
and shows how to run EpiFRIenDs to detect foci on them.

- temporal_analysis.Rmd: this Rmarkdown uses one 
artificial data set generates in the python version and shows how to 
work the temporal catalogues on them.

- detect_false_positives.Rmd: this Rmarkdown intends to guide the user on how the
approach to detect False Positives work and how it can be leveraged to obtain to
adjust the p-value of the different clusters.

- local_prevalence.Rmd: this Rmarkdown uses 2 different methodologies to account for the
local prevalence of a certain region and compares it with the base methodology, which only
considers the global prevalence.

- optimize_link_d.Rmd: this Rmarkdown shows how to optimize the linking_distance (distance used
to detect closest neighbors) through a statistical optimization of the number of real positives
detected
