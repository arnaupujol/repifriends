# repifriends
R package for Epidemiological Foci Relating Infections by Distance (EpiFRIenDs).
# Epidemiological Foci Relating Infections by Distance (EpiFRIenDs)

> Author: **Mikel Majewski Etxeberria**  
> Year: **2022**  
> Version: **1.0**  

This repository contains the EpiFRIenDs software to detect and analyse foci
(clusters, outbreaks or hotspots) of infections from a given disease.

The package is still in work.

Software requirements:
----------------------
All the packages that are required so that all the codes can run are:
- RANN
- chrone

Installation:
----------------------
To install the repository, you open your RStudio. Then in the packages window 
in the right above corner you click on install and you select the option 
"Package Archive File (.zip; .tar.gz)".

Then you select the "epifriends_0.1.0.tar.gz" or "epifriends_0.1.0.zip" on 
the package folder and epifriends will be instaled.

This version have been proven to work for R 4.1.1 version.

Structure of the repository:
----------------------------
The repository contains two main directories:
- packages/epifriends: where the code is stored
- examples: where some examples of the software implementation are shown in
Rmarkdowns.

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
this method creates a list of EpiFRIenDs catatalogues(each element of the  
list representing a time frame) and the posiblity of adding some other 
variables like a temporal id using some linking times and distances, 
and in this way assigning a temporal ID to the foci, assigning the same temporal ID
to foci form different time frames when they are close in time and space. 
Some other posible variables are those relate with the lifetimes of those foci.

#### Examples:

Examples can be found in the directory `examples`, in the following Rmarkdowns:

- epifriends_on_different_distributions.Rmd: this Rmarkdown uses three
sets of artificial data generated in the python version 
and shows how to run EpiFRIenDs to detect foci on them.

- temporal_analysis.Rmd (still in work): this Rmarkdown uses one 
artificial data set generates in the python version and show how to 
work the temporal catalogues on them.
