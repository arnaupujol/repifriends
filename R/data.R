#' Dataset with random distribution of disease
#'
#' Random distribution of positive cases
#'
#' @format ## `datarand`
#' A data frame with 7,240 rows and 60 columns:
#' \describe{
#'   \item{X, id}{Index or reference}
#'   \item{x, y}{Coordinates}
#'   \item{test}{Positive or negative tests}
#'   \item{case_count}{Number of tests with the given id}
#'   \item{dates}{Date of the test. Fixed dates}
#'   \item{geometry}{Geometrical way to express x & y coordinates}
#' }
#' @source Mock up dataset used for example used in epifriends_on_different_distributions.Rmd
"datarand"

#' Dataset with sinusoidal distribution of disease
#' 
#' Sinusoidal distribution of positive cases
#'
#' @format ## `datasin`
#' A data frame with 7,240 rows and 60 columns:
#' \describe{
#'   \item{X, id}{Index or reference}
#'   \item{x, y}{Coordinates}
#'   \item{test}{Positive or negative tests}
#'   \item{case_count}{Number of tests with the given id}
#'   \item{dates}{Date of the test. Fixed dates}
#'   \item{geometry}{Geometrical way to express x & y coordinates}
#' }
#' @source Mock up dataset used for example used in epifriends_on_different_distributions.Rmd
"datasin"

#' Dataset with clustered distribution of disease
#' 
#' Different clusters of positive cases distributed along negative ones.
#'
#' @format ## `dataclus`
#' A data frame with 7,240 rows and 60 columns:
#' \describe{
#'   \item{X, id}{Index or reference}
#'   \item{x, y}{Coordinates}
#'   \item{test}{Positive or negative tests}
#'   \item{case_count}{Number of tests with the given id}
#'   \item{dates}{Date of the test. Fixed dates}
#'   \item{geometry}{Geometrical way to express x & y coordinates}
#' }
#' @source Mock up dataset used for example used in epifriends_on_different_distributions.Rmd
"dataclus"

#' Dataset with random distribution of positive disease cases distributed against different
#' timeframes.
#' 
#' Different clusters of positive cases distributed along negative ones with different
#' temporalities to try EpiFRIenDs on different timestamps.
#'
#' @format ## `datarand_temp`
#' A data frame with 7,240 rows and 60 columns:
#' \describe{
#'   \item{X, id}{Index or reference}
#'   \item{x, y}{Coordinates}
#'   \item{test}{Positive or negative tests}
#'   \item{case_count}{Number of tests with the given id}
#'   \item{date}{Date of the test. There are different dates}
#'   \item{geometry}{Geometrical way to express x & y coordinates}
#' }
#' @source Mock up dataset used for example used in temporal_analysis.Rmd
"datarand_temp"

#' Dataset with different local-fluctuations on the prevalence of the disease
#' 
#' The twho hemispheres along the Y-axis have clusters made with different fixed 
#' prevalences. Additionally, there is a new cluster in each hemisphere with the opposite
#' prevalence with the aim to add local deviations.
#'
#' @format ## `mock_data_local`
#' A data frame with 7,240 rows and 60 columns:
#' \describe{
#'   \item{X, id}{Index or reference}
#'   \item{x, y}{Coordinates}
#'   \item{test}{Positive or negative tests}
#'   \item{case_count}{Number of tests with the given id}
#'   \item{date}{Date of the test. Fixed dates}
#'   \item{geometry}{Geometrical way to express x & y coordinates}
#' }
#' @source Mock up dataset used for example used in local_prevalence.Rmd
"mock_data_local"

#' Dataset containing clusters of population with randomly distributed positive counts.
#' 
#' Some of the generated clusters contain a higher prevalence of the disease and others 
#' have less prevalence. The aim of this dataset is to mimic the differences between urban
#' & rural regions in real data.
#' 
#' @format ## `mock_data_linkneigh`
#' A data frame with 7,240 rows and 60 columns:
#' \describe{
#'   \item{X, id}{Index or reference}
#'   \item{x, y}{Coordinates}
#'   \item{test}{Positive or negative tests}
#'   \item{case_count}{Number of tests with the given id}
#'   \item{date}{Date of the test. Fixed dates}
#'   \item{geometry}{Geometrical way to express x & y coordinates}
#' }
#' @source Mock up dataset used for example used in link_neighbours.Rmd
"mock_data_linkneigh"

#' Dataset used for the example of False Detections
#' 
#' This dataset contains differently generated clusters with random distribution of
#' positive cases & is being used for the False Detection example.
#' 
#' @format ## `mock_false_det`
#' A data frame with 7,240 rows and 60 columns:
#' \describe{
#'   \item{X, id}{Index or reference}
#'   \item{x, y}{Coordinates}
#'   \item{test}{Positive or negative tests}
#'   \item{case_count}{Number of tests with the given id}
#'   \item{date}{Date of the test. Fixed dates}
#'   \item{geometry}{Geometrical way to express x & y coordinates}
#' }
#' @source Mock up dataset used for example used in ldetect_false_positives.Rmd
"mock_false_det"