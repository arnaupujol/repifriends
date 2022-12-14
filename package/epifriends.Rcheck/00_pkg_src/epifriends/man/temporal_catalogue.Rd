% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/temporal_catalogue.R
\name{temporal_catalogue}
\alias{temporal_catalogue}
\title{This method generates a list of EpiFRIenDs catalogues representing different time frames by including only cases within a time window that moves within each time step.}
\usage{
temporal_catalogue(
  positions,
  test_result,
  dates,
  link_d,
  min_neighbours = 2,
  time_width,
  min_date = NULL,
  max_date = NULL,
  time_steps = 1,
  max_p = 1,
  min_pos = 2,
  min_total = 2,
  min_pr = 0,
  add_temporal_id = TRUE,
  linking_time,
  linking_dist,
  get_timelife = TRUE
)
}
\arguments{
\item{positions}{data.frame with the positions of parameters we want to query with shape (n,2) where n is the number of positions.}

\item{test_result}{data.frame with the test results (0 or 1).}

\item{dates}{Vector of the date times of the corresponding data in format 'y-m-d h:m:s'.}

\item{link_d}{the linking distance of the DBSCAN algorithm. Should be in the same scale as the positions.}

\item{min_neighbours}{Minium number of neighbours in the radius < link_d needed to link cases as friends.}

\item{time_width}{Number of days of the time window used to select cases in each time step.}

\item{min_date}{Initial date used in the first time step and time window selection. In format 'y-m-d h:m:s'.}

\item{max_date}{Final date to analyse, defining the last time window as the one fully overlapping the data. In format 'y-m-d h:m:s'.}

\item{time_steps}{Number of days that the time window is shifted in each time step.}

\item{max_p}{Maximum value of the p-value to consider the cluster detection.}

\item{min_pos}{Threshold of minimum number of positive cases in clusters applied.}

\item{min_total}{Threshold of minimum number of cases in clusters applied.}

\item{min_pr}{Threshold of minimum positivity rate in clusters applied.}

\item{add_temporal_id}{Boolean that indicates if we want to add a temporal id to the hotspots indicating which hotspots are near in time and space. Only posible if there 2 or more temporal windows.}

\item{linking_time}{Maximum number of timesteps of distance to link hotspots with the same temporal ID. Only necesary if add_temporal_id is TRUE.}

\item{linking_dist}{Linking distance used to link the clusters from the different time frames. Only necesary if add_temporal_id is TRUE.}

\item{get_timelife}{It specifies if the time periods and timelife of clusters are obtained. Only necesary if add_temporal_id is TRUE.}
}
\value{
List with the next objects:
temporal_catalogues: list of catalogues
  List of EpiFRIenDs catalogues, where each element contains the catalogue in each
  time step.
mean_date: vector of chorn date times.
  List of dates corresponding to the median time in each time window
}
\description{
This method generates a list of EpiFRIenDs catalogues representing different time frames by including only cases within a time window that moves within each time step.
}
\details{
The epifriends package uses the RANN package which can gives the exact nearest neighbours using the friends of friends algorithm. For more information on the RANN library please visit https://cran.r-project.org/web/packages/RANN/RANN.pdf
#'
}
\examples{
# Required packages
if(!require("RANN")) install.packages("RANN")
if(!require("chron")) install.packages("chron")
library("RANN")
library("chron")

# Creation of x vector of longitude coordinates, y vector of latitude coordinates and finaly merge them on a position data frame.
x <- c(1,2,3,4,7.5,8,8.5,9,10,13,13.1,13.2,13.3,14,15,30)
y <- c(1,2,3,4,7.5,8,8.5,9,10,13,13.1,13.2,13.3,14,15,30)
pos <- data.frame(x,y)

# Creation of test data frame with 0 for negative cases and 1 for positive clases for each position.
test <- data.frame(c(0,1,1,0,1,0,1,1,0,0,0,0,1,0,1,1))

#Creation of chron dates vector in format 'y-m-d h:m:s'.
dates <- c("2020-11-03 05:33:07","2021-05-19 10:29:59","2021-02-09 14:53:20","2021-11-21 02:35:38","2020-11-19 05:57:24",
"2021-06-09 07:50:30","2021-09-18 05:53:26","2020-03-19 17:16:56","2021-06-08 12:40:46","2020-06-26 05:01:31",
"2020-10-15 04:40:27","2021-05-28 15:23:23","2020-05-01 02:56:54","2020-08-19 22:45:35","2021-10-23 18:56:35",
"2020-10-19 00:01:25")

# Creation of temporal catalogue for this data.
tcat <- tcat <- temporal_catalogue(pos, test, dates ,link_d = 2, time_width = 305, time_steps = 305, linking_time = 3, linking_dist = 2)
}
\author{
Mikel Majewski Etxeberria based on earlier python code by Arnau Pujol.
}
