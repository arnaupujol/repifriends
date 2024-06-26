% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/false_detections.R
\name{get_false_detection}
\alias{get_false_detection}
\title{Remove false detections from the clusters identified by the catalogue}
\usage{
get_false_detection(
  positions,
  test_result,
  link_d,
  use_link_d = TRUE,
  link_neighbours = 10,
  n_simulations = 500,
  min_neighbours = 2,
  max_p = 1,
  min_pos = 2,
  min_total = 2,
  min_pr = 0,
  keep_null_tests = FALSE,
  run_seq = FALSE
)
}
\arguments{
\item{positions}{data.frame with the positions of parameters we want to query 
with shape (n,2) where n is the number of positions.}

\item{test_result}{data.frame with the test results (0 or 1).}

\item{link_d}{The linking distance to connect cases. Should be in the same 
scale as the positions.}

\item{n_simulations}{Number of random simulations to be performed.}

\item{min_neighbours}{Minium number of neighbours in the radius < link_d needed 
to link cases as friends.}

\item{max_p}{Maximum value of the p-value to consider the cluster detection.}

\item{min_pos}{Threshold of minimum number of positive cases in clusters applied.}

\item{min_total}{Threshold of minimum number of cases in clusters applied.}

\item{min_pr}{Threshold of minimum positivity rate in clusters applied.}

\item{keep_null_tests}{Whether to remove or not missings. If numeric, provide 
value to impute.}

\item{use_link_d:}{If True, use linking distance to determine the closest neighbors. 
If False, use default linking neighbours based on proximity.}

\item{min_neighbours:}{Minium number of neighbours in the radius < link_d needed 
to link cases as friends.}

\item{link_neighbours:}{Number of surrounding neighbors to link.}

\item{cluster_id}{Numeric vector with the cluster IDs of each position, with 0
for those without a cluster. Give NULL if cluster_id must be calculated.}

\item{run_seq:}{If TRUE, run sequentially instead of parallelizing.}
}
\value{
data.frame with the number of isolated positives detected per iteration and its frequency.
}
\description{
Remove false detections from the clusters identified by the catalogue
}
\examples{
# Creation of x vector of longitude coordinates, y vector of latitude coordinates 
# and finally merge them on a position data frame.
x <- c(1,2,3,4,7.5,8,8.5,9,10,13,13.1,13.2,13.3,14,15,30)
y <- c(1,2,3,4,7.5,8,8.5,9,10,13,13.1,13.2,13.3,14,15,30)
pos <- data.frame(x,y)

# Creation of test data frame with 0 for negative cases and 1 for positive classes 
# for each position.
test <- c(0,1,1,0,1,0,1,1,0,0,0,0,1,0,1,1)

# Creation of catalogue for this positions, linking distance 2 and default values.
cat <- get_false_detection(positions = pos, test_result = test, link_d = 2, 
n_simulations = 2, run_seq=TRUE)

}
\author{
Eric Matamoros
}
