# Eric Matamoros Morales (Ver 04-03-2023)

##################################################################################################
# CLEAN MISSING DATA
#################################################################################################
#' This method finds the DBSCAN clusters from a set of positions and returns their cluster IDs.
#'
#' @param positions data.frame with the positions of parameters we want to query with shape (n,2) where n is the number of positions.
#' @param keep_null_tests: Whether to remove or not missings. If numeric, provide value to impute.
#' @param verbose: If TRUE, print information of the process; else, do not print.
#'
#' @details This function cleans the position data if contains missing values. User can choose whether to remove the missings or impute them with a specific value.
#'
#' @return Data.table with cleaned/imputed coordinates based on user choose.
#' @export
#' #'
#'
#' @author Eric Matamoros Morales based on earlier python code by Arnau Pujol.
#'
#' @examples
#' # Required packages
#' if(!require("RANN")) install.packages("RANN")
#' library("RANN")
#'
#' # Creation of x vector of longitude coordinates, y vector of latitude coordinates and finaly merge them on a position data frame.
#' x <- c(1,2,3,4,7.5,8,8.5,9,10,NA,13.1,13.2,13.3,14,15,30)
#' y <- c(1,2,3,4,7.5,8,8.5,9,10,13,13.1,13.2,NA,14,15,30)
#' pos <- data.table(x,y)
#'
#' # Computation of clusters of hotspots for positions with dbscan algorithm using linking distance 2 and minimum 3 neighbours.
#' db <- clean_unknown_data(pos, FALSE, FALSE)
#'
clean_unknown_data <- function(positions, keep_null_tests,verbose){
  # This method removes all the cases with any missing value
  # in either x or y.
  #
  # Parameters:
  # -----------
  # positions: list of class data.table
  #   A list with the position parameters we want to query with shape (n,2),
  #   where n is the number of positions.
  # keep_null_tests: numeric of logical
  #   Whether to remove or not missings. If numeric, provide value to impute.
  # verbose: logical
  #   If TRUE, print information of the process; else, do not print.
  #
  # Returns:
  # --------
  # indeces: List of class data.table
  #   List with imputed or removed missing values.
  
  #Change infinites to missings
  positions[sapply(positions, is.infinite)] <- NA
  
  if(is.numeric(keep_null_tests)){
    if(verbose){print(paste0(
      "Replacing missing positions with value: ",as.character(keep_null_tests))
    )}
    positions <- data.table::setnafill(positions, fill=keep_null_tests)
  }else if(keep_null_tests == FALSE){
    positions <- positions[!(is.na(x) | is.na(y))]
  }else if(keep_null_tests == TRUE){
    if(verbose){print("Missing values in positions kept as NULL")}
  }else{
    stop("Argument keep_null_tests has an invalid argument")
  }
  
  return(positions)
  
}
