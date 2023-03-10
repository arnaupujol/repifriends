# This source code is licensed under the GNU GENERAL PUBLIC LICENSE license found in the
# LICENSE file in the root directory of this source tree.

#' This method imputes or removes missing data
#'
#' @param positions data.table with the positions of parameters we want to query with shape (n,2) where n is the number of positions.
#' @param keep_null_tests: Whether to remove or not missings. If numeric, provide value to impute.
#' @param verbose: If TRUE, print information of the process; else, do not print.
#'
#' @details This function cleans the position data if contains missing values. User can choose whether to remove the missings or impute them with a specific value.
#'
#' @return Data.table with cleaned/imputed variables based on user choose.
#' 
#' @export
#' 
#' @author Eric Matamoros Morales based on earlier python code by Arnau Pujol.
#'
#' @examples
#' # Required packages
#' if(!require("data.table")) install.packages("data.table")
#' library("data.table")
#'
#' # Creation of x vector of longitude coordinates, y vector of latitude coordinates and finaly merge them on a position data frame.
#' x <- c(1,2,3,4,7.5,8,8.5,9,10,NA,13.1,13.2,13.3,14,15,30)
#' y <- c(1,2,3,4,7.5,8,8.5,9,10,13,13.1,13.2,NA,14,15,30)
#' pos <- data.table(x,y)
#'
#' # Computation of clusters of hotspots for positions with dbscan algorithm using linking distance 2 and minimum 3 neighbours.
#' db <- clean_unknown_data(pos)
#' 
clean_unknown_data <- function(
    positions, 
    test = NULL,
    keep_null_tests = TRUE, 
    verbose = FALSE,
    cols_remove = c("x", "y")){
  
  #Change infinites to missings
  positions[sapply(positions, is.infinite)] <- NA

  if(is.null(test)){
    positions <- positions[!is.na(rowSums(positions[,..cols_remove]))]
    return(positions)
  }else{
    col_impute = "test"
    positions$test <- test
    positions <- positions[!is.na(rowSums(positions[,..cols_remove]))]
    if(is.numeric(keep_null_tests)){
      if(verbose){print(paste0(
          "Replacing missing positions with value: ",as.character(keep_null_tests))
        )}
      positions <- data.table::setnafill(positions, fill=keep_null_tests)
    }else if(keep_null_tests == FALSE){
        positions <- positions[!is.na(rowSums(positions[,..col_impute]))]
    }else if(keep_null_tests == TRUE){
        if(verbose){print("Missing values in positions kept as NULL")}
    }else{
        stop("Argument keep_null_tests has an invalid argument")
    }
  }
  
  test <- positions$test
  positions[, test := NULL]
  
  return(list("position" = positions, "test" = test))
  
}
