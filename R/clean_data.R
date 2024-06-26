# This source code is licensed under the GNU GENERAL PUBLIC LICENSE license 

#' This method removes missing data in the coordinates and imputes or removes the values
#' over the tests column
#'
#' @param positions data.table variables to be imputed.
#' @param cols_impute Vector of columns that want to be imputed based on 
#' keep_null_tests parameter.
#' @param keep_null_tests: Whether to remove or not missings. If numeric, provide 
#' value to impute.
#' @param verbose: If TRUE, print information of the process; else, do not print.
#'
#' @return Data.table with cleaned/imputed variables based on user choose.
#' 
#' @importFrom data.table setnafill
#' 
#' @export
#' 
#' @author Eric Matamoros Morales based on earlier python code by Arnau Pujol.
#'
#' @examples
#' # Creation of x vector of longitude coordinates, y vector of latitude coordinates 
#' # and finally merge them on a position data frame.
#' x <- c(1,2,3,4,7.5,8,8.5,9,10,NA,13.1,13.2,13.3,14,15,30)
#' y <- c(1,2,3,4,7.5,8,8.5,9,10,13,13.1,13.2,NA,14,15,30)
#' pos <- data.table(x,y)
#'
#' # Computation of clusters of hotspots for positions with dbscan algorithm using 
#' # linking distance 2 and minimum 3 neighbours.
#' db <- clean_unknown_data(pos)
#' 
clean_unknown_data <- function(
    positions, 
    cols_impute = NULL,
    keep_null_tests = TRUE, 
    verbose = FALSE,
    cols_remove = c("x", "y")){
  
  #Change infinites to missings
  positions[sapply(positions, is.infinite)] <- NA

  if( (is.null(cols_impute)) | (length(cols_impute) == 0)){
    positions <- positions[!is.na(rowSums(positions[,..cols_remove]))]
    if(nrow(positions) == 0){
      stop("All coordinates have either X or Y missing. Please check the data")
    }
    return(positions)
  }else{
    positions <- positions[!is.na(rowSums(positions[,..cols_remove]))]
    for(col_impute in cols_impute){
      
      if(col_impute %in% colnames(positions)){
        if(is.numeric(keep_null_tests)){
          if(verbose){print(paste0(
              "Replacing missing positions of ", col_impute," with value: ",
              as.character(keep_null_tests))
            )}
          positions <- setnafill(positions, fill=keep_null_tests, cols = col_impute)
        }else if(keep_null_tests == FALSE){
            positions <- positions[!is.na(rowSums(positions[,..col_impute]))]
        }else if(keep_null_tests == TRUE){
            if(verbose){print(paste0("Missing values for column : " , col_impute, 
                                     " kept as NULL"))}
        }else{
            stop("Argument keep_null_tests has an invalid argument")
        }
      }else{
        print(paste0("Column : ", col_impute, "does not exist."))
      }
    }
  }
  
  return(positions)
  
}

#' Method to clean an input data.table with automatically computing the columns 
#' to be imputed.
#'
#' @param positions data.table variables to be imputed.
#' @param keep_null_tests: Whether to remove or not missings. If numeric, provide 
#' value to impute.
#'
#' @return Data.table with cleaned/imputed variables based on user choose.
#' 
#' @export
#' 
#' @author Eric Matamoros Morales based on earlier python code by Arnau Pujol.
#'
#' @examples
#' # Creation of x vector of longitude coordinates, y vector of latitude coordinates 
#' # and finally merge them on a position data frame.
#' x <- c(1,2,3,4,7.5,8,8.5,9,10,NA,13.1,13.2,13.3,14,15,30)
#' y <- c(1,2,3,4,7.5,8,8.5,9,10,13,13.1,13.2,NA,14,15,30)
#' test <- c(1,0,1,1,0,NA,NA,0,1,1,1,0,0,1,1,0)
#' pos <- data.table(x,y)
#'
#' # Computation of clusters of hotspots for positions with dbscan algorithm 
#' # using linking distance 2 and minimum 3 neighbours.
#' db <- clean_data(pos, FALSE)
#' 
clean_data <- function(positions, keep_null_tests = FALSE){
  to_impute <- colnames(positions)[!(colnames(positions) %in% c("x", "y"))]
  positions = clean_unknown_data(positions,to_impute,keep_null_tests,FALSE)
  return(positions)
}
