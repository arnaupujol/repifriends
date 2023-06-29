# This source code is licensed under the GNU GENERAL PUBLIC LICENSE license 

#' This method transforms the latitude and longitud coordinates from to cartesian 
#' coordinates in meters.
#'
#' @param lon list of longitude positions.
#' @param lat List of latitude positions
#' @param to_epsg EPSG number for the projection to use
#' @param verbose It specifies if information on the process is printed
#'
#' @return  List of coordinates x & y converted to cartesian coordinates
#' 
#' @importFrom sf st_as_sf st_transform st_coordinates
#' 
#' @export
#' 
#' @author Eric Matamoros Morales based on earlier python code by Arnau Pujol.
#'
#' @examples
#' # Creation of x vector of longitude coordinates, y vector of latitude coordinates 
#' # and finally merge them on a position data frame.
#' x <- c(1,2,3,4,7.5,8,8.5,9,10,13,13.1,13.2,13.3,14,15,30)
#' y <- c(1,2,3,4,7.5,8,8.5,9,10,13,13.1,13.2,13.3,14,15,30)
#'
#'
#' # Convert to coordinates in meters
#' coords <- latlon2xy(x, y, 32733)
#' 
latlon2xy <- function(lon, lat, to_epsg = NULL, verbose = TRUE){
  
  # Create coordinates data.table
  coords <- data.table(
    "x" = lon,
    "y" = lat
  )
  
  #Define dataframe as GIS dataframe with its actual projection
  coord_df <- st_as_sf(
    x = coords ,coords = c('x','y'), crs = 4326)
  
  if(is.null(to_epsg)){
    if(verbose){print("reprojecting coordinates to EPSG: 32736")}
    coord_df <-  st_transform(x = coord_df, crs = 32736)
  }else{
    if(verbose){print(paste0("reprojecting coordinates to EPSG: ", to_epsg))}
    coord_df <- st_transform(x = coord_df, crs = to_epsg)
  }
  
  positions <- as.data.frame(st_coordinates(coord_df))
  
  return(list("x" = positions$X, "y" = positions$Y))
}


#' This method generate a 2d-vector of cartesian positions from the x, y data.
#'
#' @param x list of coordinates x
#' @param y list of coordinates y
#' @param in_latlon If True, x and y coordinates are treated as longitude and latitude 
#' respectively, otherwise they are treated as cartesian coordinates.
#' @param to_epsg If in_latlon is True, x and y are reprojected to this EPSG.
#' @param verbose It specifies if information on the process is printed.
#' 
#' @return  Dataframe with cartesian positions based on x,y coordinates
#' 
#' @export
#'
#' @author Eric Matamoros Morales based on earlier python code by Arnau Pujol.
#'
#' @examples
#' # Creation of x vector of longitude coordinates, y vector of latitude coordinates
#' # and finally merge them on a position data frame.
#' x <- c(1,2,3,4,7.5,8,8.5,9,10,13,13.1,13.2,13.3,14,15,30)
#' y <- c(1,2,3,4,7.5,8,8.5,9,10,13,13.1,13.2,13.3,14,15,30)
#'
#' # Get a data.frame with cartesian positions from x & y data 
#' positions <- get_2dpositions(x, y, TRUE)
#' 
get_2dpositions <- function(
    x, 
    y,
    in_latlon = FALSE, 
    to_epsg = NULL,
    verbose = TRUE){
  
  if(in_latlon){
    coords = latlon2xy(x, y, to_epsg = to_epsg, verbose = verbose)
    x <- coords$x
    y <- coords$y
  }
  
  positions <- data.table("x" = x, "y" = y)
  return(positions)
}