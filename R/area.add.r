#' Add areas to a Spatial Polygons Data Frame
#'
#' @description This function takes a Spatial Polygons Data Frame and calculates and adds area fields to the data frame. Areas can be calculated either treating the whole SPDF as a unit or for each polygon individually.
#' @param polygons An sf polygon object to calculate areas for
#' @param area_ha Logical. If \code{TRUE}, areas will be calculated and added in hectares. Default is \code{TRUE}.
#' @param area_sqkm Logical. If \code{TRUE}, areas will be calculated and added in square kilometers. Default is \code{TRUE}.
#' @param area_acres Logical. If \code{TRUE}, areas will be calculated and added in acres. Default is \code{TRUE}.
#' @return The original sf object with an additional field for each area unit calculated.
#' @keywords area
#' @examples
#' add_area()
#' @export

add_area <- function(polygons,
                     area_ha = TRUE,
                     area_sqkm = TRUE,
                     area_acres = TRUE
){
  if (!("sf" %in% class(polygons))) {
    stop("polygons must be an sf object.")
  }
  if (!all(sf::st_geometry_type(polygons) %in% c("POLYGON", "MULTIPOLYGON"))) {
    stop("polygons must be an sf object.")
  }
  
  ## Make sure we have the polygons in Albers equal area projection
  working_polygons <- sf::st_transform(x = polygons,
                                       crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
  
  ## Add the area in hectares, using unname() to get an unnamed vector
  # sf::st_area() will return in square meters here, so multiply by the conversion
  # factor 0.0001 to get hectares
  polygons[["AREA.HA"]] <- as.vector(sf::st_area(x = working_polygons,
                                                 by_element = TRUE) * 0.0001)
  
  if (area_sqkm) {
    # Add the area in square kilometers, converting from hectares
    polygons[["AREA.SQKM"]] <- polygons[["AREA.HA"]] * 0.01
  }
  if (area_acres) {
    # Add the area in acres, because that's how the BLM rolls
    polygons[["AREA.ACRES"]] <- polygons[["AREA.HA"]] * 2.471
  }
  
  # And remove the hectares if they're unwanted
  if (!(area_ha)) {
    polygons[["AREA.HA"]] <- NULL
  }
  
  return(polygons)
}
