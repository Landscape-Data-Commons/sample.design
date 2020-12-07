#' Add areas to a Spatial Polygons Data Frame
#'
#' @description This function takes a Spatial Polygons Data Frame and calculates and adds area fields to the data frame. Areas can be calculated either treating the whole SPDF as a unit or for each polygon individually.
#' @param polygons Spatial Polygons Data Frame to calculate areas for
#' @param area_ha Logical. If \code{TRUE}, areas will be calculated and added in hectares. Default is \code{TRUE}.
#' @param area_sqkm Logical. If \code{TRUE}, areas will be calculated and added in square kilometers. Default is \code{TRUE}.
#' @param area_acres Logical. If \code{TRUE}, areas will be calculated and added in acres. Default is \code{TRUE}.
#' @param byid Logical. If \code{TRUE}, areas will be calculated and added for each polygon by ID. If \code{FALSE} the area of the whole SPDF will be calculated and added, so every value for that field will be the same, regardless of polygon ID. Default is \code{TRUE}.
#' @return The original Spatial Polygons Data Frame with an additional field for each area unit calculated.
#' @keywords area
#' @examples
#' add_area()
#' @export

add_area <- function(polygons,
                     area_ha = TRUE,
                     area_sqkm = TRUE,
                     area_acres = TRUE,
                     byid = TRUE
){
  ## Make sure we have the polygons in Albers equal area projection
  working_polygons <- sp::spTransform(x = polygons, CRSobj = sp::CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))

  ## Add the area in hectares, using unname() to get an unnamed vector
  polygons@data[["AREA.HA"]] <- unname(rgeos::gArea(working_polygons, byid = byid) * 0.0001)

  if (area_sqkm) {
    # Add the area in square kilometers, converting from hectares
    polygons@data[["AREA.SQKM"]] <- polygons@data[["AREA.HA"]] * 0.01
  }
  if (area_acres) {
    # Add the area in acres, because that's how the BLM rolls
    polygons@data[["AREA.ACRES"]] <- polygons@data[["AREA.HA"]] * 2.471
  }

  # And remove the hectares if they're unwanted
  if (!(area_ha)) {
    polygons@data[["AREA.HA"]] <- NULL
  }

  return(polygons)
}
