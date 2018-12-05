#' Add areas to a Spatial Polygons Data Frame
#'
#' @description This function takes a Spatial Polygons Data Frame and calculates and adds area fields to the data frame. Areas can be calculated either treating the whole SPDF as a unit or for each polygon individually.
#' @param spdf Spatial Polygons Data Frame to calculate areas for
#' @param area_ha Logical. If \code{TRUE}, areas will be calculated and added in hectares. Default is \code{TRUE}.
#' @param area_sqkm Logical. If \code{TRUE}, areas will be calculated and added in square kilometers. Default is \code{TRUE}.
#' @param area_acres Logical. If \code{TRUE}, areas will be calculated and added in acres. Default is \code{TRUE}.
#' @param byid Logical. If \code{TRUE}, areas will be calculated and added for each polygon by ID. If \code{FALSE} the area of the whole SPDF will be calculated and added, so every value for that field will be the same, regardless of polygon ID. Default is \code{TRUE}.
#' @return The original Spatial Polygons Data Frame with an additional field for each area unit calculated.
#' @keywords area
#' @examples
#' add.area()
#' @export

add.area <- function(spdf,
                     area_ha = TRUE,
                     area_sqkm = TRUE,
                     area_acres = TRUE,
                     byid = TRUE
){
  # Get whatever the original projection was so we can put it back after
  original_proj <- spdf@proj4string
  ## Make sure the SPDF is in Albers equal area projection
  spdf <- sp::spTransform(x = spdf, CRSobj = sp::CRS("+proj=aea"))

  ## Add the area in hectares, using unname() to get an unnamed vector
  spdf@data$AREA.HA <- unname(rgeos::gArea(spdf, byid = byid) * 0.0001)



  if (area_sqkm) {
    # Add the area in square kilometers, converting from hectares
    spdf@data$AREA.SQKM <- spdf@data$AREA.HA * 0.01
  }
  if (area_acres) {
    # Add the area in acres, because that's how the BLM rolls
    spdf@data$AREA.ACRES <- spdf@data$AREA.HA * 2.471
  }

  # And remove the hectares if they're unwanted
  if (!(area_ha)) {
    spdf@data$AREA.HA <- NULL
  }

  return(sp::spTransform(spdf, original_proj))
}
