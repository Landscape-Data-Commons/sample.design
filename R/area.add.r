#' Add areas to a Spatial Polygons Data Frame
#'
#' @description This function takes a Spatial Polygons Data Frame and calculates and adds area fields to the data frame. Areas can be calculated either treating the whole SPDF as a unit or for each polygon individually.
#' @param spdf Spatial Polygons Data Frame to calculate areas for
#' @param area.ha Logical. If \code{TRUE}, areas will be calculated and added in hectares. Default is \code{TRUE}.
#' @param area.sqkm Logical. If \code{TRUE}, areas will be calculated and added in square kilometers. Default is \code{TRUE}.
#' @param byid Logical. If \code{TRUE}, areas will be calculated and added for each polygon by ID. If \code{FALSE} the area of the whole SPDF will be calculated and added, so every value for that field will be the same, regardless of polygon ID. Default is \code{TRUE}.
#' @return The original Spatial Polygons Data Frame with an additional field for each area unit calculated.
#' @keywords area
#' @examples
#' add.area()
#' @export

add.area <- function(spdf, ## SpatialPolygonsDataFrame to add area values to
                     area.ha = TRUE, ## Add area in hectares?
                     area.sqkm = TRUE, ## Add area in square kilometers?
                     byid = TRUE ## Do it for the whole SPDF or on a per-polygon basis? Generally don't want to toggle this
){
  original.proj <- spdf@proj4string
  ## Make sure the SPDF is in Albers equal area projection
  spdf <- sp::spTransform(x = spdf, CRSobj = sp::CRS("+proj=aea"))

  ## Add the area in hectares, stripping the IDs from gArea() output
  spdf@data$AREA.HA <- unname(rgeos::gArea(spdf, byid = byid) * 0.0001)
  ## Add the area in square kilometers, converting from hectares
  spdf@data$AREA.SQKM <- spdf@data$AREA.HA * 0.01

  if (!(area.ha)) {
    spdf@data$AREA.HA <- NULL
  }
  if (!(area.sqkm)) {
    spdf@data$AREA.SQKM <- NULL
  }
  return(sp::spTransform(spdf, original.proj))
}
