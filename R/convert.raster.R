#' Convert a raster to Spatial Polygons restriced to a frame
#' @param raster Raster. The raster to be clipped
#' @param frame_spdf Spatial polygons or spatial polygons data frame. The frame to clip the raster to.
#' @param buffer.meters Numeric. Number of meters to buffer the frame by before masking the raster. Make sure this is slightly larger than the raster cell size to make sure that no cells are lost before vectorization and clipping to the actual frame.
#' @export
raster_spdf <- function(raster,
                        frame_spdf,
                        buffer.meters = 100){
  # Reproject and buffer the frame
  frame_spdf_buffered <- raster::buffer(x = sp::spTransform(frame_spdf,
                                                            CRSobj = raster::crs(raster)),
                                        width = buffer_meters)
  # Mask the raster with the buffered frame
  raster_crop <- raster::crop(x = raster,
                              y = frame_spdf_buffered)

  # Convert the raster to vector
  raster_poly <- raster::rasterToPolygons(x = raster_crop,
                                          dissolve = TRUE)

  # Clip the vector to the original, unbuffered frame
  output <- raster::crop(x = raster_poly,
                         y = sp::spTransform(frame_spdf,
                                             CRSobj = raster::crs(raster)))

  return(output)
}
