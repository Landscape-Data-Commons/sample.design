#' Attributing a spatial data frame using another spatial data frame
#'
#' This function will take a SpatialPoints/PolygonsDataFrame and add one attribute fields from a second SpatialPoints/PolygonsDataFrame
#' @param spdf1 Spatial points or polygon data frame. The geometry that the attribute will be added to.
#' @param source_geometry Spatial polygons data frame. The polygons with the attribute to be added to \code{spdf1@@data}.
#' @param attributefield Character string. The name of the field in \code{source_geometry@@data} containing the values to add to \code{spdf1@@data}
#' @param newfield Optional character string. The name of the field in \code{spdf1@@data} to add the values from \code{source_geometry@@data$attributefield} to. If NULL, the field will use \code{attributefield}. Defaults to NULL.
#' @return The original SPDF spdf1 with the new field containing the values inherited from source_geometry.
#' @examples
#' attribute_shapefile()
#' @export

attribute_shapefile <- function(spdf1,
                                source_geometry,
                                attributefield = NULL,
                                newfield = NULL
){
  if (is.null(attributefield) | !(attributefield %in% names(source_geometry@data))) {
    stop("attributefield must be a field name found in source_geometry")
  }

  if (is.null(newfield)) {
    newfield <- attributefield
  }

  # This is so that, if something gets tripped later, we can flip to TRUE
  remove_coords <- FALSE
  coord_names <- colnames(spdf1@coords)

  if (!identical(spdf1@proj4string, source_geometry@proj4string)) {
    ## Make sure that the points also adhere to the same projection
    source_geometry <- sp::spTransform(source_geometry,
                                       CRSobj = spdf1@proj4string)
  }
  projection <- spdf1@proj4string

  ## Initialize list for attributed SPDFs
  attributed_dfs <- list()

  ## We'll check each attribute field value independently
  attributes <- unique(source_geometry@data[, attributefield])

  test <- lapply(attributes,
                 spdf1 = spdf1,
                 source_geometry = source_geometry,
                 FUN = function(X, spdf1, source_geometry){

                 })
  for (n in unique(source_geometry@data[, attributefield])) {
    ## Create a copy of the points to work with on this loop
    current_spdf <- spdf1
    ## Get the data frame from over()
    over_result <- sp::over(current_spdf,
                            source_geometry[source_geometry@data[, attributefield] == n, ])
    ## Add the values to the newfield column
    current_spdf@data[, newfield] <- over_result[, attributefield]
    if (!(coord_names[1] %in% names(current_spdf@data)) & !(coord_names[2] %in% names(current_spdf@data))){
      current_spdf@data <- cbind(current_spdf@data, current_spdf@coords)
      remove_coords <- TRUE
    }
    ## Make sure that the polygons have unique IDs
    if (class(current_spdf) == "SpatialPolygonsDataFrame") {
      current_spdf <- sp::spChFIDs(current_spdf,
                                   paste(stats::runif(n = 1, min = 0, max = 666666666),
                                         row.names(current_spdf),
                                         sep = "."))
    }
    current_df <- current_spdf@data

    ## Only if the number of coordinates is greater than 0!
    print(nrow(current_df[!is.na(current_df[, newfield]), ]))
    if (nrow(current_df[!is.na(current_df[, newfield]), ]) > 0) {
      attributed_dfs[[paste(n)]] <- current_df[!is.na(current_df[, newfield]), ]
    }
  }

  if (length(attributed_dfs) > 0) {
    if (length(attributed_dfs) == 1) {
      attributed_df <- attributed_dfs[[1]]
    } else {
      attributed_df <- do.call((attributed_dfs), rbind)
    }
    output <- sp::SpatialPointsDataFrame(data = attributed_df,
                                         coords = attributed_df[, coord_names],
                                         proj4string = projection)

    if (remove_coords) {
      output <- output[, names(output)[!(names(output) %in% coord_names)]]
    }
  } else {
    output <- NULL
  }
  return(output)
}
