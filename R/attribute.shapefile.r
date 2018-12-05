#' Attributing a spatial data frame using another spatial data frame
#'
#' This function will take a SpatialPoints/PolygonsDataFrame and add one attribute fields from a second SpatialPoints/PolygonsDataFrame
#' @param spdf1 A SpatialPoints/PolygonsDataFrame containing the geometry to add an attribute to
#' @param spdf2 A SpatialPoints/PolygonsDataFrame containing the geometry to add an attribute from
#' @param attributefield The name of the field in \code{spdf2} as a string containing the values to add to \code{spdf1}
#' @param newfield The name of the field in \code{spdf1} as a string to add the values from \code{spdf2$attributefield} to. If NULL, the field will use \code{attributefield}. Defaults to NULL.
#' @return The original SPDF spdf1 with the new field containing the values inherited from spdf2.
#' @examples
#' attribute.shapefile()
#' @export

attribute.shapefile <- function(spdf1,
                                spdf2,
                                attributefield = NULL,
                                newfield = NULL
){
  if (is.null(attributefield) | !(attributefield %in% names(spdf2@data))) {
    stop("attributefield must be a field name found in spdf2")
  }

  if (is.null(newfield)) {
    newfield <- attributefield
  }

  # This is so that, if something gets tripped later, we can flip to TRUE
  remove_coords <- FALSE
  coord_names <- colnames(spdf1@coords)

  if (spdf1@proj4string@projargs != spdf2@proj4string@projargs) {
    ## Make sure that the points also adhere to the same projection
    spdf2 <- sp::spTransform(spdf2, CRSobj = spdf1@proj4string)
  }
  projection <- spdf1@proj4string

  ## Initialize list for attributed SPDFs
  attributed_dfs <- list()

  ## We'll check each attribute field value independently
  for (n in unique(spdf2@data[, attributefield])) {
    ## Create a copy of the points to work with on this loop
    current_spdf <- spdf1
    ## Get the data frame from over()
    over_result <- sp::over(current_spdf,
                            spdf2[spdf2@data[, attributefield] == n, ])
    ## Add the values to the newfield column
    current_spdf@data[, newfield] <- over_result[, attributefield]
    if (!(coord_names[1] %in% names(current_spdf@data)) & !(coord_names[2] %in% names(current_spdf@data))){
      current_spdf@data <- cbind(current_spdf@data, current_spdf@coords)
      remove_coords <- TRUE
    }
    ## Make sure that the polygons have unique IDs
    if (class(current_spdf) == "SpatialPolygonsDataFrame") {
      current_spdf <- sp::spChFIDs(current_spdf,
                                   paste(runif(n = 1, min = 0, max = 666666666),
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
      output <- attributed_dfs[[1]] %>%
        sp::SpatialPointsDataFrame(data = .,
                                   coords = .[, coord_names],
                                   proj4string = projection)
    } else {
      output <- dplyr::bind_rows(attributed_dfs) %>%
        sp::SpatialPointsDataFrame(data = .,
                                   coords = .[, coord_names],
                                   proj4string = projection)
    }
    if (remove_coords) {
      output <- output[, names(output)[!(names(output) %in% coord_names)]]
    }
  } else {
    output <- NULL
  }
  return(output)
}
