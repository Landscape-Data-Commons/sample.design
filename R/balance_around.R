#' Select points that most closely approximate the distribution of another set of points
#' @description Sometimes you have a large collection of points which are not randomly distributed or spatially balanced and you would like a subset that more or less do. Given a template of points that are distributed the way you would like, this will return the closest existing point to each. This can be done taking into account membership in a group, either by having assigned it as a variable in both sets of points or by providing polygons that can be used to assign membership. By default, no stratification/membership is taken into account.
#' @param existing_points Spatial Points Data Frame. The points you would like to select from by comparing against \code{template_points}.
#' @param template_points  Spatial Points Data Frame. The points you would like to compare against \code{existing_points} in order to select a subset of those that most closely resemble the distribution of the template points.
#' @param strata_spdf Spatial Polygons Data Frame. Polygons assigned a variable with a name \code{stratafield} that contains the membership information (e.g. strata) to assign to \code{existing_points} and \code{template_points}. If \code{NULL} then no assignment will be attempted. Defaults to \code{NULL}.
#' @param stratafield Character string. If \code{strata_spdf} is not \code{NULL}, the name of the variable in \code{strata_spdf@@data} that contains the membership information. Otherwise, the name of the variable in both \code{template_points@@data} and \code{existing_points@@data} that contains the membership information. If \code{NULL} then the points will be considered to belong to a single group. Defaults to \code{NULL}.
#' @param projection CRS object. The projection to force all spatial objects into. Defaults to NAD83, \code{sp::CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")}.
#' @return A spatial polygons data frame made by trimming \code{existing_points} down to the points that most closely approximate the distribution of \code{template_points} while also containing the same number of points as \code{template_points}. It will be in the projection specified by \code{projection}.
#' @export
get_closest <- function(existing_points,
                        template_points,
                        strata_spdf = NULL,
                        stratafield = NULL,
                        projection = sp::CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")){

  # Reproject if necessary
  if (!(class(existing_points) %in% "SpatialPointsDataFrame")) {
    stop("existing_points must be a spatial points data frame")
  }
  if (!identical(projection, existing_points@proj4string)) {
    existing_points <- sp::spTransform(existing_points,
                                       projection)
  }
  if (!(class(template_points) %in% "SpatialPointsDataFrame")) {
    stop("template_points must be a spatial points data frame")
  }
  if (!identical(projection, template_points@proj4string)) {
    template_points <- sp::spTransform(template_points,
                                       projection)
  }

  if (nrow(template_points@data) < nrow(existing_points@data)) {
    stop("There are more points in existing_points than template_points")
  }

  # What to do about stratafield
  if (!is.null(stratafield)) {
    if (!is.character(stratafield)) {
      stop("stratafield must be a character string")
    }
    if (length(stratafield) > 1) {
      stop("stratafield must be a character string")
    }
    if (is.null(strata_spdf)) {
      if (!(stratafield %in% names(template_points@data))) {
        stop("The variable ", stratafield, " does not appear in template_points@data.")
      }
      if (!(stratafield %in% names(existing_points@data))) {
        stop("The variable ", stratafield, " does not appear in existing_points@data.")
      }
      # Sorry for the inconsistencies but "MEMBERSHIP" is trying to get away from the assumption that these will be stratifications
      # And I'm so, so tired. I'll come back to clean this up later, or so I'm currently telling myself at the end of a ten hour day
      points_spdf@data[["MEMBERSHIP"]] <- points_spdf@data[[stratafield]]
    }

  } else {
    points_spdf@data[["MEMBERSHIP"]] <- "frame"
  }

  # If there is a set of stratification polygons, use them!
  if (!is.null(strata_spdf)) {
    if (!(class(strata_spdf) %in% "SpatialPolygonsDataFrame")) {
      stop("strata_spdf must be a spatial polygons data frame.")
    }
    if (!identical(projection, strata_spdf@proj4string)) {
      strata_spdf <- sp::spTransform(strata_spdf,
                                     projection)
    }
    if (!(stratafield %in% names(strata_spdf@data))) {
      stop("The variable ", stratafield, " does not appear in strata_spdf.")
    }

    # This just puts the strata into the points
    existing_points@data[["MEMBERSHIP"]] <- sp::over(x = existing_points,
                                                     y = strata_spdf)[[stratafield]]
    template_points@data[["MEMBERSHIP"]] <- sp::over(x = template_points,
                                                     y = strata_spdf)[[stratafield]]
  }


  strata <- unique(points_spdf@data[["MEMBERSHIP"]])

  # By stratum!
  stata_selections <- lapply(X = strata,
                             template_points = template_points,
                             existing_points = existing_points,
                             FUN = function(X, template_points, existing_points){
                               stratum <- X
                               existing_points_stratum <- existing_points[existing_points@data[["MEMBERSHIP"]] == stratum, ]
                               template_points_stratum <- template_points[template_points@data[["MEMBERSHIP"]] == stratum, ]

                               # No points have been picked yet!
                               existing_points_stratum@data[["PICKED"]] <- FALSE

                               # We'll go through each new point in turn
                               # This means that each point won't even look at any existing points already flagged as "PICKED"
                               for(template_point in 1:nrow(template_points_stratum)) {
                                 # The x and y components for the current template point
                                 template_point_x <- template_points_stratum@data[template_point, "XMETERS"]
                                 template_point_y <- template_points_stratum@data[template_point, "YMETERS"]

                                 # This is the index of the existing point being used
                                 # We don't use this until every point has been tested
                                 use_index <- 0

                                 # Which points haven't been flagged as used yet?
                                 unused_indices <- (1:nrow(existing_points_stratum))[existing_points_stratum@data[["PICKED"]]]

                                 # These are going to be the squares of the distances from the current new point to every unused existing point
                                 # Basically, this is the Pythagorean theorem without the sqaure root step because why do unnecessary math?
                                 distances <- sapply(X = unused_indices,
                                                     comparison_x = template_point_x,
                                                     comparison_y = template_point_y,
                                                     points_spdf = existing_points_stratum,
                                                     FUN = function(X, comparison_x, comparison_y, points_spdf){
                                                       x_component <- (comparison_x - points_spdf@data[X, "XMETERS"])^2
                                                       y_component <- (comparison_y - points_spdf@data[X, "YMETERS"])^2
                                                       return(x_component + y_component)
                                                     })

                                 index_closest <- unused_indices[which(min(distances), distances)]

                                 existing_points_stratum$USED[index_closest] <- TRUE
                               }

                               # Only keep the points that are flagged as picked
                               output <- existing_points_stratum[existing_points_stratum@data[["PICKED"]], ]
                               output@data[["PICKED"]] <- NULL
                               return(output)
                             })

  # Each stratum's points are a separate SPDF in the list, so let's mash 'em together
  output <- do.call(rbind,
                    strata_selections)

  # If we just had to add an identity, hide that fact
  if (is.null(stratafield)) {
    output@data[["MEMBERSHIP"]] <- NULL
  } else {
    # Otherwise name that membership field according to what the user called it
    names(output@data)[names(output@data) == "MEMBERSHIP"] <- stratafield
  }

  return(output)
}


#' Find and keep the points farthest from each other
#' @description This will take a set of existing points and new points and combine them to create a set consisting of the existing points and the farthest new points. The original intended use case was to take a collection of sampling locations from one or more sample designs (\code{existing_points}) and use them as part of a new, spatially balanced sample design. The function takes a set of new, random, spatially balanced points (\code{new_points}) and determines the distance between each of them and each of the existing points. It then sequentially eliminates the new point closest to any existing point until the combined number of existing points and remaining points is equal to \code{target}.
#' @param existing_points Spatial points data frame. These are the points that will all be included in the output points and the points against which \code{new_points} will be compared against.
#' @param new_points Spatial points data frame. These are the points that may be included in the output. The number that will be is equal to \code{target - nrow(existing_points)}.
#' @param target Numeric value. The total number of points to include in the output. Defaults to \code{nrow(new_points)}.
#' @param projection CRS object. The projection to force all spatial objects into. Defaults to NAD83, \code{sp::CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")}.
#' @return A spatial points data frame containing all the points from existing points and \code{target - nrow(existing_points)} points from \code{new_points} using the CRS specified in \code{projection}. It will only have those variables that were in common between both \code{existing_points} and \code{new_points}.
#' @export
keep_farthest <- function(existing_points,
                          new_points,
                          target = NULL,
                          projection = sp::CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")){
  # TODO: Sanitize
  if (!(class(existing_points) %in% "SpatialPointsDataFrame")) {
    stop("existing_points must be a spatial points data frame")
  }
  if (!identical(projection, existing_points@proj4string)) {
    existing_points <- sp::spTransform(existing_points,
                                       projection)
  }
  if (nrow(existing_points@data) < 1) {
    stop("There are no points in existing_points")
  }
  if (!(class(new_points) %in% "SpatialPointsDataFrame")) {
    stop("new_points must be a spatial points data frame")
  }
  if (!identical(projection, new_points@proj4string)) {
    new_points <- sp::spTransform(new_points,
                                  projection)
  }
  if (nrow(new_points@data) < 1) {
    stop("There are no points in new_points")
  }

  if (is.null(target)) {
    target <- nrow(new_points@data)
  }

  if (target <= nrow(existing_points@data)) {
    stop("The target number of points is less than or equal to the number of existing points.")
  }

  common_varnames <- unique(c(names(existing_points@data)[!(names(existing_points@data) %in% names(new_points@data))],
                              names(new_points@data)[!(names(new_points@data) %in% names(existing_points@data))]))

  if (length(common_varnames) < 1) {
    stop("There are no variables in common between existing_points and new_points. There must be at least one")
  }

  if (length(names(existing_points@data)) != length(common_varnames) | length(names(new_points@data)) != length(common_varnames)) {
    message("Not all variables are in common between existing_points and new_points")
  }

  existing_points <- existing_points[, common_varnames]
  new_points <- new_points[, common_varnames]

  # How many of each point type are there? We'll use these for the loops
  n_existing <- nrow(existing_points@data)
  n_new <- nrow(new_points)

  # Get some common info added to these
  existing_points@data[["TYPE"]] <- "EXISTING"
  existing_points@data[["INDEX"]] <- 1:n_existing
  new_points@data[["TYPE"]] <- "NEW"
  new_points@data[["INDEX"]] <- 1:n_new

  # Combine the two sets of points, making sure the existing points come first!!!
  combined_points <- rbind(existing_points[, c("INDEX", "TYPE")],
                           new_points[, c("INDEX", "TYPE")])

  # Add the coordinates
  combined_points <- get_coords(combined_points,
                                x_var = "XMETERS",
                                y_var = "YMETERS",
                                projection = sp::CRS("+proj=aea"))

  # Get a distance matrix
  distance_matrix <- dist_matrix(dataframe = combined_points@data,
                                 x_var = "XMETERS",
                                 y_var = "YMETERS")

  # Remove the columns for new points
  distance_matrix <- distance_matrix[, -((n_existing + 1):(n_existing + n_new))]
  # Remove the rows for existing points
  distance_matrix <- distance_matrix[-(1:n_existing), ]

  # Convert to a data frame
  distance_df <- as.data.frame(distance_matrix,
                               stringsAsFactors = FALSE)

  # Add the new points indices
  distance_df[["INDEX"]] <- 1:n_new

  # And now we loop to remove all the closest new points
  # The looping is so that we can remove points that are already tossed
  # Each loop will:
  # 1) Check for the new point(s) closest to an existing point
  # 2) Store the index (or indices) from distance_df[["INDEX"]]
  # 3) Remove that observation(s) from distance_df
  # This will continue until the number of stored indices is equal to the number of existing points
  # If the number to remove is overshot because the finale pass through the loop identifies multiple indices, we only take the first however we need
  # The observations at the identified indices in the new points will be removed!


  # A vector to store the indices to chuck
  removal_indices <- NULL
  # Here's one we can mutilate over our iterations (always avoid violence to original data you may reference again!)
  working_distance_df <- distance_df

  while (n_removal_indices < (target - n_existing)) {
    current_min <- min(working_distance_df)
    # Get the indices from every column where that min occurs
    # Each column is an existing point, so if we check every column for the value and store that index,
    # those are the new points that are that distance from an existing point
    current_indices <- sapply(X = 1:ncol(working_distance_df),
                              value = current_min,
                              distance_df = working_distance_df,
                              FUN = function(X, value, distance_df){
                                # Grab the column as a vector
                                distances <- distance_df[[X]]
                                if (value %in% distances) {
                                  # Get the values from INDEX that correspond to the places where the value is found
                                  # We're returning from INDEX because they won't change when we remove rows from working_distance_df
                                  indices <- distance_df[which(value == distances), "INDEX"]
                                  return(indices)
                                } else {
                                  return(NULL)
                                }
                              })
    # Make sure we have the uniques, in case a new point was equidistant from multiple existing points
    current_indices <- unique(current_indices)

    # Add those to the ones we're going to remove
    removal_indices <- unique(c(removal_indices, current_indices))

    if (!is.null(removal_indices)) {
      # Make sure we don't overshoot our removal goal
      removal_indices <- removal_indices[1:min(length(removal_indices, (target - n_existing)))]

      # Remove them!
      working_distance_df <- working_distance_df[!(working_distance_df[["INDEX"]] %in% removal_indices), ]

      # Update our tracking value
      n_removal_indices <- length(removal_indices)
    }
  }

  # Now that we have our indices to remove, let's do it as we combine points
  output <- rbind(existing_points[, common_varnames],
                  new_points[-removal_indices, common_varnames])

  return(output)
}


##  Ingest inputs, call functions to access XY coords and to identify and eliminate New points, and output
##          expanded, balanced design.
#' @param existing_points_spdf Spatial ponits data frame. The existing points that will be balanced around.
#' @param new_points_spdf Spatial points data frame. The points that will be compared against the existing points and selected from to create a balanced design.
#' @param stratafield Character string. The name of the variable in common between \code{existing_points_spdf} and \code{new_points_spdf} that contains stratum identities. This is used to balance by stratum. If \code{NULL} then balancing will not take strata into account. Defaults to \code{NULL}.
#' @param projection CRS object. The projection to force on the spatial objects. Defaults to \code{sp::CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")}.
#' @return A spatial points data frame containing all the points from \code{existing_points_spdf} and the selected points from \code{new_points_spdf}. The projection will match \code{projection}.
#' @export
balance_around <- function(existing_points,
                           new_points,
                           stratafield = NULL,
                           projection = sp::CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")){
  # TODO: Sanitization (including reprojection)

  # Assign the codes that indicate if they're existing plots or freshly-drawn ones
  existing_pts_spdf@data[["TYPE"]] <- "EXISTING"
  new_points@data[["TYPE"]] <- "NEW"

  # Add in the fields that the points don't have for easy combination
  missing_vars_new <- names(existing_points@data)[!(names(existing_points@data) %in% names(new_points@data))]
  new_points@data[, missing_vars_new] <- NA
  missing_vars_existing <- names(new_points@data)[!(names(new_points@data) %in% names(existing_points@data))]
  existing_points@data[, missing_vars_existing] <- NA

  # Bind the 2 point files together
  pts <- rbind(existing_points, new_points)

  # What are the existing points' indices?
  extant_indices <- 1:nrow(existing_points@data)
  # Which are the new points' indices?
  new_indices <- (nrow(existing_points@data) + 1):(nrow(existing_points@data) + nrow(new_points@data))

  # Rename the date that the plot was sampled to "PREVDATE"
  pts@data[["PREVDATE"]] <- pts@data[["DATEVISITE"]]

  # Restrict to only relevant fields
  pts <- pts[ , c("TYPE",
                  stratafield,
                  "PLOTID",
                  "PLOTKEY",
                  "PRIMARYKEY",
                  "PANEL",
                  "PROJECTNAM",
                  "PREVDATE")]

  # Make these all NA for the new points
  new_points@data[new_indices, "PLOTID"] <- NA
  new_points@data[new_indices, "PLOTKEY"] <- NA
  new_points@data[new_indices, "PRIMARYKEY"] <- NA
  new_points@data[new_indices, "PROJECTNAM"] <- NA
  new_points@data[new_indices, "PREVDATE"] <- NA

  # Determine the number of existing points
  extant <- nrow(existing_points@data)

  # Add coordinates to the combined points for distance calculations
  pts <- get_coords(pts,
                    x_var = "XMETERS",
                    y_var = "YMETERS",
                    projection = sp::CRS("+proj=aea"))

  # This determines the number of New points to eliminate, and eliminates the points
  pts <- keep_farthest(existing_points = existing_points,
                       new_points = new_points)

  # Time to tweak the new points that were kept
  new_indices_remaining <- pts@data[["TYPE"]] == "NEW"

  # Get the existing plot ids and renumber them
  if (any(new_indices_remaining)) {
    plotids <- pts@data[new_indices_remaining, "PLOTID"]
    plotids <- gsub(plotids,
                    pattern = "\d*$",
                    replacement = "")
    pts@data[new_indices_remaining, "PLOTID"] <- paste0(plotids, 1:length(new_indices_remaining))
  }

  # TODO: Rename within strata.
  # Does this mean using over() with stratification polygons to determine new stratification assignments for old points?


  # None of the points are considered visited now!
  pts@data[["DATEVIS"]] <- ""
  pts@data[["EVALSTA"]] <- "NotEval"
  pts@data[["FINAL_DESI"]] <- ""

  # Add coordinates!
  pts <- get_coords(pts,
                    x_var = "LAT",
                    y_var = "LONG",
                    projection = projection)

  return(pts)
}
