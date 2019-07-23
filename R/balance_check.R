#' Dissolve polygons in an SPDF by an identity variable
#' @description Dissolve polygons together on an identity variable to produce an SPDF with one observation per identity.
#' @param spdf Spatial Polygons Data Frame. The polygons to dissolve. Must have an identity variable with a name patching \code{dissolve_field}.
#' @param dissolve_field Character string. The name of the variable in \code{spdf@@data} containing the polygon identities
#' @return A spatial polygons data frame with one observation per unique identity in the original \code{spdf}. The associated data frame contains only the identities in a variable named using \code{dissolve_field}.
#' @export
dissolve_spdf <- function(spdf, dissolve_field){
  if (!grepl(class(spdf), pattern = "^SpatialPolygonsDataFrame")) {
    stop("spdf must be a spatial polygons data frame")
  }
  if (!is.character(dissolve_field)) {
    stop("dissolve field must be a character string")
  }
  if (length(dissolve_field) > 1) {
    stop("dissolve field must be a character string")
  }
  if (!(dissolve_field %in% names(spdf@data))) {
    stop(dissolve_field, " does not occur in the names of spdf@data")
  }
  unique_ids <- as.character(unique(spdf@data[[dissolve_field]]))
  poly_list <- lapply(X = unique_ids,
                      spdf = spdf,
                      dissolve_field = dissolve_field,
                      FUN = function(X, spdf, dissolve_field){
                        spdf_current <- spdf[spdf@data[[dissolve_field]] == X, ]
                        spdf_current <- methods::as(sf::st_combine(sf::st_as_sf(spdf_current)), "Spatial")
                        df <- data.frame(id = X,
                                         stringsAsFactors = FALSE)
                        names(df) <- dissolve_field
                        rownames(df) <- spdf_current@polygons[[1]]@ID
                        spdf_current <- sp::SpatialPolygonsDataFrame(Sr = spdf_current,
                                                                     data = df)
                        return(spdf_current)
                      })
  output <- do.call(rbind,
                    poly_list)
  return(output)
}


#' Get the areas of polygons
#' @description Given a Spatial Polygons Data Frame, produce a data frame summarizing the areas of polygons sharing an identity
#' @param polygons Spatial Polygons or Spatial Polygons Data Frame. The polygons to be summarized.
#' @param cum_freq Logical. If \code{TRUE} then the output data frame will include the variable \code{"cum_freq"} with the cumulative frequencies of the polygon identities. Note that it sorts the polygons from smallest to largest area before doing this. Defaults to \code{TRUE}.
#' @return A data frame containing the polygon index within \code{polygons@@polygons[[1]]@@Polygons}, and if asked for, the proportional area (\code{"area_prop"}) and the cumulative frequency (\code{"cum_freq"}).
#' @export
extract_poly_area <- function(polygons,
                              cum_freq = TRUE) {
  if (!grepl(class(polygons), pattern = "^SpatialPolygons")) {
    stop("spdf must be a spatial polygons data frame")
  }

  # We need to handle what to do if the geometry is empty
  if (length(polygons@polygons) < 1) {
    stop("There's no geometry in polygons")
  }

  # And if it's not dissolved, we'll do that!
  if (length(polygons@polygons) > 1) {
    message("The polygons need to be dissolved. Dissolving now.")
    polygons <- methods::as(sf::st_combine(sf::st_as_sf(polygons)), "Spatial")
  }

  # Get the areas of the polygons
  areas_df <- data.frame(id = 1:length(polygons@polygons[[1]]@Polygons),
                         area = sapply(X = polygons@polygons[[1]]@Polygons,
                                       FUN = function(X){
                                         X@area
                                       }),
                         stringsAsFactors = FALSE)


  # Sort from largest to smallest area
  # (This can help speed up selecting from the probability distribution elsewhere)
  areas_df <- areas_df[order(-areas_df[["area"]]), ]

  # Get the total area
  total_area <- sum(areas_df[["area"]])

  # Add proportional area to the data frame
  # if (area_prop) {
  areas_df[["area_prop"]] <- areas_df[["area"]] / total_area
  # }

  # Add cumulative frequency distribution, which can be treated as a probability distribution
  if (cum_freq) {
    areas_df[["cum_freq"]] <- cumsum(areas_df[["area_prop"]])
  }

  # Return this data frame!
  return(areas_df)
}


#' Select polygon from a probability distribution
#' @description Given a data frame of IDs with associated probabilities and a number, select the ID with the smallest probability value greater than than the number given.
#' @param dataframe A data frame. Must contain an identity variable with a name matching \code{id_var} and a probability variable with a name matching \code{prob_var}. Will be sorted by ascending probability.
#' @param value Numeric. A single numeric value to compare against the probabilities in \code{dataframe}.
#' @param prob_var Character string. Must match the name of the variable in \code{dataframe} that contains the probability values. Defaults to \code{"cum_freq"}.
#' @param id_var Character string. Must match the name of the variable in \code{dataframe} that contains the identities. Defaults to \code{"id"}.
#' @return The identity value from \code{dataframe} with the smallest probabilty value greater than \code{value} OR the first identity value if none were greater than \code{value}.
#' @export
select_from_distribution <- function(dataframe,
                                     value,
                                     prob_var = "cum_freq",
                                     id_var = "id") {
  if (class(dataframe) != "data.frame") {
    stop("dataframe must be a data frame")
  }
  if (!(prob_var %in% names(dataframe))) {
    stop("prob_var must correspond to the name of a variable in dataframe")
  }
  if (!(id_var %in% names(dataframe))) {
    stop("id_var must correspond to the name of a variable in dataframe")
  }
  if (!is.numeric(dataframe[[prob_var]])) {
    stop("prob_var must correspond to the name of a numeric variable in dataframe")
  }
  if (!is.numeric(value)) {
    stop("value must be numeric")
  }
  if (length(value) > 1) {
    stop("value must be a single numeric value")
  }

  # Make sure that they're ordered!
  dataframe <- dataframe[order(dataframe[[prob_var]]), ]

  # Check to see if any of the values in prob_var are greater than value
  # If so, return the ID of the first, else return 1
  check <- value < dataframe[[prob_var]]

  if (any(check)) {
    return(dataframe[check, id_var][1])
  } else {
    return(1)
  }
}


#' Calculate the geometric mean of a numeric vector
#' @param x Numeric vector. The values to calculate a geometric mean from
#' @param na.rm Logical. If \code{TRUE} then \code{NA} values will be dropped from \code{x} before calculating the geometric mean. Defaults to \code{TRUE}
#' @return The geometric mean as a single numeric value.
#' @esport

gm_mean <- function(x,
                    na.rm = TRUE) {
  if (!is.numeric(x)) {
    stop("x must be numeric")
  }

  # Here's the math for a geometric mean, I guess
  output <- exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))

  return(output)
}


#' Calculate a distance matrix for a data frame of coordinates
#' @description Given a data frame containing X and Y Cartesian coordinates, calculate a distance matrix between the points.
#' @param dataframe A data frame. Must contain numeric variables corresponding to \code{x_var} and \code{y_var} containing the X and Y coordinate values.
#' @param x_var Character string. The name of the variable in \code{dataframe} containing the X components of the coordinates.
#' @param y_var Character string. The name of the variable in \code{dataframe} containing the Y components of the coordinates.
#' @return A matrix of the distances between the points.
#' @export
dist_matrix <- function(dataframe,
                        x_var,
                        y_var){
  if (class(dataframe) != "data.frame") {
    stop("dataframe must be a data frame")
  }
  if (!(x_var %in% names(dataframe))) {
    stop("x_var must correspond to the name of a variable in dataframe")
  }
  if (!(y_var %in% names(dataframe))) {
    stop("y_var must correspond to the name of a variable in dataframe")
  }
  if (!is.numeric(dataframe[[x_var]])) {
    stop("x_var must correspond to the name of a numeric variable in dataframe")
  }
  if (!is.numeric(dataframe[[y_var]])) {
    stop("y_var must correspond to the name of a numeric variable in dataframe")
  }

  # Calculate Euclidean distance between all points!
  a <- dataframe[[x_var]]
  b <- dataframe[[y_var]]
  a <- outer(a, a, '-')
  a_sq <- a * a
  b <- outer(b, b, '-')
  b_sq <- b * b

  # h is for hypotenuse!
  h <- a_sq + b_sq

  # Euclidean distance is the square root of a^2 + b^2
  h <- sqrt(h)

  return(h)

}


#' Find the distances to nearest neighbor for a collection of Cartesian coordinates
#' @description Given a data frame containing X and Y Cartesian coordinates, find the distance to the nearest point for each point in the data frame.
#' @param dataframe A data frame. Must contain numeric variables corresponding to \code{x_var} and \code{y_var} containing the X and Y coordinate values.
#' @param x_var Character string. The name of the variable in \code{dataframe} containing the X components of the coordinates.
#' @param y_var Character string. The name of the variable in \code{dataframe} containing the Y components of the coordinates.
#' @return A numeric vector of the distances, one for each point.
#' @export
NN <- function(dataframe,
               x_var,
               y_var) {
  if (class(dataframe) != "data.frame") {
    stop("dataframe must be a data frame")
  }
  if (!(x_var %in% names(dataframe))) {
    stop("x_var must correspond to the name of a variable in dataframe")
  }
  if (!(y_var %in% names(dataframe))) {
    stop("y_var must correspond to the name of a variable in dataframe")
  }
  if (!is.numeric(dataframe[[x_var]])) {
    stop("x_var must correspond to the name of a numeric variable in dataframe")
  }
  if (!is.numeric(dataframe[[y_var]])) {
    stop("y_var must correspond to the name of a numeric variable in dataframe")
  }

  # Get a distance matrix
  distance_matrix <- dist_matrix(dataframe = dataframe,
                                 x_var = x_var,
                                 y_var = y_var)

  # How many points are there?
  point_count <- nrow(dataframe)

  # This just "gets rid of" the difference between a point and itself
  # Enables min() function to find the true nearest neighbor of the point
  distance_matrix[distance_matrix == 0] <- Inf

  # Get the distance to the nearest neighbor for each point
  nearest_dists <- sapply(1:point_count,
                          distances = distance_matrix,
                          FUN = function(X, distances) {
                            min(distances[, X])
                          })

  return(nearest_dists)
}

#' Find the mean distance to nearest neighbor for a collection of Cartesian coordinates
#' @description Given a data frame containing X and Y Cartesian coordinates, find the arithmetic and geometric mean distance to the nearest point.
#' @param dataframe A data frame. Must contain numeric variables corresponding to \code{x_var} and \code{y_var} containing the X and Y coordinate values.
#' @param x_var Character string. The name of the variable in \code{dataframe} containing the X components of the coordinates.
#' @param y_var Character string. The name of the variable in \code{dataframe} containing the Y components of the coordinates.
#' @return A named numeric vector of the mean distances: \code{"arith_mean"} containing the arithmetic mean and \code{"geo_mean"} containing the geometric mean.
#' @export
NN_mean <- function(dataframe,
                    x_var,
                    y_var){
  if (class(dataframe) != "data.frame") {
    stop("dataframe must be a data frame")
  }
  if (!(x_var %in% names(dataframe))) {
    stop("x_var must correspond to the name of a variable in dataframe")
  }
  if (!(y_var %in% names(dataframe))) {
    stop("y_var must correspond to the name of a variable in dataframe")
  }
  if (!is.numeric(dataframe[[x_var]])) {
    stop("x_var must correspond to the name of a numeric variable in dataframe")
  }
  if (!is.numeric(dataframe[[y_var]])) {
    stop("y_var must correspond to the name of a numeric variable in dataframe")
  }

  # Get the vector of nearest neighbor distances
  nearest_dists <- NN(dataframe,
                      x_var,
                      y_var)

  # Calculate arithmetic mean
  arith_mean <- mean(nearest_dists)

  # Calculate the geometric mean
  geo_mean <- gm_mean(nearest_dists,
                      na.rm = TRUE)

  # Return the two values
  return(c(arith_mean = arith_mean, geo_mean = geo_mean))
}


#' Test points spatial balance against sets of random points
#' @description Compare a set of points to sets of random points generated from the same polygon geometry and report back the proportion of random sets which had higher mean distance to nearest neighboring point. Assuming that a higher mean distance to nearest neighbor indicates greater spatial balance, this proportion can be treated as a "probabilty that the points in \code{spdf} are not spatially balanced." Whatever is provided as \code{spdf} will be dissolved without regard for the data slot, so if you want to test subsets of the polygons, each test will need to be a separate function call provided only the relevant subset as \code{spdf}, e.g. in the case of wanting to test individual strata stored as a single SPDF you would need to call this function for each stratum (probably in a \code{lapply()} or a loop).
#' @param number Numeric. The number of sets to generate to compare. Defaults to \code{100}.
#' @param points Spatial Points Data Frame. The points that are being compared against.
#' @param polygons Spatial Polygons Data Frame. Polygons describing the boundaries of the area of interest that corresponds to \code{points}.
#' @param method Character string. Which method to use for generating random sets of points to compare against. Either \code{"sample"} to use \code{sf::st_sample()} which is the much faster option or \code{"probability"} which uses a legacy approach that depends on cumulative proportional areas of subpolygons as probability of selection for random points. Defaults to \code{"sample"}.
#' @param seed_number Numeric. The number to use in \code{set.seed()} for reproducibility. Defaults to \code{420}.
#' @param projection CRS object. The projection to force all spatial objects into to for the purpose of compatibility. Defatults to \code{sp::CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")}.
#' @return Named numeric vector. The value for \code{"p_arith"} is the proportion of comparisons that had a higher arithmetic mean nearest neighbor distance than \code{points} and \code{"p_geom"} is the proportion of comparisons that had a higher geometric mean nearest neighbor distance.
#' @export

test_points <- function(number = 100,
                        points,
                        polygons,
                        method = "sample",
                        seed_number = 420,
                        projection = sp::CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")){
  if (number < 0) {
    stop("number must be a positive integer")
  }

  if (!(method %in% c("sample", "probability"))) {
    stop("The argument method must either be 'sample' or 'probability'.")
  }

  if (!grepl(class(points), pattern = "^SpatialPoints")) {
    stop("points must be a spatial points data frame")
  }

  # We need to handle what to do if the geometry is empty
  if (nrow(points) < 1) {
    stop("There's no geometry in points")
  }

  if (!identical(projection, points@proj4string)) {
    points <- sp::spTransform(points,
                              projection)
  }

  if (!grepl(class(polygons), pattern = "^SpatialPolygons")) {
    stop("polygons must be a spatial polygons data frame")
  }

  # We need to handle what to do if the geometry is empty
  if (length(polygons@polygons) < 1) {
    stop("There's no geometry in polygons")
  }
  if (!identical(projection, polygons@proj4string)) {
    polygons <- sp::spTransform(polygons,
                                projection)
  }

  # And if it's not dissolved, we'll do that!
  if (length(polygons@polygons) > 1) {
    message("The polygons in polygons need to be dissolved. Dissolving now.")
    polygons <- methods::as(sf::st_combine(sf::st_as_sf(polygons)), "Spatial")
  }


  # NAD83 CRS for projecting
  projectionNAD83 <- sp::CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")
  # Alber's equal area CRS for projecting
  projectionAL <- sp::CRS("+proj=aea")

  # Grab the point coordinates to calculate the mean nearest neighbor
  pts_coords <- get_coords(points,
                           x_var = "XMETERS",
                           y_var = "YMETERS",
                           projection = projectionAL)@data

  # Calculate the mean nearest neighbor
  nn_means <- NN_mean(pts_coords,
                      x_var = "XMETERS",
                      y_var = "YMETERS")

  # Just splitting them out for clarity later
  design_nn_am <- nn_means[["arith_mean"]]
  design_nn_gm <- nn_means[["geo_mean"]]


  # Generate the cumulative Prob Distribution!
  if (method == "probability") {
    if (!grepl(class(polygons), pattern = "^SpatialPolygons")) {
      stop("polygons must be a spatial polygons data frame")
    }

    # We need to handle what to do if the geometry is empty
    if (length(polygons@polygons) < 1) {
      stop("There's no geometry in polygons")
    }

    # And if it's not dissolved, we'll do that!
    if (length(polygons@polygons) > 1) {
      message("The polygons in polygons need to be dissolved. Dissolving now.")
      polygons <- methods::as(sf::st_combine(sf::st_as_sf(polygons)), "Spatial")
    }

    probability_distribution <- extract_poly_area(polygons)
  }

  if (method == "sample") {
    # This is so we can generate random points with sf::st_sample()
    # An approach that stands to be like an order of magnitude faster
    probability_distribution <- NULL
    polygons <- sf::st_as_sf(x = polygons)
  }

  # How many points we want.
  point_count <- nrow(points)

  # Generate the sets of points. We're doing it over a vector of seed numbers with length = number.
  # That means that each set of generated points has their own starting seed and are therefore unique
  means <- do.call(rbind,
                   lapply(X = seed_number + 1:number,
                          point_count = point_count,
                          probability_distribution = probability_distribution,
                          nn_means = nn_means,
                          polygons = polygons,
                          projection = projection,
                          FUN = function(X, point_count, probability_distribution, nn_means, polygons, projection){
                            message("MASTER SEED NUMBER ", X)

                            # OKAY. So, we're going to implement NOT using a probability distribution
                            # This could be like a billion times more efficient
                            if (is.null(probability_distribution)) {
                              set.seed(X)
                              # Generate a random set of points within the polygons
                              # quiet = TRUE doesn't seem to have an effect here, so I'm using suppressMessages()
                              # Otherwise there are constant messages that "although coordinates are longitude/latitude, st_intersects assumes that they are planar"
                              # I'm just noting in the documentation that working in polar regions or with truly huge polygons are the only times that's relevant
                              rand_sf <- suppressMessages({sf::st_sample(x = polygons,
                                                                         size = point_count,
                                                                         type = "random",
                                                                         exact = TRUE,
                                                                         quiet = TRUE)})
                              rand_spdf <- methods::as(rand_sf, "Spatial")
                              # Translate from SpatialPoints to SpatialPointsDataFrame
                              # The ID number is just the order in which they were generated
                              rand_spdf <- sp::SpatialPointsDataFrame(coords = rand_spdf@coords,
                                                                      data = data.frame(id = 1:nrow(rand_spdf@coords)))
                              sp::proj4string(rand_spdf) <- projection
                            } else {
                              # Pick up point_count * 2 points then check them out.
                              # Because we use the polygon bounding box and not the polygon later, we may select points outside of polygon area.
                              # The extra pointcount helps to account for those non-overlapping points.
                              # We drop any excess points anyway
                              iterations <- round(point_count * 2)

                              # Get our seed numbers for these iterations
                              set.seed(X)
                              current_seeds <- sample(1:99999,
                                                      size = iterations)

                              # Generate the random points
                              # Frankly, I'm not sure why Steve did it this way instead of like sf::st_sample()
                              rand_points <- do.call(rbind,
                                                     lapply(X = current_seeds,
                                                            probability_distribution = probability_distribution,
                                                            polygons = polygons,
                                                            FUN = function(X, probability_distribution, polygons){
                                                              # We set the seed number any time it might get triggered
                                                              # Use the seed number we generated for this iteration
                                                              set.seed(X)

                                                              # Get a uniform random variate for selecting a polygon. This gets compared against the cumulative frequency
                                                              urv <- runif(1)

                                                              # Using the urv, determine the polygon number (opt) from the cumulative freq distribution
                                                              poly_index <- select_from_distribution(dataframe = probability_distribution,
                                                                                                     prob_var = "cum_freq",
                                                                                                     id_var = "id",
                                                                                                     value = urv)

                                                              # If dissolved (as it should be for polygons), then always access polygons[[1]].
                                                              poly <- polygons@polygons[[1]]@Polygons[[poly_index]]

                                                              # MAKE IT REPRODUCIBLE!!!!
                                                              set.seed(X)

                                                              # Use the bounding box of the polygon and sp:spsample() to select just 1 random point.
                                                              rand_point <- sp::spsample(poly,
                                                                                         n = 1,
                                                                                         type = "random"#,
                                                                                         #bb = sp::bbox(poly)
                                                              )

                                                              return(rand_point)
                                                            }))

                              # Translate from SpatialPoints to SpatialPointsDataFrame
                              # The ID number is just the order in which they were generated
                              rand_points_spdf <- sp::SpatialPointsDataFrame(coords = rand_points,
                                                                             data = data.frame(id = 1:nrow(rand_points@coords)))
                              sp::proj4string(rand_points_spdf) <- polygons@proj4string

                              # Find overlap!
                              overlap <- sp::over(rand_points_spdf,
                                                  polygons)

                              # Only keep the points where there was spatial overlap
                              rand_points_spdf <- rand_points_spdf[!is.na(overlap[[1]]), ]

                              # In case we have more random points than needed, just pick the first nrow(points) points.
                              # We still have a fully random sample since we effectively store the random points by accession
                              # and we elimiinate points from the bottom up.
                              rand_spdf <- rand_points_spdf[1:point_count, ]
                            }

                            # Add the x and y meters now
                            rand_spdf <- get_coords(spdf = rand_spdf,
                                                    x_var = "XMETERS",
                                                    y_var = "YMETERS",
                                                    projection = projectionAL)

                            # Get the distance to nearest neighbor for each point in the random set
                            rand_nn <- NN_mean(dataframe = rand_spdf@data,
                                               x_var = "XMETERS",
                                               y_var = "YMETERS")

                            return(data.frame("am" = rand_nn["arith_mean"],
                                              "gm" = rand_nn["geo_mean"]))
                          })
  )

  # This is the whole damn goal:
  # What proportion of our comparison draws had a HIGHER mean nearest neighbor distance than the points?
  # We assume that a higher mean nearest neighbor distance means more spatially balanced
  # So these are effectively the proportions of draws that are MORE BALANCED than the points
  # Which means that they're basically a P value for the H0 "The points provided are spatially balanced"
  nn_am_greater_prop <- sum(means[["am"]] > design_nn_am) / number
  nn_gm_greater_prop <- sum(means[["gm"]] > design_nn_gm) / number

  output <- c(p_arith = nn_am_greater_prop, p_geom = nn_gm_greater_prop)
  rownames(output) <- NULL
  return(output)
}


#' Add coordinates to the data slot of a SpatialPointsDataFrame
#' @param spdf A SpatialPointsDataFrame.
#' @param x_var Character string. The name of the variable to add the x component of the coordinates to in \code{spdf@@data}. If the variable already exists, it will be overwritten. Defaults to \code{"XMETERS"}.
#' @param y_var Character string. The name of the variable to add the y component of the coordinates to in \code{spdf@@data}. If the variable already exists, it will be overwritten. Defaults to \code{"YMETERS"}.
#' @param projection CRS object. The projection to use when determining the coordinates. Defaults to \code{sp::CRS("+proj=aea")}, Albers equal area.
#' @return The SPDF provided as \code{spdf} with the coordinates added to \code{"x_var"} and \code{"y_var"}.
#' @export
get_coords <- function(spdf,
                       x_var = "XMETERS",
                       y_var = "YMETERS",
                       projection = CRS("+proj=aea")){
  if (!grepl(class(spdf), pattern = "^SpatialPointsDataFrame")) {
    stop("spdf must be a spatial points data frame")
  }
  # Create a reprojected spdf to grab coords from
  temp_spdf <- sp::spTransform(spdf,
                               projection)

  # Write them into the original spdf
  spdf@data[[x_var]] <- temp_spdf@coords[, 1]
  spdf@data[[y_var]] <- temp_spdf@coords[, 2]

  return(spdf)
}

#' Test to see if a set of points are spatially balanced within the polygons used to draw them
#' @description Given a set of points and the polygons used to draw them, test the spatial balance of the point by comparing them to randomly located points generated within the polygons
#' @param polygons_spdf Spatial polygons data frame. The polygons that were used to draw the points. This can either be the sample frame for the points or stratification polygons. The balance check will be done for the whole frame if \code{by_frame} is \code{TRUE}. The balance check will be done by polygon identity if \code{polygons_spdf@@data} has an identity variable and that variable name is provided as the argument \code{stratafield}.
#' @param points_spdf Spatial points data frame. The points to be tested for spatial balance.
#' @param stratafield Character string. The name of the variable in \code{polygons_spdf@@data} that contains the polygon identities. If \code{NULL} then the point balance won't be checked by polygon ID. Defaults to \code{NULL}.
#' @param by_frame Logical. If \code{TRUE} then a balance check for the points will be done with the full extent of \code{polygons_spdf} ignoring polygon identities. Defaults to \code{TRUE}.
#' @param reps Numeric. The number of random draws to make to compare against \code{points}. If this is larger than \code{100} then the process can start to take more than a few minutes if \code{points} contains more than a few dozen points. Defaults to \code{100}.
#' @param seed_number Numeric. The number to supply to \code{set.seed()} for reproducibility. At multiple steps, this seed number may be used to generate additional seed numbers for function-internal use, but always reproducibly. Defaults to \code{420}.
#' @param method Character string. Which method to use for generating random sets of points to compare against. Either \code{"sample"} to use \code{sf::st_sample()} which is the much faster option or \code{"probability"} which uses a legacy approach that depends on cumulative proportional areas of subpolygons as probability of selection for random points. Defaults to \code{"sample"}.
#' @param projection CRS object. The projection to force all spatial objects into to for the purpose of compatibility. Defatults to \code{sp::CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")}.
#' @return A data frame with the variables \code{polygon} (The polygon identity, either \code{"Sample Frame"} or the ID from \code{polygons_spdf@@data$stratafield} as appropriate), \code{point_count} (Number of points from \code{points_spdf} occurring in the polygon), \code{reps} (Number of random draws compared against), \code{mean_arithmetic} (The arithmetic mean of the nearest neighbor distances for the points in \code{points_spdf}), \code{mean_geometric} (The geometric mean of the nearest neighbor distances for the points in \code{points_spdf}), \code{p_arithmetic} (The proportion of random point draws that had larger arithmetic mean neighbor distances than \code{points_spdf}), and \code{p_geometric} (The proportion of random point draws that had larger geometric mean neighbor distances than \code{points_spdf}). We treat the \code{p_geometric} as the p value for testing the H0 that \code{points_spdf} is balanced.
#' @export
check_balance <- function(polygons_spdf,
                          points_spdf,
                          reps = 100,
                          stratafield = NULL,
                          by_frame = TRUE,
                          seed_number = 420,
                          method = "sample",
                          projection = sp::CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"))

{
  # Reproject if necessary
  if (!(class(polygons_spdf) %in% "SpatialPolygonsDataFrame")) {
    stop("polygons_spdf must be a spatial polygons data frame")
  if (!(method %in% c("sample", "probability"))) {
    stop("The argument method must either be 'sample' or 'probability'.")
  }
  if (!identical(projection, polygons_spdf@proj4string)) {
    polygons_spdf <- sp::spTransform(polygons_spdf,
                                     projection)
  }
  if (!(class(points_spdf) %in% "SpatialPointsDataFrame")) {
    stop("points_spdf must be a spatial points data frame")
  }
  if (!identical(projection, points_spdf@proj4string)) {
    points_spdf <- sp::spTransform(points_spdf,
                                   projection)
  }

  if (!is.null(stratafield)) {
    if (!is.character(stratafield)) {
      stop("stratafield must be a single character string or NULL")
    }
    if (length(stratafield) > 1) {
      stop("stratafield must be a single character string or NULL")
    }
    if (!(stratafield %in% names(polygons_spdf@data))) {
      stop(stratafield, " is not a variable in polygons_spdf@data or NULL")
    }
  }

  if (reps < 1) {
    stop("reps must be a positive integer")
  }
  if (floor(reps) != ceiling(reps)) {
    stop("reps must be a positive integer")
  }
  # Add the coordinates so that we can do nearest neighbor calculations
  points_spdf <- get_coords(points_spdf)

  # The frame output is assumed to be NULL
  output_frame <- NULL

  # If requested, first analyze entire frame
  if(by_frame) {
    # Derive the arithmetic and geometric mean distance to nearest neighbor for the points
    nn_means_all <- NN_mean(points_spdf@data,
                            x_var = "XMETERS",
                            y_var = "YMETERS")

    ## Do randomization test
    proportions_frame <- test_points(number = reps,
                                     points = points_spdf,
                                     polygons = polygons_spdf,
                                     method = method,
                                     seed_number = seed_number)

    # Build the output data frame for the sample frame
    output_frame <- data.frame("polygon" = "Sample Frame",
                               "point_count" = nrow(points_spdf),
                               "reps" = reps,
                               "mean_arithmetic" = nn_means_all["arith_mean"],
                               "mean_geometric" = nn_means_all["geo_mean"],
                               "p_arithmetic" = proportions_frame["p_arith"],
                               "p_geometric" = proportions_frame["p_geom"])

  }

  ######### Analyze by strata if requested

  # Assume that there aren't strata analyses happening by default
  output_strata <- NULL

  # If there's a stratification field and it exists in the spdf, get to work
  if(!is.null(stratafield)) {
    if (!(stratafield %in% names(polygons_spdf@data))) {
      stop(paste("The variable", stratafield, "does not appear in frame_spdf@data."))
    }

    # Just to simplify things, make a STRATUM variable
    polygons_spdf@data[["STRATUM"]] <- polygons_spdf@data[[stratafield]]
    strata <- as.character(unique(polygons_spdf@data[["STRATUM"]]))

    output_strata <- do.call(rbind,
                             lapply(X = strata,
                                    strata_spdf = polygons_spdf,
                                    points = points_spdf,
                                    seed_number = seed_number,
                                    FUN = function(X, strata_spdf, points, seed_number){
                                    method = method,
                                      # For clarity
                                      stratum <- X
                                      message(stratum)

                                      # Get just the relevant polygons
                                      stratum_spdf <- strata_spdf[strata_spdf[["STRATUM"]] == stratum, ]

                                      # Get just the points that fall in this stratum
                                      current_points <- points
                                      current_points@data[["STRATUM"]] <- sp::over(points,
                                                                                   stratum_spdf)[["STRATUM"]]
                                      current_points <- points[!is.na(current_points@data[["STRATUM"]]), ]

                                      # If there are in fact points in the stratum, do the randomization test
                                      if (nrow(current_points) > 0) {

                                        # Derive the arithmetic and geometric mean distance to nearest neighbor for the points
                                        nn_means_stratum <- NN_mean(current_points@data,
                                                                    x_var = "XMETERS",
                                                                    y_var = "YMETERS")

                                        ## Do randomization test
                                        proportions_stratum <- test_points(number = reps,
                                                                           points = current_points,
                                                                           polygons = stratum_spdf,
                                                                           method = method,
                                                                           seed_number = seed_number)

                                        # Build the output data frame for the sample frame
                                        output <- data.frame("polygon" = stratum,
                                                             "point_count" = nrow(current_points),
                                                             "reps" = reps,
                                                             "mean_arithmetic" = nn_means_stratum[["arith_mean"]],
                                                             "mean_geometric" = nn_means_stratum[["geo_mean"]],
                                                             "p_arithmetic" = proportions_stratum["p_arith"],
                                                             "p_geometric" = proportions_stratum["p_geom"])

                                      } else {
                                        # If there weren't points in the stratum, just give us the empty output
                                        output <- data.frame("polygon" = stratum,
                                                             "point_count" = 0,
                                                             "reps" = NA,
                                                             "mean_arithmetic" = NA,
                                                             "mean_geometric" = NA,
                                                             "p_arithmetic" = NA,
                                                             "p_geometric" = NA)
                                      }
                                    })
    )
  }
  # Combine the outputs
  output <- rbind(output_frame, output_strata)

  return(output)
}
