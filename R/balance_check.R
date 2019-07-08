## RandomizePV4.R - 5/16/2018; SLGarman - BLM NOC, Denver.  Version 4.  Uses spsample and bb=bbox() to locate random samples in GenPts(). And, this
##                                                          version can use (buffered) NHD lines to evaluate aquatic AIM designs.

## Determines if the dispersion of points in a survey design is different from random.  Developed specifically to evaluate the spatial balance of merged (i.e., expanded) GRTS designs.
## For a specified Frame and n GRTS points, derives x number of random sets of n points, and compares
## the mean nearest neighbor (MNN) distance of the input points with that of each set of randomly selected points.

## Records the proportion of random-point sets with a MNN equal to or greater than the survey-design (GRTS) points.
## This proportion can be treated as a P value to test:
##	Ho: GRTS design is not different from random
##	Ha: GRTS design is different from random

##      E.g., a proportion of 0.03 means that the input GRTS points are different from random.
##            A proportion of 0.20 means that the input GRTS points are Likely Not different from random - that is,
##              randomly selected point sets were more dispersed than the GRTS points in 20% of the replications.

## Mean NN is based on Arithmetic mean and Geometric mean.  Arithmetic-mean results are most applicable when there is 1 highly dominant polygon per feature.
## That is, the frame or a strata is represented almost entirely by 1 polygon.  Where a frame or stratum is represented by multiple disjunct
## (and highly dispersed) polygons, the Geometric mean should be considered (better represents the central tendency of highly skewed data).  If random proportions of
##  both measures are > 0.1, then best to retain Ho:, else good bet to go with Ha: especially if the Geometric mean is <= 0.05.


## This version only works on Polygonal sample frames and strata.  Additionally, this version includes the option to perform randomization tests on a strata basis.


## I.  AIM Terrestrial analyses.  For both Frame and Strata, requires polygons to be dissolve - only 1 feature per attribute.  The strata of each point in the point file
## is derived by overlaying the points on the strata file.  Random() [see bottom of the script] is called to process AIM terrestrial designs.

## II.  AIM Aquatic analyses.  This was an add-on.  The strata and frame in this case are the same, and consists of buffered NHD lines (buffer by say 0.5 m to create a poly-line
## file - easier to process using this algorithm). Unlike AIM Terrestrial frames and strata, NHD poly-lines SHOULD NOT be dissolved; multiple polygons per feature (i.e., stream order segment)
## are EXPECTED.  The user specifies the name of the strata field in the NHD 'strata' file.  It is assumed that the Aquatic points file contains a field called STRM_ORDR that defines
## the strata (stream order), and STRM_ORDR values match the strata values in the NHD poly-line file.  STRM_ORDER is directly used to determine the strata of points - the script
## doesn't overlay points onto the NHD file to determine strata like it does with Terrestrial AIM points (aquatic point locations are field-derived coordinates and do not
## always overlap even the buffered NHD poly-line).  RandomAquatic() is called to process AIM Aquatic designs.



# How to run this script.
#1)  Copy all of section I into the R console

#2)  Modify the I/O arguments in section II, then cut and paste arguments and the call to Random() to the R console.
## Specify the working directory where the following shapefiles reside
## Specify the polygonal frame (shapefile)
## Specify the GRTS points (the input points shapefile)
## Specify the number of random replications (500 should do the trick, can try more depending on computer speed)
## Specify the output file name (summary of results - short and sweet)
## Specify the random number seed or leave as is
## Specify a strata file for an analysis based on strata, else set this to NA.  For terrestrial AIM analyses, strata file must be dissolved.
## Specify the strata field name in the strata file, else set this to NA
## Regardless if a strata file/analysis is performed, can specify an entire frame analysis that ignores strata

#3)  OR use Batch mode to run the script after you modify I/O arguments.  Your system path must include the bin subdirectory of the R version you want to use.
#     Do the following in a command-line console:   R CMD BATCH randomizepv4.r capture_console_output

#     capture_console_output is an ASCII file that records the output you'd see in the R console when running from the console.
#     You can batch multiple jobs at once (rename your modified versions of randomizepv4.r & use different file names to capture the results),
#     grab lunch, then come back and check out all of the results.......



################################################################
#' Get the areas of polygons
#' @description Given a Spatial Polygons Data Frame, produce a data frame summarizing the areas of polygons sharing an identity
#' @param polygons Spatial Polygons or Spatial Polygons Data Frame. The polygons to be summarized.
#' @param area_prop Logical. If \code{TRUE} then the output data frame will include the variable \code{"area_prop"} with the proportional areas of the polygon identities. Defaults to \code{TRUE}.
#' @param cum_freq Logical. If \code{TRUE} then the output data frame will include the variable \code{"cum_freq"} with the cumulative frequencies of the polygon identities. Note that it sorts the polygons from smallest to largest area before doing this. Defaults to \code{TRUE}.
#' @return A data frame containing the polygon index within \code{polygons@@polygons[[1]]@@Polygons}, and if asked for, the proportional area (\code{"area_prop"}) and the cumulative frequency (\code{"cum_freq"}).
#' @export
extract_poly_area <- function(polygons,
                              area_prop = TRUE,
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
  if (area_prop) {
    areas_df[["area_prop"]] <- areas_df[["area"]] / total_area
  }

  # Add cumulative frequency distribution, which can be treated as a probability distribution
  if (cum_freq) {
    areas_df[["cum_freq"]] <- cumsum(areas_df[["area_prop"]])
  }

  # Return this data frame!
  return(areas_df)
}

#########################################################################################################
#  Extract area of polygons & build a cumulative Prob. Distribution indexed by polygon number (spdf is the strata or frame)
## This only works if polygons were NOT dissolved - designed for aquatic lines buffered to form a polygonal frame
# TODO: MAKE THIS HANDLE THE SUB-POLYGONS CORRECTLY!!!! IT FUCKS UP GenPts() WHEN TRYING TO SELECT POLYGONS BECAUSE THE ID ISN'T CORRECT
# The problem being that these aren't dissolved, so there's the polygon and then subpolys as opposed to just polygons
ExtractPolyAreaAquatic <- function(spdf) {
  if (class(spdf) != "SpatialPolygonsDataFrame") {
    stop("spdf must be a spatial polygons data frame")
  }

  # Get the areas of the polygons
  areas <- sapply(X = spdf@polygons,
                  FUN = function(X) {
                    sapply(X = X@Polygons,
                           FUN = function(X) {
                             X@area
                           })
                  })

  # Make a data frame with the areas and the within-polygon ID
  areas_df <- data.frame(area = areas,
                         id = 1:length(areas))

  # Sort from largest to smallest area
  # (This can help speed up selecting from the probability distribution)
  areas_df <- areas_df[order(-areas_df[["area"]]), ]

  # Get the total area
  total_area <- sum(areas_df[["area"]])

  # Add proportional area to the data frame
  areas_df[["area_prop"]] <- areas_df[["area"]] / total_area

  # Add cumulative frequency distribution, which can be treated as a probability distribution
  areas_df[["cum_freq"]] <- cumsum(areas_df[["area_prop"]])

  # Return this data frame!
  return(areas_df)
}

###############################################################################
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


###########################################################################GenPts(number,thepts,tempS,StrataNN,StrataBox)
## Generate sets of random points
# derive Mean NN of each set
# compare with the specified NN (nearest_dists) of the GRTS (input) pts
# record the proportion of random replicates where the random NN >= NN of the input pts.
#' @description Compare a set of points to sets of random points generated from the same polygon geometry and report back the proportion of random sets which had higher mean distance to nearest neighboring point. Assuming that a higher mean distance to nearest neighbor indicates greater spatial balance, this proportion can be treated as a "probabilty that the points in \code{spdf} are not spatially balanced." Whatever is provided as \code{spdf} will be dissolved without regard for the data slot, so if you want to test subsets of the polygons, each test will need to be a separate function call provided only the relevant subset as \code{spdf}, e.g. in the case of wanting to test individual strata stored as a single SPDF you would need to call this function for each stratum (probably in a \code{lapply()} or a loop).
#' @param number Numeric. The number of sets to generate to compare. Defaults to \code{3}.
#' @param pts_spdf Spatial Points Data Frame. The points that are being compared against.
#' @param aoi_spdf Spatial Polygons Data Frame. Polygons describing the boundaries of the area of interest that corresponds to \code{pts_spdf}.
#' @param type Numeric. Use \code{1} when aoi_spdf is a dissolved set of polygons and \code{2} when aoi_spdf is an undissolved set of polylines. Defaults to \code{1}
#' @param seed_number Numeric. The number to use in \code{set.seed()} for reproducibility. Defaults to \code{1}.
#' @return Named numeric vector. The value for \code{"p_arith"} is the proportion of comparisons that had a higher arithmetic mean nearest neighbor distance than \code{pts_spdf} and \code{"p_geom"} is the proportion of comparisons that had a higher geometric mean nearest neighbor distance.
#' @export

test_points <- function(number,
                        pts_spdf,
                        aoi_spdf,
                        type = 1, ## type ==2 for aquatic analysis where poly-lines are not dissolved, else 1 TOD
                        seed_number = 1){
  if (number < 0) {
    stop("number must be a positive integer")
  }
  if (!grepl(class(pts_spdf), pattern = "^SpatialPoints")) {
    stop("pts_spdf must be a spatial points data frame")
  }

  # We need to handle what to do if the geometry is empty
  if (nrow(pts_spdf) < 1) {
    stop("There's no geometry in pts_spdf")
  }

  if (!grepl(class(aoi_spdf), pattern = "^SpatialPolygons")) {
    stop("aoi_spdf must be a spatial polygons data frame")
  }

  # We need to handle what to do if the geometry is empty
  if (length(aoi_spdf@polygons) < 1) {
    stop("There's no geometry in aoi_spdf")
  }

  # And if it's not dissolved, we'll do that!
  if (length(aoi_spdf@polygons) > 1) {
    message("The polygons in aoi_spdf need to be dissolved. Dissolving now.")
    aoi_spdf <- methods::as(sf::st_combine(sf::st_as_sf(aoi_spdf)), "Spatial")
  }


  # NAD83 CRS for projecting
  projectionNAD83 <- sp::CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")
  # Alber's equal area CRS for projecting
  projectionAL <- sp::CRS("+proj=aea")

  # Grab the point coordinates to calculate the mean nearest neighbor
  pts_coords <- pts_spdf@coords
  names(pts_coords) <- c("x", "y")

  # Calculate the mean nearest neighbor
  nn_means <- NN_mean(pts_coords,
                      x_var = "x",
                      y_var = "y")

  # Just splitting them out for clarity later
  design_nn_am <- nn_means[["arith_mean"]]
  design_nn_gm <- nn_means[["geo_mean"]]

  # Count is just used to trigger telling the user how things are going at major intervals
  count <- 0

  # The requested number of random simulations
  samplesize <- number

  # The proportion of random sets with a mean NN => than MNN of the input pts
  nn_am_greater_count <- 0
  nn_gm_greater_count <- 0


  # Generate the cumulative Prob Distribution!
  if (type == 1) {
    if (!grepl(class(aoi_spdf), pattern = "^SpatialPolygons")) {
      stop("aoi_spdf must be a spatial polygons data frame")
    }

    # We need to handle what to do if the geometry is empty
    if (length(aoi_spdf@polygons) < 1) {
      stop("There's no geometry in aoi_spdf")
    }

    # And if it's not dissolved, we'll do that!
    if (length(aoi_spdf@polygons) > 1) {
      message("The polygons in aoi_spdf need to be dissolved. Dissolving now.")
      aoi_spdf <- methods::as(sf::st_combine(sf::st_as_sf(aoi_spdf)), "Spatial")
    }

    probability_distribution <- extract_poly_area(aoi_spdf)
  }
  if (type == 2) {
    probability_distribution <- ExtractPolyAreaAquatic(aoi_spdf)
  }

  while(number > 0) {
    # How many points we want. On each pass, counter is reduced by however many points were just grabbed
    # I think that should always be 1, but just to be safe I'm leaving Steve's incrementer
    counter <- point_count <- nrow(pts_spdf)

    while(counter > 0) {
      # Pick up counter * 1.25 draws then check them out.
      # Because we use the polygon bounding box and not the polygon later, we may select points outside of polygon area.
      # The 25% extra helps to account for non-overlapping points.
      # By adding perhaps more than we need, we at least cut down on repeat conversion from Spatial Points to SpatialPointsDataFrame
      # We drop any excess points anyway

      iterations <- round(counter * 1.25)

      # Get our seed numbers for these iterations
      # We add the value of number and counter so that we don't get identical seeds on different passes
      set.seed(seed_number + number + counter)
      current_seeds <- sample(1:99999,
                              size = iterations)
      for(iteration in 1:interations) {
        # We set the seed number any time it might get triggered
        # Use the seed number we generated for this iteration
        set.seed(current_seeds[iteration])

        # Get a uniform random variate for selecting a polygon. This gets compared against the cumulative frequency
        urv <- runif(1)

        # Using the urv, determine the polygon number (opt) from the cumulative freq distribution
        poly_index <- select_from_distribution(dataframe = probability_distribution,
                                               prob_var = "cum_freq",
                                               id_var = "id",
                                               value = urv)

        # If dissolved (as it should be for polygons), then always access polygons[[1]].
        if (type == 1) {
          poly <- aoi_spdf@polygons[[1]]@Polygons[[poly_index]]
        }
        # If not dissolved (a polylines situation) then things are different
        if (type == 2) {
          poly <- aoi_spdf@polygons[[poly_index]]@Polygons[[1]]
        }

        # MAKE IT REPRODUCIBLE!!!!
        set.seed(current_seeds[iteration])

        # Use the bounding box of the polygon and sp:spsample() to select just 1 random point.
        rand_point <- sp::spsample(poly,
                                   n = 1,
                                   type = "random",
                                   bb = sp::bbox(poly))

        # Build a data frame of all the randomly selected points as we loop
        if(is.null(rand_coords_df)) {
          rand_points <- rand_point
        } else {
          rand_points <- rbind(rand_points, rand_point)
        }
      }


      # Translate from SpatialPoints to SpatialPointsDataFrame
      sp::proj4string(rand_points) <- projectionNAD83
      # The ID number is just the order in which they were generated
      rand_points_spdf <- sp::SpatialPointsDataFrame(coords = rand_points,
                                                     data = data.frame(id = 1:nrow(rand_points)))

      # Find overlap!
      overlap <- sp::over(rand_points_spdf,
                          aoi_spdf)

      # Only keep the points where there was spatial overlap
      rand_points_spdf <- rand_points_spdf[!is.na(overlap[[1]]), ]

      ## Bind the points that overlap aoi_spdf, and adjust counter.
      if(nrow(rand_points_spdf) > 0) {
        if(is.null(store_spdf)) {
          store_spdf <- rand_points_spdf
        } else {
          store_spdf <- rbind(store_spdf, rand_points_spdf)
        }
        counter <- counter - nrow(temp_spdf)
      }
    }

    # In case we have more random points than needed, just pick the first nrow(pts_spdf) points.
    # We still have a fully random sample since we effectively store the random points by accession
    # and we elimiinate points from the bottom up.
    rand_spdf <- store_spdf[1:point_count, ]

    # Decrement rep counter
    number <- number - 1

    # Add the latitude and longitude to the data frame
    rand_spdf <- get_coords(spdf = rand_spdf,
                            x_var = "LONG",
                            y_var = "LAT",
                            projection = projectionNAD83)

    # Add the x and y meters now
    rand_spdf <- get_coords(spdf = rand_spdf,
                            x_var = "XMETERS",
                            y_var = "YMETERS",
                            projection = projectionAL)

    # Get the distance to nearest neighbor for each point in the random set
    rand_nn <- NN_mean(dataframe = rand_spdf@data,
                       x_var = "XMETERS",
                       y_var = "YMETERS")

    # Split the arithmetic and geometric means out for clarity
    rand_nn_am <- rand_nn[["arith_mean"]]
    rand_nn_gm <- rand_nn[["geo_mean"]]

    # If the mean nearest neighbor distance is >= that of the overall design
    # then increment the count for that stat up by 1
    if(rand_nn_am >= design_nn_am) {
      nn_am_greater_count <- nn_am_greater_count + 1
    }
    if(rand_nn_gm >= design_nn_gm) {
      nn_gm_greater_count <- nn_gm_greater_count + 1
    }

    # Just for show: Let's you see the evolution of the P value as number of reps increase.  Output to console every 100 reps.
    #print(number)
    count <- samplesize - number
    a <- count / 100
    b <- round(a, digits = 0)
    if(count == round(count / 100, digits = 0) * 100 & count != 0){
      message(paste0("Working on rep # = ",
                     count,
                     "\nEvolving Arithmetic P value = ",
                     nn_am_greater_count / count))
      message(paste0("Working on rep # = ",
                     count,
                     "\nEvolving Geometric P value = ",
                     nn_gm_greater_count / count))
    }
  }

  # Convert count to proportion
  nn_am_greater_prop <- nn_am_greater_count / samplesize
  nn_gm_greater_prop <- nn_gm_greater_count / samplesize

  # s.l. P value based on means
  output <- c(p_arith = nn_am_greater_prop, p_geom = nn_gm_greater_prop)
  return(output)
}

###############################################################################
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
  spdf@data[[x_var]] <- temp_spdf@coords[[1]]
  spdf@data[[y_var]] <- temp_spdf@coords[[2]]

  return(spdf)
}

test_point_balance <- function(aoi_spdf,
                               aoi_name,
                               pts_spdf,
                               reps,
                               type,
                               seed_number) {
  if (!grepl(class(aoi_spdf), pattern = "^SpatialPolygons")) {
    stop("aoi_spdf must be a spatial polygons data frame")
  }

  # We need to handle what to do if the geometry is empty
  if (length(aoi_spdf@polygons) < 1) {
    stop("There's no geometry in aoi_spdf")
  }

  # And if it's not dissolved, we'll do that!
  if (length(aoi_spdf@polygons) > 1) {
    message("The polygons in aoi_spdf need to be dissolved. Dissolving now.")
    aoi_spdf <- methods::as(sf::st_combine(sf::st_as_sf(aoi_spdf)), "Spatial")
  }

  # Add the coordinates so that we can do nearest neighbor calculations
  pts_spdf <- get_coords(pts_spdf)

  # Derive the arithmetic and geometric mean distance to nearest neighbor for the points
  nn_means <- NN_mean(pts_spdf@data,
                      x_var = "XMETERS",
                      y_var = "YMETERS")

  ## Do randomization test
  proportions_frame <- test_points(number = reps,
                                   pts_spdf = pts_spdf,
                                   aoi_spdf = aoi_spdf,
                                   type = type,
                                   seed_number = seed_number)

  # Build the output data frame for the sample frame
  output_frame <- data.frame("polygon" = aoi_name,
                             "point_count" = nrow(pts_spdf),
                             "reps" = reps,
                             "mean_arithmetic" = nn_means_all[["arith_mean"]],
                             "mean_geometric" = nn_means_all[["geo_mean"]],
                             "p_arithmetic" = proportions_frame["p_arith"],
                             "p_geometric" = proportions_frame["p_geom"])

  return(output)
}

######################################################################
#' Test to see if a set of points are spatially balanced within the polygons used to draw them
#' @description Given a
check_balance <- function(frame_spdf = NULL,		## the sample frame as a spdf (shapefile)
                          pts_spdf,		## the GRTS points as a spdf  (shapefile)
                          reps,		## Number of sets of random points
                          strata_spdf = NULL,	## Name of strata file (analyses will be conducted at the strata level) or NA
                          stratafield,	## Strata-field name in strata or NA
                          doFrame = TRUE,
                          seed_number = 1)

{

  # Add the coordinates so that we can do nearest neighbor calculations
  pts_spdf <- get_coords(pts_spdf)

  # If requested, first analyze entire frame
  if(doFrame) {
    # Derive the arithmetic and geometric mean distance to nearest neighbor for the points
    nn_means_all <- NN_mean(pts_spdf@data,
                            x_var = "XMETERS",
                            y_var = "YMETERS")

    ## Do randomization test
    proportions_frame <- test_points(number = reps,
                                     pts_spdf = pts_spdf,
                                     aoi_spdf = frame_spdf,
                                     type = 1,
                                     seed_number = seed_number)

    # Build the output data frame for the sample frame
    output_frame <- data.frame("polygon" = "Sample Frame",
                               "point_count" = nrow(pts_spdf),
                               "reps" = reps,
                               "mean_arithmetic" = nn_means_all["arith_mean"],
                               "mean_geometric" = nn_means_all["geo_mean"],
                               "p_arithmetic" = proportions_frame["p_arith"],
                               "p_geometric" = proportions_frame["p_geom"])

  } else {
    output_frame <- NULL
  }

  ######### Analyze by strata if requested
  # If there's a stratification field and it exists in the spdf, get to work
  if(!is.null(stratafield)) {
    if (!(stratafield %in% names(frame_spdf@data))) {
      stop(paste("The variable", stratfield, "does not appear in frame_spdf@data."))
    }

    # Just to simplify things, make a STRATUM variable
    frame_spdf@data[["STRATUM"]] <- frame_spdf@data[["STRATUM"]]
    strata <- unique(frame_spdf@data[["STRATUM"]])

    output_strata <- lapply(X = strata,
                            strata_spdf = frame_spdf,
                            pts_spdf = pts_spdf,
                            FUN = function(X, strata_spdf, pts_spdf){
                              # For clarity
                              stratum <- X

                              # Get just the relevant polygons
                              stratum_spdf <- strata_spdf[strata_spdf[["STRATUM"]] == stratum, ]

                              # Get just the points that fall in this stratum
                              current_pts_spdf <- pts_spdf
                              current_pts_spdf@data[["STRATUM"]] <- sp::over(pts_spdf,
                                                                             stratum_spdf)[["STRATUM"]]
                              current_pts_spdf <- pts_spdf[!is.na(current_pts_spdf@data[["STRATUM"]]), ]

                              # If there are in fact points in the stratum, do the randomization test
                              if (nrow(current_pts_spdf) > 1) {

                                # Add the coordinates so that we can do nearest neighbor calculations
                                current_pts_spdf <- get_coords(current_pts_spdf)

                                # Derive the arithmetic and geometric mean distance to nearest neighbor for the points
                                nn_means_stratum <- NN_mean(current_pts_spdf@data,
                                                            x_var = "XMETERS",
                                                            y_var = "YMETERS")

                                ## Do randomization test
                                proportions_frame <- test_points(number = reps,
                                                                 pts_spdf = current_pts_spdf,
                                                                 aoi_spdf = stratum_spdf,
                                                                 type = 1,
                                                                 seed_number = seed_number)

                                # Build the output data frame for the sample frame
                                output <- data.frame("polygon" = stratum,
                                                     "point_count" = nrow(current_pts_spdf),
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
  }
  # Combine the outputs
  output <- rbind(output_frame, output_strata)

  return(output)
}



######################################################################
##  This more-or-less dups Random() but is customized to work on buffered NHD polygons for aquatic analyses.
RandomAquatic<-function(WD,		##working directory
                        frame,		## the sample frame as a spdf (shapefile)
                        pts,		## the GRTS points as a spdf  (shapefile)
                        reps,		## Number of sets of random points
                        output,	## Name of file to store results
                        strata,	## Name of strata file (analyses will be conducted at the strata level) or NA
                        stratafield,	## Strata-field name in strata or NA
                        doFrame) 	## Set to T for analysis across the entire FRAME (ignores strata if specified), else F

{
  projection = CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")
  projectionAL <- CRS("+proj=aea")

  setwd(WD)
  numberofreps<-reps

  blank<-"Randomization Results"
  write.table(blank,output,row.names=F,col.names=F,quote=F,append=F)


  ## Ingest all of the points corresponding to the entire Frame.  Even if doFrame==F, we'll use the full suite of pts for the strata-based analysis
  Framepts.spdf<-Ingest(pts)
  a<-GetXY(Framepts.spdf)	## Access X & Y
  a$STRATA<-Framepts.spdf$STRM_ORDR			## Assumes pts have field STRM_ORDR
  Framepts.spdf<-a

  # quickie for experimenting with N
  # Framepts.spdf$DELETE<-0
  #for(i in 201:nrow(Framepts.spdf)) {
  #   Framepts.spdf$DELETE[i]<-1
  #}
  #Framepts.spdf<-Framepts.spdf[Framepts.spdf$DELETE==0 ,]



  if(doFrame==T) {		## If requested, first analyze entire frame
    ## ingest the frame
    frame.spdf<-Ingest(frame)
    frame.spdf$CODE<-1
    frame.spdf<-frame.spdf[, c("CODE")]			## Render Frame to just CODE

    # Derive Mean NN of the Frame points
    FrameNN<-NN(Framepts.spdf)		## Mean NN; [1]=arithmetic, [2]=geometric
    print(paste("Arithmetic Mean nearest neighbor of the input points= ",FrameNN[1]," for N= ",nrow(Framepts.spdf)," points",sep=" "))
    print(paste("Geometric Mean nearest neighbor of the input points= ",FrameNN[2]," for N= ",nrow(Framepts.spdf)," points",sep=" "))

    ## Do randomization test for the entire frame
    number<-reps
    thepts<-Framepts.spdf
    proportion<-GenPts(number,thepts,frame.spdf,FrameNN,2)	## Number of requested sim reps, the pts file, the aoi file (frame or a stratum file),
    ##   the previously calculated vector (arithmetic and geometric) of MNN distances of the input pts

    print(paste("Proportion of random draws with Arithmetic NN >=  NN of the input points = ",proportion[1],sep=" "))
    a<-data.frame(COL="Arithmetic Mean NN(meters) of the ",COL2="",COL3="",N=nrow(Framepts.spdf),COL4="input points = ",VAL=FrameNN[1])
    b<-data.frame(COL="Proportion of",COL2=as.character(numberofreps),COL3=" random draws with an Arithmetic mean NN >= mean NN of the",N=nrow(Framepts.spdf),COL4=" input points = ",VAL=proportion[1])
    ans<-rbind(a,b)

    print(paste("Proportion of random draws with Geometric NN >=  NN of the input points = ",proportion[2],sep=" "))
    a<-data.frame(COL="Geometric Mean NN(meters) of the ",COL2="",COL3="",N=nrow(Framepts.spdf),COL4="input points = ",VAL=FrameNN[2])
    b<-data.frame(COL="Proportion of",COL2=as.character(numberofreps),COL3=" random draws with a Geometric mean NN >= mean NN of the",N=nrow(Framepts.spdf),COL4=" input points = ",VAL=proportion[2])
    ans<-rbind(ans,a)
    ans<-rbind(ans,b)
    write.table(ans,output,row.names=F,col.names=F,quote=F,append=T)

  }else {
    frame.spdf<-NULL
  }



  ######### Analyze by strata if requested
  strata.spdf<-NULL
  # Ingest the strata file if specified, and analyze on a stratum basis
  if(!is.na(strata)) {
    strata.spdf<-Ingest(strata)
  }else {
    q()
  }

  if(!is.null(strata.spdf)) {
    names(strata.spdf)[names(strata.spdf)==stratafield]<-"STRATA"
    strata.spdf<-strata.spdf[, c("STRATA")]
    slist<-unique(strata.spdf$STRATA)

    for(i in 1:length(slist)) {
      tempS<-strata.spdf[strata.spdf$STRATA %in% slist[i] ,]
      tempS$CODE<-1
      tempPTS<-Framepts.spdf[Framepts.spdf$STRATA==slist[i] ,]
      if(nrow(tempPTS)>1) {			## If we have >1 points, do a randomization
        # Derive Mean NN of the Stratum points
        StrataNN<-NN(tempPTS)		## Mean NN; [1]=arithmetic, [2]=geometric
        print(paste("For Stratum = ",slist[i],":",sep=""))
        print(paste("Arithmetic Mean nearest neighbor of the input points= ",StrataNN[1]," for N= ",nrow(tempPTS)," points",sep=" "))
        print(paste("Geometric Mean nearest neighbor of the input points= ",StrataNN[2]," for N= ",nrow(tempPTS)," points",sep=" "))

        ## Do randomization test for this stratum
        number<-reps
        thepts<-tempPTS
        proportion<-GenPts(number,thepts,tempS,StrataNN,2)	## Number of requested sim reps, the pts file, the aoi file (frame or a stratum file),
        ##   the previously calculated vector (arithmetic and geometric) of NN distances of the input pts
        print(paste("Proportion of random draws with Arithmetic NN >=  NN of the input points = ",proportion[1],sep=" "))
        a<-data.frame(COL="Arithmetic Mean NN(meters) of the ",COL2="",COL3="",N=nrow(tempPTS),COL4="input points = ",VAL=StrataNN[1])
        b<-data.frame(COL="Proportion of",COL2=as.character(numberofreps),COL3=" random draws with an Arithmetic mean NN >= mean NN of the",N=nrow(tempPTS),COL4=" input points = ",VAL=proportion[1])
        ans<-rbind(a,b)

        print(paste("Proportion of random draws with Geometric NN >=  NN of the input points = ",proportion[2],sep=" "))
        a<-data.frame(COL="Geometric Mean NN(meters) of the ",COL2="",COL3="",N=nrow(tempPTS),COL4="input points = ",VAL=StrataNN[2])
        b<-data.frame(COL="Proportion of",COL2=as.character(numberofreps),COL3=" random draws with a Geometric mean NN >= mean NN of the",N=nrow(tempPTS),COL4=" input points = ",VAL=proportion[2])
        ans<-rbind(ans,a)
        ans<-rbind(ans,b)
        blank<-c(" ")
        write.table(blank,output,row.names=F,col.names=F,quote=F,append=T)
        blank<-paste("For Stratum= ",slist[i],":",sep="")
        write.table(blank,output,row.names=F,col.names=F,quote=F,append=T)
        write.table(ans,output,row.names=F,col.names=F,quote=F,append=T)
      }else {
        print(paste("For Stratum = ",slist[i],":",sep=""))
        print(paste("WARNING, number of points only = ",nrow(tempPTS),sep=""))
        blank<-" "
        write.table(blank,output,row.names=F,col.names=F,quote=F,append=T)
        blank<-paste("For Stratum = ",slist[i],":",sep="")
        write.table(blank,output,row.names=F,col.names=F,quote=F,append=T)
        blank<-paste("WARNING, number of points only = ",nrow(tempPTS),sep="")
        write.table(blank,output,row.names=F,col.names=F,quote=F,append=T)
      }
    }
  }
} ## End of RandomAquatic()
##  Paste all of the above into the R console.

##################################################################################################################
## II. Modify the following arguments, then call Random for Terrestrial AIM analysis, else RandomAquatic for Aquatic AIM analysis.
##     For Terrestrial AIM analyses, all shapefiles must be dissolved (1 row per feature in the shapefile).
##     For Aquatic AIM analyses, the strata/frame file (1 in the same) should not be dissolved (pertains to the standard way strata are generated in the aquatic world).


WD<-c("c:/projects/mastergrts")	# Working directory - all input files must reside in this directory
frame<-c("Frame")			# Name of the frame file in the working directory - this will be the same as strata for AIM Aquatic analyses
pts<-c("grtssample")			# Name of the GRTS pts file in the working directory.  A field called USE is created in the terrestrial overlay function, so
##    the incoming pts file shouldn't have a field named USE.
reps<-500				# No. of randomly selected sets of points
output<-c("test")		        # Name of file to record the results
set.seed(1)				# Set the random number seed or not!

############## Strata info
strata<-NA
#strata<-c("eagleSTRDISS")		# Set to NA if you don't want to analyze by strata (e.g., strata<-NA) , else enter the name of the strata shapefile.
#stratafield<-c("GRIDCODE")		# If strata is set, then specify the strata field in the strata shapefile (upper case), else stratafield is ignored.
# stratafield is the strata attribute within file=strata.
## NOTE:  The strata field is renamed STRATA in prcessing, so the strata file can't already have this field name.
## For Terrestrial AIM processing, the over() function is used to derive the 'actual' strata of points; thus,
## any strata designation within the pts file is not used!  THUS, the strata file used to generate in input GRTS pts should be same as used here.

doFrame<-T				# Set to T if you want to analyze across the enter Frame without regard to strata, else set to F.  Even if you set strata, setting this to T
#     will generate results where strata are ignored; in essence you do 2 analyses at once - frame-based and strata-based.  If strata  == NA,
#     make sure this is set to T, otherwise no analyses are performed!

## Nothing is returned from Random() - the ASCII-formatted results are output during processing to the file specified by output (set above)
Random(WD,frame,pts,reps,output,strata,stratafield,doFrame)		## Call this if doing Terrestrial AIM
#RandomAquatic(WD,frame,pts,reps,output,strata,stratafield,doFrame)	## Call this if doing Aquatic AIM
q()									## Use this to terminate when running a batch file

