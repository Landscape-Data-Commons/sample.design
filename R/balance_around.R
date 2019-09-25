#' Determine the preference of two sets of points for each other based on distance
#' @description Given two spatial points data frames, compare the members of each against all the members of the other, ranking their "preference" for the compared points based on distance. By default this assigns higher preference for points that are closer.
#' @param template_points Spatial points data frame. The first set of points to be used in the comparison.
#' @param comparison_points Spatial points data frame. The second set of points to be used in the comparison.
#' @return A named list of two data frames: \code{"template"} and \code{"comparison"}. Both data frames contain \code{template_index} (the indices of points in \code{template_points}) and \code{comparison_index} (the indices of points in \code{template_points}). Both data frames contain all unique combinations of the two sets of points. \code{"template"} contains a variable named \code{"rank_by_template"} which is the relative preference of the template point at that index for the comparison point at that index (lower values representing higher preference). \code{"comparison"} likewise has \code{"rank_by_comparison"} which is the preference of the comparison point for the template point.
#' @export
find_preferences <- function(template_points,
                             comparison_points){
  if (!(class(template_points) %in% "SpatialPointsDataFrame")) {
    stop("template_points must be a data frame.")
  }
  if (!(class(comparison_points) %in% "SpatialPointsDataFrame")) {
    stop("comparison_points must be a data frame.")
  }

  n_template <- nrow(template_points@data)
  n_comparison <- nrow(comparison_points@data)

  if (n_template < 1) {
    stop("There are no points in template_points")
  }
  if (n_comparison < 1) {
    stop("There are no points in comparison_points")
  }

  if (!identical(template_points@proj4string, comparison_points@proj4string)) {
    warning("comparison_points and template_points aren't in the same projection. Reprojecting comparison_points to match.")
    comparison_points <- sp::spTransform(comparison_points,
                                         CRSobj = template_points@proj4string)
  }

  template_points@data[["source"]] <- "template"
  template_points <- get_coords(points = template_points,
                                x_var = "XMETERS",
                                y_var = "YMETERS",
                                projection = sp::CRS("+proj=aea"))
  comparison_points@data[["source"]] <- "comparison"
  comparison_points <- get_coords(points = comparison_points,
                                  x_var = "XMETERS",
                                  y_var = "YMETERS",
                                  projection = sp::CRS("+proj=aea"))

  distance_matrix <- dist_matrix(dataframe = rbind(template_points@data,
                                      comparison_points@data),
                     x_var = "XMETERS",
                     y_var = "YMETERS")

  template_df <- do.call(rbind,
                         lapply(X = 1:n_template,
                                distance_matrix = distance_matrix,
                                comparison_indices = (n_template + 1):(n_template + n_comparison),
                                n_comparison = n_comparison,
                                FUN = function(X, distance_matrix, comparison_indices, n_comparison){
                                  distances <- distance_matrix[comparison_indices, X]
                                  output <- data.frame(template_index = X,
                                                       comparison_index = 1:n_comparison,
                                                       rank_by_template = order(distance_matrix[comparison_indices, X]))
                                  return(output)
                                }))

  comparison_df <- do.call(rbind,
                           lapply(X = (n_template + 1):(n_template + n_comparison),
                                  distance_matrix = distance_matrix,
                                  template_indices = 1:n_template,
                                  n_template = n_template,
                                  FUN = function(X, distance_matrix, template_indices, n_template){
                                    distances <- distance_matrix[template_indices, X]
                                    output <- data.frame(comparison_index = X - n_template,
                                                         template_index = 1:n_template,
                                                         rank_by_comparison = order(distance_matrix[template_indices, X]))
                                    return(output)
                                  }))

  output <- list(template = template_df,
                 comparison = comparison_df)
  return(output)
}

#' Find optimal pairings from two sets based on their rankings
#' @description For two sets of unique identifiers, given each member's ranked preference for every member of the other set, find optimal pairings. The two sets are referred to as \code{match_to} and \code{match_from}. All the members of \code{match_to} will be paired to members of \code{match_from}, but if there are more members of \code{match_from} than \code{match_to} there will be unpaired, discarded members of \code{match_from}. The inputs are two data frames, one where each observation is a unique combination of each member of \code{match_to} and \code{match_from} and the preference of that \code{match_to} member for \code{match_from}, numeric with lower values treated as more preferred. The other data frame contains the preferences of the \code{match_from} members for \code{match_to} members.
#' @param match_to Data frame. The preferences of the members of the set to match to for the members of the set to match from. There must be one observation for each unique pairing between the two sets. The variable \code{match_to_idvar} must contain the identities of the members of the set to match to; \code{match_from_idvar}, the set to match from; and \code{match_to_rankvar} the rank preference of the MATCH FROM member for the MATCH FROM member.
#' @param match_from Data frame. The preferences of the members of the set to match from for the members of the set to match to. There must be one observation for each unique pairing between the two sets. The variable \code{match_to_idvar} must contain the identities of the members of the set to match to; \code{match_from_idvar}, the set to match from; and \code{match_from_rankvar} the rank preference of the MATCH FROM member for the MATCH TO member.
#' @param match_to_idvar Character string. The name of the variable in both \code{match_to} and \code{match_from} that contains the identities of the members in the set to match to. Defaults to \code{"match_to_id"}.
#' @param match_from_idvar Character string. The name of the variable in both \code{match_to} and \code{match_from} that contains the identities of the members in the set to match from. Defaults to \code{"match_from_id"}.
#' @param match_to_rankvar Character string. The name of the variable in \code{match_to} that contains the rank preference of the match to members for the match from members as a relative numeric value (lower values are more preferred). Defaults to \code{"match_to_rank"}.
#' @param match_from_rankvar Character string. The name of the variable in \code{match_from} that contains the rank preference of the match from members for the match to members as a relative numeric value (lower values are more preferred). Defaults to \code{"match_from_rank"}.
#' @return A data frame with the variables \code{match_to_id} and \code{match_from_id} containing the paired members of the two sets.
#' @export
ranked_sort <- function(match_to,
                        match_from,
                        match_to_idvar = "match_to_id",
                        match_from_idvar = "match_from_id",
                        match_to_rankvar = "match_to_rank",
                        match_from_rankvar = "match_from_rank",
                        iteration_limit = NULL){
  if (class(match_to) != "data.frame") {
    stop("match_to must be a data frame.")
  }
  match_to_names <- names(match_to)
  if (class(match_from) != "data.frame") {
    stop("match_from must be a data frame.")
  }
  match_from_names <- names(match_from)

  if (!all(match_to_idvar %in% match_to_names, match_to_idvar %in% match_from_names)) {
    stop("The match_to_idvar, ", match_to_idvar, ", must appear in both match_to and match_from.")
  }
  if (!all(match_from_idvar %in% match_to_names, match_from_idvar %in% match_from_names)) {
    stop("The match_from_idvar, ", match_from_idvar, ", must appear in both match_to and match_from.")
  }
  if (!(match_to_rankvar %in% match_to_names)) {
    stop("The match_to_rankvar, ", match_to_rankvar, ", must appear in match_to.")
  }
  if (!(match_from_rankvar %in% match_from_names)) {
    stop("The match_from_rankvar, ", match_from_rankvar, ", must appear in match_from.")
  }

  # Add the assigned variable, but set them all to FALSE for now because nothing's been assigned
  match_to[["assigned"]] <- FALSE
  match_from[["assigned"]] <- FALSE

  # How many to match from?
  n_matchfrom <- length(unique(match_from[[match_from_idvar]]))

  # Make a tracker for whether the existing points have been confirmed
  matchfrom_confirmed <- data.frame(matchfrom_id = unique(match_from[[match_from_idvar]]),
                                    confirmed = FALSE)

  # Have all the points been confirmed to be checked and assigned or rejected?
  all_confirmed <- FALSE
  if (is.null(iteration_limit)) {
    iteration_limit <- max(n_matchfrom, n_matchto)
  } else {
    if (!is.numeric(iteration_limit)) {
      stop("The maximum number of iterations given as iteration_limit must be numeric")
    }
    if (iteration_limit < 0) {
      stop("The maximum number of iterations given as iteration_limit must be positive")
    }
  }

  # Limit how many loops can be made. This counter gets compared against iteration_limit
  iterations <- 0

  # If we have unconfirmed points then let's work on assigning them
  # This may take a couple of passes
  while (!all_confirmed & iterations < iteration_limit) {
    # Get the still-unconfirmed existing points' indices
    unconfirmed_existing_ids <- matchfrom_confirmed[!matchfrom_confirmed[["confirmed"]], "matchfrom_id"]

    # Now work through those points!
    for (current_id in unconfirmed_existing_ids) {
      # This is used to flag that it's still unconfirmed
      # It'll flip to TRUE if it gets assigned to or fails to qualify for assignment to any match_to point
      current_id_unconfirmed <- TRUE

      # Which rows in match_from correspond to this current point?
      match_from_current_vector <- match_from[[match_from_idvar]] == current_id
      # And in match_to?
      match_to_current_vector <- match_to[[match_from_idvar]] == current_id

      # Now to get the order of indices to work through to try to assign the point to
      match_to_info <- match_from[match_from_current_vector, c(match_from_rankvar, match_to_idvar)]

      # Work through the comparisons in rank order from most preferred match_to ID to least preferred.
      # As soon as it qualifies, this while() will trip and it'll stop looking
      # It'll also trip if it tries every ID and can't find anywhere for it
      while (current_id_unconfirmed) {
        for (rank_by_current in sort(match_to_info[[match_from_rankvar]])) {
          # Get the match_to ID
          current_match_to_id <- match_to_info[match_to_info[[match_from_rankvar]] == rank_by_current, match_to_idvar]

          # Which rows in match_to correspond to this match_to ID?
          match_to_match_to_vector <- match_to[[match_to_idvar]] == current_match_to_id
          # And in match_from?
          match_from_match_to_vector <- match_from[[match_to_idvar]] == current_match_to_id

          # This is the row where the point and match_to point meet in match_to!
          current_row_match_to <- match_to_current_vector & match_to_match_to_vector
          # And in match_from
          current_row_match_from <- match_from_current_vector & match_from_match_to_vector

          # Is the match_to point paired off yet?
          # We'll do this by treating T/F as 1/0 and finding the sum so that we can confirm that only one is TRUE
          match_to_assigned_count <- sum(match_to[match_to_match_to_vector, "assigned"])

          # If it's already has another point assigned
          # There are also two separate checks to prevent multiple assignments from happening
          if (match_to_assigned_count == 1 & current_id_unconfirmed & !any(match_to[match_to[[match_to_idvar]] == current_id, "assigned"])) {
            # We're going to compare the preference of the MATCH_TO ID for the current and assigned MATCH_FROM IDS

            # This is the rank for the index of the unconfirmed point
            current_rank_by_match_to <- match_to[current_row_match_to, match_to_rankvar]

            # This is the rank of the currently assigned point
            comparison_row_match_to <- match_to[["assigned"]] & match_to_match_to_vector
            comparison_rank_by_match_to <- match_to[comparison_row_match_to, match_to_rankvar]

            # And the index for it in existing points
            comparison_id <- match_to[comparison_row_match_to, match_from_idvar]

            # Which rows in match_from correspond to this comparison point?
            match_from_comparison_vector <- match_from[[match_from_idvar]] == comparison_id
            # And in match_to?
            match_to_comparison_vector <- match_to[[match_from_idvar]] == comparison_id

            # This is the row where the comparison point and match_to point meet in match_to!
            comparison_row_match_to <- match_to_comparison_vector & match_to_match_to_vector
            # And in match_from
            comparison_row_match_from <- match_from_comparison_vector & match_from_match_to_vector

            # Compare the two
            # If the current point is more highly ranked, assign it and unassign the point that's already attached
            if (current_rank_by_match_to < comparison_rank_by_match_to) {
              # For the already assigned point, set values to FALSE in both match_from and match_to
              matchfrom_confirmed[matchfrom_confirmed[["matchfrom_id"]] == comparison_id, "confirmed"] <- FALSE
              match_from[comparison_row_match_from, "assigned"] <- FALSE
              match_to[comparison_row_match_to, "assigned"] <- FALSE

              # For the current point that's being assigned, set those values to TRUE
              matchfrom_confirmed[matchfrom_confirmed[["matchfrom_id"]] == current_id, "confirmed"] <- TRUE
              match_from[current_row_match_from, "assigned"] <- TRUE
              match_to[current_row_match_to, "assigned"] <- TRUE

              # It's not unconfirmed anymore
              current_id_unconfirmed <- FALSE
            }
          } else if (match_to_assigned_count == 0 & current_id_unconfirmed & !any(match_to[match_to[[match_to_idvar]] == current_id, "assigned"])) {
            # There are TWO separate checks to prevent multiple assignments from happening
            # If there isn't a point assigned yet, this point gets assigned
            matchfrom_confirmed[matchfrom_confirmed[["matchfrom_id"]] == current_id, "confirmed"] <- TRUE
            match_from[current_row_match_from, "assigned"] <- TRUE
            match_to[current_row_match_to, "assigned"] <- TRUE

            # It's not unconfirmed anymore
            current_id_unconfirmed <- FALSE
          }
        }

        # If it doesn't qualify for any, it gets confirmed but not assigned
        matchfrom_confirmed[matchfrom_confirmed[["matchfrom_id"]] == current_id, "confirmed"] <- TRUE

        # It's not unconfirmed anymore
        current_id_unconfirmed <- FALSE
      }
    }

    # Update whether every point is confirmed!
    all_confirmed <- all(matchfrom_confirmed[["confirmed"]])

    # Increment the iteration counter
    iterations <- iterations + 1
  }

  # If we hit the iteration limit, something is wrong!!!!!!!!!!!
  if (iterations >= iteration_limit) {
    message("The iteration limit of ", iteration_limit, " has been reached without reaching an optimal solution. Check your rankings or use another method.")
  }

  # Time to return just the pairings!
  output <- match_to[match_to[["assigned"]], c(match_to_idvar, match_from_idvar)]

  return(output)
}

#' Select points that most closely approximate the distribution of another set of points
#' @description Sometimes you have a large collection of points which are not randomly distributed or spatially balanced and you would like a subset that more or less do. Given a template of points that are distributed the way you would like, this will return the closest existing point to each. This can be done taking into account membership in a group, either by having assigned it as a variable in both sets of points or by providing polygons that can be used to assign membership. By default, no stratification/membership is taken into account. NOTA BENE: As currently implemented, the evaluation to select points occurs in the order of the points in \code{template_points} which means that the solution found may not be the optimal one, just the best for that order. If that's important to you, the "easiest" workaround is to run this several times with randomized \code{template_points} orders and to keep the result with the smallest mean distance between result points and the template points.
#' @param existing_points Spatial Points Data Frame. The points you would like to select from by comparing against \code{template_points}.
#' @param template_points  Spatial Points Data Frame. The points you would like to compare against \code{existing_points} in order to select a subset of those that most closely resemble the distribution of the template points.
#' @param strata_polygons Spatial Polygons Data Frame. Polygons assigned a variable with a name \code{stratafield} that contains the membership information (e.g. strata) to assign to \code{existing_points} and \code{template_points}. If \code{NULL} then no assignment will be attempted. Defaults to \code{NULL}.
#' @param stratafield Character string. If \code{strata_polygons} is not \code{NULL}, the name of the variable in \code{strata_polygons@@data} that contains the membership information. Otherwise, the name of the variable in both \code{template_points@@data} and \code{existing_points@@data} that contains the membership information. If \code{NULL} then the points will be considered to belong to a single group. Defaults to \code{NULL}.
#' @param projection CRS object. The projection to force all spatial objects into. Defaults to NAD83, \code{sp::CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")}.
#' @return A spatial polygons data frame made by trimming \code{existing_points} down to the points that most closely approximate the distribution of \code{template_points} while also containing the same number of points as \code{template_points}. It will be in the projection specified by \code{projection}.
#' @export
get_closest <- function(existing_points,
                        template_points,
                        strata_polygons = NULL,
                        stratafield = NULL,
                        projection = sp::CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")){


  if (!(class(existing_points) %in% "SpatialPointsDataFrame")) {
    stop("existing_points must be a spatial points data frame")
  }

  if (!(class(template_points) %in% "SpatialPointsDataFrame")) {
    stop("template_points must be a spatial points data frame")
  }
  if (nrow(template_points@data) < 1) {
    stop("There are no points in template_points.")
  }
  if (nrow(existing_points@data) < 1) {
    stop("There are no points in existing_points.")
  }

  # Reproject if necessary
  if (!identical(projection, template_points@proj4string)) {
    template_points <- sp::spTransform(template_points,
                                       projection)
  }
  if (!identical(projection, existing_points@proj4string)) {
    existing_points <- sp::spTransform(existing_points,
                                       projection)
  }


  # What to do about stratafield
  if (!is.null(stratafield)) {
    if (!is.character(stratafield)) {
      stop("stratafield must be a character string")
    }
    if (length(stratafield) > 1) {
      stop("stratafield must be a character string")
    }
    if (is.null(strata_polygons)) {
      if (!(stratafield %in% names(template_points@data))) {
        stop("The variable ", stratafield, " does not appear in template_points@data.")
      }
      if (!(stratafield %in% names(existing_points@data))) {
        stop("The variable ", stratafield, " does not appear in existing_points@data.")
      }
      # Sorry for the inconsistencies but "MEMBERSHIP" is trying to get away from the assumption that these will be stratifications
      # And I'm so, so tired. I'll come back to clean this up later, or so I'm currently telling myself at the end of a ten hour day
      template_points@data[["MEMBERSHIP"]] <- template_points@data[[stratafield]]
      existing_points@data[["MEMBERSHIP"]] <- existing_points@data[[stratafield]]
    }

  } else {
    template_points@data[["MEMBERSHIP"]] <- "frame"
    existing_points@data[["MEMBERSHIP"]] <- "frame"
  }

  # If there is a set of stratification polygons, use them!
  if (!is.null(strata_polygons)) {
    if (!(class(strata_polygons) %in% "SpatialPolygonsDataFrame")) {
      stop("strata_polygons must be a spatial polygons data frame.")
    }
    if (!identical(projection, strata_polygons@proj4string)) {
      strata_polygons <- sp::spTransform(strata_polygons,
                                     projection)
    }
    if (!(stratafield %in% names(strata_polygons@data))) {
      stop("The variable ", stratafield, " does not appear in strata_polygons.")
    }

    # This just puts the strata into the points
    existing_points@data[["MEMBERSHIP"]] <- sp::over(x = existing_points,
                                                     y = strata_polygons)[[stratafield]]
    template_points@data[["MEMBERSHIP"]] <- sp::over(x = template_points,
                                                     y = strata_polygons)[[stratafield]]
  }


  # What strata are there?
  strata <- unique(c(existing_points@data[["MEMBERSHIP"]], template_points@data[["MEMBERSHIP"]]))

  # By stratum!
  stata_selections <- lapply(X = strata,
                             template_points = template_points,
                             existing_points = existing_points,
                             FUN = function(X, template_points, existing_points){
                               # Narrow it down to the points in the current stratum
                               stratum <- X
                               existing_points_stratum <- existing_points[existing_points@data[["MEMBERSHIP"]] == stratum, ]
                               template_points_stratum <- template_points[template_points@data[["MEMBERSHIP"]] == stratum, ]

                               # What are their ranking of each other between template and comparison based on distance?
                               preferences <- find_preferences(template_points = template_points_stratum,
                                                               comparison_points = existing_points_stratum)

                               # What's the optimal solution for minimizing distances for pairing?
                               sorting <- ranked_sort(match_to = preferences[["template"]],
                                                      match_from = preferences[["comparison"]],
                                                      match_to_idvar = "template_index",
                                                      match_from_idvar = "comparison_index",
                                                      match_to_rankvar = "rank_by_template",
                                                      match_from_rankvar = "rank_by_comparison")

                               # Only keep the points that got paired with the template points
                               output <- existing_points_stratum[sorting[["comparison_index"]], ]

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

  if (!(class(existing_points) %in% "SpatialPointsDataFrame")) {
    stop("existing_points must be a spatial points data frame")
  }

  if (nrow(existing_points@data) < 1) {
    stop("There are no points in existing_points")
  }
  if (!(class(new_points) %in% "SpatialPointsDataFrame")) {
    stop("new_points must be a spatial points data frame")
  }
  if (nrow(new_points@data) < 1) {
    stop("There are no points in new_points")
  }

  if (!identical(projection, existing_points@proj4string)) {
    existing_points <- sp::spTransform(existing_points,
                                       projection)
  }
  if (!identical(projection, new_points@proj4string)) {
    new_points <- sp::spTransform(new_points,
                                  projection)
  }

  if (is.null(target)) {
    target <- nrow(new_points@data)
  }

  if (target <= nrow(existing_points@data)) {
    stop("The target number of points is less than or equal to the number of existing points.")
  }

  common_varnames <- unique(c(names(existing_points@data)[(names(existing_points@data) %in% names(new_points@data))],
                              names(new_points@data)[(names(new_points@data) %in% names(existing_points@data))]))

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
  n_removal_indices <- 0
  # Here's one we can mutilate over our iterations (always avoid violence to original data you may reference again!)
  working_distance_df <- distance_df
  # We'll be ignoring the index column because of course sequential ordinals starting at 1 are the smallest numbers in the matrix!
  # The indices screw up min() results without removing this for the finding minimum evaluation step
  index_colnum <- grep(names(distance_df), pattern = "INDEX")

  while (n_removal_indices < (n_existing + n_new - target)) {
    current_min <- min(working_distance_df[, -index_colnum])
    # Get the indices from every column where that min occurs
    # Each column is an existing point, so if we check every column for the value and store that index,
    # those are the new points that are that distance from an existing point
    current_indices <- sapply(X = 1:(ncol(working_distance_df) - 1),
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
    # Drop any NULLs
    removal_indices <- removal_indices[!sapply(removal_indices, is.null)]

    if (!is.null(removal_indices)) {
      # Make sure we don't overshoot our removal goal
      removal_indices <- removal_indices[1:min(length(removal_indices), (target - n_existing))]

      # Remove them!
      working_distance_df <- working_distance_df[!(working_distance_df[["INDEX"]] %in% removal_indices), ]

      # Update our tracking value
      n_removal_indices <- length(removal_indices)
    }
  }


  # removal_indices <- removal_indices[!sapply(removal_indices, is.null)]

  # Now that we have our indices to remove, let's do it as we combine points
  # The as.numeric() is because due to the NULL that's in there from the pre-loop setup removal_indices is a list, not a vector
  output <- rbind(existing_points[, common_varnames],
                  new_points[-as.numeric(removal_indices), common_varnames])

  return(output)
}


#' Combine existing and new points to create a spatially balanced design
#' @param existing_points Spatial points data frame. The existing points that will be balanced around.
#' @param new_points Spatial points data frame. The points that will be compared against the existing points and selected from to create a balanced design.
#' @param stratafield Character string. The name of the variable in common between \code{existing_points} and \code{new_points} that contains stratum identities. This is used to balance by stratum. If \code{NULL} then balancing will not take strata into account. Defaults to \code{NULL}.
#' @param projection CRS object. The projection to force on the spatial objects. Defaults to \code{sp::CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")}.
#' @return A spatial points data frame containing all the points from \code{existing_points} and the selected points from \code{new_points}. The projection will match \code{projection}.
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
                    pattern = "\\d*$",
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


#' Select points that most closely approximate the distribution of another set of points and substitute them
#' @description When creating a design that incorporates points from another design to revisit, it can be important to approximate the same spatial distribution, especially if the design is spatially balanced. This will substitute in revisit points that will cause the least shift in the spatial distribution.
#' @param existing_points Spatial Points Data Frame. The points you would like to select from by comparing against \code{template_points}.
#' @param template_points  Spatial Points Data Frame. The points you would like to compare against \code{existing_points} in order to select a subset of those that most closely resemble the distribution of the template points.
#' @param existing_idvar Character string. The name of the variable in \code{existing_points@@data} that contains the unique identifiers. Defaults to \code{"plotid"}.
#' @param template_idvar Character string. The name of the variable in \code{template_points@@data} that contains the unique identifiers. Defaults to \code{"plotid"}.
#' @param replacement_count Optional integer. The number of template points to replace with the closest approximation from the existing points. At least one point in each stratum will be replaced and the rest of the remain replacements will be made proportionally to the weight of the strata (e.g. a stratum that makes up 1/2 of the area or contains 1/2 of the template points will also contain 1/2 of the replacements). If \code{NULL} then points will be allocated according to the vector \code{keep_counts}. Defaults to \code{NULL}.
#' @param keep_counts Optional named numeric vector. A vector of the number of template points to replace with the closest approximation from the existing points in each stratum. Only use if there is more than one stratum. The names of the vector must exactly match the strata of the points. If \code{NULL} then points will be allocated proportionally by area using the value in \code{replacement_count}. Defaults to \code{NULL}.
#' @param strata_polygons Optional Spatial Polygons Data Frame. Polygons assigned a variable with a name \code{polygons_stratavar} that contains the membership information (e.g. strata) to assign to \code{existing_points} and \code{template_points}. If \code{NULL} then no assignment will be attempted. Defaults to \code{NULL}.
#' @param polygons_stratavar Optional character string. If \code{strata_polygons} is not \code{NULL}, the name of the variable in \code{strata_polygons@@data} that contains the membership information. Otherwise, the name of the variable in both \code{template_points@@data} and \code{existing_points@@data} that contains the membership information. If \code{NULL} then the points will be considered to belong to a single group. Defaults to \code{"stratum"}.
#' @param existing_stratavar Optional character string. The name of the variable in \code{existing_points@@data} that contains the membership information. If both \code{existing_stratavar} and \code{template_stratavar} are \code{NULL} then all the points in both will be considered to belong to a single group. Defaults to \code{"stratum"}.
#' @param template_stratavar Optional character string. The name of the variable in \code{template_points@@data} that contains the membership information. If both \code{existing_stratavar} and \code{template_stratavar} are \code{NULL} then all the points in both will be considered to belong to a single group. Defaults to \code{"stratum"}.
#' @param projection CRS object. The projection to force all spatial objects into. If \code{NULL} then the projection from \code{existing_points@@proj4string} will be used. Defaults to NAD83, \code{sp::CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")}.
#' @return A spatial points data frame of the \code{template_points} with points from \code{existing_points} substituted in. It will be in the projection specified by \code{projection}.
#' @export
revisit_design <- function(existing_points,
                           template_points,
                           existing_idvar = "plotid",
                           template_idvar = "plotid",
                           replacement_count = NULL,
                           keep_counts = NULL,
                           strata_polygons = NULL,
                           polygons_stratavar = "stratum",
                           existing_stratavar = "stratum",
                           template_stratavar = "stratum",
                           projection = sp::CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")){
  # SO MUCH SANITIZATION
  if (is.null(replacement_count) & is.null(keep_counts)) {
    stop("You must supply either a number of points to replace (allocated to strata by proportional area) or a named vector of counts per stratum as keep_counts")
  }

  if (!(class(existing_points) %in% "SpatialPointsDataFrame")) {
    stop("existing_points must be a spatial points data frame")
  }
  if (!(existing_idvar %in% names(existing_points@data))) {
    stop("The variable ", existing_idvar, " does not appear in existing_points@data")
  } else {
    if (length(unique(existing_points@data[[existing_idvar]])) != nrow(existing_points@data)) {
      stop("The variable ", existing_idvar, " does not contain unique identifiers in existing_points")
    }
  }

  if (!(class(template_points) %in% "SpatialPointsDataFrame")) {
    stop("template_points must be a spatial points data frame")
  }
  if (!(template_idvar %in% names(template_points@data))) {
    stop("The variable ", template_idvar, " does not appear in template_points@data")
  } else {
    if (length(unique(template_points@data[[template_idvar]])) != nrow(template_points@data)) {
      stop("The variable ", template_idvar, " does not contain unique identifiers in template_points")
    }
  }

  if (is.null(projection)) {
    projection <- existing_points@proj4string
  }

  if (!identical(existing_points@proj4string, projection)) {
    existing_points <- sp::spTransform(existing_points, projection)
  }
  if (!identical(template_points@proj4string, projection)) {
    template_points <- sp::spTransform(template_points, projection)
  }

  if (is.null(strata_polygons)) {
    if (xor(is.null(existing_stratavar), is.null(template_stratavar))) {
      stop("You must either provide the stratification variables for both existing_points and template_points (for a stratified approach without strata polygons) or neither (for a non-stratified approach).")
    }
    if (is.null(existing_stratavar) & is.null(template_stratavar)) {
      existing_points@data[["MEMBERSHIP"]] <- "frame"
      template_points@data[["MEMBERSHIP"]] <- "frame"
      existing_stratavar <- template_stratavar <- "MEMBERSHIP"
    } else {
      if (!(existing_stratavar %in% names(existing_points@data))) {
        stop("The variable", existing_stratavar, "does not occur in existing_points@data")
      }
      if (!(template_stratavar %in% names(template_points@data))) {
        stop("The variable", template_stratavar, "does not occur in template_points@data")
      }
      existing_points@data[["MEMBERSHIP"]] <- existing_points@data[[existing_stratavar]]
      template_points@data[["MEMBERSHIP"]] <- template_points@data[[template_stratavar]]
    }
  } else {
    if (!(class(strata_polygons) %in% "SpatialPolygonsDataFrame")) {
      stop("strata must be a spatial polygons data frame")
    }
    if (!(polygons_stratavar %in% names(strata_polygons@data))) {
      stop("The variable ", polygons_stratavar, " does not occur in strata_polygons@data")
    }
    if (!identical(strata_polygons@proj4string, projection)) {
      strata_polygons <- sp::spTransform(strata_polygons, projection)
    }

    if (is.null(existing_stratavar) & is.null(template_stratavar)) {
      warning("Ignoring existing_stratavar and template_stratavar in favor of strata_polygons")
    }

    existing_points[["MEMBERSHIP"]] <- sp::over(x = existing_points,
                                                y = strata_polygons)[[polygons_stratavar]]
    template_points[["MEMBERSHIP"]] <- sp::over(x = template_points,
                                                y = strata_polygons)[[polygons_stratavar]]
  }

  if (any(is.na(existing_points@data[["MEMBERSHIP"]])) | any(is.na(template_points@data[["MEMBERSHIP"]]))) {
    warning("Excluding points that have no stratification membership")
    existing_points <- existing_points[is.na(existing_points[["MEMBERSHIP"]]), ]
    template_points <- template_points[is.na(template_points[["MEMBERSHIP"]]), ]
  }

  # What strata are there?
  strata <- unique(c(existing_points@data[["MEMBERSHIP"]], template_points@data[["MEMBERSHIP"]]))

  if (!is.null(keep_counts)) {
    if (!all(names(keep_counts) %in% strata)) {
      stop("Not all strata named in keep_counts are associated with the points")
    }
    if (!all(strata %in% names(keep_counts))) {
      stop("Not all strata associated with the points are found in keep_counts")
    }
    if (!is.null(replacement_count)) {
      warning("Provided both replacement_count and keep_counts. Using keep_counts")
      replacement_count <- NULL
    }
  }


  # Make keep_counts, a named vector of how many points per stratum to pull
  # This is based on proportional allocation
  if (!is.null(replacement_count)) {
    allocation <- allocate_panels(spdf = template_points,
                                  stratum_field = "MEMBERSHIP",
                                  panel_names = c("revisit"),
                                  panel_sample_size = replacement_count,
                                  points_min = 1,
                                  oversample_proportion = 0,
                                  oversample_min = 0)

    keep_counts <- sapply(X = strata,
                          allocation = allocation,
                          FUN = function(X, allocation){
                            allocation[[X]][["panel"]]["revisit"]
                          })

    names(keep_counts) <- strata
  }


  # By stratum!
  strata_selections <- lapply(X = strata,
                              template_points = template_points,
                              existing_points = existing_points,
                              keep_counts = keep_counts,
                              FUN = function(X, template_points, existing_points, keep_counts){
                                # Narrow it down to the points in the current stratum
                                stratum <- X
                                existing_points_stratum <- existing_points[existing_points@data[["MEMBERSHIP"]] == stratum, ]
                                template_points_stratum <- template_points[template_points@data[["MEMBERSHIP"]] == stratum, ]

                                # What are their ranking of each other between template and comparison based on distance?
                                preferences <- find_preferences(template_points = template_points_stratum,
                                                                comparison_points = existing_points_stratum)

                                # What's the optimal solution for minimizing distances for pairing?
                                sorting <- ranked_sort(match_to = preferences[["template"]],
                                                       match_from = preferences[["comparison"]],
                                                       match_to_idvar = "template_index",
                                                       match_from_idvar = "comparison_index",
                                                       match_to_rankvar = "rank_by_template",
                                                       match_from_rankvar = "rank_by_comparison")

                                # Get a data frame of the paired IDs
                                pairs <- data.frame(existing_id = existing_points_stratum[sorting[["comparison_index"]], existing_idvar],
                                                    template_id = template_points_stratum[sorting[["template_index"]], template_idvar],
                                                    stringsAsFactors = FALSE)

                                # Just some quick renaming
                                names(pairs) <- gsub(names(pairs),
                                                     pattern = "_id\\.",
                                                     replacement = "_")
                                names(pairs) <- gsub(names(pairs),
                                                     pattern = "_coords\\.x1$",
                                                     replacement = "_xcoord")
                                names(pairs) <- gsub(names(pairs),
                                                     pattern = "_coords\\.x2$",
                                                     replacement = "_ycoord")

                                # And remove optional variables
                                pairs <- pairs[, !grepl(names(pairs), pattern = "_optional$")]

                                # Now we'll keep the n closest where n is the value in keep_counts for this stratum
                                # First up is the distance between these pairs (apply() was acting up, so sapply() it is)
                                pairs[["distance"]] <- sapply(X = 1:nrow(pairs),
                                                              df = pairs,
                                                              FUN = function(X, df){
                                                                (df[X, "existing_xcoord"] - df[X, "template_xcoord"])^2 + (df[X, "existing_ycoord"] - df[X, "template_ycoord"])^2
                                                              })

                                pairs <- pairs[order(pairs[["distance"]]), ]

                                n_keep <- keep_counts[[stratum]]

                                if (n_keep > nrow(pairs)) {
                                  stop("Attempting to keep ", n_keep, " points in ", stratum, " but only ", nrow(pairs), " plot pairs are available.")
                                }

                                return(pairs[1:n_keep, ])
                              })

  strata_selection_df <- do.call(rbind, strata_selections)

  # This part builds a new spatial points data frame that subs in the existing points that were identified
  # IT ASSUMES THAT THE TEMPLATE POINTS ARE IN ORDER!!!!!!!!!!!!!!!
  template_df <- cbind(template_points@data, template_points@coords)
  template_df[["order_number"]] <- 1:nrow(template_df)
  template_df <- template_df[, c("order_number", template_idvar, "xcoord", "ycoord")]
  names(template_df)[names(template_df) == template_idvar] <- "template_plotid"

  merged_points <- merge(x = template_df,
                         y = strata_selection_df[, c("existing_plotid", "template_plotid", "existing_xcoord", "existing_ycoord")],
                         by.x = "template_plotid",
                         by.y = "template_plotid",
                         all.x = TRUE)

  output_df <- merged_points
  output_df[["plotid"]] <- output_df[["template_plotid"]]
  output_df[["revisit"]] <- !is.na(output_df[["existing_plotid"]])
  output_df[output_df[["revisit"]], c("plotid", "xcoord", "ycoord")] <- output_df[output_df[["revisit"]], c("existing_plotid", "existing_xcoord", "existing_ycoord")]

  output_spdf <- sp::SpatialPointsDataFrame(coords = output_df[, c("xcoord", "ycoord")],
                                            data = output_df[, c("order_number", "plotid", "revisit")],
                                            proj4string = projection)

  return(output_spdf)
}
