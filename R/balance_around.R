#' Determine the preference of two sets of points for each other based on distance
#' @description Given two sf point objects, compare the members of each against all the members of the other, ranking their "preference" for the compared points based on distance. By default this assigns higher preference for points that are closer.
#' @param template_points An sf point object. The first set of points to be used in the comparison.
#' @param comparison_points An sf point object. The second set of points to be used in the comparison.
#' @return A named list of two data frames: \code{"template"} and \code{"comparison"}. Both data frames contain \code{template_index} (the indices of points in \code{template_points}) and \code{comparison_index} (the indices of points in \code{template_points}). Both data frames contain all unique combinations of the two sets of points. \code{"template"} contains a variable named \code{"rank_by_template"} which is the relative preference of the template point at that index for the comparison point at that index (lower values representing higher preference). \code{"comparison"} likewise has \code{"rank_by_comparison"} which is the preference of the comparison point for the template point.
#' @export
find_preferences <- function(template_points,
                             comparison_points){
  if (!("sf" %in% class(template_points))) {
    stop("template_points must be an sf point object.")
  }
  if (!all(sf::st_geometry_type(template_points) %in% c("POINT"))) {
    stop("template_points must be an sf point object.")
  }
  if (!("sf" %in% class(comparison_points))) {
    stop("comparison_points must be an sf point object.")
  }
  if (!all(sf::st_geometry_type(comparison_points) %in% c("POINT"))) {
    stop("comparison_points must be an sf point object.")
  }
  
  n_template <- nrow(template_points)
  n_comparison <- nrow(comparison_points)
  
  if (n_template < 1) {
    stop("There are no points in template_points")
  }
  if (n_comparison < 1) {
    stop("There are no points in comparison_points")
  }
  
  if (!identical(sf::st_crs(template_points), sf::st_crs(comparison_points))) {
    warning("comparison_points and template_points aren't in the same projection. Reprojecting comparison_points to match.")
    comparison_points <- sf::st_transform(x = comparison_points,
                                          crs = sf::st_crs(template_points))
  }
  
  template_points[["source"]] <- "template"
  template_points <- get_coords(points = template_points,
                                x_var = "XMETERS",
                                y_var = "YMETERS",
                                projection = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
  comparison_points[["source"]] <- "comparison"
  comparison_points <- get_coords(points = comparison_points,
                                  x_var = "XMETERS",
                                  y_var = "YMETERS",
                                  projection = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
  
  distance_matrix <- dist_matrix(dataframe = rbind(sf::st_drop_geometry(template_points[, c("XMETERS", "YMETERS")]),
                                                   sf::st_drop_geometry(comparison_points[, c("XMETERS", "YMETERS")])),
                                 x_var = "XMETERS",
                                 y_var = "YMETERS")
  
  template_indices <- 1:n_template
  comparison_indices <- (n_template + 1):(n_template + n_comparison)
  
  template_df <- do.call(rbind,
                         lapply(X = template_indices,
                                distance_matrix = distance_matrix,
                                comparison_indices = comparison_indices,
                                n_comparison = n_comparison,
                                FUN = function(X, distance_matrix, comparison_indices, n_comparison){
                                  distances <- distance_matrix[comparison_indices, X]
                                  output <- data.frame(template_index = X,
                                                       comparison_index = 1:n_comparison)
                                  output <- output[order(distances), ]
                                  output[["rank_by_template"]] <- 1:nrow(output)
                                  return(output)
                                }))
  
  comparison_df <- do.call(rbind,
                           lapply(X = comparison_indices,
                                  distance_matrix = distance_matrix,
                                  template_indices = template_indices,
                                  n_template = n_template,
                                  FUN = function(X, distance_matrix, template_indices, n_template){
                                    distances <- distance_matrix[template_indices, X]
                                    output <- data.frame(comparison_index = X - n_template,
                                                         template_index = 1:n_template)
                                    output <- output[order(distances), ]
                                    output[["rank_by_comparison"]] <- 1:nrow(output)
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
#' @param iteration_limit Numeric. The maximum number of iterations to attempt to sort before giving up. Defaults to \code{5000}.
#' @return A data frame with the variables \code{match_to_id} and \code{match_from_id} containing the paired members of the two sets.
#' @export
ranked_sort <- function(match_to,
                        match_from,
                        match_to_idvar = "match_to_id",
                        match_from_idvar = "match_from_id",
                        match_to_rankvar = "match_to_rank",
                        match_from_rankvar = "match_from_rank",
                        iteration_limit = 5000){
  # Just as a heads up, I am not proud of the variable names in this function, particularly inside the while().
  # They're definitely on the confusing side, but I promise this is an improvement over
  # what I called them at first.
  
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
  
  # How many to match from and to?
  # We care because we'll reference to see if the looping is done
  n_matchfrom <- length(unique(match_from[[match_from_idvar]]))
  n_matchto <- length(unique(match_to[[match_to_idvar]]))
  
  # For if we get in an endless loop, we'll have a way to hit the brakes
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
  
  # Here are the two data frames of pairs
  # Since nothing is paired yet, their partners are NA rather than id values
  matchfrom_pairs <- data.frame(id = unique(match_from[[match_from_idvar]]),
                                paired = NA,
                                stringsAsFactors = FALSE)
  # For some reason it's carrying over the variable name from match_from_idvar
  # so this corrects the variable names
  names(matchfrom_pairs) <- c("id", "paired")
  matchto_pairs <- data.frame(id = unique(match_to[[match_to_idvar]]),
                              paired = NA,
                              stringsAsFactors = FALSE)
  names(matchto_pairs) <- c("id", "paired")
  
  # Initializing these. Each pass through the while() will update these two values
  # It stops once the number of points matched from one of these groups is the same as the number of points in that group
  matchfrom_matched_n <- sum(!is.na(matchfrom_pairs[["paired"]]))
  matchto_matched_n <- sum(!is.na(matchto_pairs[["paired"]]))
  
  # This is just a panic option. If somehow the while() wouldn't exit otherwise, we set a limit on iterations
  iterations <- 0
  
  # This while() will keep the process running as long as not all the points have been paired
  # and we haven't hit the iteration limit.
  # while ((n_matchfrom > matchfrom_matched_n & n_matchto > matchto_matched_n) & iterations < iteration_limit) {
  while (any(is.na(matchto_pairs$paired)) & iterations < iteration_limit) {
    # Get the match_from ID that we'll work on. It's the first one that isn't already paired
    current_matchto_id <- matchto_pairs[is.na(matchto_pairs[["paired"]]), "id"][1]
    # It's unpaired. We use this so that once it finds a suitable partner the loop will stop searching
    current_unassigned <- TRUE
    
    # Find out how it feels about the other set of points
    current_preferences <- match_to[match_to[[match_to_idvar]] == current_matchto_id, ]
    # And get those other point IDs in order of preference to work through
    current_preferred_order <- current_preferences[current_preferences[[match_to_rankvar]], match_from_idvar]
    
    # Work through the potential pairings
    for (current_matchfrom_id in current_preferred_order) {
      # Only bother to even compare if it doesn't have a more preferred partner yet
      if (current_unassigned) {
        # Does the potential partner already have a point?
        current_matchfrom_id_pair <- matchfrom_pairs[matchfrom_pairs[["id"]] == current_matchfrom_id, "paired"]
        # If no, pair them off
        if (is.na(current_matchfrom_id_pair)) {
          matchfrom_pairs[matchfrom_pairs[["id"]] == current_matchfrom_id, "paired"] <- current_matchto_id
          matchto_pairs[matchto_pairs[["id"]] == current_matchto_id, "paired"] <- current_matchfrom_id
          # It's no longer unassigned
          current_unassigned <- FALSE
        } else {
          # If yes, find out how much it prefers the pair it's a part of
          matchfrom_preference_already <- match_from[match_from[[match_from_idvar]] == current_matchfrom_id & match_from[[match_to_idvar]] == current_matchfrom_id_pair, match_from_rankvar]
          # And how much it prefers the point we're working on
          matchfrom_preference_current <- match_from[match_from[[match_from_idvar]] == current_matchfrom_id & match_from[[match_to_idvar]] == current_matchto_id, match_from_rankvar]
          # If it likes the point we're asking about more, then pair them
          if (matchfrom_preference_current < matchfrom_preference_already) {
            matchfrom_pairs[matchfrom_pairs[["id"]] == current_matchfrom_id, "paired"] <- current_matchto_id
            matchto_pairs[matchto_pairs[["id"]] == current_matchto_id, "paired"] <- current_matchfrom_id
            # Make sure that its previous partner is flagged as unpaired!
            matchto_pairs[matchto_pairs[["id"]] == current_matchfrom_id_pair, "paired"] <- NA
            # It's no longer unassigned
            current_unassigned <- FALSE
          }
        }
      }
    }
    
    # Increment the iterations
    iterations <- iterations + 1
    # Update the numbers of matches made
    matchfrom_matched_n <- sum(!is.na(matchfrom_pairs[["paired"]]))
    matchto_matched_n <- sum(!is.na(matchto_pairs[["paired"]]))
  }
  
  
  # If we hit the iteration limit, something is wrong!!!!!!!!!!!
  if (iterations >= iteration_limit) {
    message("The iteration limit of ", iteration_limit, " has been reached without reaching an optimal solution. Check your rankings to make sure they make sense or use another method.")
  }
  
  names(matchfrom_pairs) <- c("match_from_id", "match_to_id")
  names(matchto_pairs) <- c("match_to_id", "match_from_id")
  
  # Return the complete set of pairs (so, not the one that had more points)
  if (!any(is.na(matchfrom_pairs[["match_to_id"]]))) {
    output <- matchfrom_pairs[, c("match_from_id", "match_to_id")]
  } else if (!any(is.na(matchto_pairs[["match_from_id"]]))) {
    output <- matchto_pairs[, c("match_from_id", "match_to_id")]
  } else {
    stop("Something went wrong and there are still points in both match_to and match_from without partners. I don't know what to tell you?")
  }
  
  names(output) <- c(match_from_idvar, match_to_idvar)
  
  return(output)
}

#' Select points that most closely approximate the distribution of another set of points
#' @description Sometimes you have a large collection of points which are not randomly distributed or spatially balanced and you would like a subset that more or less do. Given a template of points that are distributed the way you would like, this will return the closest existing point to each. This can be done taking into account membership in a group, either by having assigned it as a variable in both sets of points or by providing polygons that can be used to assign membership. By default, no stratification/membership is taken into account.
#' @param existing_points Point sf object. The points you would like to select from by comparing against \code{template_points}.
#' @param template_points  Point sf object. The points you would like to compare against \code{existing_points} in order to select a subset of those that most closely resemble the distribution of the template points.
#' @param strata_polygons Optional polygon sf object. Polygons assigned a variable with a name \code{stratafield} that contains the membership information (e.g. strata) to assign to \code{existing_points} and \code{template_points}. If \code{NULL} then no assignment will be attempted. Defaults to \code{NULL}.
#' @param stratafield Character string. If \code{strata_polygons} is not \code{NULL}, the name of the variable in \code{strata_polygons} that contains the membership information. Otherwise, the name of the variable in both \code{template_points} and \code{existing_points} that contains the membership information. If \code{NULL} then the points will be considered to belong to a single group. Defaults to \code{NULL}.
#' @param projection Character string. The projection to force all spatial objects into. Defaults to NAD83, \code{sp::CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")}.
#' @param iteration_limit Numeric. The maximum number of iterations to attempt to sort before giving up. Defaults to \code{5000}.
#' @return A spatial points data frame made by trimming \code{existing_points} down to the points that most closely approximate the distribution of \code{template_points} while also containing the same number of points as \code{template_points}. It will be in the projection specified by \code{projection}.
#' @export
get_closest <- function(existing_points,
                        template_points,
                        strata_polygons = NULL,
                        stratafield = NULL,
                        projection = "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0",
                        iteration_limit = 5000){
  
  if (!("sf" %in% class(existing_points))) {
    stop("existing_points must be a point sf object")
  }
  if (!all(sf::st_geometry_type(existing_points) %in% c("POINT"))) {
    stop("existing_points must be a point sf object")
  }
  if (!("sf" %in% class(template_points))) {
    stop("template_points must be a point sf object")
  }
  if (!all(sf::st_geometry_type(template_points) %in% c("POINT"))) {
    stop("template_points must be a point sf object")
  }
  if (nrow(template_points) < 1) {
    stop("There are no points in template_points.")
  }
  if (nrow(existing_points) < 1) {
    stop("There are no points in existing_points.")
  }
  
  # Reproject if necessary
  if (!identical(sf::st_crs(projection),
                 sf::st_crs(template_points))) {
    template_points <- sf::st_transform(x = template_points,
                                        crs = projection)
  }
  if (!identical(sf::st_crs(projection),
                 sf::st_crs(existing_points))) {
    existing_points <- sf::st_transform(x = existing_points,
                                        crs = projection)
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
      if (!(stratafield %in% names(template_points))) {
        stop("The variable ", stratafield, " does not appear in template_points.")
      }
      if (!(stratafield %in% names(existing_points))) {
        stop("The variable ", stratafield, " does not appear in existing_points.")
      }
      # Sorry for the inconsistencies but "MEMBERSHIP" is trying to get away from the assumption that these will be stratifications
      # And I'm so, so tired. I'll come back to clean this up later, or so I'm currently telling myself at the end of a ten hour day
      template_points[["MEMBERSHIP"]] <- template_points[[stratafield]]
      existing_points[["MEMBERSHIP"]] <- existing_points[[stratafield]]
    }
    
  } else {
    template_points[["MEMBERSHIP"]] <- "frame"
    existing_points[["MEMBERSHIP"]] <- "frame"
  }
  
  # If there is a set of stratification polygons, use them!
  if (!is.null(strata_polygons)) {
    if (!("sf" %in% class(strata_polygons))) {
      stop("strata_polygons must be a polygon sf object.")
    }
    if (!all(sf::st_geometry_type(strata_polygons) %in% c("POLYGON", "MULTIPOLYGON"))) {
      stop("strata_polygons must be a polygon sf object.")
    }
    if (!identical(sf::st_crs(projection), sf::st_crs(strata_polygons))) {
      strata_polygons <- sf::st_transform(x = strata_polygons,
                                          crs = projection)
    }
    if (!(stratafield %in% names(strata_polygons))) {
      stop("The variable ", stratafield, " does not appear in strata_polygons.")
    }
    
    # This just puts the strata into the points
    existing_points[["MEMBERSHIP"]] <- sf::st_intersection(x = existing_points,
                                                           y = strata_polygons)[[stratafield]]
    template_points[["MEMBERSHIP"]] <- sf::st_intersection(x = template_points,
                                                           y = strata_polygons)[[stratafield]]
  }
  
  
  # What strata are there?
  strata <- unique(c(existing_points[["MEMBERSHIP"]], template_points[["MEMBERSHIP"]]))
  
  # By stratum!
  strata_selections <- lapply(X = strata,
                              template_points = template_points,
                              existing_points = existing_points,
                              iteration_limit = iteration_limit,
                              FUN = function(X, template_points, existing_points, iteration_limit){
                                # Narrow it down to the points in the current stratum
                                stratum <- X
                                existing_points_stratum <- existing_points[existing_points[["MEMBERSHIP"]] == stratum, ]
                                template_points_stratum <- template_points[template_points[["MEMBERSHIP"]] == stratum, ]
                                
                                # What are their ranking of each other between template and comparison based on distance?
                                preferences <- find_preferences(template_points = template_points_stratum,
                                                                comparison_points = existing_points_stratum)
                                
                                # What's the optimal solution for minimizing distances for pairing?
                                sorting <- ranked_sort(match_to = preferences[["template"]],
                                                       match_from = preferences[["comparison"]],
                                                       match_to_idvar = "template_index",
                                                       match_from_idvar = "comparison_index",
                                                       match_to_rankvar = "rank_by_template",
                                                       match_from_rankvar = "rank_by_comparison",
                                                       iteration_limit = iteration_limit)
                                
                                # Only keep the points that got paired with the template points
                                output <- existing_points_stratum[sorting[["comparison_index"]], ]
                                
                                return(output)
                              })
  
  # Each stratum's points are a separate SPDF in the list, so let's mash 'em together
  output <- do.call(rbind,
                    strata_selections)
  
  # If we just had to add an identity, hide that fact
  if (is.null(stratafield)) {
    output[["MEMBERSHIP"]] <- NULL
  } else {
    # Otherwise name that membership field according to what the user called it
    names(output)[names(output) == "MEMBERSHIP"] <- stratafield
  }
  
  return(output)
}


#' Find and keep the points farthest from each other
#' @description This will take a set of existing points and new points and combine them to create a set consisting of the existing points and the farthest new points. The original intended use case was to take a collection of sampling locations from one or more sample designs (\code{existing_points}) and use them as part of a new, spatially balanced sample design. The function takes a set of new, random, spatially balanced points (\code{new_points}) and determines the distance between each of them and each of the existing points. It then sequentially eliminates the new point closest to any existing point until the combined number of existing points and remaining points is equal to \code{target}.
#' @param existing_points Point sf object. These are the points that will all be included in the output points and the points against which \code{new_points} will be compared against.
#' @param new_points Point sf object. These are the points that may be included in the output. The number that will be is equal to \code{target - nrow(existing_points)}.
#' @param target Numeric value. The total number of points to include in the output. Defaults to \code{nrow(new_points) - nrow(existing_points)}.
#' @param projection Character string. The projection to force all spatial objects into. Defaults to NAD83, \code{"+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"}.
#' @param verbose Logical value. If \code{TRUE} then additional diagnostic messages will be produced while the function runs. Defaults to \code{FALSE}.
#' @return A point sf object containing all the points from \code{existing_points} and \code{target - nrow(existing_points)} points from \code{new_points} using the CRS specified in \code{projection}. It will only have those variables that were in common between both \code{existing_points} and \code{new_points}.
#' @export
keep_farthest <- function(existing_points,
                          new_points,
                          target = NULL,
                          projection = "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0",
                          verbose = FALSE){
  
  if (!("sf" %in% class(existing_points))) {
    stop("existing_points must be a point sf object.")
  }
  if (!all(sf::st_geometry_type(existing_points) %in% c("POINT"))) {
    stop("existing_points must be a point sf object.")
  }
  if (nrow(existing_points) < 1) {
    stop("There are no points in existing_points")
  }
  if (!("sf" %in% class(new_points))) {
    stop("existing_points must be a point sf object.")
  }
  if (!all(sf::st_geometry_type(new_points) %in% c("POINT"))) {
    stop("existing_points must be a point sf object.")
  }
  if (nrow(new_points) < 1) {
    stop("There are no points in new_points")
  }
  
  if (!identical(sf::st_crs(projection), sf::st_crs(existing_points))) {
    existing_points <- sf::st_transform(x = existing_points,
                                        crs = projection)
  }
  
  # Putting aside the input new_points so we can
  if (!identical(sf::st_crs(projection), sf::st_crs(new_points))) {
    new_points <- sf::st_transform(x = new_points,
                                   crs = projection)
  }
  
  # How many of each point type are there? We'll use these for the loops
  n_existing <- nrow(existing_points)
  n_new <- nrow(new_points)
  
  # Just in cae we're deciding the target count right now
  count_difference <- n_new - n_existing
  if (is.null(target)) {
    target <- count_difference
  }
  if (target > nrow(new_points)) {
    stop("The target number of points to return is greater than the number of points available")
  }
  if (target < 1) {
    stop("The target number of points to return is less than 1")
  }
  
  # I don't think this matters????
  # if (target <= nrow(existing_points@data)) {
  #   stop("The target number of points is less than or equal to the number of existing points.")
  # }
  
  existing_point_vars <- names(existing_points)
  new_point_vars <- names(new_points)
  if (!all(new_point_vars %in% existing_point_vars) | !all(existing_point_vars %in% new_point_vars)) {
    message("existing_points and new_points must have all the same variables as each other")
  } else if (ncol(existing_points) > 1) {
    existing_points <- existing_points[, new_point_vars]
  }
  
  
  
  # Get some common info added to these
  existing_points[["TYPE"]] <- "EXISTING"
  existing_points[["INDEX"]] <- 1:n_existing
  new_points[["TYPE"]] <- "NEW"
  new_points[["INDEX"]] <- 1:n_new
  
  # Combine the two sets of points, making sure the existing points come first!!!
  combined_points <- rbind(existing_points[, c("INDEX", "TYPE")],
                           new_points[, c("INDEX", "TYPE")])
  
  # Add the coordinates
  combined_points <- get_coords(combined_points,
                                x_var = "XMETERS",
                                y_var = "YMETERS",
                                projection = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
  
  # Get a distance matrix
  distance_matrix <- dist_matrix(dataframe = as.data.frame(sf::st_drop_geometry(combined_points)),
                                 x_var = "XMETERS",
                                 y_var = "YMETERS")
  
  # Remove the columns for new points
  distance_matrix <- distance_matrix[, -((n_existing + 1):(n_existing + n_new))]
  # Remove the rows for existing points
  distance_matrix <- distance_matrix[-(1:n_existing), ]
  
  # Convert to a data frame
  # SO!!!! Each row should represent a new point and each column represents an existing point
  # That means that to find the new point closest to an existing point, you look for the minimum value in a column
  distance_df <- as.data.frame(distance_matrix,
                               stringsAsFactors = FALSE)
  
  # Add the new points indices
  distance_df[["INDEX"]] <- 1:n_new
  
  # And now we loop to remove all the closest new points
  # The looping is so that we don't get hung up on the same minimum distance
  # and can remove progressively more distant points until we reach our goal
  # Each loop will:
  # 1) Check for the new point(s) closest to an existing point
  # 2) Store the index (or indices) from distance_df[["INDEX"]]
  # 3) Remove that observation(s) from distance_df
  # This will continue until the number of stored indices is equal to the number of existing points
  # If the number to remove is overshot because the final pass through the loop identifies multiple indices,
  # we only take the first however-many-we-need
  # The observations at the identified indices in the new points will be removed!
  
  
  # A vector to store the indices to chuck
  removal_indices <- NULL
  n_removal_indices <- 0
  # Here's one we can mutilate over our iterations (always avoid violence to original data you may reference again!)
  working_distance_df <- distance_df
  # We'll be ignoring the index column because of course sequential ordinals starting at 1 are the smallest numbers in the matrix!
  # The indices screw up min() results without removing this for the finding minimum evaluation step
  index_colnum <- grep(names(distance_df), pattern = "INDEX")
  # So, given the number of points we have and the target
  n_indices_to_remove <- (n_new - target)
  if (verbose) {
    message("Aiming to drop ", n_indices_to_remove, " points")
  }
  
  while (n_removal_indices < n_indices_to_remove) {
    if (verbose) {
      message("Starting while() iteration with ", n_removal_indices, " indices to remove identified")
    }
    current_min <- min(working_distance_df[, -index_colnum])
    if (verbose) {
      message("Current minimum distance is ", current_min)
    }
    # Get the indices from every column where that min occurs
    # Each column is an existing point, so if we check every column for the value and store that index,
    # those are the new points that are that distance from an existing point
    current_indices <- unlist(sapply(X = 1:(ncol(working_distance_df) - 1),
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
                                     }))
    # Make sure we have the uniques, in case a new point was equidistant from multiple existing points
    current_indices <- unique(current_indices)
    if (verbose) {
      message("The current indices are: ", paste(current_indices, collapse = ", "))
    }
    # Remove the NULLs
    # current_indices <- current_indices[!is.na(current_indices)]
    
    # Add those to the ones we're going to remove
    removal_indices <- unique(c(removal_indices, current_indices))
    # Drop any NULLs
    # removal_indices <- removal_indices[!sapply(removal_indices, is.null)]
    
    # Out of paranoia
    removal_indices <- unlist(removal_indices)
    if (verbose) {
      message("The current full set of removal indices is: ", paste(removal_indices, collapse = ", "))
    }
    
    if (!is.null(removal_indices)) {
      # Make sure we don't overshoot our removal goal
      removal_indices <- removal_indices[1:min(length(removal_indices), n_indices_to_remove)]
      
      # Removing that point so that we can identify a new minimum in the next pass
      working_distance_df <- working_distance_df[!(working_distance_df[[index_colnum]] %in% removal_indices), ]
      
      # Update our tracking value
      n_removal_indices <- length(removal_indices)
    }
  }
  
  
  # removal_indices <- removal_indices[!sapply(removal_indices, is.null)]
  
  # Now that we have our indices to remove, let's do it
  # The as.numeric() is because due to the NULL that's in there from the pre-loop setup removal_indices is a list, not a vector
  output <- new_points[-as.numeric(removal_indices), new_point_vars]
  
  return(output)
}


#' Combine existing and new points to create a spatially balanced design
#' @param existing_points Point sf object. The existing points that will be balanced around.
#' @param new_points Point sf object. The points that will be compared against the existing points and selected from to create a balanced design.
#' @param stratafield Character string. The name of the variable in common between \code{existing_points} and \code{new_points} that contains stratum identities. This is used to balance by stratum. If \code{NULL} then balancing will not take strata into account. Defaults to \code{NULL}.
#' @param id_var Character string. The name of the variable common between \code{existing_points} and \code{new_points} that contains point identities. This assumes that the format of plot IDs is "character-number".
#' @param projection CRS object or character string. The projection to force on the spatial objects. Defaults to \code{"+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"}.
#' @return A point sf object containing all the points from \code{existing_points} and the selected points from \code{new_points}. The projection will match \code{projection}.
#' @export
balance_around <- function(existing_points,
                           new_points,
                           stratafield = NULL,
                           id_var,
                           projection = "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"){
  if (!("sf" %in% class(existing_points))) {
    stop("existing_points must be a point sf object.")
  }
  if (!all(sf::st_geometry_type(existing_points) %in% c("POINT"))) {
    stop("existing_points must be a point sf object.")
  }
  if (!(id_var %in% names(existing_points))) {
    stop("id_var must correspond to a variable in existing_points.")
  }
  if (!identical(sf::st_crs(projection), sf::st_crs(existing_points))) {
    existing_points <- sf::st_transform(x = existing_points,
                                        crs = projection)
  }
  
  if (!("sf" %in% class(new_points))) {
    stop("new_points must be a point sf object.")
  }
  if (!all(sf::st_geometry_type(new_points) %in% c("POINT"))) {
    stop("new_points must be a point sf object.")
  }
  if (!(id_var %in% names(existing_points))) {
    stop("id_var must correspond to a variable in new_points.")
  }
  if (!identical(sf::st_crs(projection), sf::st_crs(new_points))) {
    new_points <- sf::st_transform(x = new_points,
                                   crs = projection)
  }
  
  
  # Assign the codes that indicate if they're existing plots or freshly-drawn ones
  existing_points[["TYPE"]] <- "EXISTING"
  new_points[["TYPE"]] <- "NEW"
  
  # Add in the fields that the points don't have for easy combination
  missing_vars_new <- names(existing_points)[!(names(existing_points) %in% names(new_points))]
  new_points[, missing_vars_new] <- NA
  missing_vars_existing <- names(new_points)[!(names(new_points) %in% names(existing_points))]
  existing_points[, missing_vars_existing] <- NA
  
  # Bind the 2 point files together
  pts <- rbind(existing_points, new_points)
  
  # What are the existing points' indices?
  extant_indices <- 1:nrow(existing_points)
  # Which are the new points' indices?
  new_indices <- (nrow(existing_points) + 1):(nrow(existing_points) + nrow(new_points))
  
  # Rename the date that the plot was sampled to "PREVDATE"
  # pts@data[["PREVDATE"]] <- pts@data[["DATEVISITE"]]
  
  # Restrict to only relevant fields
  # pts <- pts[ , c("TYPE",
  #                 stratafield,
  #                 "PLOTID",
  #                 "PLOTKEY",
  #                 "PRIMARYKEY",
  #                 "PANEL",
  #                 "PROJECTNAM")]
  
  # Make these all NA for the new points
  # new_points@data[new_indices, "PLOTID"] <- NA
  # new_points@data[new_indices, "PLOTKEY"] <- NA
  # new_points@data[new_indices, "PRIMARYKEY"] <- NA
  # new_points@data[new_indices, "PROJECTNAM"] <- NA
  # new_points@data[new_indices, "PREVDATE"] <- NA
  
  # Determine the number of existing points
  extant <- nrow(existing_points)
  
  # Add coordinates to the combined points for distance calculations
  pts <- get_coords(pts,
                    x_var = "XMETERS",
                    y_var = "YMETERS",
                    projection = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
  
  # This determines the number of New points to eliminate, and eliminates the points
  pts <- keep_farthest(existing_points = existing_points,
                       new_points = new_points)
  
  # Time to tweak the new points that were kept
  new_indices_remaining <- pts[["TYPE"]] == "NEW"
  
  # Get the existing plot ids and renumber them
  if (any(new_indices_remaining)) {
    plotids <- pts[new_indices_remaining, id_var]
    plotids <- gsub(plotids,
                    pattern = "\\d*$",
                    replacement = "")
    pts[new_indices_remaining, id_var] <- paste0(plotids, 1:length(new_indices_remaining))
  }
  
  
  # None of the points are considered visited now!
  # pts@data[["DATEVIS"]] <- ""
  # pts@data[["EVALSTA"]] <- "NotEval"
  # pts@data[["FINAL_DESI"]] <- ""
  
  # Add coordinates!
  pts <- get_coords(pts,
                    x_var = "LAT",
                    y_var = "LONG",
                    projection = projection)
  
  return(pts)
}


#' Select points that most closely approximate the distribution of another set of points and substitute them
#' @description When creating a design that incorporates points from another design to revisit, it can be important to approximate the same spatial distribution, especially if the design is spatially balanced. This will substitute in revisit points that will cause the least shift in the spatial distribution.
#' @param sub_points Point sf object. The points you would like to substitute from by comparing against \code{template_points}.
#' @param template_points  Point sf object. The points you would like to compare against \code{sub_points} in order to select a subset of those that most closely resemble the distribution of the template points.
#' @param sub_idvar Character string. The name of the variable in \code{sub_points} that contains the unique identifiers. Defaults to \code{"plotid"}.
#' @param template_idvar Character string. The name of the variable in \code{template_points} that contains the unique identifiers. Defaults to \code{"plotid"}.
#' @param sub_counts Optional named numeric vector. A vector of the number of template points to replace with the closest approximation from the existing points in each stratum. Only use if there is more than one stratum. The names of the vector must exactly match the strata of the points. If \code{NULL} then points will be allocated proportionally by area using the value in \code{replacement_count}. Defaults to \code{NULL}.
#' @param strata_polygons Optional polygon sf object. Polygons assigned a variable with a name \code{polygons_stratavar} that contains the membership information (e.g. strata) to assign to \code{sub_points} and \code{template_points}. If \code{NULL} then no assignment will be attempted. Defaults to \code{NULL}.
#' @param polygons_stratavar Optional character string. If \code{strata_polygons} is not \code{NULL}, the name of the variable in \code{strata_polygons} that contains the membership information. Otherwise, the name of the variable in both \code{template_points} and \code{sub_points} that contains the membership information. If \code{NULL} then the points will be considered to belong to a single group. Defaults to \code{"stratum"}.
#' @param sub_stratavar Optional character string. The name of the variable in \code{sub_points} that contains the membership information. If both \code{sub_stratavar} and \code{template_stratavar} are \code{NULL} then all the points in both will be considered to belong to a single group. Defaults to \code{"stratum"}.
#' @param template_stratavar Optional character string. The name of the variable in \code{template_points} that contains the membership information. If both \code{sub_stratavar} and \code{template_stratavar} are \code{NULL} then all the points in both will be considered to belong to a single group. Defaults to \code{"stratum"}.
#' @param projection Optional character string or CRS object. The projection to force all spatial objects into, e.g. \code{"+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"}. If \code{NULL} then the projection from \code{template_points@@proj4string} will be used. Defaults to \code{NULL}.
#' @param iteration_limit Numeric. The upper limit on the number of iterations to try to match points in \code{sub_points} to their best-fitting mates in \code{template_points}. Defaults to \code{5000}.
#' @return A spatial points data frame of the \code{template_points} with points from \code{sub_points} substituted in. It will be in the projection specified by \code{projection}.
#' @export
combine_designs <- function(sub_points,
                            template_points,
                            sub_idvar = "plotid",
                            template_idvar = "plotid",
                            sub_counts = NULL,
                            strata_polygons = NULL,
                            polygons_stratavar = "stratum",
                            sub_stratavar = "stratum",
                            template_stratavar = "stratum",
                            projection = NULL,
                            iteration_limit = 5000){
  # SO MUCH SANITIZATION
  if (is.null(sub_counts)) {
    stop("You must supply a named vector of counts per stratum as sub_counts")
  }
  
  if (!("sf" %in% class(sub_points))) {
    stop("sub_points must be a point sf object")
  }
  if (!all(sf::st_geometry_type(sub_points) %in% c("POINT"))) {
    stop("sub_points must be a point sf object")
  }
  
  if (!(sub_idvar %in% names(sub_points))) {
    stop("The variable ", sub_idvar, " does not appear in sub_points")
  } else {
    if (length(unique(sub_points[[sub_idvar]])) != nrow(sub_points)) {
      stop("The variable ", sub_idvar, " does not contain unique identifiers in sub_points")
    }
  }
  
  if (!("sf" %in% class(template_points))) {
    stop("template_points must be a point sf object")
  }
  if (!all(sf::st_geometry_type(template_points) %in% c("POINT"))) {
    stop("template_points must be a point sf object")
  }
  
  if (!(template_idvar %in% names(template_points))) {
    stop("The variable ", template_idvar, " does not appear in template_points")
  } else {
    if (length(unique(template_points[[template_idvar]])) != nrow(template_points)) {
      stop("The variable ", template_idvar, " does not contain unique identifiers in template_points")
    }
  }
  
  
  
  if (is.null(projection)) {
    projection <- sf::st_crs(template_points)
  }
  
  if (!is.null(iteration_limit)) {
    if (class(iteration_limit) != "numeric" | length(iteration_limit) > 1) {
      stop("iteration_limit must be a single numeric value and determines the number of iterations to run while trying to find pairs. It defaults to the number of potential substitution points.")
    }
  }
  
  if (!identical(sf::st_crs(sub_points), sf::st_crs(projection))) {
    sub_points <- sf::st_transform(x = sub_points,
                                   crs = projection)
  }
  if (!identical(sf::st_crs(template_points), sf::st_crs(projection))) {
    template_points <- sf::st_transform(x = template_points,
                                        crs = projection)
  }
  
  if (is.null(strata_polygons)) {
    if (xor(is.null(sub_stratavar), is.null(template_stratavar))) {
      stop("You must either provide the stratification variables for both sub_points and template_points (for a stratified approach without strata polygons) or neither (for a non-stratified approach).")
    }
    if (is.null(sub_stratavar) & is.null(template_stratavar)) {
      sub_points[["MEMBERSHIP"]] <- "frame"
      template_points[["MEMBERSHIP"]] <- "frame"
      sub_stratavar <- template_stratavar <- "MEMBERSHIP"
    } else {
      if (!(sub_stratavar %in% names(sub_points))) {
        stop("The variable ", sub_stratavar, " does not occur in sub_points")
      }
      if (!(template_stratavar %in% names(template_points))) {
        stop("The variable ", template_stratavar, " does not occur in template_points")
      }
      sub_points[["MEMBERSHIP"]] <- sub_points[[sub_stratavar]]
      template_points[["MEMBERSHIP"]] <- template_points[[template_stratavar]]
    }
  } else {
    if (!("sf" %in% class(strata_polygons))) {
      stop("strata_polygons must be a polygon sf object")
    }
    if (!all(sf::st_geometry_type(strata_polygons) %in% c("POLYGON", "MULTIPOLYGON"))) {
      stop("strata_polygons must be a polygon sf object")
    }
    if (!(polygons_stratavar %in% names(strata_polygons))) {
      stop("The variable ", polygons_stratavar, " does not occur in strata_polygons")
    }
    if (!identical(sf::st_crs(strata_polygons), sf::st_crs(projection))) {
      strata_polygons <- sf::st_transform(x = strata_polygons,
                                          crs = projection)
    }
    
    if (is.null(sub_stratavar) & is.null(template_stratavar)) {
      warning("Ignoring sub_stratavar and template_stratavar in favor of strata_polygons")
    }
    
    sub_points[["MEMBERSHIP"]] <- sf::st_intersection(x = sub_points,
                                                      y = strata_polygons)[[polygons_stratavar]]
    template_points[["MEMBERSHIP"]] <- sf::st_intersection(x = template_points,
                                                           y = strata_polygons)[[polygons_stratavar]]
  }
  
  if (any(is.na(sub_points[["MEMBERSHIP"]])) | any(is.na(template_points[["MEMBERSHIP"]]))) {
    warning("Excluding points that have no stratification membership")
    sub_points <- sub_points[is.na(sub_points[["MEMBERSHIP"]]), ]
    template_points <- template_points[is.na(template_points[["MEMBERSHIP"]]), ]
  }
  
  # The as.character() is in case there are factors
  sub_points[["MEMBERSHIP"]] <- as.character(sub_points[["MEMBERSHIP"]])
  template_points[["MEMBERSHIP"]] <- as.character(template_points[["MEMBERSHIP"]])
  
  
  # What strata are there?
  strata <- unique(c(sub_points[["MEMBERSHIP"]], template_points[["MEMBERSHIP"]]))
  
  if (!is.null(sub_counts)) {
    if (!all(names(sub_counts) %in% strata)) {
      stop("Not all strata named in sub_counts are associated with the points")
    }
    if (!all(strata %in% names(sub_counts))) {
      stop("Not all strata associated with the points are found in sub_counts")
    }
  }
  
  
  # By stratum!
  strata_selections <- lapply(X = strata,
                              template_points = template_points,
                              sub_points = sub_points,
                              sub_counts = sub_counts,
                              iteration_limit = iteration_limit,
                              FUN = function(X, template_points, sub_points, sub_counts, iteration_limit){
                                # Narrow it down to the points in the current stratum
                                stratum <- X
                                n_keep <- sub_counts[[stratum]]
                                sub_points_stratum <- sub_points[sub_points[["MEMBERSHIP"]] == stratum, ]
                                template_points_stratum <- template_points[template_points[["MEMBERSHIP"]] == stratum, ]
                                
                                # What if there are no points to sub?????
                                if (n_keep < 1) {
                                  return(NULL)
                                }
                                if (nrow(sub_points_stratum) < 1) {
                                  warning("No substitution points available in ", stratum,
                                          ". Making zero substitutions instead of ", n_keep, ".")
                                  return(NULL)
                                }
                                if (nrow(template_points_stratum) < 1) {
                                  warning("No template points available in ", stratum,
                                          ". Making zero substitutions instead of ", n_keep, ".")
                                  return(NULL)
                                }
                                
                                # Add the coordinates
                                sub_points_stratum_coords <- as.data.frame(sf::st_coordinates(sub_points_stratum))
                                names(sub_points_stratum_coords) <- c("xcoord", "ycoord")
                                sub_points_stratum <- cbind(sub_points_stratum,
                                                            sub_points_stratum_coords)
                                template_points_stratum_coords <- as.data.frame(sf::st_coordinates(template_points_stratum))
                                names(template_points_stratum_coords) <- c("xcoord", "ycoord")
                                # We specifically remove any existing xcoord and ycoord variables
                                # in case they're from a different projection
                                template_points_stratum <- cbind(dplyr::select(template_points_stratum,
                                                                               -dplyr::matches("xcoord"),
                                                                               -dplyr::matches("ycoord")),
                                                                 template_points_stratum_coords)
                                
                                # What are their ranking of each other between template and comparison based on distance?
                                preferences <- find_preferences(template_points = template_points_stratum,
                                                                comparison_points = sub_points_stratum)
                                
                                # What's the optimal solution for minimizing distances for pairing?
                                sorting <- ranked_sort(match_to = preferences[["template"]],
                                                       match_from = preferences[["comparison"]],
                                                       match_to_idvar = "template_index",
                                                       match_from_idvar = "comparison_index",
                                                       match_to_rankvar = "rank_by_template",
                                                       match_from_rankvar = "rank_by_comparison",
                                                       iteration_limit = iteration_limit)
                                
                                # Get a data frame of the paired IDs
                                pairs <- data.frame(sub_id = sub_points_stratum[sorting[["comparison_index"]], c(sub_idvar, "xcoord", "ycoord")],
                                                    template_id = template_points_stratum[sorting[["template_index"]], c(template_idvar, "xcoord", "ycoord")],
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
                                # I don't know why the template var names changed, but this catches them if it happens
                                names(pairs) <- gsub(names(pairs),
                                                     pattern = "_X$",
                                                     replacement = "_xcoord")
                                names(pairs) <- gsub(names(pairs),
                                                     pattern = "_Y$",
                                                     replacement = "_ycoord")
                                
                                # And remove optional variables
                                pairs <- pairs[, !grepl(names(pairs), pattern = "_optional$")]
                                
                                # Now we'll keep the n closest where n is the value in sub_counts for this stratum
                                # First up is the distance between these pairs (apply() was acting up, so sapply() it is)
                                pairs[["distance"]] <- unlist(sapply(X = 1:nrow(pairs),
                                                                     df = pairs,
                                                                     FUN = function(X, df){
                                                                       (df[X, "sub_xcoord"] - df[X, "template_xcoord"])^2 + (df[X, "sub_ycoord"] - df[X, "template_ycoord"])^2
                                                                     }))
                                
                                pairs <- pairs[order(pairs[["distance"]]), ]
                                
                                if (n_keep > nrow(pairs)) {
                                  stop("Attempting to keep ", n_keep, " points in ", stratum, " but only ", nrow(pairs), " plot pairs are available.")
                                }
                                
                                return(pairs[1:n_keep, ])
                              })
  
  strata_selection_df <- do.call(rbind,
                                 strata_selections)
  
  # This part builds a new spatial points data frame that subs in the existing points that were identified
  # IT ASSUMES THAT THE TEMPLATE POINTS ARE IN ORDER!!!!!!!!!!!!!!!
  template_points_coords <- as.data.frame(sf::st_coordinates(template_points))
  names(template_points_coords) <- c("xcoord", "ycoord")
  template_df <- cbind(dplyr::select(sf::st_drop_geometry(template_points),
                                     -dplyr::matches("xcoord"),
                                     -dplyr::matches("ycoord")),
                       template_points_coords)
  template_df[["order_number"]] <- 1:nrow(template_df)
  template_df <- template_df[, c("order_number", template_idvar, "xcoord", "ycoord")]
  names(template_df)[names(template_df) == template_idvar] <- "template_plotid"
  
  # Because they now have prefixes
  current_sub_idvar <- paste0("sub_", sub_idvar)
  current_template_idvar <- paste0("template_", template_idvar)
  
  merged_points <- merge(x = template_df,
                         y = strata_selection_df[, c(current_sub_idvar, current_template_idvar, "sub_xcoord", "sub_ycoord")],
                         by.x = "template_plotid",
                         by.y = current_template_idvar,
                         all.x = TRUE)
  
  output_df <- merged_points
  output_df[["plotid"]] <- output_df[["template_plotid"]]
  output_df[["revisit"]] <- !is.na(output_df[[current_sub_idvar]])
  output_df[output_df[["revisit"]], c("plotid", "xcoord", "ycoord")] <- output_df[output_df[["revisit"]], c(current_sub_idvar, "sub_xcoord", "sub_ycoord")]
  
  output_sf <- sf::st_as_sf(x = output_df[, c("order_number", "plotid", "revisit", "xcoord", "ycoord")],
                            coords = c("xcoord", "ycoord"),
                            crs = projection)
  
  return(output_sf)
}
