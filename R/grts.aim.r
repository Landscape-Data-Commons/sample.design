#' Draw spatially balanced points with \code{spsurvey::grts()}, output in AIM format
#'
#' @description A wrapper for \code{spsurvey::grts()} that modifies the output SPDF to match the format expected by other parts of the AIM workflow, including distribution to field crews and populating design databases.
#' @param design_object Optional if not providing \code{n_base} and \code{n_over}. A list of lists structured and named according to the requirements of \code{spsurvey::grts()} prior to version 5.3.0. Can be automatically constructed with \code{allocate.panels()}.
#' @param n_base Optional if not providing \code{design_object}. A named numeric value for the numbers of base points to draw in each stratum in the design. Names must match the values in \code{frame} in the variable specified by \code{stratum_field}.
#' @param n_over Optional if not providing \code{design_object}. A named numeric value for the numbers of oversample (spare) points to draw  in each stratum in the design. Names must match the values in \code{frame} in the variable specified by \code{stratum_field}.
#' @param project_name A character string of the name to use for the design, often the name of the project.
#' @param frame Polygon sf object with a variable name matching \code{stratum_field} which contains values matching the strata names in \code{design_object} OR in \code{n_base} and \code{n_over}.
#' @param stratum_var A character string representing the name of the variable in \code{frame} containing the strata names matching those in \code{design_object} OR in \code{n_base} and \code{n_over}. Defaults to \code{"STRATUM"}.
#' @param seed_number An optional numeric value to be passed to \code{set.seed()} to make the output reproducible. If \code{NULL} then a random seed will be used. Defaults to \code{NULL}.
#' @param projection Optional character string or CRS object. The projection to reproject the frame into if no \code{NULL}. Defaults to \code{NULL}.
#' @return A Spatial Points Data Frame of the sampling locations with the fields \code{PLOTID}, \code{STRATUM}, \code{PANEL}, \code{IntPtWt} (initial point weight), \code{xcoord}, and \code{ycoord}
#' @export
grts_aim <- function(design_object = NULL,
                     n_base = NULL,
                     n_over = NULL,
                     project_name = "Design name",
                     frame,
                     stratum_var = "STRATUM",
                     seed_number = NULL,
                     projection = NULL
){
  if (!("sf" %in% class(frame))) {
    stop("frame must be a polygon sf object")
  }
  if (!all(sf::st_geometry_type(frame) %in% c("POLYGON", "MULTIPOLYGON"))) {
    stop("frame must be a polygon sf object")
  }
  
  ## Reproject the frame
  if (!is.null(projection)) {
    frame <- sf::st_transform(x = frame,
                              crs = projection)
  }
  
  if (!(stratum_var %in% names(frame))) {
    stop(paste0("The variable ", stratum_var, " does not appear in frame."))
  }
  
  if (is.null(design_object) & (is.null(n_base) | is.null(n_over))) {
    stop("You must either provide a stuctured list as design_object or numeric vectors as n_base and n_over.")
  }
  
  if (!is.null(design_object)) {
    message("Using provided design object.")
    
    n_base <- sapply(X = design_object,
                     FUN = function(X){
                       sum(X$panel)
                     })
    names(n_base) <- names(design_object)
    n_over <- sapply(X = design_object,
                     FUN = function(X){
                       X$over
                     })
    names(n_over) <- names(design_object)
    
    panel_vector <- unname(unlist(sapply(X = design_object,
                                         FUN = function(X){
                                           base_panels <- sapply(X = names(X$panel),
                                                                 base_vector = X$panel,
                                                                 FUN = function(X, base_vector){
                                                                   rep(x = X,
                                                                       times = base_vector[X])
                                                                 })
                                           over_panels <- rep(x = "Oversample",
                                                              times = X$over)
                                           
                                           c(base_panels,
                                             over_panels)
                                         })))
    
  } else if (is.null(n_base) & is.null(n_over)) {
    stop("If not providing design_object, you must provide both n_base and n_over")
  } else {
    message("Using n_base and n_over.")
  }
  
  if(!any(names(n_base) %in% frame[[stratum_var]])) {
    stop("There are no matches between the names of the vector n_base and the strata in frame.")
  }
  if(!any(names(n_over) %in% frame[[stratum_var]])) {
    stop("There are no matches between the names of the vector n_over and the strata in frame.")
  }
  
  if (!all(names(n_base) %in% names(n_over)) | !all(names(n_over) %in% names(n_base))) {
    stop("n_base and n_over do not share all their strata")
  }
  
  if (!is.null(seed_number)) {
    if(class(seed_number) == "numeric") {
      set.seed(seed_number)
    } else {
      stop("seed_number is not a numeric value.")
    }
  }
  
  ## Invoke spsurvey::grts() first
  sample_design <- spsurvey::grts(sframe = frame,
                                  n_base = n_base,
                                  stratum_var = stratum_var,
                                  seltype = "equal",
                                  n_over = n_over,
                                  DesignID = "Site",
                                  sep = "-")
  
  sample_points_list <- lapply(X = names(n_base),
                               base_points = sample_design$sites_base,
                               over_points = sample_design$sites_over,
                               FUN = function(X, base_points, over_points) {
                                 rbind(base_points[base_points$stratum == X, ],
                                       over_points[over_points$stratum == X, ])
                               })
  sample_points <- do.call(rbind,
                           sample_points_list)
  
  # Change the variable name
  names(sample_points)[names(sample_points) == "siteuse"] <- "panel"
  
  # If there are panels to worry about from a design object, put them in here!
  if (!is.null(design_object)) {
    sample_points[["panel"]] <- panel_vector
  }
  
  # Add in the design name
  sample_points[["ProjectName"]] <- project_name
  
  # Add in the design type
  sample_points[["Design"]] <- "Random"
  
  # Add in the visit number (these are new, so V1)
  sample_points[["Visit"]] <- "V1"
  
  # Add in the evaluation status (these are new, so 'NotEval')
  sample_points[["EvalStatus"]] <- "NotEval"
  
  # Change the plot ID variable name
  names(sample_points)[names(sample_points) == "siteID"] <- "PlotID"
  
  # Change the coordinate variable names
  names(sample_points)[names(sample_points) == "lon_WGS84"] <- "DesignLongWGS"
  names(sample_points)[names(sample_points) == "lat_WGS84"] <- "DesignLatWGS"
  
  ## Reproject the sample sites to WGS84
  sample_points <- sf::st_transform(x = sample_points,
                                    crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

  ## Dropping the extra fields and renaming the remaining fields
  fields_relevant <- c("stratum", "panel", "EvalStatus", "PlotID", "Visit", "ProjectName", "Design", "DesignLongWGS", "DesignLatWGS")

  return(sample_points[, fields_relevant])
}

