#' A quick wrapper to execute simple designs
#' @description Take some combination of strata, sample frame, and points with a design object or values for constructing one and produce a sample design as a Spatial Points Data Frame. Automatic generation of a design object is done using allocate_panels().
#' @param design_name A character string of the name to use for the design, often the name of the project.
#' @param strata_spdf An optional Spatial Polygons Data Frame. Must either contain an attribute which exactly matches the string \code{stratum_field} and holds the stratum identities OR, if using \code{strata_lut}, an attribute which exactly matches \code{strata_lut_field} and contains values found in \code{strata_lut$strata_lut_field}. If \code{points_spdf} is provided, this will be used to attribute it with strata identities.
#' @param sampleframe_spdf An optional Spatial Polygons Data Frame. Necessary only to restrict \code{strata_spdf} or \code{points_spdf} if they are not already limited to the extent of the intended sample frame. If it is the only SPDF provided because there is no stratification scheme, it still must contain an attribute which exactly matches the string \code{stratum_field}, although there will only be one unique value in it.
#' @param points_spdf An optional Spatial Points Data Frame. This is only necessary if the draw should be from within a set of discrete points and not from polygon areas. If neither \code{strata_spdf} nor \code{sampleframe_spdf} with strata identity information is provided, this must either contain an attribute which exactly matches the string \code{stratum_field} and holds the stratum identities OR, if using \code{strata_lut}, an attribute which exactly matches \code{strata_lut_field} and contains values found in \code{strata_lut$strata_lut_field}.
#' @param strata_lut An optional data frame. If provided, this must contain a variable with a name that exactly matches the string provided in \code{stratum_field} and a variable with a name that exactly matches the string provided in \code{strata_lut_field}. There must be unique relationships between the values in the two fields, i.e. a value in \code{strata_lut$strata_lut_field} CANNOT correspond with more than one value in \code{strata_lut$stratum_field}. The values in \code{strata_lut$strata_lut_field} must match values in the \code{data$strata_lut_field} slot of the provided SPDFs.
#' @param strata_lut_field An optional string specifying a column/attribute/variable name. If \code{strata_lut} is provided, then a variable with a name exactly matching this string must exist in both \code{strata_lut} and one of the SPDFs provided.
#' @param design_object An optional specially-constructed list to use instead of an automatically-generated design object. This will be passed to \code{spsurvey::grts()} as the argument \code{design} and must exactly follow the structure required for that. If this is not provided, a design object will be automatically constructed with points allocated proportionally by stratum size.
#' @param panel_names A character vector of the names to assign to the panels for automatically generating a proportional-allocation design object. All values must be unique. This must be the same length and in the same order as \code{panel_sample_size}. Defaults to \code{c("One", "Two", "Three", "Four", "Five")}.
#' @param panel_sample_size A numeric vector of the number of base points to be drawn in each panel for automatically generating a proportional-allocation design object. This must be the same length and in the same order as \code{panel_names} UNLESS all panels have the same number of base points in which case \code{panel_sample_size} optionally may be a single numeric value despite \code{panel_names} containing more than one value. Defaults to \code{c(50)}.
#' @param points_min A numeric value of the minimum number of base points to allocate to a stratum within a panel regardless of its relative size for automatically generating a proportional-allocation design object. Defaults to \code{3}.
#' @param oversample_proportion A numeric value between 0 and 1 representing the minimum relative proportion of oversample points to allocate per stratum per panel using the formula \code{panel_sample_size * min.oversample_proportion} for automatically generating a proportional-allocation design object. Defaults to \code{0.25}.
#' @param oversample_min A numeric value representing the minimum number of oversample points to allocate per stratum per panel for automatically generating a proportional-allocation design object. This is only used if it is greater than \code{panel_sample_size * min.oversample_proportion}. Defaults to \code{3}.
#' @param seed_number An optional numeric value to be passed to \code{set.seed()} to make the output reproducible. Defaults to \code{NULL}.
#' @return  A Spatial Points Data Frame of the sampling locations with the fields \code{PLOTID}, \code{STRATUM}, \code{PANEL}, \code{IntPtWt} (initial point weight), \code{xcoord}, and \code{ycoord}
#' @export

draw <- function(design_name = "design",
                 strata_spdf = NULL,
                 stratum_field = "STRATUM",
                 sampleframe_spdf = NULL,
                 points_spdf = NULL,
                 strata_lut = NULL,
                 strata_lut_field = NULL,
                 design_object = NULL,
                 panel_names = c("One", "Two", "Three", "Four", "Five"),
                 panel_sample_size = 50,
                 points_min = 3,
                 oversample_proportion = 0.25,
                 oversample_min = 3,
                 seed_number = NULL
){
  ## Validity checks
  # Enough data
  if (is.null(strata_spdf) & is.null(sampleframe_spdf) & is.null(points_spdf)) {
    stop("Must provide at least a valid SPDF for one of strata_spdf, sampleframe_spdf, and points_spdf.")
  }

  # Correct format data checks
  if (!is.null(points_spdf)) {
    if (class(points_spdf) != "SpatialPointsDataFrame") {
      stop("points_spdf must be a Spatial Points Data Frame")
    }
  }
  if (!is.null(strata_spdf)) {
    if (class(strata_spdf) != "SpatialPolygonsDataFrame") {
      stop("strata_spdf must be a Spatial Polygons Data Frame")
    }
  }
  if (!is.null(sampleframe_spdf)) {
    if (class(sampleframe_spdf) != "SpatialPolygonsDataFrame") {
      stop("sampleframe_spdf must be a Spatial Polygons Data Frame")
    }
  }
  if (!is.null(strata_lut)) {
    if (class(strata_lut) != "data.frame") {
      stop("strata_lut must be a data frame")
    }
  }
  if (!is.null(strata_lut_field)) {
    if (!is.character(strata_lut_field) | length(strata_lut_field) != 1) {
      stop("strata_lut_field must be be a single character string")
    }
  }
  if (!is.null(seed_number)) {
    if (!is.numeric(seed_number)) {
      stop("seed_number must be a single numeric value")
    }
  }

  # Making sure that there's correlation
  if (!(stratum_field %in% c(names(strata_spdf),
                             names(sampleframe_spdf),
                             names(points_spdf)))) {
    stop("The string stratum_field must exactly match the name of a variable in at least one of the provided SPDFs.")
  }
  if (!(stratum_field %in% names(strata_lut)) & !is.null(strata_lut)) {
    stop("A variable must exist in strata_lut matching stratum_field.")
  }
  if (is.character(strata_lut_field) & !(strata_lut_field %in% names(strata_lut))) {
    if (is.null(strata_lut)) {
      message(paste0("There is no strata_lut, so the value ", strata_lut_field, " as strata_lut_field will be ignored"))
    } else {
      stop("The string strata_lut_field must exactly match the name of a variable in strata_lut.")
    }
  }
  if (is.character(strata_lut_field) & !(strata_lut_field %in% c(names(strata_spdf),
                                                                 names(sampleframe_spdf),
                                                                 names(points_spdf)))) {
    stop("The string strata_lut_field must exactly match the name of a variable in at least one of the provided SPDFs.")
  }

  ## Make sure the lookup table is a proper lookup table
  if (!is.null(strata_lut)) {
    strata_lut <- dplyr::unique(strata_lut[, c(stratum_field, strata_lut_field)])
    if (length(unique(strata_lut_field) != nrow(strata_lut))) {
      stop(paste0("The lookup table has non-unique relationships between one or more of the values in ",
                  stratum_field, " and ", strata_lut_field, "."))
    }
  }

  ## Spatial stuff
  # Get a population SPDF pulled together
  if (is.null(strata_spdf)) {
    population_spdf <- strata_spdf
    if (is.null(sampleframe_spdf)) {
      population_spdf <- restrict(spdf1 = frame.spdf,
                                  spdf2 = sampleframe_spdf)
    }
  } else if (is.null(sampleframe_spdf)) {
    population_spdf <- sampleframe_spdf
  } else {
    population_spdf <- NULL
  }

  # If at least one of strata_spdf and sampleframe_spdf existed, try to merge it with strata_lut
  if (!is.null(strata_lut) & (strata_lut_field %in% names(population@data))) {
    population_spdf <- merge(x = population_spdf,
                             y = strata_lut,
                             by = strata_lut_field)
  }

  # If the points exist
  if (is.null(points_spdf)) {
    # And the can be merged with strata_lut, then do it
    if (!is.null(strata_lut) & (strata_lut_field %in% names(points@data))) {
      points_spdf <- merge(x = points_spdf,
                           y = strata_lut,
                           by = strata_lut_field)
    }
    # If there weren't any other SPDFs provided
    if (is.null(population_spdf)) {
      # And the strata are already in the points
      if (stratum_field %in% names(points_spdf@data)) {
        population_spdf <- points_spdf
      } else {
        # Otherwise we've made it to here without strata, so we have to abort
        stop("If no strata or sample frame SPDFs are provided, then points_spdf must contain a variable exactly matching stratum_field or a variable exactly matching strata_lut_field.")
      }
      # Otherwise, if population_spdf already exists
    } else {
      # If population_spdf already has the stratum field, restrict by it and inherit that attribute
      if (stratum_field %in% names(population_spdf)) {
        population_spdf <- restrict(spdf1 = points_spdf,
                                    spdf2 = population,
                                    inherit = TRUE,
                                    inherit_field = stratum_field,
                                    bookend_inherit_field = TRUE)
        # If the points already have the stratum field, just restrict by population_spdf
      } else if (stratum_field %in% names(points_spdf@data)) {
        population_spdf <- restrict(spdf1 = points_spdf,
                                    spdf2 = population)
        # If there's STILL no stratum field, abort. I don't think this is possible with the error checking in place?
      } else {
        stop("Somehow between all provided SPDFs and lookup table, no strata were assigned. This is an ambiguous error for now, sorry!")
      }
    }

  }

  ## And now time to build the design object, if it isn't already provided
  if (is.null(design_object)) {
    design_object <- allocate_panels(spdf = population_spdf,
                                     stratum_field = stratum_field,
                                     panel_names = panel_names,
                                     points_min = points_min,
                                     oversample_proportion = oversample_proportion,
                                     oversample_min = oversample_min)
  }

  ## The use the design object to draw points
  output <- grts.aim(design_object = design_object,
                     design_name = design_name,
                     sp.object = population_spdf,
                     stratum_field = stratum_field,
                     seed_number = seed_number)

  return(output)
}
