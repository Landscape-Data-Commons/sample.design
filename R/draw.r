#' A quick wrapper to execute simple designs
#' @description Take some combination of strata, sample frame, and points with a design object or values for constructing one and produce a sample design as a Spatial Points Data Frame. Automatic generation of a design object is done using allocate.panels().
#' @param design.name A character string of the name to use for the design, often the name of the project.
#' @param strata.spdf An optional Spatial Polygons Data Frame. Must either contain an attribute which exactly matches the string \code{stratum.field} and holds the stratum identities OR, if using \code{strata.lut}, an attribute which exactly matches \code{strata.lut.field} and contains values found in \code{strata.lut$strata.lut.field}. If \code{points.spdf} is provided, this will be used to attribute it with strata identities.
#' @param sampleframe.spdf An optional Spatial Polygons Data Frame. Necessary only to restrict \code{strata.spdf} or \code{points.spdf} if they are not already limited to the extent of the intended sample frame. If it is the only SPDF provided because there is no stratification scheme, it still must contain an attribute which exactly matches the string \code{stratum.field}, although there will only be one unique value in it.
#' @param points.spdf An optional Spatial Points Data Frame. This is only necessary if the draw should be from within a set of discrete points and not from polygon areas. If neither \code{strata.spdf} nor \code{sampleframe.spdf} with strata identity information is provided, this must either contain an attribute which exactly matches the string \code{stratum.field} and holds the stratum identities OR, if using \code{strata.lut}, an attribute which exactly matches \code{strata.lut.field} and contains values found in \code{strata.lut$strata.lut.field}.
#' @param strata.lut An optional data frame. If provided, this must contain a variable with a name that exactly matches the string provided in \code{stratum.field} and a variable with a name that exactly matches the string provided in \code{strata.lut.field}. There must be unique relationships between the values in the two fields, i.e. a value in \code{strata.lut$strata.lut.field} CANNOT correspond with more than one value in \code{strata.lut$stratum.field}. The values in \code{strata.lut$strata.lut.field} must match values in the \code{data$strata.lut.field} slot of the provided SPDFs.
#' @param strata.lut.field An optional string specifying a column/attribute/variable name. If \code{strata.lut} is provided, then a variable with a name exactly matching this string must exist in both \code{strata.lut} and one of the SPDFs provided.
#' @param design.object An optional specially-constructed list to use instead of an automatically-generated design object. This will be passed to \code{spsurvey::grts()} as the argument \code{design} and must exactly follow the structure required for that. If this is not provided, a design object will be automatically constructed with points allocated proportionally by stratum size.
#' @param panel.names A character vector of the names to assign to the panels for automatically generating a proportional-allocation design object. All values must be unique. This must be the same length and in the same order as \code{panel.sample.size}. Defaults to \code{c("One", "Two", "Three", "Four", "Five")}.
#' @param panel.sample.size A numeric vector of the number of base points to be drawn in each panel for automatically generating a proportional-allocation design object. This must be the same length and in the same order as \code{panel.names} UNLESS all panels have the same number of base points in which case \code{panel.sample.size} optionally may be a single numeric value despite \code{panel.names} containing more than one value. Defaults to \code{c(50)}.
#' @param points.min A numeric value of the minimum number of base points to allocate to a stratum within a panel regardless of its relative size for automatically generating a proportional-allocation design object. Defaults to \code{3}.
#' @param oversample.proportion A numeric value between 0 and 1 representing the minimum relative proportion of oversample points to allocate per stratum per panel using the formula \code{panel.sample.size * min.oversample.proportion} for automatically generating a proportional-allocation design object. Defaults to \code{0.25}.
#' @param oversample.min A numeric value representing the minimum number of oversample points to allocate per stratum per panel for automatically generating a proportional-allocation design object. This is only used if it is greater than \code{panel.sample.size * min.oversample.proportion}. Defaults to \code{3}.
#' @param seed.number An optional numeric value to be passed to \code{set.seed()} to make the output reproducible. Defaults to \code{NULL}.
#' @return  A Spatial Points Data Frame of the sampling locations with the fields \code{PLOTID}, \code{STRATUM}, \code{PANEL}, \code{IntPtWt} (initial point weight), \code{xcoord}, and \code{ycoord}
#' @export

draw <- function(design.name = "design",
                 strata.spdf = NULL,
                 stratum.field = "STRATUM",
                 sampleframe.spdf = NULL,
                 points.spdf = NULL,
                 strata.lut = NULL,
                 strata.lut.field = NULL,
                 design.object = NULL,
                 panel.names = c("One", "Two", "Three", "Four", "Five"),
                 panel.sample.size = 50,
                 points.min = 3,
                 oversample.proportion = 0.25,
                 oversample.min = 3,
                 seed.number = NULL
){
  ## Validity checks
  # Enough data
  if (is.null(strata.spdf) & is.null(sampleframe.spdf) & is.null(points.spdf)) {
    stop("Must provide at least a valid SPDF for one of strata.spdf, sampleframe.spdf, and points.spdf.")
  }

  # Correct format data checks
  if (!is.null(points.spdf) & class(points.spdf) != "SpatialPointsDataFrame") {
    stop("points.spdf must be a Spatial Points Data Frame")
  }
  if (!is.null(strata.spdf) & class(strata.spdf) != "SpatialPolygonsDataFrame") {
    stop("strata.spdf must be a Spatial Polygons Data Frame")
  }
  if (!is.null(sampleframe.spdf) & class(sampleframe.spdf) != "SpatialPolygonsDataFrame") {
    stop("sampleframe.spdf must be a Spatial Polygons Data Frame")
  }
  if (!is.null(strata.lut) & class(strata.lut) != "data.frame") {
    stop("strata.lut must be a data frame")
  }
  if (!is.null(strata.lut.field) & (!is.character(strata.lut.field) | length(strata.lut.field) != 1)) {
    stop("strata.lut.field must be be a single character string")
  }
  if (!is.null(seed.number) & !is.numeric(seed.number)) {
    stop("seed.number must be a single numeric value")
  }

  # Making sure that there's correlation
  if (!(stratum.field %in% c(names(strata.spdf), names(sampleframe.spdf), names(points.spdf)))) {
    stop("The string stratum.field must exactly match the name of a variable in at least one of the provided SPDFs.")
  }
  if (!(stratum.field %in% names(strata.lut)) & !is.null(strata.lut)) {
    stop("A variable must exist in strata.lut matching stratum.field.")
  }
  if(is.character(strata.lut.field) & !(strata.lut.field %in% names(strata.lut))) {
    if (is.null(strata.lut)) {
      message(paste0("There is no strata.lut, so the value ", strata.lut.field, " as strata.lut.field will be ignored"))
    } else {
      stop("The string strata.lut.field must exactly match the name of a variable in strata.lut.")
    }
  }
  if(is.character(strata.lut.field) & !(strata.lut.field %in% c(names(strata.spdf), names(sampleframe.spdf), names(points.spdf)))) {
    stop("The string strata.lut.field must exactly match the name of a variable in at least one of the provided SPDFs.")
  }

  ## Make sure the lookup table is a proper lookup table
  if (!is.null(strata.lut)) {
    strata.lut <- dplyr::unique(strata.lut[, c(stratum.field, strata.lut.field)])
    if (length(unique(strata.lut.field) != nrow(strata.lut))) {
      stop(paste0("The lookup table has non-unique relationships between one or more of the values in ", stratum.field, " and ", strata.lut.field, "."))
    }
  }

  ## Spatial stuff
  # Get a population SPDF pulled together
  if (is.null(strata.spdf)) {
    population.spdf <- strata.spdf
    if (is.null(sampleframe.spdf)) {
      population.spdf <- restrict(spdf1 = frame.spdf,
                                  spdf2 = sampleframe.spdf)
    }
  } else if (is.null(sampleframe.spdf)) {
    population.spdf <- sampleframe.spdf
  } else {
    population.spdf <- NULL
  }

  # If at least one of strata.spdf and sampleframe.spdf existed, try to merge it with strata.lut
  if (!is.null(strata.lut) & (strata.lut.field %in% names(population@data))) {
    population.spdf <- merge(x = population.spdf,
                             y = strata.lut,
                             by = strata.lut.field)
  }

  # If the points exist
  if (is.null(points.spdf)) {
    # And the can be merged with strata.lut, then do it
    if (!is.null(strata.lut) & (strata.lut.field %in% names(points@data))) {
      points.spdf <- merge(x = points.spdf,
                           y = strata.lut,
                           by = strata.lut.field)
    }
    # If there weren't any other SPDFs provided
    if (is.null(population.spdf)) {
      # And the strata are already in the points
      if (stratum.field %in% names(points.spdf@data)) {
        population.spdf <- points.spdf
      } else {
        # Otherwise we've made it to here without strata, so we have to abort
        stop("If no strata or sample frame SPDFs are provided, then points.spdf must contain a variable exactly matching stratum.field or a variable exactly matching strata.lut.field.")
      }
    # Otherwise, if population.spdf already exists
    } else {
      # If population.spdf already has the stratum field, restrict by it and inherit that attribute
      if (stratum.field %in% names(population.spdf)) {
        population.spdf <- restrict(spdf1 = points.spdf,
                                    spdf2 = population,
                                    inherit = T,
                                    inherit.field = stratum.field,
                                    bookend.inherit.field = T)
      # If the points already have the stratum field, just restrict by population.spdf
      } else if (stratum.field %in% names(points.spdf@data)) {
        population.spdf <- restrict(spdf1 = points.spdf,
                                    spdf2 = population)
      # If there's STILL no stratum field, abort. I don't think this is possible with the error checking in place?
      } else {
        stop("Somehow between all procided SPDFs and lookup table, no strata were assigned. This is an ambiguous error for now, sorry!")
      }
    }

  }

  ## And now time to build the design object, if it isn't already provided
  if (is.null(design.object)) {
    design.object <- allocate.panels(spdf = population.spdf,
                                     stratum.field = stratum.field,
                                     panel.names = panel.names,
                                     points.min = points.min,
                                     oversample.proportion = oversample.proportion,
                                     oversample.min = oversample.min)
  }

  ## The use the design object to draw points
  output <- grts.aim(design.object = design.object,
                     design.name = design.name,
                     sp.object = population.spdf,
                     stratum.field = stratum.field,
                     seed.number = seed.number)

  return(output)
}
