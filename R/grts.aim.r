#' Draw spatially balanced points with \code{spsurvey::grts()}, output in AIM format
#'
#' @description A wrapper for \code{spsurvey::grts()} that modifies the output SPDF to match the format expected by other parts of the AIM workflow, including distribution to field crews and populating design databases.
#' @param design_object A list of lists structured and named according to the requirements of \code{spsurvey::grts()}. Can be automatically constructed with \code{allocate.panels()}.
#' @param design_name A character string of the name to use for the design, often the name of the project.
#' @param source_frame Character string. To be passed to \code{spsurvey::grts()} as the argument \code{src.frame}; see that function's documentation for further details. Defaults to \code{"sp.object"}.
#' @param sp_object If using an SPDF instead of a shapefile, an SPDF with a variable name matching \code{stratum_field} which contains values matching the strata names in \code{design_object}. Defaults to \code{NULL}.
#' @param in_shape If using a shapefile instead of an SPDF, a string representing the filepath to a shapefile (without file extension) to pass to \code{spsurvey::grts()} as its argument \code{in_shape}. Defaults to \code{NULL}.
#' @param stratum_field A character string representing the name of the variable in either \code{sp_object} or \code{in_shape} containing the strata names matching those in \code{design_object}. Defaults to \code{"STRATUM"}.
#' @param seed_number An optional numeric value to be passed to \code{set.seed()} to make the output reproducible. Defaults to \code{NULL}.
#' @return A Spatial Points Data Frame of the sampling locations with the fields \code{PLOTID}, \code{STRATUM}, \code{PANEL}, \code{IntPtWt} (initial point weight), \code{xcoord}, and \code{ycoord}
#' @export
grts_aim <- function(design_object,
                     design_name = "Design name",
                     source_frame = "sp.object",
                     sp_object = NULL,
                     in_shape = NULL,
                     stratum_field = "STRATUM",
                     seed_number = NULL
){
  if (!is.null(seed_number)) {
    set.seed(seed_number)
  }
  if (!is.null(sp_object)) {
    source_frame <- "sp.object"
    if (!(stratum_field %in% names(sp_object))) {
      stop("The variable stratum_field was not found in sp_object. Check case and spelling.")
    }
  } else if (!is.null(in_shape)) {
    source_frame <- "shapefile"
    in_shape <- gsub(in_shape, pattern = "\\.(shp)|(shp)$", replacement = "")
  } else {
    stop("Provide either an SPDF as sp_object or a filepath to a shapefile as in_shape.")
  }
  if (source_frame == "sp.object" & is.null(sp_object)) {
    stop("Please provide an SPDF as sp_object.")
  }
  if (source_frame == "shapefile" & is.null(in_shape)) {
    stop("Please provide a filepath to a shapefile as in_shape.")
  }

  ## Invoke spsurvey::grts() first
  sample_sites <- spsurvey::grts(design = design_object,
                                 DesignID = design_name,
                                 type.frame = "area",
                                 src.frame = source_frame,
                                 sp.object = sp_object,
                                 in.shape = in_shape,
                                 stratum = stratum_field,
                                 shapefile = FALSE
  )

  ## Assign projection info to the sample sites SPDF
  if (!is.null(sp_object)) {
    sp::proj4string(sample_sites) <- sp_object@proj4string
  } else {
    sp::proj4string(sample_sites) <- rgdal::readOGR(dsn = gsub(in_shape,
                                                               pattern = "/([A-z]|[0-9])+$",
                                                               replacement = ""),
                                                    layer = gsub(in_shape,
                                                                 pattern = "/([A-z]|[0-9])+$"))@proj4string
  }
  ## Reproject the sample sites to Geographic DD NAD83
  sample_sites <- sp::spTransform(sample_sites, sp::CRS("+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs"))
  ## update the X and Y coordinate values
  sample_sites[["xcoord"]] <- sp::coordinates(sample_sites)[, 1]
  sample_sites[["ycoord"]] <- sp::coordinates(sample_sites)[, 2]

  ## Dropping the extra fields from the master sample, just keeping the ones specific to the draw plus the master sample ID
  fields_relevant <- c("siteID", "stratum", "panel", "wgt", "xcoord", "ycoord")
  sample_sites@data <- sample_sites@data[, fields_relevant]
  names(sample_sites@data) <- c("PLOTID", "STRATUM", "PANEL", "IntPtWt", "xcoord", "ycoord")

  # Change "OverSamp" to "OverSample + [STRATUM]"
  oversample_panel_names <- paste(sample_sites@data[["PANEL"]][sample_sites@data[["PANEL"]] == "OverSamp"], "Oversample")
  sample_sites@data[["PANEL"]][sample_sites@data[["PANEL"]] == "OverSamp"] <- oversample_panel_names

  ## Rename the plots with the strata
  sample_sites@data[["PLOTID"]] <- paste0(sample_sites@data[["STRATUM"]], stringr::str_extract(string = sample_sites@data[["PLOTID"]], pattern = "-[0-9]{1,4}$"))

  return(sample_sites)
}

