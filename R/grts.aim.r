#' Draw spatially balanced points with \code{spsurvey::grts()}, output in AIM format
#'
#' @description A wrapper for \code{spsurvey::grts()} that modifies the output SPDF to match the format expected by other parts of the AIM workflow, including distribution to field crews and populating design databases.
#' @param design.object A list of lists structured and named according to the requirements of \code{spsurvey::grts()}. Can be automatically constructed with \code{allocate.panels()}.
#' @param design.name A character string of the name to use for the design, often the name of the project.
#' @param sp.object If using an SPDF instead of a shapefile, an SPDF with a variable name matching \code{stratum.field} which contains values matching the strata names in \code{design.object}. Defaults to \code{NULL}.
#' @param in.shape If using a shapefile instead of an SPDF, a string representing the filepath to a shapefile (without file extension) to pass to \code{spsurvey::grts()} as its argument \code{in.shape}. Defaults to \code{NULL}.
#' @param stratum.field A character string representing the name of the variable in either \code{sp.object} or \code{in.shape} containing the strata names matching those in \code{design.object}. Defaults to \code{"STRATUM"}.
#' @param seed.number An optional numeric value to be passed to \code{set.seed()} to make the output reproducible. Defaults to \code{NULL}.
#' @return A Spatial Points Data Frame of the sampling locations with the fields \code{PLOTID}, \code{STRATUM}, \code{PANEL}, \code{IntPtWt} (initial point weight), \code{xcoord}, and \code{ycoord}
#' @export
grts.aim <- function(design.object,
                     design.name = "Design name",
                     src.frame = "sp.object",
                     sp.object = NULL,
                     in.shape = NULL,
                     stratum.field = "STRATUM",
                     seed.number = NULL
){
  if (!is.null(seed.number)) {
    set.seed(seed.number)
  }
  if (!is.null(sp.object)) {
    src.frame <- "sp.object"
    if (!(stratum.field %in% names(sp.object))) {
      stop("The variable stratum.field was not found in sp.object. Check case and spelling.")
    }
  } else if (!is.null(in.shape)) {
    src.frame <- "shapefile"
    in.shape <- in.shape %>% stringr::str_replace(pattern = "\\.(shp)|(shp)$", replacement = "")
  } else {
    stop("Provide either an SPDF as sp.object or a filepath to a shapefile as in.shape.")
  }
  if (src.frame == "sp.object" & is.null(sp.object)) {
    stop("Please provide an SPDF as sp.object.")
  }
  if (src.frame == "shapefile" & is.null(in.shape)) {
    stop("Please provide a filepath to a shapefile as in.shape.")
  }

  ## Invoke spsurvey::grts() first
  sample.sites <- spsurvey::grts(design = design.object,
                                 DesignID = design.name,
                                 type.frame = "area",
                                 src.frame = src.frame,
                                 sp.object = sp.object,
                                 in.shape = in.shape,
                                 stratum = stratum.field,
                                 shapefile = FALSE
  )

  ## Assign projection info to the sample sites SPDF
  if (is.null(sp.object)) {
    proj4string(sample.sites) <- sp.object@proj4string
  } else {
    proj4string(sample.sites) <- rgdal::readOGR(dsn = stringr::str_replace(in.shape, pattern = "/([a-Z]|[0-9]){1,256}$", replacement = ""),
                                                layer = stringr::str_extract(in.shape, pattern = "/([a-Z]|[0-9]){1,256}$")) %>% .@proj4string
    # proj4string(sample.sites) <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
  }
  ## Reproject the sample sites to Geographic DD NAD83
  sample.sites <- spTransform(sample.sites, CRS("+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs"))
  ## update the X and Y coordinate values
  sample.sites$xcoord <- coordinates(sample.sites)[, 1]
  sample.sites$ycoord <- coordinates(sample.sites)[, 2]

  ## Dropping the extra fields from the master sample, just keeping the ones specific to the draw plus the master sample ID
  fields.relevant <- c("siteID", "stratum", "panel", "wgt", "xcoord", "ycoord")
  sample.sites@data <- sample.sites@data[, fields.relevant]
  names(sample.sites@data) <- c("PLOTID", "STRATUM", "PANEL", "IntPtWt", "xcoord", "ycoord")

  ## Last step is to rename the oversample panels so that they're broken up into the years instead of being an oversample chunk
  panel.names <- unique(sample.sites@data$PANEL[!(sample.sites@data$PANEL %in% "OverSamp")])
  oversample.df <- sample.sites@data[sample.sites@data$PANEL == "OverSamp",] %>% group_by(STRATUM) %>% summarize(oversample.pts.per.panel = floor(n()/length(panel.names)), total.oversample.drawn = n())
  oversample.panels <- list()
  for (s in oversample.df$STRATUM) {
    oversample.count.per.panel <- oversample.df$oversample.pts.per.panel[oversample.df$STRATUM == s]
    oversample.drawn <- oversample.df$total.oversample.drawn[oversample.df$STRATUM == s]
    oversample.panels.current <- c()
    for (n in panel.names) {
      oversample.panels.current <- c(oversample.panels.current,
                                     rep(x = paste("OverSample", n),
                                         times = oversample.count.per.panel))
    }
    ## If there was a difference introduced by rounding, here's the place to put the extra points into the last panel
    if (length(oversample.panels.current) < oversample.drawn) {
      oversample.panels.current <- c(oversample.panels.current,
                                     rep(x = paste("OverSample", last(panel.names)),
                                         times = oversample.drawn - length(oversample.panels.current)))
    }
    oversample.panels[[s]] <- oversample.panels.current
  }

  for (s in names(oversample.panels)) {
    sample.sites@data$PANEL[sample.sites@data$PANEL == "OverSamp" & sample.sites@data$STRATUM == s] <- oversample.panels[[s]]
  }

  ## Rename the plots with the strata
  sample.sites@data$PLOTID <- paste0(sample.sites@data$STRATUM, sample.sites@data$PLOTID %>% stringr::str_extract(pattern = "-[0-9]{1,4}$"))

  return(sample.sites)
}

