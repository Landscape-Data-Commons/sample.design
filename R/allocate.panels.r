#' Create a design object programmatically to use with \code{spsurvey::grts()}.
#'
#' @description Creates a structured list that meets the requirements for a design object for \code{spsurvey::grts()} distributing points proportionally to strata based on their relative areas/abundances.
#' @param spdf SpatialPointsDataFrame or data frame containing the possible sampling points being used or a SpatialPolygonsDataFrame of polygons describing the strata boundaries.
#' @param stratum.field A character string of the name of the variable in \code{points} containing the stratum identities of the points. Defaults to \code{"stratum"}.
#' @param panel.names A character vector of the names to assign to the panels. All values must be unique. This must be the same length and in the same order as \code{panel.sample.size}. Defaults to \code{c("Year1", "Year2", "Year3", "Year4", "Year5")}.
#' @param panel.sample.size A numeric vector of the number of base points to be drawn in each panel. This must be the same length and in the same order as \code{panel.names} UNLESS all panels have the same number of base points in which case \code{panel.sample.size} optionally may be a single numeric value despite \code{panel.names} containing more than one value. Defaults to \code{c(50)}.
#' @param points.min A numeric value of the minimum number of base points to allocate to a stratum within a panel regardless of its relative size. Defaults to \code{3}.
#' @param oversample.proportion A numeric value between 0 and 1 representing the minimum relative proportion of oversample points to allocate per stratum per panel using the formula \code{panel.sample.size * min.oversample.proportion}. Defaults to \code{0.25}.
#' @param oversample.min A numeric value representing the minimum number of oversample points to allocate per stratum per panel. This is only used if it is greater than \code{panel.sample.size * min.oversample.proportion}. Defaults to \code{3}.
#' @return A named list of named lists conforming to the requirements for the design object for \code{spsurvey::grts()}.
#' @export

allocate.panels <- function(spdf,
                            stratum.field = "stratum",
                            panel.names = c("Year1", "Year2", "Year3", "Year4", "Year5"),
                            panel.sample.size = 50,
                            points.min = 3,
                            oversample.proportion = 0.25,
                            oversample.min = 3
){
  ## Error checking
  if (length(panel.sample.size) != (length(panel.names)) & length(panel.sample.size) != 1) {
    stop("Error: panel.sample.size either needs to be a single numeric value or a vector of numeric values with a length equal to panel.number")
  }

  ## Sanitization
  if (class(spdf) == "SpatialPolygonsDataFrame") {
    df <- dplyr::mutate(add.area(spdf)@data, AREA = AREA.HA)
  }

  if (class(spdf) == "SpatialPointsDataFrame") {
    df <- spdf@data
  }

  if (!(stratum.field %in% names(df))) {
    stop("Error: Couldn't find the specified stratum field in the supplied points or polygons. Check spelling and case.")
  } else {
    df$STRATUM <- df[, stratum.field]
  }

  ## Remove all points or areas without strata assigned
  df <- dplyr::filter(df, !is.na(STRATUM))

  ## Create a data frame of strata and "area"
  if ("AREA" %in% names(df)) {
    workingframe <- dplyr::summarize(dplyr::group_by(df, STRATUM), AREA = sum(AREA))
  } else (
    workingframe <- dplyr::summarize(dplyr::group_by(df, STRATUM), AREA = n())
  )

  ## Prep the working frame for the mutates
  workingframe <- stratum.sizes
  workingframe$min <- points.min
  workingframe$remainder <- panel.sample.size - nrow(workingframe)*points.min
  workingframe$oversample.proportion <- oversample.proportion
  workingframe$oversample.min <- oversample.min
  workingframe$panel.number <- panel.number
  print(workingframe)

  ## Create all the support values then the list that goes into the design object for each stratum
  workingframe <- workingframe %>% dplyr::mutate(PROPORTION = AREA/sum(AREA)) %>%
    dplyr::mutate(PER.PANEL.BASE = round(PROPORTION*remainder) + min) %>%
    dplyr::mutate(PER.PANEL.OVERSAMPLE = pmax(PER.PANEL.BASE*oversample.proportion, oversample.min) %>% ceiling()) %>%
    dplyr::mutate(TOTAL.OVERSAMPLE = panel.number*PER.PANEL.OVERSAMPLE)

  ## Create the output design object list.
  output <- lapply(workingframe$STRATUM,
                   function(X, workingframe, panel.names) {
                     list(panel = rep(workingframe$PER.PANEL.BASE[workingframe$STRATUM == X],
                                      times = workingframe$panel.number[workingframe$STRATUM == X]) %>%
                            setNames(panel.names),
                          seltype = "Equal",
                          over = workingframe$TOTAL.OVERSAMPLE[workingframe$STRATUM == X])
                   },
                   workingframe = workingframe,
                   panel.names = panel.names) %>%
    setNames(workingframe$STRATUM)

  return(output)
}

#' Use a data frame of strata and point counts to build a design object for \code{spsurvey::grts()}
#' @param dataframe Data frame. This must have at least a variable for the strata identities and one variable for each panel, e.g. "Stratum", "Year1", "Year2", "Year3", where each row is a stratum and the number of points desired in each panel for that stratum.
#' @param stratum.field Character string. This must exactly match the name of the variable in \code{dataframe} that contains the stratum identities. Defaults to \code{"stratum"}.
#' @param panel.names Optional vector of character strings. Necessary if the data frame given in \code{dataframe} has additional fields beyond the ones for the strata, panels, and oversample. This character vector specifies which variables correspond to panels and must match those names exactly.
#' @param oversample.field Optional character string. If used, this must exactly match the name of the variable in \code{dataframe} that contains the TOTAL number of oversample points desired for the stratum in the design. If not provided or specified as \code{NULL}, oversample point counts will be calculated using \code{oversample.proportion} and \code{oversample.min}. Defaults to \code{NULL}.
#' @param oversample.proportion Optional numeric value. If not providing an oversample point count in a variable specified by \code{oversample.field}, this must be between 0 and 1, representing the minimum relative proportion of oversample points to allocate per stratum per panel using the formula \code{panel point count * min.oversample.proportion}. Defaults to \code{0.25}.
#' @param oversample.min Optional numeric value. If not providing an oversample point count in a variable specified by \code{oversample.field}, this must be a positive integer, representing the minimum number of oversample points to allocate per stratum per panel. This is only used if it is greater than \code{ppanel point count * min.oversample.proportion}. Defaults to \code{3}.
#' @export
read.panels <- function(dataframe,
                        stratum.field = "stratum",
                        panel.names = NULL,
                        oversample.field = NULL,
                        oversample.proportion = 0.25,
                        oversample.min = 3) {

  # For each stratum, make the design list specifying number of sample points per panel and the oversample count
  design <- lapply(X = dataframe[[stratum.field]],
                   FUN =   function(X,
                                    dataframe,
                                    stratum.field,
                                    panel.names,
                                    oversample.field,
                                    oversample.proportion,
                                    oversample.min){
                     # Create a vector of panel names from the variable names if there isn't one yet
                     if (is.null(panel.names)) {
                       panel.names <- names(dataframe)[!(names(dataframe) %in% c(stratum.field, oversample.field))]
                     }
                     # Get just the relevant variables in the current data frame for the stratum
                     dataframe.current <- dataframe[dataframe[[stratum.field]] == X, names(dataframe)[names(dataframe) %in% c(stratum.field, panel.names, oversample.field)]]
                     # Pull the panel values and create a named vector from them
                     panel <- setNames(unlist(lapply(X = panel.names,
                                                     FUN = function(X, df){return(df[,X])},
                                                     df = dataframe.current)),
                                       panel.names)

                     # E pull the oversample value or calculate it
                     if (is.null(oversample.field)) {
                       over <- dataframe.current[[oversample.field]]
                     } else {
                       # For each panel, get the number of oversample points, then sum the vector
                       over <- sum(unlist(lapply(X = panel,
                                                 FUN = function(X, oversample.proportion, oversample.min){
                                                   return(max(round(X * oversample.proportion), oversample.min))
                                                 },
                                                 oversample.proportion = oversample.proportion,
                                                 oversample.min = oversample.min)))
                     }

                     # Return the list for this stratum
                     return(list(panel = panel, seltype = "Equal", over = over))
                   },
                   dataframe = dataframe,
                   stratum.field = stratum.field,
                   panel.names = panel.names,
                   oversample.field = oversample.field,
                   oversample.proportion = oversample.proportion,
                   oversample.min = oversample.min)

  # Name the lists for each stratum with the stratum name
  design <- setNames(design, dataframe[[stratum.field]])

  return(design)
}

#' Convert a GRTS design list to a data frame
#' @export
design.dataframe <- function(design){
  strata <- names(design)

  lapply(X = strata,
         FUN = function(X, design){
           df <- data.frame(stratum = X,
                            panel = names(design[[X]][["panel"]]), count = unname(design[[X]][["panel"]]))
           df.wide <- tidyr::spread(data = df,
                                    key = panel,
                                    value = count)
           df.wide$oversample <- design[[X]][["over"]]
         }, design = design)
}
