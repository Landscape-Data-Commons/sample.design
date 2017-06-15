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
    df <- spdf %>% area.add() %>% dplyr::mutate(AREA = AREA.HA) %>% .@data
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
    workingframe <- df %>% dplyr::group_by(STRATUM) %>% dplyr::summarize(AREA = sum(AREA))
  } else (
    workingframe <- df %>% dplyr::group_by(STRATUM) %>% dplyr::summarize(AREA = n())
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
