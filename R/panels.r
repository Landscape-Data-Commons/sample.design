#' Create a design object programmatically to use with \code{spsurvey::grts()}.
#'
#' @description Creates a structured list that meets the requirements for a design object for \code{spsurvey::grts()} distributing points proportionally to strata based on their relative areas/abundances.
#' @param polygons Polygon sf object. The polygons describing the strata boundaries. Must contain a variable containing the stratum identities.
#' @param stratum_field A character string of the name of the variable in \code{points} containing the stratum identities of the polygons. Defaults to \code{"stratum"}.
#' @param panel_names A character vector of the names to assign to the panels. All values must be unique. This must be the same length and in the same order as \code{panel_sample_size}. Defaults to \code{c("1", "2", "3", "4", "5")}.
#' @param panel_sample_size A numeric vector of the number of base points to be drawn in each panel. This must be the same length and in the same order as \code{panel_names} UNLESS all panels have the same number of base points in which case \code{panel_sample_size} optionally may be a single numeric value despite \code{panel_names} containing more than one value. Defaults to \code{c(50)}.
#' @param points_min A numeric value of the minimum number of base points to allocate to a stratum within a panel regardless of its relative size. Defaults to \code{3}.
#' @param oversample_proportion A numeric value between 0 and 1 representing the minimum relative proportion of oversample points to allocate per stratum per panel using the formula \code{panel_sample_size * min.oversample_proportion}. Defaults to \code{0.25}.
#' @param oversample_min A numeric value representing the minimum number of oversample points to allocate per stratum per panel. This is only used if it is greater than \code{panel_sample_size * min.oversample_proportion}. Defaults to \code{3}.
#' @return A named list of named lists conforming to the requirements for the design object for \code{spsurvey::grts()}.
#' @export

allocate_panels <- function(polygons,
                            stratum_field = "stratum",
                            panel_names = c("1", "2", "3", "4", "5"),
                            panel_sample_size = 50,
                            points_min = 3,
                            oversample_proportion = 0.25,
                            oversample_min = 3
){
  ## Error checking
  if (length(panel_sample_size) != (length(panel_names)) & length(panel_sample_size) != 1) {
    stop("Error: panel_sample_size either needs to be a single numeric value or a vector of numeric values with a length equal to the number of panels.")
  }

  ## Sanitization
  if (!("sf") %in% class(polygons)) {
    stop("polygons must be a polygon sf object.")
  }
  if (!all(sf::st_geometry_type(polygons) %in% c("POLYGON", "MULTIPOLYGON"))) {
    stop("polygons must be a polygon sf object.")
  }
  
  # Make sure that there's area information!
  polygons <- add_area(polygons)
  polygons[["AREA"]] <- polygons[["AREA.HA"]]
  
  df <- sf::st_drop_geometry(polygons)

  if (!(stratum_field %in% names(df))) {
    stop("Error: Couldn't find the specified stratum field in the supplied points or polygons. Check spelling and case.")
  } else {
    df[["STRATUM"]] <- df[, stratum_field]
  }

  ## Remove all polygons without strata assigned
  df <- df[!is.na(df[["STRATUM"]]), ]

  # How many points minimum will be needed?
  required_base <- points_min * length(unique(df[["STRATUM"]]))
  required_oversample <- round(max(required_base * oversample_proportion, oversample_min * length(unique(df[["STRATUM"]]))))

  ## Create a data frame of strata and "area"
  workingframe <- dplyr::summarize(dplyr::group_by(df, STRATUM),
                                   AREA = sum(AREA))

  # After the minimum points are allocated, how many remain to be allocated?
  remainder <- panel_sample_size - nrow(workingframe) * points_min
  # How many panels are there?
  panel_count <- length(panel_names)

  ## Create all the support values then the list that goes into the design object for each stratum
  workingframe[["PROPORTION"]] <- workingframe[["AREA"]] / sum(workingframe[["AREA"]])
  workingframe[["PER.PANEL.BASE"]] <- round(workingframe[["PROPORTION"]] * remainder) + points_min
  workingframe[["PER.PANEL.OVERSAMPLE"]] <- ceiling(pmax(workingframe[["PER.PANEL.BASE"]] * oversample_proportion, oversample_min))
  workingframe[["TOTAL.OVERSAMPLE"]] <- workingframe[["PER.PANEL.OVERSAMPLE"]] * panel_count

  if (any(workingframe[["PER.PANEL.BASE"]] < 0)) {
    stop("One or more strata ended up with a negative number of base points allocated. This usually happened if the minimum number of points per stratum is high enough that the minimum number of points exceeds the total points per panel. Check to make sure you aren't asking for too many points.")
  }
  if (any(workingframe[["TOTAL.OVERSAMPLE"]] < 0)) {
    stop("One or more strata ended up with a negative number of oversample points allocated. Check to make sure you aren't asking for too many points.")
  }

  ## Create the output design object list.
  output <- lapply(split(workingframe, workingframe[["STRATUM"]]),
                   panel_names = panel_names,
                   panel_count = panel_count,
                   function(X, panel_names, panel_count) {
                     # Just for clarity because X isn't obvious
                     df <- X
                     # Make the list. It's made of a named vector of panel sizes in base point count
                     list(panel = unlist(stats::setNames(rep(df[1, "PER.PANEL.BASE"],
                                                             times = panel_count),
                                                         panel_names)),
                          # The selection type (always equal here)
                          seltype = "Equal",
                          # And total oversample points
                          over = unname(unlist(df[1, "TOTAL.OVERSAMPLE"])))
                   })

  # The list needs to be named by stratum
  output <- stats::setNames(output,
                            workingframe[["STRATUM"]])

  return(output)
}

#' Use a data frame of strata and point counts to build a design object for \code{spsurvey::grts()}
#' @param dataframe Data frame. This must have at least a variable for the strata identities and one variable for each panel, e.g. "Stratum", "Year1", "Year2", "Year3", where each row is a stratum and the number of points desired in each panel for that stratum.
#' @param stratum_field Character string. This must exactly match the name of the variable in \code{dataframe} that contains the stratum identities. Defaults to \code{"stratum"}.
#' @param panel_names Optional vector of character strings. Necessary if the data frame given in \code{dataframe} has additional fields beyond the ones for the strata, panels, and oversample. This character vector specifies which variables correspond to panels and must match those names exactly.
#' @param oversample_field Optional character string. If used, this must exactly match the name of the variable in \code{dataframe} that contains the TOTAL number of oversample points desired for the stratum in the design. If not provided or specified as \code{NULL}, oversample point counts will be calculated using \code{oversample_proportion} and \code{oversample_min}. Defaults to \code{NULL}.
#' @param oversample_proportion Optional numeric value. If not providing an oversample point count in a variable specified by \code{oversample_field}, this must be between 0 and 1, representing the minimum relative proportion of oversample points to allocate per stratum per panel using the formula \code{panel point count * min.oversample_proportion}. Defaults to \code{0.25}.
#' @param oversample_min Optional numeric value. If not providing an oversample point count in a variable specified by \code{oversample_field}, this must be a positive integer, representing the minimum number of oversample points to allocate per stratum per panel. This is only used if it is greater than \code{ppanel point count * min.oversample_proportion}. Defaults to \code{3}.
#' @export
read_panels <- function(dataframe,
                        stratum_field = "stratum",
                        panel_names = NULL,
                        oversample_field = NULL,
                        oversample_proportion = 0.25,
                        oversample_min = 3) {

  if (class(dataframe) == "character") {
    if (!file.exists(dataframe) | grepl(dataframe, pattern = "\\.csv$", ignore.case = TRUE)) {
      stop("dataframe must either be a data frame or a filepath to a .CSV file")
    }
    dataframe <- read.csv(dataframe,
                          stringsAsFactors = FALSE)
  }
  # Create a vector of panel names from the variable names if there isn't one yet
  # This assumes that there are no variables other than stratum, oversample, and panels!!!
  if (is.null(panel_names)) {
    message("No panel_names vector was provided. Guessing panel names from data frame variables.")
    panel_names <- names(dataframe)[!(names(dataframe) %in% c(stratum_field, oversample_field))]
  }

  # For each stratum, make the design list specifying number of sample points per panel and the oversample count
  design <- lapply(X = dataframe[[stratum_field]],
                   dataframe = dataframe,
                   stratum_field = stratum_field,
                   panel_names = panel_names,
                   oversample_field = oversample_field,
                   oversample_proportion = oversample_proportion,
                   oversample_min = oversample_min,
                   FUN =   function(X,
                                    dataframe,
                                    stratum_field,
                                    panel_names,
                                    oversample_field,
                                    oversample_proportion,
                                    oversample_min){
                     # In case your data frame is from a .csv from an excel workbook that added an X to variable names that were entirely numeric
                     names_damaged_indices <- grepl(names(dataframe), pattern = "^X\\d+$")
                     names_damaged <- names(dataframe)[names_damaged_indices]
                     names_repaired <- gsub(names_damaged,
                                            pattern = "^X",
                                            replacement = "")
                     names(dataframe)[names_damaged_indices] <- names_repaired

                     # Get just the relevant variables in the current data frame for the stratum
                     variables_relevant <- names(dataframe)[names(dataframe) %in% c(stratum_field, panel_names, oversample_field)]
                     df_current <- dataframe[dataframe[[stratum_field]] == X, variables_relevant]
                     # Pull the panel values and create a named vector from them
                     panel <- stats::setNames(sapply(X = panel_names,
                                                     FUN = function(X, df){
                                                       return(df[, X])
                                                     },
                                                     df = df_current),
                                              panel_names)

                     # E pull the oversample value or calculate it
                     if (!is.null(oversample_field)) {
                       over <- df_current[[oversample_field]]
                     } else {
                       # For each panel, get the number of oversample points, then sum the vector
                       over <- sum(sapply(X = panel,
                                          oversample_proportion = oversample_proportion,
                                          oversample_min = oversample_min,
                                          FUN = function(X, oversample_proportion, oversample_min){
                                            return(max(round(X * oversample_proportion), oversample_min))
                                          }))
                     }

                     # Return the list for this stratum
                     return(list(panel = panel, seltype = "Equal", over = over))
                   })

  # Name the lists for each stratum with the stratum name
  design <- stats::setNames(design, dataframe[[stratum_field]])

  return(design)
}

#' Convert a GRTS design list to a data frame
#' @param design A list in the format taken by \code{spsurvey::grts()}.
#' @export
design_dataframe <- function(design){
  # Grab the strata names
  strata <- names(design)

  # Work through the strata, extracting the information and converting it into a data frame
  output <- do.call(rbind, lapply(X = strata,
                                  design = design,
                                  FUN = function(X, design){
                                    # Make a tidy data frame
                                    df <- data.frame(stratum = X,
                                                     panel = names(design[[X]][["panel"]]),
                                                     count = unname(design[[X]][["panel"]]))
                                    # Spread it into a wide format because users expect the variables to correspond to panels
                                    df_wide <- tidyr::spread(data = df,
                                                             key = panel,
                                                             value = count)
                                    df_wide[["oversample"]] <- design[[X]][["over"]]

                                    return(df_wide)
                                  }))
  return(output)
}
