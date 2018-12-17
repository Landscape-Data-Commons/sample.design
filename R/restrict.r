#' Restrict an SPDF with the option to inherit values from a Spatial Polygons Data Frame
#' @description Restrict an SPDF to geometry where specified values are found in a given attribute in the data frame, to geometry overlapping a Spatial Polygons Data Frame, or both. If using a Spatial Polygons Data Frame as \code{spdf2}, that can also be filtered by values in an attribute before being used to restrict \code{spdf1}. Also, if using \code{spdf2} then the output can inherit the values in one attribute from that SPDF.
#' @param spdf1 An SPDF to restrict.
#' @param spdf2 An optional Spatial Polygons Data Frame to restrict \code{spdf1} by. Defaults to \code{NULL}.
#' @param inherit Logical. If \code{TRUE} then the output will be the restricted \code{spdf1} with an attribute called by the name that matched \code{inherit_field} in \code{spdf2} and the values from that attribute in \code{spdf2}. Defaults to \code{FALSE}.
#' @param inherit_field An optional character string to be used as a regular expression to find a single matching attribute name in \code{spdf2}. The associated values will be inherited by the output where there is spatial overlap and added to an attribute with the name that matched the regular expression. Required if \code{inherit} is \code{TRUE}.
#' @param ignore_case_inherit_field Logical. If \code{TRUE} then finding the attribute name in \code{spdf2} with \code{inherit_field} will be case insensitive. Defaults to \code{FALSE}.
#' @param bookend_inherit_field Logical. If \code{TRUE} then \code{^} and \code{$} will be added to the ends of the regular expression passed to \code{grepl()} when searching using \code{inherit_field}. Defaults to \code{FALSE}.
#' @param filter_field_spdf1 An optional character string to be used as a regular expression to find a single matching attribute name in \code{spdf1}. The associated values will be used to filter geometry to only observations where the matching attribute values match \code{filter_values_spdf1}.
#' @param filter_values_spdf1  An optional list or vector of values to filter \code{spdf1} by. Only geometry where the values in the attribute matching \code{filter_field_spdf1} which match \code{filter_values_spdf1} will be retained. By default the comparison is done with \code{match()} but uses \code{grepl()} \code{ignore_case_values_spdf1 == T} or \code{bookend_values_spdf1 == T}. If \code{grepl()} is used then these values are used to create a regular expression.
#' @param use_grep_values_spdf1 Logical. If \code{TRUE} then \code{filter_values_spdf1} will be used to create a regular expression and passed to \code{grepl()} instead of just being compared using \code{match}. Defaults to \code{FALSE}.
#' @param ignore_case_field_spdf1 Logical. If \code{TRUE} then finding the attribute name in \code{spdf1} using \code{filter_field_spdf1} will use \code{grepl()} instead of \code{match()} and be case insensitive. Defaults to \code{FALSE}.
#' @param ignore_case_values_spdf1 Logical. If \code{TRUE} then finding the values matching \code{filter_values_spdf1} in the attribute name in \code{spdf1} matching \code{filter_field_spdf1} will use \code{grepl()} instead of \code{match()} and be case insensitive. Defaults to \code{FALSE}.
#' @param bookend_field_spdf1 Logical. If \code{TRUE} then finding the attribute name in \code{spdf1} using \code{filter_field_spdf1} will use \code{grepl()} instead of \code{match()} and add \code{"^"} and \code{"$"} to the ends of the regular expression to force exact matches. Defaults to \code{FALSE}.
#' @param bookend_values_spdf1 Logical. If \code{TRUE} then finding the values matching \code{filter_values_spdf1} in the attribute name in \code{spdf1} matching \code{filter_field_spdf1} will use \code{grepl()} instead of \code{match()} and add \code{"^"} and \code{"$"} to the ends of the regular expression to force exact matches. Defaults to \code{FALSE}.
#' @param filter_field_spdf2 An optional character string to be used as a regular expression to find a single matching attribute name in \code{spdf2}. The associated values will be used to filter geometry to only observations where the matching attribute values match \code{filter_values_spdf2}.
#' @param filter_values_spdf2  An optional list or vector of values to filter \code{spdf2} by. Only geometry where the values in the attribute matching \code{filter_field_spdf2} which match \code{filter_values_spdf2} will be retained. By default the comparison is done with \code{match()} but uses \code{grepl()} \code{ignore_case_values_spdf2 == T} or \code{bookend_values_spdf2 == T}. If \code{grepl()} is used then these values are used to create a regular expression.
#' @param use_grep_values_spdf1 Logical. If \code{TRUE} then \code{filter_values_spdf2} will be used to create a regular expression and passed to \code{grepl()} instead of just being compared using \code{match}. Defaults to \code{FALSE}.
#' @param ignore_case_field_spdf2 Logical. If \code{TRUE} then finding the attribute name in \code{spdf2} using \code{filter_field_spdf2} will use \code{grepl()} instead of \code{match()} and be case insensitive. Defaults to \code{FALSE}.
#' @param ignore_case_values_spdf2 Logical. If \code{TRUE} then finding the values matching \code{filter_values_spdf2} in the attribute name in \code{spdf2} matching \code{filter_field_spdf2} will use \code{grepl()} instead of \code{match()} and be case insensitive. Defaults to \code{FALSE}.
#' @param bookend_field_spdf2Logical. If \code{TRUE} then finding the attribute name in \code{spdf2} using \code{filter_field_spdf2} will use \code{grepl()} instead of \code{match()} and add \code{"^"} and \code{"$"} to the ends of the regular expression to force exact matches. Defaults to \code{FALSE}.
#' @param bookend_values_spdf2 Logical. If \code{TRUE} then finding the values matching \code{filter_values_spdf2} in the attribute name in \code{spdf2} matching \code{filter_field_spdf2} will use \code{grepl()} instead of \code{match()} and add \code{"^"} and \code{"$"} to the ends of the regular expression to force exact matches. Defaults to \code{FALSE}.
#' @return An SPDF of geometry and values from \code{spdf1} where the filtering criteria were met. If \code{inherit == T} then there will be an additional attribute from \code{spdf2}.
#' @examples
#' restrict(spdf1 = wyoming.spdf,
#'          filter_field_spdf1 = "COUNTY",
#'          filter_values_spdf1 = "Teton")
#' restrict(spdf1 = wyoming.spdf,
#'          spdf2 = grand.tetons.np.spdf)
#' restrict(spdf1 = wyoming.spdf,
#'          spdf2 = grand.tetons.np.spdf,
#'          filter_field_spdf1 = "COUNTY",
#'          filter_values_spdf1 = "Teton",
#'          inherit = TRUE,
#'          inherit_field = "OWNERSHIP")
#' restrict(spdf1 = wyoming.spdf,
#'          spdf2 = national.ownership,
#'          filter_field_spdf2 = "OWNERSHIP",
#'          filter_values_spdf2 = c("Bureau of Land Management", "BLM"))
#' restrict(spdf1 = wyoming.spdf,
#'          spdf2 = national.ownership,
#'          filter_field_spdf1 = "COUNTY",
#'          filter_values_spdf1 = "Teton",
#'          filter_field_spdf2 = "OWNERSHIP",
#'          # Because I can't confidently spell "bureau"
#'          filter_values_spdf2 = "Bur[(e)|(a)|(u)]{2,5} of Land Management",
#'          use_grep_values_spdf2 = TRUE)
#' @export

restrict <- function(spdf1 = NULL,
                     spdf2 = NULL,
                     inherit = FALSE,
                     inherit_field = NULL,
                     ignore_case_inherit_field = FALSE,
                     bookend_inherit_field = FALSE,
                     filter_field_spdf1 = NULL,
                     filter_values_spdf1 = NULL,
                     use_grep_values_spdf1 = FALSE,
                     ignore_case_field_spdf1 = FALSE,
                     ignore_case_values_spdf1 = FALSE,
                     bookend_field_spdf1 = FALSE,
                     bookend_values_spdf1 = FALSE,
                     filter_field_spdf2 = NULL,
                     filter_values_spdf2 = NULL,
                     use_grep_values_spdf2 = FALSE,
                     ignore_case_field_spdf2 = FALSE,
                     ignore_case_values_spdf2 = FALSE,
                     bookend_field_spdf2 = FALSE,
                     bookend_values_spdf2 = FALSE
){
  ## Validity checks
  # Make sure spdf1 is a proper SPDF
  if (!grepl(class(spdf1), pattern = "^(Spatial).{1,100}(DataFrame)$")) {
    stop("An SPDF must be provided as spdf1")
  }
  # Make sure that spdf2 is a proper SPDF if provided
  if (!is.null(spdf2) & class(spdf2) != "SpatialPolygonsDataFrame") {
    stop("spdf2 must be a Spatial Polygons Data Frame")
  }

  # Making sure that values are character strings if provided
  if (!is.null(filter_field_spdf1) & !is.character(filter_field_spdf1)) {
    stop("filter_field_spdf1 must be a character string")
  }
  if (!is.null(filter_field_spdf2) & !is.character(filter_field_spdf2)) {
    stop("filter_field_spdf2 must be a character string")
  }
  if (!is.null(filter_values_spdf1) & !is.character(filter_values_spdf1)) {
    stop("filter_values_spdf1 must be a character string")
  }
  if (!is.null(filter_values_spdf2) & !is.character(filter_values_spdf2)) {
    stop("filter_values_spdf2 must be a character string")
  }
  if (inherit & (is.null(inherit_field) | !is.character(inherit_field))) {
    stop("inherit_field must be a character string")
  }
  if (!is.null(inherit_field) & !is.character(inherit_field)) {
    stop("inherit_field must be a character string")
  }

  # Making sure that if a filter field is provided, there are also values, and vice versa
  if (xor(!is.null(filter_field_spdf1), !is.null(filter_values_spdf1))) {
    stop("If filter_field_spdf1 or filter_values_spdf1 is provided, the other one must be as well.")
  }
  if (xor(!is.null(filter_field_spdf2), !is.null(filter_values_spdf2))) {
    stop("If filter_field_spdf2 or filter_values_spdf2 is provided, the other one must be as well.")
  }

  # Making sure that the fields actually exist
  if (is.character(filter_field_spdf1)) {
    filter_field_spdf1 <- find_name(obj = spdf1@data,
                                    name = filter_field_spdf1,
                                    ignore_case_name = ignore_case_field_spdf1,
                                    bookend = bookend_field_spdf1,
                                    multiple = FALSE)
  }
  if (is.character(filter_field_spdf2)) {
    filter_field_spdf2 <- find_name(obj = spdf2@data,
                                    name = filter_field_spdf2,
                                    ignore_case_name = ignore_case_field_spdf2,
                                    bookend = bookend_field_spdf2,
                                    multiple = FALSE)
  }
  if (is.character(inherit_field)) {
    inherit_field <- find_name(obj = spdf2@data,
                               name = inherit_field,
                               ignore_case_name = ignore_case_inherit_field,
                               bookend = bookend_inherit_field,
                               multiple = FALSE)
    if (inherit_field %in% names(spdf1)) {
      message(paste0("The column/variable ", inherit_field, " already exists in spdf1 and will be overwritten."))
    }
  }

  # Provide warning messages if arguments will be unused
  if (is.null(spdf2) & inherit) {
    message("inherit is TRUE but ignored because there is no spdf2")
  }
  if (!inherit & is.character(inherit_field)) {
    message("inherit_field is valid but ignored because inherit is FALSE")
  }
  if (inherit & !is.character(inherit_field)) {
    stop("A valid inherit_field must be provided when inherit is TRUE")
  }
  if (!inherit & ignore_case_inherit_field) {
    message("ignore_case_inherit_field is TRUE but ignored because inherit is FALSE")
  }
  if (!inherit & bookend_inherit_field) {
    message("bookend_inherit_field is TRUE but ignored because inherit is FALSE")
  }
  if (is.null(filter_field_spdf1) & bookend_field_spdf1) {
    message("bookend_field_spdf1 is TRUE but ignored because no filter_field_spdf1 was provided.")
  }
  if (is.null(filter_values_spdf1) & bookend_values_spdf1) {
    message("bookend_values_spdf1 is TRUE but ignored because no filter_values_spdf1 was provided.")
  }
  if (is.null(filter_field_spdf2) & bookend_field_spdf2) {
    message("bookend_field_spdf2 is TRUE but ignored because no filter_field_spdf2 was provided.")
  }
  if (is.null(filter_values_spdf2) & bookend_values_spdf2) {
    message("bookend_values_spdf2 is TRUE but ignored because no filter_values_spdf2 was provided.")
  }

  ## If filtering the SPDFs, do it here
  if (!is.null(filter_field_spdf1)) {
    spdf1 <- spdf1[spdf1@data[[filter_field_spdf1]] %in% search(df = spdf1,
                                                                values = filter_values_spdf1,
                                                                namestring = filter_field_spdf1,
                                                                ignore_case_namestring = FALSE,
                                                                bookend_namestring = TRUE,
                                                                use_grepl = use_grep_values_spdf1,
                                                                ignore_case_values = ignore_case_values_spdf1,
                                                                bookend_values = bookend_values_spdf1),
                   ]
  }

  if (!is.null(filter_field_spdf2)) {
    spdf2 <- spdf2[spdf2@data[[filter_field_spdf2]] %in% search(df = spdf2,
                                                                values = filter_values_spdf2,
                                                                namestring = filter_field_spdf2,
                                                                ignore_case_namestring = FALSE,
                                                                bookend_namestring = TRUE,
                                                                use_grepl = use_grep_values_spdf2,
                                                                ignore_case_values = ignore_case_values_spdf2,
                                                                bookend_values = bookend_values_spdf2),
                   ]
  }

  if (is.null(spdf2)) {
    if (class(spdf1) == class(spdf2)) {
      ## Get the intersection as a SpatialPolygons
      current.drop <- rgeos::get_RGEOS_dropSlivers()
      current.warn <- rgeos::get_RGEOS_warnSlivers()
      current.tol <- rgeos::get_RGEOS_polyThreshold()
      rgeos::set_RGEOS_dropSlivers(TRUE)
      rgeos::set_RGEOS_warnSlivers(TRUE)
      rgeos::set_RGEOS_polyThreshold(0.01)
      intersection <- rgeos::gIntersection(spgeom1 = spdf1,
                                           spgeom2 = spdf2,
                                           drop_lower_td = TRUE)
      rgeos::set_RGEOS_dropSlivers(current.drop)
      rgeos::set_RGEOS_warnSlivers(current.warn)
      rgeos::set_RGEOS_polyThreshold(current.tol)

      ## Now we need to build the data frame that goes back into this. It's a pain
      ## Get the rownames from the polygons. This will consist of the two row names from spdf1 and spdf2 separated by a " "
      intersection_rownames <- strsplit(row.names(intersection), split = " ")

      ## Create an empty data frame that we can add the constructed rows to
      intersection_dataframe <- data.frame()

      ## For each of the intersection polygons, create a row with the attributes from the source polygons
      for (row in seq_along(intersection_rownames)) {
        intersection_dataframe <- rbind(intersection_dataframe,
                                        cbind(
                                          spdf1@data[intersection_rownames[[row]][1], ],
                                          spdf2@data[intersection_rownames[[row]][2], ]
                                        ))
      }
      rownames(intersection_dataframe) <- row.names(intersection)

      ## Create the output SPDF
      output <- sp::SpatialPolygonsDataFrame(Sr = intersection,
                                             data = intersection_dataframe)

      ## Return only the spdf1 fields and the inherit field, for consistency
      output@data <- output@data[, names(output@data)[names(output@data) %in% c(names(spdf1@data), inherit_field)]]
    } else {
      ## Attribute spdf1 with spdf2
      output <- attribute_shapefile(spdf1 = spdf1,
                                    spdf2 = spdf2,
                                    newfield = if (inherit) {
                                      inherit_field
                                      } else {
                                        NULL
                                      },
                                    attributefield = if (inherit) {
                                      inherit_field
                                      } else {
                                        NULL
                                      })
    }
  } else {
    output <- spdf1
  }


  return(output)
}
