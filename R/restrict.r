#' Restrict an SPDF with the option to inherit values from a Spatial Polygons Data Frame
#' @description Restrict an SPDF to geometry where specified values are found in a given attribute in the data frame, to geometry overlapping a Spatial Polygons Data Frame, or both. If using a Spatial Polygons Data Frame as \code{spdf2}, that can also be filtered by values in an attribute before being used to restrict \code{spdf1}. Also, if using \code{spdf2} then the output can inherit the values in one attribute from that SPDF.
#' @param spdf1 An SPDF to restrict.
#' @param spdf2 An optional Spatial Polygons Data Frame to restrict \code{spdf1} by. Defaults to \code{NULL}.
#' @param inherit Logical. If \code{T} then the output will be the restricted \code{spdf1} with an attribute called by the name that matched \code{inherit.field} in \code{spdf2} and the values from that attribute in \code{spdf2}. Defaults to \code{F}.
#' @param inherit.field An optional character string to be used as a regular expression to find a single matching attribute name in \code{spdf2}. The associated values will be inherited by the output where there is spatial overlap and added to an attribute with the name that matched the regular expression. Required if \code{inherit} is \code{T}.
#' @param ignore.case.inherit.field Logical. If \code{T} then finding the attribute name in \code{spdf2} with \code{inherit.field} will be case insensitive. Defaults to \code{F}.
#' @param bookend.inherit.field Logical. If \code{T} then \code{^} and \code{$} will be added to the ends of the regular expression passed to \code{grepl()} when searching using \code{inherit.field}. Defaults to \code{F}.
#' @param filter.field.spdf1 An optional character string to be used as a regular expression to find a single matching attribute name in \code{spdf1}. The associated values will be used to filter geometry to only observations where the matching attribute values match \code{filter.values.spdf1}.
#' @param filter.values.spdf1  An optional list or vector of values to filter \code{spdf1} by. Only geometry where the values in the attribute matching \code{filter.field.spdf1} which match \code{filter.values.spdf1} will be retained. By default the comparison is done with \code{match()} but uses \code{grepl()} \code{ignore.case.values.spdf1 == T} or \code{bookend.values.spdf1 == T}. If \code{grepl()} is used then these values are used to create a regular expression.
#' @param use.grep.values.spdf1 Logical. If \code{T} then \code{filter.values.spdf1} will be used to create a regular expression and passed to \code{grepl()} instead of just being compared using \code{match}. Defaults to \code{F}.
#' @param ignore.case.field.spdf1 Logical. If \code{T} then finding the attribute name in \code{spdf1} using \code{filter.field.spdf1} will use \code{grepl()} instead of \code{match()} and be case insensitive. Defaults to \code{F}.
#' @param ignore.case.values.spdf1 Logical. If \code{T} then finding the values matching \code{filter.values.spdf1} in the attribute name in \code{spdf1} matching \code{filter.field.spdf1} will use \code{grepl()} instead of \code{match()} and be case insensitive. Defaults to \code{F}.
#' @param bookend.field.spdf1 Logical. If \code{T} then finding the attribute name in \code{spdf1} using \code{filter.field.spdf1} will use \code{grepl()} instead of \code{match()} and add \code{"^"} and \code{"$"} to the ends of the regular expression to force exact matches. Defaults to \code{F}.
#' @param bookend.values.spdf1 Logical. If \code{T} then finding the values matching \code{filter.values.spdf1} in the attribute name in \code{spdf1} matching \code{filter.field.spdf1} will use \code{grepl()} instead of \code{match()} and add \code{"^"} and \code{"$"} to the ends of the regular expression to force exact matches. Defaults to \code{F}.
#' @param filter.field.spdf2 An optional character string to be used as a regular expression to find a single matching attribute name in \code{spdf2}. The associated values will be used to filter geometry to only observations where the matching attribute values match \code{filter.values.spdf2}.
#' @param filter.values.spdf2  An optional list or vector of values to filter \code{spdf2} by. Only geometry where the values in the attribute matching \code{filter.field.spdf2} which match \code{filter.values.spdf2} will be retained. By default the comparison is done with \code{match()} but uses \code{grepl()} \code{ignore.case.values.spdf2 == T} or \code{bookend.values.spdf2 == T}. If \code{grepl()} is used then these values are used to create a regular expression.
#' @param use.grep.values.spdf1 Logical. If \code{T} then \code{filter.values.spdf2} will be used to create a regular expression and passed to \code{grepl()} instead of just being compared using \code{match}. Defaults to \code{F}.
#' @param ignore.case.field.spdf2 Logical. If \code{T} then finding the attribute name in \code{spdf2} using \code{filter.field.spdf2} will use \code{grepl()} instead of \code{match()} and be case insensitive. Defaults to \code{F}.
#' @param ignore.case.values.spdf2 Logical. If \code{T} then finding the values matching \code{filter.values.spdf2} in the attribute name in \code{spdf2} matching \code{filter.field.spdf2} will use \code{grepl()} instead of \code{match()} and be case insensitive. Defaults to \code{F}.
#' @param bookend.field.spdf2Logical. If \code{T} then finding the attribute name in \code{spdf2} using \code{filter.field.spdf2} will use \code{grepl()} instead of \code{match()} and add \code{"^"} and \code{"$"} to the ends of the regular expression to force exact matches. Defaults to \code{F}.
#' @param bookend.values.spdf2 Logical. If \code{T} then finding the values matching \code{filter.values.spdf2} in the attribute name in \code{spdf2} matching \code{filter.field.spdf2} will use \code{grepl()} instead of \code{match()} and add \code{"^"} and \code{"$"} to the ends of the regular expression to force exact matches. Defaults to \code{F}.
#' @return An SPDF of geometry and values from \code{spdf1} where the filtering criteria were met. If \code{inherit == T} then there will be an additional attribute from \code{spdf2}.
#' @examples
#' restrict(spdf1 = wyoming.spdf,
#'          filter.field.spdf1 = "COUNTY",
#'          filter.values.spdf1 = "Teton")
#' restrict(spdf1 = wyoming.spdf,
#'          spdf2 = grand.tetons.np.spdf)
#' restrict(spdf1 = wyoming.spdf,
#'          spdf2 = grand.tetons.np.spdf,
#'          filter.field.spdf1 = "COUNTY",
#'          filter.values.spdf1 = "Teton",
#'          inherit = T,
#'          inherit.field = "OWNERSHIP")
#' restrict(spdf1 = wyoming.spdf,
#'          spdf2 = national.ownership,
#'          filter.field.spdf2 = "OWNERSHIP",
#'          filter.values.spdf2 = "(Bureau of Land Management)|(BLM)",
#'          bookend.values.spdf2 = T)
#' restrict(spdf1 = wyoming.spdf,
#'          spdf2 = national.ownership,
#'          filter.field.spdf1 = "COUNTY",
#'          filter.values.spdf1 = "Teton",
#'          filter.field.spdf2 = "OWNERSHIP",
#'          filter.values.spdf2 = "(Bureau of Land Management)|(BLM)",
#'          use.grep.values.spdf2 = T)
#' @export

restrict <- function(spdf1 = NULL,
                     spdf2 = NULL,
                     inherit = F,
                     inherit.field = NULL,
                     ignore.case.inherit.field = F,
                     bookend.inherit.field = F,
                     filter.field.spdf1 = NULL,
                     filter.values.spdf1 = NULL,
                     use.grep.values.spdf1 = F,
                     ignore.case.field.spdf1 = F,
                     ignore.case.values.spdf1 = F,
                     bookend.field.spdf1 = F,
                     bookend.values.spdf1 = F,
                     filter.field.spdf2 = NULL,
                     filter.values.spdf2 = NULL,
                     use.grep.values.spdf2 = F,
                     ignore.case.field.spdf2 = F,
                     ignore.case.values.spdf2 = F,
                     bookend.field.spdf2 = F,
                     bookend.values.spdf2 = F
){
  ## Validity checks
  # Make sure spdf1 is a proper SPDF
  if (!grepl(class(spdf1), pattern = "(^Spatial).{1,100}(DataFrame$)")) {
    stop("An SPDF must be provided as spdf1")
  }
  # Make sure that spdf2 is a proper SPDF if provided
  if (!is.null(spdf2) & class(spdf2) != "SpatialPolygonsDataFrame") {
    stop("spdf2 must be a Spatial Polygons Data Frame")
  }

  # Making sure that values are character strings if provided
  if (!is.null(filter.field.spdf1) & !is.character(filter.field.spdf1)) {
    stop("filter.field.spdf1 must be a character string")
  }
  if (!is.null(filter.field.spdf2) & !is.character(filter.field.spdf2)) {
    stop("filter.field.spdf2 must be a character string")
  }
  if (!is.null(filter.values.spdf1) & !is.character(filter.values.spdf1)) {
    stop("filter.values.spdf1 must be a character string")
  }
  if (!is.null(filter.values.spdf2) & !is.character(filter.values.spdf2)) {
    stop("filter.values.spdf2 must be a character string")
  }
  if (inherit & (is.null(inherit.field) | !is.character(inherit.field))) {
    stop("inherit.field must be a character string")
  }
  if (!is.null(inherit.field) & !is.character(inherit.field)) {
    stop("inherit.field must be a character string")
  }

  # Making sure that if a filter field is provided, there are also values, and vice versa
  if (xor(!is.null(filter.field.spdf1), !is.null(filter.values.spdf1))) {
    stop("If filter.field.spdf1 or filter.values.spdf1 is provided, the other one must be as well.")
  }
  if (xor(!is.null(filter.field.spdf2), !is.null(filter.values.spdf2))) {
    stop("If filter.field.spdf2 or filter.values.spdf2 is provided, the other one must be as well.")
  }

  # Making sure that the fields actually exist
  if (is.character(filter.field.spdf1)) {
    filter.field.spdf1 <- find.name(obj = spdf1@data,
                                    name = filter.field.spdf1,
                                    ignore.case = ignore.case.field.spdf1,
                                    bookend = bookend.field.spdf1,
                                    multiple = F)
  }
  if (is.character(filter.field.spdf2)) {
    filter.field.spdf2 <- find.name(obj = spdf2@data,
                                    name = filter.field.spdf2,
                                    ignore.case = ignore.case.field.spdf2,
                                    bookend = bookend.field.spdf2,
                                    multiple = F)
  }
  if (is.character(inherit.field)) {
    inherit.field <- find.name(obj = spdf2@data,
                               name = inherit.field,
                               ignore.case = ignore.case.inherit.field,
                               bookend = bookend.inherit.field,
                               multiple = F)
    if (inherit.field %in% names(spdf1)) {
      message(paste0("The column/variable ", inherit.field, " already exists in spdf1 and will be overwritten."))
    }
  }

  # Provide warning messages if arguments will be unused
  if (is.null(spdf2) & inherit) {
    message("inherit is T but ignored because there is no spdf2")
  }
  if (!inherit & is.character(inherit.field)) {
    message("inherit.field is valid but ignored because inherit is F")
  }
  if (inherit & !is.character(inherit.field)) {
    stop("A valid inherit.field must be provided when inherit is T")
  }
  if (!inherit & ignore.case.inherit.field) {
    message("ignore.case.inherit.field is T but ignored because inherit is F")
  }
  if (!inherit & bookend.inherit.field) {
    message("bookend.inherit.field is T but ignored because inherit is F")
  }
  if (is.null(filter.field.spdf1) & bookend.field.spdf1) {
    message("bookend.field.spdf1 is T but ignored because no filter.field.spdf1 was provided.")
  }
  if (is.null(filter.values.spdf1) & bookend.values.spdf1) {
    message("bookend.values.spdf1 is T but ignored because no filter.values.spdf1 was provided.")
  }
  if (is.null(filter.field.spdf2) & bookend.field.spdf2) {
    message("bookend.field.spdf2 is T but ignored because no filter.field.spdf2 was provided.")
  }
  if (is.null(filter.values.spdf2) & bookend.values.spdf2) {
    message("bookend.values.spdf2 is T but ignored because no filter.values.spdf2 was provided.")
  }

  ## If filtering the SPDFs, do it here
  if (!is.null(filter.field.spdf1)) {
    spdf1 <- spdf1[spdf1@data[[filter.field.spdf1]] %in% search(df = spdf1,
                                                                values = filter.values.spdf1,
                                                                namestring = filter.field.spdf1,
                                                                ignore.case.namestring = F,
                                                                bookend.namestring = T,
                                                                use.grepl = use.grep.values.spdf1,
                                                                ignore.case.values = ignore.case.values.spdf1,
                                                                bookend.values = bookend.values.spdf1),
                   ]
  }

  if (!is.null(filter.field.spdf2)) {
    spdf2 <- spdf2[spdf2@data[[filter.field.spdf2]] %in% search(df = spdf2,
                                                                values = filter.values.spdf2,
                                                                namestring = filter.field.spdf2,
                                                                ignore.case.namestring = F,
                                                                bookend.namestring = T,
                                                                use.grepl = use.grep.values.spdf2,
                                                                ignore.case.values = ignore.case.values.spdf2,
                                                                bookend.values = bookend.values.spdf2),
                   ]
  }

  if (is.null(spdf2)) {
    if (class(spdf1) == class(spdf2)) {
      ## Get the intersection as a SpatialPolygons
      current.drop <- rgeos::get_RGEOS_dropSlivers()
      current.warn <- rgeos::get_RGEOS_warnSlivers()
      current.tol <- rgeos::get_RGEOS_polyThreshold()
      rgeos::set_RGEOS_dropSlivers(T)
      rgeos::set_RGEOS_warnSlivers(T)
      rgeos::set_RGEOS_polyThreshold(0.01)
      intersection <- rgeos::gIntersection(spgeom1 = spdf1,
                                           spgeom2 = spdf2,
                                           drop_lower_td = T)
      rgeos::set_RGEOS_dropSlivers(current.drop)
      rgeos::set_RGEOS_warnSlivers(current.warn)
      rgeos::set_RGEOS_polyThreshold(current.tol)

      ## Now we need to build the data frame that goes back into this. It's a pain
      ## Get the rownames from the polygons. This will consist of the two row names from spdf1 and spdf2 separated by a " "
      intersection.rownames <- strsplit(row.names(intersection), split = " ")

      ## Create an empty data frame that we can add the constructed rows to
      intersection.dataframe <- data.frame()

      ## For each of the intersection polygons, create a row with the attributes from the source polygons
      for (row in 1:length(intersection.rownames)) {
        intersection.dataframe <- rbind(intersection.dataframe,
                                        cbind(
                                          spdf1@data[intersection.rownames[[row]][1],],
                                          spdf2@data[intersection.rownames[[row]][2],]
                                        ))
      }
      rownames(intersection.dataframe) <- row.names(intersection)

      ## Create the output SPDF
      output <- sp::SpatialPolygonsDataFrame(Sr = intersection,
                                             data = intersection.dataframe)

      ## Return only the spdf1 fields and the inherit field, for consistency
      output@data <- output@data[, names(output@data)[names(output@data) %in% c(names(spdf1@data), inherit.field)]]
    } else {
      ## Attribute spdf1 with spdf2
      output <- attribute.shapefile(spdf1 = spdf1,
                                    spdf2 = spdf2,
                                    newfield = if (inherit) {inherit.field} else {NULL},
                                    attributefield = if (inherit) {inherit.field} else {NULL})
    }
  } else {
    output <- spdf1
  }


  return(output)
}
