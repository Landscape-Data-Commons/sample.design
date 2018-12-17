#' Find names in an object matching a pattern
#' @param obj The named object (a data frame, Spatial Data Frame, list, etc.) to search for matching names.
#' @param name A string to search the names of \code{obj} using. This will be used as a regular expression in \code{grepl()}.
#' @param ignore_case_name Logical. Passed to \code{grepl()} to decide if the regular expression should be case insensitive. Defaults to \code{FALSE}.
#' @param bookend Logical. If \code{TRUE} then \code{^} and \code{$} will be added to the ends of the regular expression. Defaults to \code{TRUE}.
#' @param multiple Logical. If \code{TRUE} then the function can return a vector of multiple names instead of just one. If \code{FALSE} then it will produce an error if multiple matches are found. Defaults to \code{FALSE}.
#' @return A vector of character strings.
#' @export
find_name <- function(obj,
                      name = "",
                      ignore_case_name = FALSE,
                      bookend,
                      multiple = FALSE
){
  ## Get the data frame if it's an SPDF
  if (grepl(class(obj), pattern = "^(Spatial).{1,100}(DataFrame)$")) {
    obj <- obj@data
  }

  if (is.null(names(obj))) {
    stop("The provided object has no names.")
  }

  ## Get the fieldname as it occurs in the data frame
  fieldname <- names(obj)[grepl(x = names(obj),
                                pattern = if (bookend) {
                                  paste("^", name, "$")
                                  } else {
                                    name
                                  },
                                ignore.case = ignore_case_name)]

  ## Make sure that the fieldname exists
  if (length(fieldname) < 1) {
    stop(paste0("No name matching the regular expression '",
                if (bookend) {
                  paste("^", name, "$")
                  } else {
                    name
                  },
                "' exists in the object with ignore_case_name = ",
                ignore_case_name,
                " and bookend = ", bookend, "."))
  }
  if (length(fieldname) > 1 & !multiple) {
    message(paste0("The regular expression '",
                   if (bookend) {
                     paste("^", name, "$")
                   } else {
                     name
                   },
                   "' with ignore_case_name = ", ignore_case_name,
                   " and bookend = ", bookend,
                   " matched more than one name:"))
    stop(paste(fieldname, collapse = ", "))
  }

  return(fieldname)
}
