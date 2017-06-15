#' Find names in an object matching a pattern
#' @param obj The named object (a data frame, Spatial Data Frame, list, etc.) to search for matching names.
#' @param name A string to search the names of \code{obj} using. This will be used as a regular expression in \code{grepl()}.
#' @param ignore.case Logical. Passed to \code{grepl()} to decide if the regular expression should be case insensitive. Defaults to \code{F}.
#' @param bookend Logical. If \code{T} then \code{^} and \code{$} will be added to the ends of the regular expression. Defaults to \code{T}.
#' @param multiple Logical. If \code{T} then the function can return a vector of multiple names instead of just one. If \code{F} then it will produce an error if multiple matches are found. Defaults to \code{F}.
#' @return A vector of character strings.
#' @export
find.name <- function(obj,
                      name = "",
                      ignore.case.name = F,
                      bookend,
                      multiple = F
){
  ## Get the data frame if it's an SPDF
  if (grepl(class(obj), pattern = "(^Spatial).{1,100}(DataFrame$)")) {
    obj <- obj@data
  }

  if(is.null(names(obj))) {
    stop("The provided object has no names.")
  }

  ## Get the fieldname as it occurs in the data frame
  fieldname <- names(obj)[grepl(x = names(obj),
                                pattern = if (bookend) {paste("^", name, "$")} else {name},
                                ignore.case = ignore.case.name)]

  ## Make sure that the fieldname exists
  if (length(fieldname) < 1) {
    stop("The fieldname does not exist in the data frame. Did you intend to set ignore.case.fieldname = T?")
  }
  if (length(fieldname) > 1 & !multiple) {
    message("The fieldname matched more than one variable/column:")
    stop(paste(fieldname, collapse = ", "))
  }

  return(fieldname)
}
