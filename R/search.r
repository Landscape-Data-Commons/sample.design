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
                               ignore.case = ignore.case)]

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

#' Find values in a data frame column/variable that match a list or vector
#' @description Return all unique values from a column/variable in a data frame which match values in a list or vector.
#' @param df The data frame to search through
#' @param name A string to search the names of \code{df} using. This will be used as a regular expression in \code{grepl()}.
#' @param values A list or vector of values to search for. By default the comparison is done with \code{match()} but uses \code{grepl()} if \code(use.grepl = T) or \code(ignore.case.values = T). If \code{grepl()} is used then these values are used to create a regular expression.
#' @param use.grepl Logical. If \code{T} then the search will be done using \code{grepl()} instead of \code{match()}. Defaults to \code{F}.
#' @param ignore.case.name Logical. If \code{T} then finding the column/variable name will be case insensitive. Defaults to \code{F}.
#' @param ignore.case.value Logical. If \code{T} then finding the values will be case insensitive and will use \code{grepl()} instead of \code{match()}. Defaults to \code{F}.
#' @param bookend.name Logical. If \code{T} then \code{^} and \code{$} will be added to the ends of the regular expression passed to \code{grepl()} whean searching using \code{name}. Defaults to \code{F}.
#' @param bookend.values Logical. If \code{T} then \code{^} and \code{$} will be added to the ends of the regular expression passed to \code{grepl()} when searching using \code{values}. Defaults to \code{F}.
#' @return A vector of unique values.
#' @export
search <- function(df = NULL,
                   values,
                   name = "",
                   use.grepl = F,
                   ignore.case.name = F,
                   ignore.case.values = F,
                   bookend.name = F,
                   bookend.values = F
){
  if (is.null(df)) {
    stop("Missing a data frame as df.")
  }
  ## Get the data frame from df if it's an SPDF
  if (grepl(class(df), pattern = "(^Spatial).{1,100}(DataFrame$)")) {
    df <- df@data
  }

  values <- unique(values)

  ## Get the fieldname as it occurs in the data frame
  fieldname <- find.name(obj = df,
                         name = name,
                         ignore.case = ignore.case.name,
                         bookend = bookend.name,
                         multiple = F)

  ## If grepl() is needed
  if (use.grepl | ignore.case.values) {
    if (any(sapply(values, function (X) {grepl(class(X), pattern = "logical")}))) {
      message("At least one of the values is logical and using grepl() may produce unexpected results depending on if the logical values are stored as T/F or TRUE/FALSE.")
    }
    ## Create a regular expression looking for the values
    output <- df[grepl(df[,fieldname],
                       pattern = if (bookend.values) {paste0("^(", paste(values, collapse = ")|("), ")$")} else {paste0("(", paste(values, collapse = ")|("), ")")},
                       ignore.case = ignore.case.values)]
  } else {
    output <- df[df$fieldname %in% values,]
  }

  return(unique(output[,fieldname]))
}



