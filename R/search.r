#' Find values in a data frame column/variable that match a list or vector
#' @description Return all unique values from a column/variable in a data frame which match values in a list or vector.
#' @param df The data frame to search through
#' @param namestring A string to search the names of \code{df} using. This will be used as a regular expression in \code{grepl()}.
#' @param values A list or vector of values to search for. By default the comparison is done with \code{match()} but uses \code{grepl()} if \code{use.grepl = T} or \code{ignore.case.values = T}. If \code{grepl()} is used then these values are used to create a regular expression.
#' @param use.grepl Logical. If \code{T} then the search will be done using \code{grepl()} instead of \code{match()}. Defaults to \code{F}.
#' @param ignore.case.namestring Logical. If \code{T} then finding the column/variable name will be case insensitive. Defaults to \code{F}.
#' @param ignore.case.values Logical. If \code{T} then finding the values will be case insensitive and will use \code{grepl()} instead of \code{match()}. Defaults to \code{F}.
#' @param bookend.namestring Logical. If \code{T} then \code{^} and \code{$} will be added to the ends of the regular expression passed to \code{grepl()} when searching using \code{name}. Defaults to \code{F}.
#' @param bookend.values Logical. If \code{T} then \code{^} and \code{$} will be added to the ends of the regular expression passed to \code{grepl()} when searching using \code{values}. Defaults to \code{F}.
#' @return A vector of unique values.
#' @export
search <- function(df = NULL,
                   values,
                   namestring = "",
                   use.grepl = F,
                   ignore.case.namestring = F,
                   ignore.case.values = F,
                   bookend.namestring = F,
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
                         name = namestring,
                         ignore.case = ignore.case.namestring,
                         bookend = bookend.namestring,
                         multiple = F)

  ## If grepl() is needed
  if (use.grepl | ignore.case.values) {
    if (any(sapply(values, function (X) {grepl(class(X), pattern = "logical")}))) {
      message("At least one of the values is logical and using grepl() may produce unexpected results depending on if the logical values are stored as T/F or TRUE/FALSE.")
    }
    ## Create a regular expression looking for the values
    output <- df[grepl(df[,fieldname],
                       pattern = if (bookend.values) {paste0("^(", paste(values, collapse = ")|("), ")$")} else {paste0("(", paste(values, collapse = ")|("), ")")},
                       ignore.case = ignore.case.values), ]
  } else {
    output <- df[df$fieldname %in% values,]
  }

  return(unique(output[,fieldname]))
}



