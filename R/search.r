#' Find values in a data frame column/variable that match a list or vector
#' @description Return all unique values from a column/variable in a data frame which match values in a list or vector.
#' @param df The data frame to search through
#' @param namestring A string to search the names of \code{df} using. This will be used as a regular expression in \code{grepl()}.
#' @param values A list or vector of values to search for. By default the comparison is done with \code{match()} but uses \code{grepl()} if \code{use_grepl = TRUE} or \code{ignore_case_values = TRUE}. If \code{grepl()} is used then these values are used to create a regular expression.
#' @param use_grepl Logical. If \code{TRUE} then the search will be done using \code{grepl()} instead of \code{match()}. Defaults to \code{FALSE}.
#' @param ignore_case_namestring Logical. If \code{TRUE} then finding the column/variable name will be case insensitive. Defaults to \code{FALSE}.
#' @param ignore_case_values Logical. If \code{TRUE} then finding the values will be case insensitive and will use \code{grepl()} instead of \code{match()}. Defaults to \code{FALSE}.
#' @param bookend_namestring Logical. If \code{TRUE} then \code{^} and \code{$} will be added to the ends of the regular expression passed to \code{grepl()} when searching using \code{name}. Defaults to \code{FALSE}.
#' @param bookend_values Logical. If \code{TRUE} then \code{^} and \code{$} will be added to the ends of the regular expression passed to \code{grepl()} when searching using \code{values}. Defaults to \code{FALSE}.
#' @return A vector of unique values.
#' @export
search <- function(df,
                   values,
                   namestring = "",
                   use_grepl = FALSE,
                   ignore_case_namestring = FALSE,
                   ignore_case_values = FALSE,
                   bookend_namestring = FALSE,
                   bookend_values = FALSE
){
  ## Get the data frame from df if it's an SPDF
  if (grepl(class(df), pattern = "^(Spatial).{1,100}(DataFrame)$")) {
    df <- df@data
  }

  values <- unique(values)

  ## Get the fieldname as it occurs in the data frame
  fieldname <- find.name(obj = df,
                         name = namestring,
                         ignore_case = ignore_case_namestring,
                         bookend = bookend_namestring,
                         multiple = FALSE)

  ## If grepl() is needed
  if (use_grepl | ignore_case_values) {
    if (any(sapply(values,
                   function (X) {
                     grepl(class(X), pattern = "logical")
                     }
                   ))
        ) {
      message("At least one of the values is logical and using grepl() may produce unexpected results depending on if the logical values are stored as T/F or TRUE/FALSE.")
    }
    ## Create a regular expression looking for the values
    output <- df[grepl(df[, fieldname],
                       pattern = if (bookend_values) {
                         paste0("^(", paste(values, collapse = ")|("), ")$")
                       } else {
                         paste0("(", paste(values, collapse = ")|("), ")")
                       },
                       ignore.case = ignore_case_values), ]
  } else {
    output <- df[df[[fieldname]] %in% values, ]
  }

  return(unique(output[, fieldname]))
}
