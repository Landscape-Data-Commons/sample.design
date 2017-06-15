#' Search a variable in an SPDF for values
#' @description A wrapper for \code{[grepl(),]} to help quickly find values matching a regular expression in a given variable in an SPDF.
#' @param df The data frame or SPDF to search through.
#' @param field A character string specifying the field name to search within.
#' @param ignore.case.field Logical. If \code{T} then search for \code{field} in the variable names in \code{spdf} is case insensitive. Defaults to \code{T}.
#' @param regex A character string regular expression to be passed to \code{grepl()}.
#' @param ignore.case.regex  Logical. If \code{T} then search for \code{regex} in \code{field} is case insensitive. Defaults to \code{T}.
#' @return A vector of unique values found in \code{field} matching \code{regex}.
#' @export
df.search <- function(df,
                        field = "",
                        ignore.case.field = T,
                        regex = "",
                        ignore.case.regex = T
){
  if (grepl(x = class(df), pattern = "^Spatial")) {
    df <- df@data
  }
  field <- names(df)[grepl(x = names(df),
                                  pattern = paste0("^", field, "$"),
                                  ignore.case = ignore.case.field)]
  values <- df[grepl(x = df[, field], pattern = regex, ignore.case = ignore.case.regex), field] %>% unique()
  return(values)
}
