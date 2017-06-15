#' Read in a lookup table from .CSV or .XLSX
#' @description Reads in either a .CSV using \code{read.csv()} or a sheet from a .XLSX file using \code{readxl::read_excel()}, limits it to specified columns/variables, and returns only distinct rows.
#' @param filepath A string specifying either the path to the folder containing the filename provided in \code{filename} OR the full filepath to a file, including the filename and extension. If unprovided, the current working directory will be used.
#' @param filename An optional string specifying the filename and extension of the file to open from the location \code{filepath}. Only use this if the filename and extension aren't included in \code{filepath}. Defaults to \code{NULL}.
#' @param sheet An optional string or numeric index value to be passed to \code{readxl::read_excel()} specifying the sheet from the .XLSX workbook to read in. Do not use if reading in a .CSV. Defaults to \code{NULL}.
#' @param ... Strings corresponding to column/variable names in the source file. Must provide at least one string.
#' @return A data frame of the unique rows/observations in the source file for the given columns/variables.
#' @examples
#' # Read in a table of data about ecological sites and soils and create a lookup table with the columns "Ecological.Site" and "Soil.Map.Unit.Component"
#' ecosite.lookup <- read.lut(filename = "ecosite_soil_lookup.csv",
#'                            "Ecological.Site",
#'                            "Soil.Map.Unit.Component")
#' @export
read.lut <- function(filepath = NULL,
                     ...,
                     filename = NULL,
                     sheet = NULL) {
  fieldnames <- list(...)
  if (length(fieldnames) < 1) {
    stop("Please provide at least one string column/variable name to restrict the lookup table to.")
  }

  if (is.null(filepath)) {
    filepath <- getwd()
  }

  if (!is.null(filename)) {
    filepath <- paste0(filepath, "/", filename)
  }

  if (grepl(x = filepath, pattern = "\\.csv$", ignore.case = T)) {
    lut.raw <- read.csv(filepath, stringsAsFactors = F)
  } else if (grepl(x = filepath, pattern = "\\.xlsx$", ignore.case = T)) {
    lut.raw <- readxl::read_excel(path = filepath,
                                  sheet = sheet)
  } else {
    stop("Must include either .csv/.CSV or .xlsx/.XLSX as a file extension")
  }

  if (length(fieldnames) < length(fieldnames %in% names(lut.raw))) {
    message("The following fieldnames were not found in the lookup table. Check spelling and case.")
    stop(paste(fieldnames[fieldnames %in% names(lut.raw)], collapse = ", "))
  }

  lut <- lut.raw[, unlist(fieldnames)]

  return(dplyr::distinct(lut))
}
