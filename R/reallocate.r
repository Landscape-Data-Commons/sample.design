#' Reallocate points in an existing design_object.
#' @description Allows the surgical reallocation of points in an existing design_object output from allocate.panels(). This can only be done one stratum at a time, so if multiple strata need to be changed, they'll have to be done in separate calls of this function. Likewise, changes are applied to all the panels provided in \code{panels}, so if different changes are required in different panels for a stratum, each of those will need to be a separate call.
#' @param design_object A list structured according to the requirements of \code{spsurvey::grts()}. May be the output from \code{allocate.panels()}.
#' @param stratum A character string specifying the name of the stratum to change in \code{design_object}. If \code{grep_stratum} is \code{TRUE}, then this will be treated as a regular expression.
#' @param grep_stratum Logical. If \code{TRUE} then the value in \code{stratum} will be used as a regular expression to match the names of strata in \code{design_object}. Defaults to \code{FALSE}.
#' @param ignore_case_stratum Logical. If \code{TRUE} and \code{grep_stratum} is \code{TRUE} then the \code{grepl()} will be case insensitive. Defaults to \code{FALSE}.
#' @param bookend_stratum Logical. If \code{TRUE} and \code{grep_stratum} is \code{TRUE} then the regular expression passed to \code{grepl()} will start with \code{"^"} and end with \code{"$"} to prevent partial matches with the value in \code{stratum}. Defaults to \code{FALSE}.
#' @param panels A character string or vector of character strings specifying the panels in the design_object to alter. If \code{grep_panels} is \code{TRUE}, then these will be treated as regular expressions.
#' @param grep_panels Logical. If \code{TRUE} then the values in \code{panels} will be used as regular expressions to match the names of panels in \code{design_object}. Defaults to \code{FALSE}.
#' @param ignore_case_panels Logical. If \code{TRUE} and \code{grep_panels} is \code{TRUE} then the \code{grepl()} will be case insensitive. Defaults to \code{FALSE}.
#' @param bookend_panels Logical. If \code{TRUE} and \code{grep_panels} is \code{TRUE} then the regular expression passed to \code{grepl()} will start with \code{"^"} and end with \code{"$"} to prevent partial matches with any of the values in \code{panels}. Defaults to \code{FALSE}.
#' @param base_point_change An optional numeric value representing net change in base point count for the stratum and panels specified, e.g. 3 or -4. Defaults to \code{NULL}.
#' @param base_point_set An optional numeric value representing absolute the base point count to replace the current value for the stratum and panels specified. This will override \code{base_point_change} if both are provided. Defaults to \code{NULL}.
#' @param recalc_oversample Logical. If \code{TRUE} then the oversample point counts will be recalculated for the specified stratum and panels. Defaults to \code{FALSE}.
#' @param oversample_proportion A numeric value representing the proportion of base points that should be drawn as oversample points per panel within the specified stratum. Only used if \code{recalc_oversample} is \code{TRUE}. Defaults to \code{0.25}.
#' @param oversample_min A numeric value representing the minimum number of points per stratum per panel. Only used if \code{recalc_oversample} is \code{TRUE} and it is larger than the value calculated using \code{oversample_proportion}. Defaults to \code{3}.
#' @return The input \code{design_object} with the relevant base and oversample point values modified.
#' @export
reallocate <- function(design_object = list(),
                       stratum = character(),
                       grep_stratum = FALSE,
                       ignore_case_stratum = FALSE,
                       bookend_stratum = FALSE,
                       panels = c("Year1", "Year2", "Year3", "Year4", "Year5"),
                       grep_panels = FALSE,
                       ignore_case_panels = FALSE,
                       bookend_panels = FALSE,
                       base_point_change = NULL,
                       base_point_set = NULL,
                       recalc_oversample = FALSE,
                       oversample_proportion = 0.25,
                       oversample_min = 3
){
  if (!is.numeric(base_point_change) & !is.numeric(base_point_set)) {
    stop("At least one of either base_point_change and base_point_set must be specified as a numeric value.")
  }

  if (grep_stratum) {
    stratum <- find_name(obj = design_object,
                         name = stratum,
                         ignore.case.name = ignore_case_stratum,
                         bookend = bookend_stratum,
                         multiple = FALSE)
  }


  ## Get the values that we'll be tweaking
  current_stratum_object <- design_object[[stratum]]

  if (grep_panels) {
    panels <- find.name(obj = current_stratum_object,
                        name = panels,
                        ignore.case.name = ignore_case_panels,
                        bookend = bookend_panels,
                        multiple = TRUE)
  }


  ## Get the particular panels
  for (panel in panels) {
    current.panel <- current_stratum_object$panel[[panel]]

    ## Modify the panels
    current.panel <- current.panel + base_point_change

    if (!is.null(base_point_set)) {
      current.panel <- base_point_set
    }

    ## Write the panels back in
    current_stratum_object$panel[[panel]] <- current.panel
  }


  ## If oversample is being recalculated, recalculate it
  if (recalc_oversample) {
    current_stratum_object$over <- 0
    for (panel in names(current_stratum_object$panel)) {
      current_stratum_object$over <- ceiling(pmax(oversample_min, current_stratum_object$panel[[panel]] * oversample_proportion)) + current_stratum_object$over
    }
  }

  ## Write the altered piece of the design_object back in
  design_object[[stratum]] <- current_stratum_object

  ## Return the altered design_object
  return(design_object)
}
