#' Reallocate points in an existing design object.
#' @description Allows the surgical reallocation of points in an existing design object output from allocate.panels(). This can only be done one stratum at a time, so if multiple strata need to be changed, they'll have to be done in separate calls of this function. Likewise, changes are applied to all the panels provided in \code{panels}, so if different changes are required in different panels for a stratum, each of those will need to be a separate call.
#' @param design.object A list structured according to the requirements of \code{spsurvey::grts()}. May be the output from \code{allocate.panels()}.
#' @param stratum A character string specifying the name of the stratum to change in \code{design.object}. If \code{grep.stratum} is \code{T}, then this will be treated as a regular expression.
#' @param grep.stratum Logical. If \code{T} then the value in \code{stratum} will be used as a regular expression to match the names of strata in \code{design.object}. Defaults to \code{F}.
#' @param ignore.case.stratum Logical. If \code{T} and \code{grep.stratum} is \code{T} then the \code{grepl()} will be case insensitive. Defaults to \code{F}.
#' @param bookend.stratum Logical. If \code{T} and \code{grep.stratum} is \code{T} then the regular expression passed to \code{grepl()} will start with \code{"^"} and end with \code{"$"} to prevent partial matches with the value in \code{stratum}. Defaults to \code{F}.
#' @param panels A character string or vector of character strings specifying the panels in the design object to alter. If \code{grep.panels} is \code{T}, then these will be treated as regular expressions.
#' @param grep.panels Logical. If \code{T} then the values in \code{panels} will be used as regular expressions to match the names of panels in \code{design.object}. Defaults to \code{F}.
#' @param ignore.case.panels Logical. If \code{T} and \code{grep.panels} is \code{T} then the \code{grepl()} will be case insensitive. Defaults to \code{F}.
#' @param bookend.panels Logical. If \code{T} and \code{grep.panels} is \code{T} then the regular expression passed to \code{grepl()} will start with \code{"^"} and end with \code{"$"} to prevent partial matches with any of the values in \code{panels}. Defaults to \code{F}.
#' @param base.point.change An optional numeric value representing net change in base point count for the stratum and panels specified, e.g. 3 or -4. Defaults to \code{NULL}.
#' @param base.point.set An optional numeric value representing absolute the base point count to replace the current value for the stratum and panels specified. This will override \code{base.point.change} if both are provided. Defaults to \code{NULL}.
#' @param recalc.oversample Logical. If \code{T} then the oversample point counts will be recalculated for the specified stratum and panels. Defaults to \code{F}.
#' @param oversample.proportion A numeric value representing the proportion of base points that should be drawn as oversample points per panel within the specified stratum. Only used if \code{recalc.oversample} is \code{T}. Defaults to \code{0.25}.
#' @param oversample.min A numeric value representing the minimum number of points per stratum per panel. Only used if \code{recalc.oversample} is \code{T} and it is larger than the value calculated using \code{oversample.proportion}. Defaults to \code{3}.
#' @return The input \code{design.object} with the relevant base and oversample point values modified.
#' @export
reallocate <- function(design.object = list(),
                       stratum = character(),
                       grep.stratum = F,
                       ignore.case.stratum = F,
                       bookend.stratum = F,
                       panels = c("Year1", "Year2", "Year3", "Year4", "Year5"),
                       grep.panels = F,
                       ignore.case.panels = F,
                       bookend.panels = F,
                       base.point.change = NULL,
                       base.point.set = NULL,
                       recalc.oversample = F,
                       oversample.proportion = 0.25,
                       oversample.min = 3
){
  if (!is.numeric(base.point.change) & !is.numeric(base.point.set)) {
    stop("At least one of either base.point.change and base.point.set must be specified as a numeric value.")
  }

  if (grep.stratum) {
    stratum <- find.name(obj = design.object,
                         name = stratum,
                         ignore.case.name = ignore.case.stratum,
                         bookend = bookend.stratum,
                         multiple = F)
  }


  ## Get the values that we'll be tweaking
  current.stratum.object <- design.object[[stratum]]

  if (grep.panels) {
    panels <- find.name(obj = current.stratum.object,
                        name = panels,
                        ignore.case.name = ignore.case.panels,
                        bookend = bookend.panels,
                        multiple = T)
  }


  ## Get the particular panels
  for (panel in panels) {
    current.panel <- current.stratum.object$panel[[panel]]

    ## Modify the panels
    current.panel <- current.panel + base.point.change

    if (!is.null(base.point.set)) {
      current.panel <- base.point.set
    }

    ## Write the panels back in
    current.stratum.object$panel[[panel]] <- current.panel
  }


  ## If oversample is being recalculated, recalculate it
  if (recalc.oversample) {
    current.stratum.object$over <- 0
    for (panel in names(current.stratum.object$panel)) {
      current.stratum.object$over <- ceiling(pmax(oversample.min, current.stratum.object$panel[[panel]] * oversample.proportion)) + current.stratum.object$over
    }
  }

  ## Write the altered piece of the design object back in
  design.object[[stratum]] <- current.stratum.object

  ## Return the altered design object
  return(design.object)
}
