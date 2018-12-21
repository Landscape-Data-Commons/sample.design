#' Comparing two design objects
#' @description Takes two different design objects formatted for \code{spsurvy::grts()} and returns a data frame detailing the differences between them as errors
#' @param design1 List (see: \code{spsurvey::grts()}) or point allocation data frame (see: \code{sample.design::read_panels()}) describing a design expected to be identical to \code{design2}.
#' @param design2 List (see: \code{spsurvey::grts()}) or point allocation data frame (see: \code{sample.design::read_panels()}) describing a design expected to be identical to \code{design1}.
#' @return Data frame with three variables: panel, stratum, and error
#' @export
compare_designs <- function(design1,
                            design2) {
  if (class(design1) == "list") {
    design1 <- sample.design::design_dataframe(design1)
  }
  if (class(design2) == "list") {
    design2 <- sample.design::design_dataframe(design2)
  }

  names(design1) <- tolower(names(design1))
  names(design2) <- tolower(names(design2))

  if (!"stratum" %in% names(design1)) {
    stop("design1 must contain a variable called 'stratum'")
  }
  if (!"stratum" %in% names(design2)) {
    stop("design2 must contain a variable called 'stratum'")
  }

  errors <- data.frame("panel" = NA,
                       "stratum" = NA,
                       "error" = NA,
                       stringsAsFactors = FALSE)

  # Are there strata that straight-up don't appear?
  if (any(!design1[["stratum"]] %in% design2[["stratum"]])) {
    strata_only_design1 <- design1[["stratum"]][!design1[["stratum"]] %in% design2[["stratum"]]]

    errors <- rbind(errors,
                    data.frame("panel" = NA,
                               "stratum" = strata_only_design1,
                               "error" = "Stratum appears in design1 but not design2"))
  }

  if (any(!design2[["stratum"]] %in% design1[["stratum"]])) {
    strata_only_design2 <- design2[["stratum"]][!design2[["stratum"]] %in% design1[["stratum"]]]

    errors <- rbind(errors,
                    data.frame("panel" = NA,
                               "stratum" = strata_only_design2,
                               "error" = "Stratum appears in design2 but not design1"))
  }

  # What about panels?
  if (any(!names(design1) %in% names(design2))) {
    panels_only_design1 <- names(design1)[!names(design1) %in% names(design2)]
    errors <- rbind(errors,
                    data.frame("panel" = panels_only_design1,
                               "stratum" = NA,
                               "error" = "Panel appears in design1 but not design2"))
  }
  if (any(!names(design2) %in% names(design1))) {
    panels_only_design2 <- names(design2)[!names(design2) %in% names(design1)]
    errors <- rbind(errors,
                    data.frame("panel" = panels_only_design2,
                               "stratum" = NA,
                               "error" = "Panel appears in design2 but not design1"))
  }

  # And now we'll get them merged so we can compare points-per-stratum-per-panel
  # This'll reduce them to only strata and panels they share in common

  comparison <- merge(x = tidyr::gather(design1,
                                        key = "panel",
                                        value = "point_count_design1",
                                        -stratum),
                      y = tidyr::gather(design2,
                                        key = "panel",
                                        value = "point_count_design2",
                                        -stratum))

  comparison[["error"]][comparison[["point_count_design1"]] > comparison[["point_count_design2"]]] <- "Design1 has more points than design2"
  comparison[["error"]][comparison[["point_count_design2"]] > comparison[["point_count_design1"]]] <- "Design2 has more points than design1"

  errors <- rbind(errors,
                  comparison[, c("panel", "stratum", "error")])

  return(errors[!is.na(errors[["error"]]),])
}
