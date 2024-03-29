% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/balance_around.R
\name{get_closest}
\alias{get_closest}
\title{Select points that most closely approximate the distribution of another set of points}
\usage{
get_closest(
  existing_points,
  template_points,
  strata_polygons = NULL,
  stratafield = NULL,
  projection = "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0",
  iteration_limit = 5000
)
}
\arguments{
\item{existing_points}{Point sf object. The points you would like to select from by comparing against \code{template_points}.}

\item{template_points}{Point sf object. The points you would like to compare against \code{existing_points} in order to select a subset of those that most closely resemble the distribution of the template points.}

\item{strata_polygons}{Optional polygon sf object. Polygons assigned a variable with a name \code{stratafield} that contains the membership information (e.g. strata) to assign to \code{existing_points} and \code{template_points}. If \code{NULL} then no assignment will be attempted. Defaults to \code{NULL}.}

\item{stratafield}{Character string. If \code{strata_polygons} is not \code{NULL}, the name of the variable in \code{strata_polygons} that contains the membership information. Otherwise, the name of the variable in both \code{template_points} and \code{existing_points} that contains the membership information. If \code{NULL} then the points will be considered to belong to a single group. Defaults to \code{NULL}.}

\item{projection}{Character string. The projection to force all spatial objects into. Defaults to NAD83, \code{sp::CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")}.}

\item{iteration_limit}{Numeric. The maximum number of iterations to attempt to sort before giving up. Defaults to \code{5000}.}
}
\value{
A spatial points data frame made by trimming \code{existing_points} down to the points that most closely approximate the distribution of \code{template_points} while also containing the same number of points as \code{template_points}. It will be in the projection specified by \code{projection}.
}
\description{
Sometimes you have a large collection of points which are not randomly distributed or spatially balanced and you would like a subset that more or less do. Given a template of points that are distributed the way you would like, this will return the closest existing point to each. This can be done taking into account membership in a group, either by having assigned it as a variable in both sets of points or by providing polygons that can be used to assign membership. By default, no stratification/membership is taken into account.
}
