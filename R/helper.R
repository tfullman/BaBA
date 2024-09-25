#### Helper functions used by BaBA()

#' Helper function for calculating Euclidean distance between locations
#'
#' \code{calc_dist} is a helper function used by the \code{\link{BaBA}} function
#' to simplify calculation of Euclidean distance between subsequent telemetry
#' locations.
#'
#' @param x.start x-coordinate of the starting location.
#' @param x.end x-coordinate of the ending location.
#' @param y.start y-coordinate of the starting location.
#' @param y.end y-coordinate of the ending location.
#'
#' @details The coordinate reference system for the arguments is assumed to be
#'   the same as that for all location data in the \code{animal} argument of
#'   \code{\link{BaBA_caribou}}, matching the \code{crs} argument. Euclidean
#'   distance is simply calculated using the Pythagorean theorem and does not
#'   take into account the curvature of the earth. It is thus most appropriate
#'   for use at local scales and may not be valid if the time difference between
#'   locations is excessively large and the animal exhibits large displacement
#'   between locations.
#'
#' @return Numeric value indicating the Euclidean distance (in units of the
#'   telemetry locations, often meters or degrees) between the input locations.
#' 
calc_dist <- function(x.start, x.end, y.start, y.end){
  sqrt((x.end - x.start)^2 + (y.end - y.start)^2)
}


## Helper function used by BaBA_default() to increase movement segment by one
## points before and one point after the focused encounter
movement.segment.b <- function(animal, pt1, pt2) {
  pts_tmp <- animal[animal$ptsID >= pt1 - 1 & animal$ptsID <= pt2 + 1, ]
  pts_comb <- dplyr::summarize(pts_tmp, do_union = FALSE)
  segments_out <- sf::st_cast(pts_comb, to = 'LINESTRING')
  return(segments_out)
}


