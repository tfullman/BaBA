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

#' Calculate straightness of movement segment
#'
#' Helper function used by the \code{\link{BaBA_caribou}} function to calculate
#' straightness of movement during an encounter as the ratio of Euclidean to
#' path distance.
#'
#' @param mov_seg \code{\href{https://cran.r-project.org/package=sf}{sf}}
#'   \code{POINT} object representing locations that are part of the encounter
#'   being analyzed.
#'
#' @return Numeric value indicating the straightness of movement during the
#'   encounter. Values range between 0-1, with values closer to 0 indicating
#'   more sinuous movement.
#' 
strtns <- function(mov_seg) {
  locs_tmp <- sf::st_coordinates(mov_seg)
  if (sum(duplicated(mov_seg$date)) > 0 ) {
    straightness <- NA
  } else if(nrow(locs_tmp) == 1){
    straightness <- NA
  } else {
    ## Calculate Euclidean distance between the first and last point
    euc_dist <- 
      as.numeric(
        calc_dist(x.start = locs_tmp[1,1], x.end = locs_tmp[nrow(locs_tmp),1],
                  y.start = locs_tmp[1,2], y.end = locs_tmp[nrow(locs_tmp),2]))
    
    ## Calculate path distance as the sum of all step lengths
    mov_seg$dist <- NA
    for(j in 2:nrow(mov_seg)){
      mov_seg$dist[j] <- calc_dist(x.start = locs_tmp[j - 1, 1],
                                   x.end = locs_tmp[j, 1],
                                   y.start = locs_tmp[j - 1, 2],
                                   y.end = locs_tmp[j, 2])
    }
    path_dist <- sum(mov_seg$dist, na.rm = TRUE)
    
    ## Calculate straightness as the ratio of Euclidean to path distance.
    ## Straightness ranges from 0 to 1 with values closer to 0 being more
    ## sinuous.
    straightness <- euc_dist/path_dist
  }
  
  return(straightness)
}


## Helper function used by BaBA_default() to increase movement segment by one
## points before and one point after the focused encounter
movement.segment.b <- function(animal, pt1, pt2) {
  pts_tmp <- animal[animal$ptsID >= pt1 - 1 & animal$ptsID <= pt2 + 1, ]
  pts_comb <- dplyr::summarize(pts_tmp, do_union = FALSE)
  segments_out <- sf::st_cast(pts_comb, to = 'LINESTRING')
  return(segments_out)
}


