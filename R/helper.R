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
#' @param mov_seg \href{https://cran.r-project.org/package=sf}{\code{sf}}
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


#' Calculate headings of one or more linear steps
#'
#' Helper function used by \code{\link{BaBA_caribou}} to calculate absolute
#' angles between a movement step or barrier segment and the x-axis to provide
#' an indication of heading. This is adapted from code in the \code{as.ltraj}
#' function of
#' \href{https://cran.r-project.org/package=adehabitatLT}{\code{adehabitatLT}}.
#'
#' @param x \href{https://cran.r-project.org/package=sf}{\code{sf}}
#'   \code{POINT} object representing points of movement locations
#'   during an encounter or segment points along a linear barrier. This should
#'   have columns labelled \code{x} and \code{y} indicating the xy coordinates
#'   of feature \code{x}.
#'
#' @return Numeric vector with length equal to \code{nrow(x)} indicating the
#'   absolute angles in radians between each pair of points and the x-axis. The
#'   final value will always be NA as an angle cannot be calculated for the last
#'   point in a sequence.
#' 
calc_angle <- function(x){
  x1 <- x[-1, ]
  x2 <- x[-nrow(x), ]
  dist <- c(sqrt((x1$x - x2$x)^2 + (x1$y - x2$y)^2), NA)
  dx <- c(x1$x - x2$x, NA)
  dy <- c(x1$y - x2$y, NA)
  abs.angle <- ifelse(dist < 1e-07, NA, atan2(dy, dx))
  return(abs.angle)
}


#' Check whether an angle is within a specified range
#'
#' Helper function used by \code{\link{BaBA_caribou}} to check whether a given
#' angle is within the specified range of values. Code adapted from
#' \href{https://stackoverflow.com/questions/66799475/how-to-elegantly-find-if-an-angle-is-between-a-range}{StackOverflow}.
#'
#' @param x Numeric value indicating the angle to be checked, in units of
#'   degrees on a 360 degree scale. May be of class \code{circular} or
#'   \code{numeric}.
#' @param lower Numeric value indicating the lower end of the angular range to
#'   be evaluated. Units, range, and class as with \code{x}.
#' @param upper Numeric value indicating the upper end of the angular range to
#'   be evaluated. Units, range, and class as with \code{x}.
#'
#' @return Logical value indicating whether \code{x} is within the interval
#'   \code{c(lower, upper)}.
#' 
angle_in_range <- function(x, lower, upper){
  (x - lower) %% 360 <= (upper - lower) %% 360
}


#' Project a point along an angle and distance
#'
#' Helper function used by \code{\link{BaBA_caribou}} to project a point
#' elsewhere in space given an initial point, angle for projection, and distance
#' to move the point. This is adapted from the \code{CreateSegmentAngle}
#' function in the
#' \href{https://cran.r-project.org/package=LearnGeom}{\code{LearnGeom}}
#' package.
#'
#' @param pt \code{data.frame} with one row and columns \code{x} and \code{y}
#'   indicating the location of the point to be projected into space.
#' @param angle Numerical value indicating the angle, in units of degrees,
#'   indicating the direction along which \code{pt} should be projected.
#' @param len Numerical value indicating the distance which {pt} should be
#'   moved, in the same units as \code{pt}.
#'
#' @return \code{data.frame} with two rows and two columns. Columns indicate the
#'   x and y coordinates of the points. The first row represents the location of
#'   \code{pt}. The second row indicates the location of the newly projected
#'   point \code{len} units away from \code{pt} at a bearing of \code{angle}.
#'
line_extend <- function (pt, angle, len){
  P1 = pt
  angle = pi * angle / 180
  P2 = pt + c(len * cos(angle), len * sin(angle))
  Segment = rbind(P1, P2)
  return(Segment)
}


## Helper function used by BaBA_default() to increase movement segment by one
## points before and one point after the focused encounter
movement.segment.b <- function(animal, pt1, pt2) {
  pts_tmp <- animal[animal$ptsID >= pt1 - 1 & animal$ptsID <= pt2 + 1, ]
  pts_comb <- dplyr::summarize(pts_tmp, do_union = FALSE)
  segments_out <- sf::st_cast(pts_comb, to = 'LINESTRING')
  return(segments_out)
}


