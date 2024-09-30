## Wrapper function - to be developed at the end
# BaBA <- function(caribou = FALSE, ...){
#   if(caribou){
#     BaBA_caribou()
#   } else{
#     BaBA_default()
#   }
# }


#'Barrier Behavior Analysis (BaBA)
#'
#'This function classifies animal movement near linear barriers into 6 different
#'behavioral categories: quick cross, average movement, bounce, back-and-forth,
#'trace, and trapped. The classification unit of BaBA first identifies
#'\emph{'encounter events'}, defined by a series of continuous locations that
#'fall within the barrier buffers. Each event will then be classified into one
#'of the 6 barrier behaviors. Even though many arguments have default settings,
#'we recommend users adjust these values to best fit the ecology and management
#'goals of the studied species. This is the default version of the \code{BaBA}
#'function, as developed by Xu et al. (2021).
#'
#'@param animal An \code{sf POINT} object representing animal telemetry
#'  locations with a column named \code{"Animal.ID"} that identifies each
#'  individual; a column named \code{"date"} in \code{"POSIXct"} format.
#'  \code{CRS} has to be the same as \code{barrier}, preferably in a projected
#'  UTM coordinate system. More see Details.
#'@param barrier An \code{sf LINESTRING} or \code{MULTILINESTRING} object
#'  showing barrier locations in the area overlapped with \code{animal} movement
#'  data. In the same spatial projection as \code{animal}.
#'@param d Barrier buffer size in meters if \code{barrier} has a projected
#'  coordinate system CRS, in the units of \code{barrier} otherwise.
#'@param interval Time interval of the movement data (unit specified in
#'  \code{units}). If not specified, BaBA will use the most frequent time
#'  difference between steps as \code{interval}, which might affect result
#'  accuracy if the data has an irregular time interval or a lot of missing
#'  data.
#'@param b_time Maximum duration, in the same unit as \code{interval}, that an
#'  encounter event would be considered as a short event (\emph{'bounce'} or
#'  \emph{'quick cross'}). Must be evenly divisible by \code{interval}.
#'@param p_time Minimum duration, in the same unit as \code{interval}, that an
#'  encounter event would be considered as a \emph{'trapped'} condition. Must be
#'  evenly divisible by \code{interval}.
#'@param w The length of time, in the same unit as \code{interval}, to include
#'  around the encounter event to calculate average movement straightness using
#'  a moving window. Locations included are all locations within \code{w/2}
#'  before the first location of the encounter event and \code{w/2} after the
#'  last location of the event. More see Details.
#'@param tolerance The maximum duration, in the same unit as \code{interval}, to
#'  allow points that are outside of barrier buffer but between 2 sets of points
#'  within the buffer to be included, so that all points are considered as one
#'  continous \emph{encounter event}. Useful when movement data has a very high
#'  temporal resolution. Must be evenly divisible by \code{interval}.
#'@param units The temporal units of \code{interval}, \code{b_time},
#'  \code{p_time}, \code{w}, and \code{tolerance}. One of "secs", "mins",
#'  "hours", "days", and "weeks".
#'@param max_cross The maximum number of crosses in an encounter event allowed
#'  for in tracing and back-and-forth behavior. More see Details.
#'@param sd_multiplier A \code{numeric} value to determine the numbers of
#'  standard deviations used that defines the normal range of movement
#'  straightness. Default is 1. More in Details.
#'@param exclude_buffer \code{Logical}. Whether to consider movement locations
#'  within barrier buffers when calculating average movement straightness. More
#'  see Details.
#'@param round_fixes \code{Logical}. Indicates whether elements of the
#'  \code{"date"} column of \code{animal} should be rounded to the nearest
#'  \code{units}. This is intended to account for minor variations in fix
#'  acquisition time that may lead to some intervals being slightly longer than
#'  \code{interval}. Note that this is not a replacement for cleaning the
#'  \code{animal} data.frame and removing repetitive timesteps as described in
#'  Details.
#'@param export_images \code{Logical}. If \code{TRUE}, will export snapshot of
#'  event locations as \code{.png} files named with the classification of the
#'  event, the \code{Animal.ID} of the animal, the \code{burstID} of the event,
#'  and potential suffix specified by argument \code{img_suffix}.
#'@param img_path \code{character} string indicating the name of the folder to
#'  which images should be exported when \code{export_images == TRUE}. The
#'  folder will be created in the working directory if it does not exist.
#'  Defaults to "event_imgs".
#'@param img_suffix \code{character} string to be added at the end of image file
#'  names written when \code{export_images == TRUE}.
#'
#'@details \code{BaBA} works better with cleaned \code{animal} data. Make sure
#'  your \code{date} column is cleaned so that the time difference between any
#'  step is evenly divisible by \code{interval}. It is OK to have missing steps
#'  but remove repetitive timesteps and bursts of locations smaller than
#'  \code{interval} before running \code{BaBA}. Minor fix acquisition
#'  differences can be addressed by setting \code{round_fixes} to \code{TRUE}.
#'
#'  The barrier buffer distance \code{d} represents the distance on both sides
#'  of the barrier lines within which animal movement locations are considered
#'  'encounters' with the barrier, and continuous locations form an
#'  \emph{'encounter event'}. This might be the most important parameter to set
#'  the BaBA. The buffer distance will affect numbers and durations of
#'  trajectories identified as 'encounter events'. For species with different
#'  movement capacities or ecological attributes, barrier effect distance might
#'  be different. You can decide the distance either by \emph{a priori}
#'  knowledge or testing a range of distance and comparing the results.
#'
#'  The average straightness is compared to the straightness of a
#'  to-be-classified encounter event in order to determine whether the event is
#'  \emph{'trace'} (more straight than average), \emph{'back-n-forth'} (less
#'  straight than average), or \emph{'normal movement'} (similar to average).
#'  The reason to apply a moving window method with a width of \code{w} to
#'  calculate the average straightness is because some animal movements, such as
#'  migrations, show great seasonal variations. The moving window method
#'  calculates localized average straightness around the time when the
#'  to-be-classified encounter event occurs.
#'
#'  The default average straightness calculation considers movement segments
#'  \code{w/2} days before and after the focal encounter event, which can
#'  include movement locations within and outside of the barrier buffers. One
#'  might want to control potential impacts of barriers on movement and only
#'  calculate average movement straightness based on locations outside of fence
#'  buffers. Use \code{exclude_buffer} to indicate whether or not to include
#'  movement locations inside barrier buffers. The default is \code{FALSE}
#'  because in the sample study fence density is relatively high and animals
#'  fall in barrier buffers much of the time. If excluded, not enough continuous
#'  locations within the time window would be outside of the buffer.
#'
#'  The average and standard deviation of straightness measurements calculated
#'  by the moving window method are used to define the normal range of movement
#'  straightness. Any encounter event with a straightness < (average
#'  straightness + \code{sd_multiplier} * sd straightness) would be classified
#'  as \emph{'trace'} (more straight than normal), < (average straightness -
#'  \code{sd_multiplier} * sd straightness) would be \emph{'back-and-forth'}
#'  (more tortuous than normal), and in between would be \emph{'average
#'  movement'}. To make the standard more strict for an event to be classified
#'  as not 'normal' (i.e. \emph{'back-and-forth'} or \emph{'trace'}), use a
#'  larger number for the \code{sd_multiplier}.
#'
#'  When the barriers are curvy or the temporal resolution of the movement data
#'  is coarse, straight lines between movement locations (movement segments)
#'  might appear to be crossing the barrier even if the locations are on the
#'  same side of the barrier. \code{max_cross} allows some numbers of
#'  intersections between movement segments and barriers to be included in the
#'  \emph{'trace'} and \emph{'back-and-forth'} behavior. When the intersections
#'  are larger than \code{max_cross}, the encounter event will be classified as
#'  \emph{'unknown'}.

#'@return \code{BaBA} returns a list with two components:
#'
#'  \code{$classification} is a \code{data.frame}. Each row represents one
#'  encounter event. \code{AnimalID} is the ID of the individual; \code{burstID}
#'  is the ID of the encounter event and it matches the ID of the exported
#'  images; \code{easting} and \code{northing} are the coordinates of the
#'  starting location of the encounter event; \code{start_time},
#'  \code{end_time}, and \code{duration} describe the temporal characteristics
#'  of the encounter event; \code{cross} is the number of intersections between
#'  the encounter event trajectory and the barrier; \code{straightness} is the
#'  straightness of the encounter event trajectory and is \code{NaN} when the
#'  encounter event is \emph{'bounce} or \emph{'quick cross'}; finally,
#'  \code{eventTYPE} is the barrier behavior classification of the encounter
#'  event.
#'
#'  \code{$encounters} is a \code{sf POINT} object showing the locations of
#'  classified encounter events represented by the starting location.

#'@export
#'
#'@author Wenjing Xu \email{wenjing.xu@berkeley.edu} and Valentine Herrmann
#'  \email{HerrmannV@si.edu}
#'
#'@references Xu W, Dejid N, Herrmann V, Sawyer H, Middleton AD. 2021. Barrier
#'  Behaviour Analysis (BaBA) reveals extensive effects of fencing on
#'  wide-ranging ungulates. Journal of Applied Ecology 58: 690-698.
#'  https://doi.org/10.1111/1365-2664.13806.
#'
#' @examples
#' \dontrun{
#' library(BaBA)
#' 
#' # load sample data.
#' data("pronghorn")
#' data("fences")
#' 
#' # Animal data must be an "sf POINT" object created by the "sf" package. 
#' # individual ID column should be named "Animal.ID" and timestamp column should be named "date"
#' class(pronghorn)
#' names(pronghorn)
#' 
#' # run BaBA on the pronghorn data
#' results_prong <- BaBA(animal = pronghorn, barrier = fences, d = 110, max_cross = 4)
#' 
#' # View BaBA results
#' head(results_prong$classification)
#' 
#' # plot encounter event locations
#' plot(fences)
#' plot(results_prong$encounters, add = TRUE)
#' 
#' # write the encounters as shapefile, changing the date column to character to
#' # avoid the shapefile dropping the time
#' sf::st_write(dplyr::mutate(results_prong$encounters, date = as.character(date)),
#'              ".", "encountersPRONG", driver = "ESRI Shapefile")
#' 
#' # export event images to visually check event classifications. Using mule deer data as an example
#' data("muleDeer")
#' results_deer <- BaBA(animal = muleDeer, barrier = fences, d = 90,
#'                      export_images = TRUE, img_suffix = "DEER")
#' }
BaBA_default <-
  function(animal, barrier, d, interval = NULL, b_time = 4, p_time = 36, w = 168,
           tolerance = 0, units = "hours", max_cross = 0,  sd_multiplier = 1,
           exclude_buffer = F, round_fixes = F, export_images = F,
           img_path = "event_imgs", img_suffix = NULL) {
    
    # initial checks ----------------------------------------------------------
    if(export_images) {
      if(!dir.exists(img_path)) dir.create(img_path)
    }
    
    ## prepare parameters and check input
    if (class(animal)[1] != "sf") stop("animal needs to be an sf object")
    if (sf::st_geometry_type(animal)[1] != 'POINT') stop("animal needs to have a POINT geometry")
    if (class(barrier)[1] != "sf") stop("barrier needs to be an sf object")
    if (!(sf::st_geometry_type(barrier)[1] %in% c('MULTILINESTRING', 'LINESTRING'))) stop("barrier needs to have either a LINESTRING or MULTILINESTRING geometry")
    if (!"date" %in% names(animal)) stop("please rename the date column to 'date'")
    if (!"Animal.ID" %in% names(animal)) stop("please rename the individual ID column to 'Animal.ID'")
    if (!(inherits(animal$date, "POSIXct"))) stop("date needs to be 'POSIXct' format")
    if (sum(is.na(animal$date)) > 0) stop("please exclude rows where date is NA")
    
    if(round_fixes){
      interval_per_individual <- tapply(animal$date, animal$Animal.ID, function(x) names(which.max(table(round(as.numeric(diff(x), units = units),0)))))
    } else {
      interval_per_individual <- tapply(animal$date, animal$Animal.ID, function(x) names(which.max(table(as.numeric(diff(x), units = units)))))
    }
    if(is.null(interval)) { ## figure out interval (as the most frequent difference in timestamp) if not provided but give an error if not the same for all individuals
      if(all(interval_per_individual == interval_per_individual[1])) interval <- as.numeric(interval_per_individual[1]) else stop("Not all individuals have been sampled at the same frequency. Run individuals with different intervals seperately, or double-check whether your date column is cleaned.")
    } else {
      if (any(as.numeric(interval_per_individual) > interval, na.rm = T)) stop("BaBA interval needs to be no smaller than the actual data interval. Also double-check whether your date column is cleaned.") 
    }
    
    b <- b_time / interval
    if(b < 1) stop("interval needs to be set no bigger than b_time")
    if (round(b) != b) stop("b_time must be divisible by interval")
    p <- p_time / interval
    if (round(p) != p) stop("p_time must be divisible by interval")
    
    
    # classification step 1: generate encounter event data.frame -------------------------------------------------
    
    ## create point ID by individual
    animal <-
      animal %>% 
      dplyr::arrange(Animal.ID, date) %>% 
      dplyr::group_by(Animal.ID) %>% 
      dplyr::mutate(ptsID = 1:dplyr::n()) %>% 
      dplyr::ungroup()
    
    ## explicitly suppress constant geometry assumption warning by confirming attribute is constant throughout the geometry. See https://github.com/r-spatial/sf/issues/406 for details.
    sf::st_agr(animal) <- 'constant'   
    sf::st_agr(barrier) <- 'constant'
    
    ## create buffer around barrier
    print("locating encounter events...")
    barrier_buffer <- 
      barrier %>% 
      sf::st_buffer(dist = d, nQuadSegs = 5) %>%   ## Note that nQuadSegs is set to 5 as this was the default value for rgeos::gBuffer in previous versions of BaBA
      sf::st_union()
    
    ## extract points that fall inside the buffer
    encounter <- sf::st_intersection(animal, barrier_buffer)
    
    if (nrow(encounter) == 0) stop("no barrier encounter detected.")
    
    ## create unique burstIDs
    for(i in unique(encounter$Animal.ID)){
      if (nrow(encounter %>% dplyr::filter(Animal.ID == i)) == 0) {
        warning(paste0 ("Individual ", i, " has no locations overlapped with the barrier buffer and is eliminated from analysis." ))
        next()
      }
      encounter_i <-
        encounter %>% 
        dplyr::filter(Animal.ID == i) %>% 
        ## add time difference
        dplyr::mutate(
          ## time difference between all points in the buffer
          timediff = c(interval, as.numeric(diff(date), units = units)),
          ## remove the interval from that so when there is no missing point, timediff2 should be 0. If <0, replicated timestamp; if >0, missing timestamp
          timediff2 = round(timediff - interval, digits = 1))
      
      ## if any timediff2 is >= interval but <= tolerance, bring in the missing points from outside the buffer
      if(any(encounter_i$timediff2 >= interval & encounter_i$timediff2 <= tolerance, na.rm = T )) {
        idx_pts_of_interest <- which(encounter_i$timediff2 >= interval & encounter_i$timediff2 <= tolerance)
        
        for(pt in idx_pts_of_interest) {
          ## find out what pts to fetch
          ptsID_of_interest_B <- encounter_i$ptsID[pt]
          ptsID_of_interest_A <- encounter_i$ptsID[pt-1]
          
          ## fetch the points outside of the buffer and placehold timediff as NA and timediff2 as 0
          fetched_pt <- 
            animal %>% 
            dplyr::filter(Animal.ID == i & 
                            ptsID > ptsID_of_interest_A & 
                            ptsID < ptsID_of_interest_B)
          
          if (nrow(fetched_pt) == 0) {  ## if there's no point outside of the buffer between the timestamp that means there's missing data
            ## since the missing data is still within the tolerance, we consider timediff2=0 so the points before and after will be in the same event
            encounter_i$timediff2[pt] <- 0
            next() } 
          else {
            fetched_pt$timediff <- NA
            fetched_pt$timediff2 <- 0 
            ## reset timediff2 of pts_of_interests to 0
            encounter_i$timediff2[pt] <- 0 
            ## append fetched points to each other 
            if(pt == idx_pts_of_interest[1]) {fetched_pts <- fetched_pt} else if (exists("fetched_pts")) { fetched_pts <- rbind(fetched_pts, fetched_pt) } else {fetched_pts <- fetched_pt}
          }
        }
        
        ## append fetched pts
        encounter_i <- rbind(encounter_i, fetched_pts)
        ## recorder animal i's encounter event data.frame
        encounter_i <- encounter_i[order(encounter_i$ptsID), ]
      }
      
      ## do the cumulative sum of the new data.frame based on timediff2, using that as the updated unique burst ID (with animalID) 
      encounter_i$burstID <- paste(i, cumsum(encounter_i$timediff2), sep = "_")
      
      ## save into encounter_complete
      if(i == unique(encounter$Animal.ID[1])) encounter_complete <- encounter_i else encounter_complete <- rbind(encounter_complete, encounter_i)
    }
    
    encounter <- encounter_complete ## save back as encounter (encounter_complete is bigger as it includes extra points that are within tolerance)
    
    
    # classification step 2: classify events -------------------------------------------------
    
    print("classifying behaviors...") 
    ## open progress bar
    pb <- utils::txtProgressBar(style = 3)
    
    ## create empty object that will hold results
    event_df <- NULL
    
    ## run classification procedure for each encounter
    for(i in unique(encounter$burstID)) {
      ## update progressbar
      utils::setTxtProgressBar(pb, which(unique(encounter$burstID) == i)/length(unique(encounter$burstID)))
      
      ## get what we need from the encounter
      encounter_i <- encounter[encounter$burstID == i, ]
      animal_i <- animal[animal$Animal.ID == encounter_i$Animal.ID[1],]
      start_time <- encounter_i$date[1]
      end_time <- encounter_i$date[nrow(encounter_i)]
      duration <-  difftime (end_time, start_time, units = units)
      
      ## calculating straightness of the encounter event
      ## this will be used for median duration events but is output for reference for other events
      straightness_i <- strtns(encounter_i)
      
      
      ### classify short events (bounce and quick cross) ---------------------------------------------------
      
      ## if no more than b*interval, only spend small amount of time in this burst
      if (duration <= b_time) {
        pt.first <- encounter_i$ptsID[1] ## first point in the burst
        pt.last <- encounter_i$ptsID[nrow(encounter_i)]
        
        ## extract movement segment with one point before and one point after the segmentation
        mov_seg_i <- movement.segment.b(animal_i, pt.first, pt.last)
        
        ## count the number of crossings
        int.num <-
          mov_seg_i %>% 
          sf::st_intersection(barrier) %>% 
          sf::st_cast(to = 'MULTIPOINT') %>% 
          sf::st_coordinates() %>% 
          nrow()
        
        ## if no crossing is indicated and both before and after points were missing then we cannot tell if the animal crossed
        if (int.num == 0 & nrow(sf::st_coordinates(mov_seg_i)) != (nrow(encounter_i)+2)) {
          classification <- "unknown"
        } else {
          ## if there was not a crossing, classify as bounce, otherwise quick cross
          classification <- ifelse(int.num == 0, "Bounce", "Quick_Cross")
        }
      }
      
      ## dummy variable to ensure desired plotting and output
      tbd.plot <- 0
      
      if (duration > b_time) {
        
        ### classify trapped events -------------------------------------------------
        
        ## first calculate number of crossings (without looking at extra points like we did for short encounter)
        mov_seg_i <- 
          encounter_i %>% 
          dplyr::summarize(do_union = FALSE) %>% 
          sf::st_cast(to = 'LINESTRING')
        
        int.num <-
          mov_seg_i %>% 
          sf::st_intersection(barrier) %>% 
          sf::st_cast(to = 'MULTIPOINT') %>% 
          sf::st_coordinates() %>% 
          nrow()
        
        ## check if duration is smaller of bigger than p and classify accordingly
        if(duration > p_time) {
          classification <- "Trapped"
        } else {
          classification <- "TBD" ## these will be further classified in the next loop
        }
        
        
        ### classify back-n-forth, trace, and average movement -----------------------------------------
        
        ## process the "TBD" event types as back-n-forth, trace, or average movement
        ## back-n-forth and trace are based on comparing average straightness around the encounter event
        if(classification == 'TBD'){
          
          tbd.plot <- 1
          
          ## remove points that are inside the buffer if user said so
          if (exclude_buffer) {
            animal_i <- animal_i[!animal_i$ptsID %in% encounter$ptsID[encounter$Animal.ID == animal_i$Animal.ID[1]], ]
          }
          
          ## keep only data w/2 units before and w/2 after event
          animal_i <- animal_i[animal_i$date >= start_time - as.difftime(w/2, units = units) & animal_i$date <= end_time +  as.difftime(w/2, units = units), ]
          
          ## identify continuous sections in the remaining movement data
          animal_i$continuousID <- cumsum(abs(c(interval, round(as.numeric(diff(animal_i$date), units = units), digits = 1) - interval))) # abs() is to accommodate potential data points with smaller time intervals
          ## for each continuous sections, calculate straightness of all movements lasting the duration of our event (moving window of the size of the encounter)
          straightnesses_i <- NULL
          for(ii in unique(animal_i$continuousID)) {
            animal_ii <- animal_i[animal_i$continuousID == ii, ]
            
            ## duration of period
            duration_ii <- difftime(animal_ii$date[nrow(animal_ii)], animal_ii$date[1], units = units)
            
            ## calculate straightness only if at least as long as encounter event
            if(duration_ii >= duration) {
              for(iii in 1:(which(animal_ii$date > (animal_ii$date[nrow(animal_ii)] - as.difftime(duration, units = units)))[1] -1)) {
                mov_seg <- animal_ii[iii:(iii + duration/interval), ]
                straightnesses_i <- c(straightnesses_i, strtns(mov_seg))
              }
            }
          }
          
          ## make sure there are enough data to calculate average straightness before and after the encounter event
          ## (w/interval + 1) is the total possible segments if all data are present. 
          ## We define "enough" as at least 1/4 of the total possible segments are present to calculate average straightness.
          if (length(straightnesses_i) >= (w/interval + 1)/4) {
            ## minimum max number possible/2 to calculate sd
            upper <- mean(straightnesses_i) + sd_multiplier * stats::sd(straightnesses_i)
            lower <- mean(straightnesses_i) - sd_multiplier * stats::sd(straightnesses_i)
            if(straightness_i < lower) classification <- ifelse(int.num <= max_cross, "Back_n_forth", "unknown")
            if (straightness_i > upper) classification <- ifelse(int.num <= max_cross, "Trace", "unknown")
            if(straightness_i >= lower & straightness_i <= upper) classification <- "Average_Movement"
          } else {
            classification <- "unknown"
            if(is.null(straightnesses_i)) {straightnesses_i <- NA} ## add this to avoid warning message when plotting
          }
        }
      }
      
      
      ### Consolidate outputs -----------------------------------------------------
      
      ## plot the encounters to check later, if desired
      if (export_images) {
        grDevices::png(paste0(img_path, "/", classification, "_", i, "_", img_suffix, ".png"), width = 6, height = 6, units = "in", res = 90)
        if(tbd.plot == 0){
          plot(mov_seg_i, main = classification, sub = paste("cross =", int.num, ", duration =", round(duration, 0), units))
          plot(barrier_buffer, border = scales::alpha("red", 0.5), lty = "dashed", add = T)
          plot(sf::st_geometry(barrier), col = 'red', lwd = 2, add = TRUE)
          plot(sf::st_geometry(encounter_i), pch = 20, col = "cyan3", type = "o", lwd = 2, add = TRUE)
        } else{
          A <-
            animal_i %>% 
            dplyr::group_by(continuousID) %>% 
            ## Check sample size per group
            dplyr::mutate(n = n()) %>%
            ## Exclude groups with only a single point to avoid errors
            dplyr::filter(n > 1) %>% 
            ## Convert to a line
            dplyr::summarize(do_union = FALSE) %>% 
            sf::st_cast(to = 'LINESTRING')
          plot(sf::st_geometry(A), main = classification,
               sub = paste0("cross = ", int.num, ", duration =", round(duration, 0), ", stri =", round(straightness_i, 2), ", str_mean = ",  round(mean(straightnesses_i), 2), ", str_sd = ",  round(stats::sd(straightnesses_i), 2)))
          plot(sf::st_geometry(barrier_buffer), border = scales::alpha("red", 0.5), lty = "dashed", add = TRUE)
          plot(sf::st_geometry(barrier), col = "red", lwd = 2, add = TRUE)
          plot(sf::st_geometry(encounter_i), pch = 20, col = "cyan3", type = "o", lwd = 2, add = TRUE)
        }
        grDevices::dev.off()
      }
      
      ## combine output
      event_df <- rbind(event_df, data.frame(
        AnimalID = encounter_i$Animal.ID[1],
        burstID = i,
        easting = sf::st_coordinates(encounter_i)[1, 1],
        northing = sf::st_coordinates(encounter_i)[1, 2],
        start_time,
        end_time,
        duration,
        cross = int.num,
        str_i = straightness_i,
        str_mean = ifelse(tbd.plot == 0, NA, mean(straightnesses_i)),
        str_sd = ifelse(tbd.plot == 0, NA, stats::sd(straightnesses_i)),
        eventTYPE = classification,
        stringsAsFactors = F
      ))
    }
    
    ## close progress bar
    close(pb)
    
    
    # finalize data -----------------------------------------------------------
    
    print("creating dataframe...")
    ## clean the encounter data
    encounter_final <- 
      encounter %>% 
      dplyr::filter(!duplicated(burstID)) %>% 
      dplyr::left_join(event_df %>% 
                         dplyr::select(burstID, eventTYPE),
                       by = 'burstID') %>% 
      dplyr::select(Animal.ID, burstID, date, eventTYPE)
    
    ## return output as a list
    return(list(encounters = encounter_final,
                classification = event_df))
  }



#'Barrier Behavior Analysis (caribou-specific)
#'
#'This is a heavily modified adaptation of \code{\link{BaBA_default}}, tailored
#'to barrier analysis for Western Arctic Herd (WAH) caribou in northwestern
#'Alaska. The analysis workflow of Xu et al. (2021), on which
#'\code{BaBA_default} is built was reworked to improve accuracy for the WAH. See
#'details below.
#'
#'@param animal An \code{sf POINT} object representing animal telemetry
#'  locations with a column named \code{"Animal.ID"} that identifies each
#'  individual and a column named \code{"date"} in \code{"POSIXct"} format. The
#'  coordinate reference system of \code{animal} should match \code{barrier},
#'  preferably in a projected coordinate system.
#'@param barrier An \code{sf LINESTRING} or \code{MULTILINESTRING} object
#'  showing barrier locations to be evaluated against \code{animal} movement
#'  data. This should have the same spatial projection as \code{animal}. It also
#'  is assumed to have a column labelled \code{"Name"} that uniquely identifies
#'  each barrier segment, ideally with a single term (without spaces).
#'@param d Barrier buffer size in meters if \code{barrier} has a projected
#'  coordinate system CRS, in the units of \code{barrier} otherwise.
#'@param interval Optional. Numeric value specifying the time interval of the
#'  movement data (in units specified in \code{units}). If not specified, BaBA
#'  will use the most frequent time difference between steps as \code{interval},
#'  which might affect result accuracy if the data has an irregular time
#'  interval or a large amount of missing data.
#'@param tolerance Numeric value indicating the maximum number of points that
#'  can occur outside of the barrier buffer but still be counted part of the
#'  same encounter. Note that this is a different definition of tolerance than
#'  that used in \code{\link{BaBA_default}}, which defined tolerance as the
#'  maximum duration in which points could be outside of barrier buffer and
#'  still count as a single encounter. The two definitions are conceptually
#'  similar, but have different units. Defaults to 0.
#'@param units Character string indicating the temporal units of
#'  \code{interval}. Values should be one of \code{"secs"}, \code{"mins"},
#'  \code{"hours"}, \code{"days"}, or \code{"weeks"}. Defaults to
#'  \code{"hours"}.
#'@param sd_multiplier Numeric value indicating the number of standard
#'  deviations used to define the normal range of movement straightness.
#'  Increasing this value results in widening the confidence interval used to
#'  identify normal behavior. Defaults to 1.
#'@param round_fixes Logical indicator of whether elements of the \code{"date"}
#'  column of \code{animal} should be rounded to the nearest \code{units}. This
#'  is intended to account for minor variations in fix acquisition time that may
#'  lead to some intervals being slightly longer than \code{interval}. Note that
#'  this is not a replacement for cleaning the \code{animal} data.frame and
#'  removing repetitive timesteps as described in Details. Defaults to
#'  \code{FALSE}.
#'@param crs Optional character string indicating the coordinate reference
#'  system used for all spatial data in the analysis. If not specified, this
#'  will be assumed to match the coordinate reference system of \code{animal},
#'  with a warning. Note that this is a convenience variable used in defining
#'  intermediate products in the analysis and is not a replacement for data
#'  cleaning, preparation, and standardizing projections before the analysis.
#'@param export_images Logical indicator of whether to export images of burst
#'  events as \code{.png} files named with the classification of the event, the
#'  \code{Animal.ID} of the animal, the \code{burstID} of the event, and
#'  potential prefix and suffix values specified by the \code{img_prefix} and
#'  \code{img_suffix} arguments. Defaults to \code{FALSE}.
#'@param img_path Character string indicating the name of the folder to which
#'  images should be exported when \code{export_images == TRUE}. The folder will
#'  be created in the working directory if it does not exist. Defaults to
#'  \code{"event_imgs"}.
#'@param img_prefix Optional character string to be added at the beginning of
#'  image file names written when \code{export_images == TRUE}.
#'@param img_suffix Optional character string to be added at the end of image
#'  file names written when \code{export_images == TRUE}.
#'@param img_background Optional \code{list} object containing one or more
#'  \code{sf} spatial objects to be included in the background of plots when
#'  \code{export_images == TRUE}. For example, this could include spatial
#'  objects representing coastline, rivers, cities, or other features of
#'  interest.
#'
#'@details Terminology, behavioral classifications, and analysis approach in
#'  \code{BaBA_caribou} differ from that used in \code{BaBA_default} and by Xu
#'  et al. (2021).
#'
#'  \cr
#'  \cr
#'  \strong{Terminology}
#'  \itemize{
#'    \item{\emph{Encounter}}{
#'      - One interaction of an animal with one or more barrier(s), consisting of the time from which an animal enters within a barrier buffer (\code{d}) until it leaves that buffer, with brief steps outside the buffer possible as long as they are within a pre-specified \code{tolerance}.
#'    }
#'    \item{\emph{Burst}}{
#'      - An \emph{encounter} may be split into multiple \emph{bursts}, with bursts distinguished if a barrier is crossed or the animal moves from proximity of one barrier to another while remaining within the barrier buffer. This allows identification of multiple behavioral responses to a single barrier, or to multiple nearby barriers with overlapping buffers.
#'    }
#'  }
#'
#'  \cr
#'  \cr
#'  \strong{Behavioral classifications}
#'  
#'  Potential values include:
#'  \itemize{
#'    \item{Normal movement}
#'    \item{Quick cross}
#'    \item{Back-and-forth}
#'    \item{Bounce}
#'    \item{Trace}
#'  }
#'
#'  "Unaltered movements" are those in which movement is indistinguishable from
#'  that typical for a given season in the absence of potential barriers, or
#'  movement in which crossing does not appear hindered. These consist of
#'  “normal movement” and “quick cross” behaviors. Note that in Xu et al. (2021)
#'  and \code{\link{BaBA_default}} “normal movement” is referred to as “average
#'  movement”. This has been changed here to avoid confusion between a
#'  behavioral classification and the model approach of comparing of movement
#'  statistics from each burst with the season-specific average and standard
#'  deviation of movements outside of barrier buffers. Also note that one could
#'  argue that “quick cross” behavior is altered movement because the animal may
#'  increase its movement speed and thus expend additional energy. Nevertheless,
#'  we followed Xu et al. (2021) in labeling this as unaltered due to the
#'  barrier not conspicuously reducing an animal’s mobility.
#'
#'  "Altered movements" are behavior states in which movement behavior appears
#'  to be altered by the presence of a barrier. These consist of
#'  “back-and-forth,” “bounce,” and “trace” behaviors.
#'
#'  \cr
#'  \cr
#'  \strong{Analysis approach}
#'  
#'  Extensive changes to the analysis approach were made in \code{BaBA_caribou}
#'  compared to \code{BaBA_default} to better represent behavioral responses for
#'  caribou of the Western Arctic Herd. These changes may not be appropriate for
#'  all study systems and species but certain aspects may be highly applicable
#'  in other contexts. For details of changes see comments in the code below.
#'  Major changes include:
#'  \itemize{
#'    \item{Always excluding locations within buffers when calculating average movement and the standard deviation (sd) of movement.}
#'    \item{Using season-specific average/sd movement calculations, rather than a user-specific moving window.}
#'    \item{Dropped user-specified time thresholds for short-duration and long-duration events.}
#'    \item{Dropped the "trapped" movement behavior, instead relying on duration and behavior to indicating when such an event occurred.}
#'    \item{Added analysis of turning angles to identify behaviors such as "bounce" based not just on their duration, but also changes in direction and consistency of angles in a burst.}
#'    \item{Refined how "trace" behavior is reflected to look at the direction and persistence of movement relative to the angle of the nearest barrier.}
#'    \item{Distinguished between true and false crossings and tested for this in the dataset.}
#'    \item{Divided encounters into bursts based on barrier crossing and proximity to nearest barrier.}
#'  }
#'
#'@return A \code{list} consisting of three objects:
#'  \itemize{
#'    \item{\code{$encounters}} {
#'      An \code{sf POINT} object with one record per burst in the dataset. These indicate the date and spatial location of the first point in each burst, along with burstID and classification information. Retained for historical reasons as a holdover from \code{\link{BaBA_default}}. The informatino here is redundant with that in \code{encounter_is}.
#'    }
#'    \item{\code{$encounter_is}} {
#'      A named \code{list} object with a length equal to the number of bursts in the dataset. Names correspond to the burstID for each entry. Each entry consists of a \code{sf POINT} object with all locations for each burst in the dataset, along with corresponding information about season, barrier proximity, and crossing. Useful for visualizing specific bursts after an analysis is completed.
#'    }
#'    \item{\code{$classification}} {
#'      A \code{tibble} containing the BaBA results for each burst analyzed in the dataset. Number of rows corresponds to the number of bursts. Columns include:
#'      \itemize{
#'        \item{\code{AnimalID}} {Unique character indicator of each collared animal in the dataset. There will likely be multiple records for the same animal in the dataset.}
#'        \item{\code{encounter}} {Unique character indicator of each encounter in the dataset, defined as one interaction of an animal with one or more barrier(s). This consists of the time from which an animal enters within a barrier buffer (\code{d}) until it leaves that buffer, with brief steps outside the buffer possible as long as they are within a pre-specified \code{tolerance}. There will likely be multiple records for the same encounter in the same dataset. Encounters are specified with the \code{AnimalID} and a uniuqe numeric indicator representing the specific encounter, separated by an underscore (e.g., '0903_1833').}
#'        \item{\code{burstID}} {Unique character indicator of the specific burst evaluated. An \emph{encounter} may be split into multiple \emph{bursts}, with bursts distinguished if a barrier is crossed or the animal moves from proximity of one barrier to another while remaining within the barrier buffer (\code{d}). This allows identification of multiple behavioral responses to a single barrier, or to multiple nearby barriers. Bursts are specified with the \code{encounter} label and a uniuqe integer indicator representing the specific burst within the encounter, separated by an underscore (e.g., '0903_1833_2').}
#'        \item{\code{season}} {Character indicator of the predominant season in which each burst takes place. Based on the dates of the locations in the burst and the season boundaries identified for the WAH by Joly and Cameron (2023).}
#'        \item{\code{barrier}} {Character indicator of the barrier(s) with which the animal comes into proximity (i.e., within the barrier buffer distance, \code{d}) during a given burst.}
#'        \item{\code{barrier_n}} {Character indicator of the number of locations in closest proximity to each barrier in the burst. For example, if a caribou had three locations closer to the Kivalina road and seven closer to the Red Dog Road in a particular burst, the value for that burst would be "Kivalina-3-RedDog-7".}
#'        \item{\code{barrier_min_dist}} {Numeric value indicating the minimum distance, in km, at which a point was observed near a barrier during the given burst. Note that because observations are periodic, not continuous, the animal may have come closer to the barrier, or even crossed the barrier, but this distance is recorded based on the nearest observed point.}
#'        \item{\code{closest_bar}} {Character indicator of the single barrier that is closest to most points within the given burst.}
#'        \item{\code{closest_dist}} {Numeric value indicating the closest distance (in km) of any point in the burst to the closest barrier. As with \code{barrier_min_dist} the animal may have come closer to the barrier, or even crossed the barrier, at an unobserved time, but this distance reflects the nearest observed location.}
#'        \item{\code{start_time}} {Datetime object indicating the date and time in \code{POSIXct} format of the first location for the given burst.}
#'        \item{\code{end_time}} {Datetime object indicating the date and time in \code{POSIXct} format of the last location for the given burst.}
#'        \item{\code{duration}} {Duration of the given burst, in units specified by the \code{units} argument of \code{BaBA}.}
#'        \item{\code{cross_any}} {Integer indicator of the number of barrier crossings apparent in the given burst. Note, these are crossings based on straight lines connecting subsequent locations and may not be representative of real crossings.}
#'        \item{\code{cross_true}} {Integer indicator of whether a crossing actually is expected to have occurred (1) or not (0). If subsequent locations are on the same side of the barrier, the animal is assumed not to have crossed the barrier.}
#'        \item{\code{cross_bar}} {Character indicator of which barrier was crossed. Takes a value of \code{NA} unless \code{cross_true == 1}.}
#'        \item{\code{cross_x}} {Numeric indicator of the x-coordinate of the location at which the straight line between subsequent points crosses the barrier indicated in \code{cross_bar}. Note that crossing may have occurred at a different location, as the crossing location is nearly always unobserved, but this is the best available indication of where the barrier was crossed. Takes a value of \code{NA} unless \code{cross_true == 1}.}
#'        \item{\code{cross_y}} {Numeric indicator of the y-coordinate of the location at which the straight line between subsequent points crosses the barrier indicated in \code{cross_bar}. Note that crossing may have occurred at a different location, as the crossing location is nearly always unobserved, but this is the best available indication of where the barrier was crossed. Takes a value of \code{NA} unless \code{cross_true == 1}.}
#'        \item{\code{class}} {Character indicator of the movement classification assigned to the given burst. Primary classes include "Average_Movement", "Quick_Cross", "Back_n_forth", "Bounce", "Trace", and "Unknown". Additional temporary indicators include "Avg_lclNA", "Unknown_insufficient_n", "Average_sd0", and "Unknown_cross_bounce". These tier to the previous classifications but provide additional information about why that categorization was reached, for diagnostic purposes.}
#'        \item{\code{str_i}} {Numeric value indicating the straightness of movement in the given burst, as calculated by the \code{\link{strtns}} function. Values range between 0-1, with with values closer to 0 indicating more sinuous movement.}
#'        \item{\code{str_mean}} {Numeric value indicating the average straightness of movement for other movements of the same duration as the given burst, that occur outside of barrier buffers in the same season for the given individual across all years of observation. Values range between 0-1, with with values closer to 0 indicating more sinuous movement.}
#'        \item{\code{str_sd}} {Numeric value indicating the standard deviation of straightness of movement values for other movements of the same duration as the given burst, that occur outside of barrier buffers in the same season for the given individual across all years of observation.}
#'        \item{\code{ang_i}} {Numeric value indicating the average encounter angle, as a general sense of the predominant heading of movement during the encounter. Values are given in units of degrees, ranging between 0-360.}
#'        \item{\code{ang_mean}} {Numeric value indicating the average barrier angle, as a general sense of the predominant heading of the nearest barrier in space. Values are given in units of degrees, ranging between 0-360.}
#'        \item{\code{ang_sd}} {Numeric value indicating the standard deviation of barrier angles along the nearest barrier, giving a general indication of how concentrated the direction of the barrier is in space. Values are given in units of degrees, ranging between 0-360.}
#'        \item{\code{mrl_1}} {Numeric value indicating the mean resultant length of angles between the burst points occurring \emph{before} the closest point to the given barrier. The mean resultant length indicates the concentration of data points around a circle. Used for determining whether a change in direction occurred before/after the nearest point to a barrier to help distinguish between \emph{bounce} and \emph{back-and-forth} movement.}
#'        \item{\code{mrl_2}} {Numeric value indicating the mean resultant length of angles between the burst points occurring \emph{after} the closest point to the given barrier. The mean resultant length indicates the concentration of data points around a circle. Used for determining whether a change in direction occurred before/after the nearest point to a barrier to help distinguish between \emph{bounce} and \emph{back-and-forth} movement.}
#'        \item{\code{mrl_mean}} {Numeric value indicating the average value of \code{mrl_1} and \code{mrl_2}. Used in distinguishing \emph{bounce} and \emph{back-and-forth} movements, based on preliminary testing.}
#'        \item{\code{lcl_prop}} {Numeric value indicating the proportion of movement steps in the given burst that are parallel to the closest barrier. Used to distinguish \emph{trace} behavior, based on preliminary testing.}
#'        \item{\code{easting}} {Numeric value indicating the x-coordinate of the first location of a given burst. Retained for legacy purposes as this occurs in \code{BaBA_default}.}
#'        \item{\code{northing}} {Numeric value indicating the y-coordinate of the first location of a given burst. Retained for legacy purposes as this occurs in \code{BaBA_default}.}
#'      }  
#'    }
#'  }
#'
#'@export
#'
#'@references Joly K, Cameron MD. 2023. Caribou vital sign annual report for the
#'  Arctic Network Inventory and Monitoring Program: September 2022–August 2023.
#'  Natural Resource Report NPS/ARCN/NRR—2023/2612. National Park Service, Fort
#'  Collins, Colorado.
#'
#'  Xu W, Dejid N, Herrmann V, Sawyer H, Middleton AD. 2021. Barrier Behaviour
#'  Analysis (BaBA) reveals extensive effects of fencing on wide-ranging
#'  ungulates. Journal of Applied Ecology 58: 690-698.
#'
#'@examples
#'\dontrun{
#'wah.out <- BaBA_caribou(animal = wah_all, barrier = rd_final,
#'                        d = c(20000, 20000, 5000, 20000, 20000),
#'                        interval = 8, tolerance = 0, units = 'hours',
#'                        round_fixes = TRUE, crs = 'EPSG:6393',
#'                        export_images = TRUE,
#'                        img_path = 'C:/BaBA_runs'),
#'                        img_prefix = 'wah_all', img_suffix = Sys.Date(),
#'                        img_background = list(ak_border))
#'}
BaBA_caribou <-
  function(animal, barrier, d, interval = NULL, tolerance = 0, units = "hours",
           sd_multiplier = 1, round_fixes = FALSE, crs = NULL, 
           export_images = FALSE, img_path = "event_imgs", img_prefix = NULL,
           img_suffix = NULL, img_background = NULL) {
    
    if(export_images) {
      if(!dir.exists(img_path)) dir.create(img_path)
    }
    
    ## Prepare parameters and check input
    if (class(animal)[1] != "sf") stop("animal needs to be an sf object")
    if (sf::st_geometry_type(animal)[1] != 'POINT') stop("animal needs to have a POINT geometry")
    if (class(barrier)[1] != "sf") stop("barrier needs to be an sf object")
    if (!(sf::st_geometry_type(barrier)[1] %in% c('MULTILINESTRING', 'LINESTRING'))) stop("barrier needs to have either a LINESTRING or MULTILINESTRING geometry")
    if (!"date" %in% names(animal)) stop("Please rename the date column to 'date'")
    if (!"Animal.ID" %in% names(animal)) stop("Please rename the individual ID column to 'Animal.ID'")
    if (!(inherits(animal$date, "POSIXct"))) stop("Date needs to be 'POSIXct' format")
    if (sum(is.na(animal$date)) > 0) stop("Please exclude rows where date is NA")
    if(is.null(crs)) warning("No crs specified. assuming crs for animal applies to all spatial data") 
    if(is.null(crs)) crs <- sf::st_crs(animal)
    if(sf::st_crs(animal) != sf::st_crs(crs)) stop("Coordinate reference system of animal must match crs")
    if(sf::st_crs(barrier) != sf::st_crs(crs)) stop("Coordinate reference system of barrier must match crs")
    
    ## Round fixes, if desired. NOTE: This is only intended to account for minor
    ## variation in fix acquisition time, not to clean irregular data. Please
    ## clean data to remove bursts of different fix intervals prior to
    ## conducting this analysis.
    if(round_fixes){
      animal$date <- lubridate::round_date(animal$date, unit = units)
      interval_per_individual <- tapply(animal$date, animal$Animal.ID, function(x) names(which.max(table(round(as.numeric(diff(x), units = units),0)))))
    } else {
      interval_per_individual <- tapply(animal$date, animal$Animal.ID, function(x) names(which.max(table(as.numeric(diff(x), units = units)))))
    }
    if(is.null(interval)) { ## figure out interval (as the most frequent difference in timestamp) if not provided but give an error if not the same for all individuals
      if(all(interval_per_individual == interval_per_individual[1])) interval <- as.numeric(interval_per_individual[1]) else stop("Not all individuals have been sampled at the same frequency. Run individuals with different intervals seperately, or double-check whether your date column is cleaned.")
    } else {
      if (any(as.numeric(interval_per_individual) > interval, na.rm = TRUE)) stop("BaBA interval needs to be no smaller than the actual data interval. Also double-check whether your date column is cleaned.") 
    }
    
    
    # Classification step 1: generate encounter event data.frame --------------

    ## Create point ID by individual and add season indicators
    animal <-
      animal %>% 
      dplyr::arrange(Animal.ID, date) %>% 
      dplyr::group_by(Animal.ID) %>% 
      dplyr::mutate(ptsID = 1:dplyr::n()) %>% 
      dplyr::ungroup() %>% 
      ## Add season indicator here to allow separation of points by season for
      ## seasonal buffer exclusion below. Season breaks derived from Joly and
      ## Cameron (2023):
      ## Spring migration: Apr 1 - May 27
      ## Calving: May 28 - Jun 14
      ## Insect relief: Jun 15 - Jul 14
      ## Late summer: Jul 15 - Aug 31
      ## Fall migration: Sep 1 - Nov 30
      ## Winter: Dec 1 - Mar 31
      mutate(season = dplyr::case_when(date >= paste0(lubridate::year(date), '-04-01') & date < paste0(lubridate::year(date), '-05-28') ~ 'spring_mig',
                                date >= paste0(lubridate::year(date), '-05-28') & date < paste0(lubridate::year(date), '-06-15') ~ 'calving',
                                date >= paste0(lubridate::year(date), '-06-15') & date < paste0(lubridate::year(date), '-07-15') ~ 'insect',
                                date >= paste0(lubridate::year(date), '-07-15') & date < paste0(lubridate::year(date), '-09-01') ~ 'late_summer',
                                date >= paste0(lubridate::year(date), '-09-01') & date < paste0(lubridate::year(date), '-12-01') ~ 'fall_mig',
                                TRUE ~ 'winter'))
    
    ## Calculate the distance to nearest barrier for each point in the animal
    ## dataset, as well as the identity of that barrier. This assumes that the
    ## barrier file has a "Name" column that uniquely identifies each barrier
    ## segment. Add the nearest barrier and minimum distance to barrier to the
    ## animal object.
    bar_dist <- sf::st_distance(x = animal, y = barrier)
    animal$bar_min <- barrier$Name[apply(bar_dist, 1, which.min)]
    animal$bar_dist_km <- apply(bar_dist, 1, function(x) x[which.min(x)])/1000
    
    ## Explicitly suppress constant geometry assumption warning by confirming
    ## attribute is constant throughout the geometry. See
    ## https://github.com/r-spatial/sf/issues/406 for details.
    sf::st_agr(animal) <- 'constant'   
    sf::st_agr(barrier) <- 'constant'
    
    ## Create a point version of the barrier lines
    barrier_pts_all <-
      barrier %>% 
      ## Convert the line to a point geometry
      sf::st_cast(to = 'MULTIPOINT', warn = FALSE) %>% ## Do this first to retain all pieces. Turning off warnings since we are not worried about attributes, just locations.
      sf::st_cast(to = 'POINT', warn = FALSE) %>%   ## Do this so that each point gets its own distance. Turning off warnings since we are not worried about attributes, just locations.
      ## Add helpful information used in later steps
      dplyr::mutate(
        ## Assign a unique id (barID) to each point
        barID = 1:nrow(.),
        ## Add in the xy coords as columns, as expected by calc_angle()
        x = sf::st_coordinates(.)[,1],
        y = sf::st_coordinates(.)[,2])
    
    ## Iterate through each barrier, creating a polygon on each side of the
    ## barrier that will be used below to identify true crossings, storing the
    ## results in a pre-made list with length and names corresponding to the
    ## unique barriers
    bar_list <- vector('list', length(unique(barrier_pts_all$Name)))
    names(bar_list) <- unique(barrier_pts_all$Name)
    for(z in 1:length(unique(barrier_pts_all$Name))){
      ## Identify the barrier of interest
      bar_tmp <- 
        barrier_pts_all %>% 
        dplyr::filter(Name == unique(barrier_pts_all$Name)[z])
      
      ## Calculate the mean direction of the barrier
      bar_ang_mean <-
        bar_tmp %>% 
        ## Calculate angles
        calc_angle() %>% 
        circular::circular(zero = 0) %>% 
        ## Calculate the mean value
        circular::mean.circular(na.rm = TRUE) %>% 
        ## Convert units to degrees in the 360 deg range
        circular::conversion.circular(units = 'degrees') %% 360
      
      ## Identify the distance that should be extended around the barrier. For
      ## most barriers this should be 2d, using the barrier-specific d buffer
      ## distance. However, the Dalton Highway is so large that this needs to be
      ## increased. Preliminary testing found that 10d is suitable.
      d_tmp <- 
        ifelse(length(d) > 1,
               d[which(barrier$Name == bar_tmp$Name[1])],
               d)
      d_target <- ifelse(bar_tmp$Name[1] == 'Dalton', 10*d, 2*d)
      
      ## First extend the barrier on either end in the predominant barrier
      ## direction
      pt0a <-
        line_extend(pt = bar_tmp[1, c('x','y')] %>% 
                      sf::st_drop_geometry(),
                    angle = bar_ang_mean + 180 %% 360,
                    len = d_target) %>% 
        sf::st_as_sf(coords = c('x', 'y'), crs = crs, agr = 'constant') %>% 
        dplyr::slice(2)
      pt0b <-
        line_extend(pt = bar_tmp[nrow(bar_tmp), c('x','y')] %>%
                      sf::st_drop_geometry(),
                    angle = bar_ang_mean,
                    len = d_target) %>% 
        sf::st_as_sf(coords = c('x', 'y'), crs = crs, agr = 'constant') %>% 
        dplyr::slice(2)
      
      ## Then extend these points 90 degrees in either direction to get the
      ## outer polygon edges
      pt1a <-  
        line_extend(pt = pt0a %>%
                      sf::st_coordinates(),
                    angle = bar_ang_mean + 90 %% 360,
                    len = d_target) %>% 
        as.data.frame() %>% 
        sf::st_as_sf(coords = c('X', 'Y'), crs = crs, agr = 'constant') %>% 
        dplyr::slice(2)
      pt2a <-  
        line_extend(pt = pt0a %>%
                      sf::st_coordinates(),
                    angle = bar_ang_mean - 90 %% 360,
                    len = d_target) %>% 
        as.data.frame() %>% 
        sf::st_as_sf(coords = c('X', 'Y'), crs = crs, agr = 'constant') %>% 
        dplyr::slice(2)
      pt1b <-  
        line_extend(pt = pt0b %>%
                      sf::st_coordinates(),
                    angle = bar_ang_mean + 90 %% 360,
                    len = d_target) %>% 
        as.data.frame() %>% 
        sf::st_as_sf(coords = c('X', 'Y'), crs = crs, agr = 'constant') %>% 
        dplyr::slice(2)
      pt2b <-  
        line_extend(pt = pt0b %>%
                      sf::st_coordinates(),
                    angle = bar_ang_mean - 90 %% 360,
                    len = d_target) %>% 
        as.data.frame() %>% 
        sf::st_as_sf(coords = c('X', 'Y'), crs = crs, agr = 'constant') %>% 
        dplyr::slice(2)
      
      ## Create a polygon on each side of the barrier using these points and the
      ## barrier points
      poly1 <-
        dplyr::bind_rows(pt0a,
                         bar_tmp %>% 
                           dplyr::select(geometry),
                         pt0b,
                         pt1b,
                         pt1a,
                         pt0a) %>% 
        dplyr::summarize(geometry = sf::st_combine(geometry)) %>% 
        sf::st_cast(to = 'POLYGON')
      poly2 <-
        dplyr::bind_rows(pt0a,
                         bar_tmp %>% 
                           dplyr::select(geometry),
                         pt0b,
                         pt2b,
                         pt2a,
                         pt0a) %>% 
        dplyr::summarize(geometry = sf::st_combine(geometry)) %>% 
        sf::st_cast(to = 'POLYGON')
      
      ## Combine these into a single two-feature polygon object and store it in
      ## the appropriate place in the barrier list
      bar_list[z] <- rbind(poly1, poly2)
    }
    
    ## Create buffer around barrier for identifying encounters
    print("locating encounter events...")
    barrier_buffer <- 
      barrier %>% 
      sf::st_buffer(dist = d, nQuadSegs = 5) %>%   ## Note that nQuadSegs is set to 5 as this was the default value for rgeos::gBuffer in previous versions of BaBA
      sf::st_union()
    
    ## Extract points that fall inside the buffer
    encounter <- sf::st_intersection(animal, barrier_buffer)
    
    if (nrow(encounter) == 0) stop("no barrier encounter detected.")
    
    ## Create an object to hold crossing coordinates
    cross_coords <- NULL
    
    ## Create unique burstIDs
    for(i in unique(encounter$Animal.ID)){
      if (nrow(encounter %>% dplyr::filter(Animal.ID == i)) == 0) {
        warning(paste0 ("Individual ", i, " has no locations overlapped with the barrier buffer and is eliminated from analysis." ))
        next()
      }
      
      ## Prep encounter data for the current animal
      encounter_i <-
        encounter %>% 
        dplyr::filter(Animal.ID == i) %>% 
        ## Add indicator of point difference. Subtract one from each so that
        ## subsequent ptsIDs get a value of zero, which will be useful for
        ## calculating the burstIDs below.
        dplyr::mutate(ptdiff = ptsID - dplyr::lag(ptsID) - 1)
      ## Set the first ptdiff value to zero
      encounter_i$ptdiff[1] <- 0
      
      ## Include points where the animal stepped outside the buffer during an
      ## encounter but within the tolerance threshold. These points will still
      ## be counted as within the same encounter. Note this is a different use
      ## of tolerance from the original BaBA() code. Instead of the amount of
      ## time in which the animal can be outside the buffer and still be part of
      ## the same encounter, it is now the number of steps (i.e., ptIDs) that an
      ## animal can be out of the buffer and still count. This is conceptually
      ## the same, but with different units.
      if(any(encounter_i$ptdiff > 0 & encounter_i$ptdiff <= tolerance, na.rm = TRUE)){
        idx_pts_of_interest <- which(encounter_i$ptdiff > 0 & encounter_i$ptdiff <= tolerance)
        for(pt in idx_pts_of_interest) {
          ## Identify pts to fetch
          ptsID_of_interest_B <- encounter_i$ptsID[pt]
          ptsID_of_interest_A <- encounter_i$ptsID[pt-1]
          
          ## Fetch points outside of the buffer within tolerance
          fetched_pt <- 
            animal %>% 
            dplyr::filter(Animal.ID == i & 
                            ptsID > ptsID_of_interest_A & 
                            ptsID < ptsID_of_interest_B)
          
          ## If there are no points outside of the buffer that means there is
          ## missing data. Since the missing data are still within the
          ## tolerance, we consider ptdiff = 0 so the points before and after
          ## will be in the same event.
          if (nrow(fetched_pt) == 0) {  
            encounter_i$ptdiff[pt] <- 0
            next() } 
          else {
            fetched_pt$ptdiff <- 0 
            ## Reset ptdiff  to 0
            encounter_i$ptdiff[pt] <- 0 
            ## Append fetched points 
            if(pt == idx_pts_of_interest[1]) {fetched_pts <- fetched_pt} else if (exists("fetched_pts")) { fetched_pts <- rbind(fetched_pts, fetched_pt) } else {fetched_pts <- fetched_pt}
          }
        }
        
        ## Add fetched pts to the encounter
        encounter_i <- rbind(encounter_i, fetched_pts)
        ## Reorder the encounter data
        encounter_i <- encounter_i[order(encounter_i$ptsID), ]
      }
      
      ## Identify unique burstIDs using the cumulative sum of the ptdiff column to identify unique burst IDs (with animalID) 
      encounter_i$burstID <- paste(i, cumsum(encounter_i$ptdiff), sep = "_")
      
      ## Check for barrier crossing
      
      ## Convert the encounter points with more than one record for a given
      ## burstID into a movement line
      mov_seg_i <-
        encounter_i %>% 
        dplyr::add_count(burstID) %>% 
        dplyr:: filter(n > 1) %>% 
        dplyr::group_by(burstID) %>%
        dplyr::summarize(do_union = FALSE) %>%
        sf::st_cast(to = 'LINESTRING')
      ## Identify the intersection coordinates
      sf::st_agr(mov_seg_i) <- 'constant'   ## Suppress warning
      int_pts <-
        mov_seg_i %>% 
        sf::st_intersection(barrier) %>% 
        sf::st_cast(to = 'MULTIPOINT')
      if(nrow(int_pts) > 0){
        int_coords <-
          int_pts %>% 
          sf::st_coordinates() %>% 
          dplyr::as_tibble() %>% 
          dplyr::mutate(burstID = int_pts$burstID[.$L1]) %>% 
          dplyr::select(-L1)  ## For cleanliness remove L1
        ## Add the coordinates to the output
        cross_coords <- rbind(cross_coords, int_coords)
      }
      
      ## Add into encounter_complete
      if(i == unique(encounter$Animal.ID[1])) encounter_complete <- encounter_i else encounter_complete <- rbind(encounter_complete, encounter_i)
    }
    
    ## Add indicators to encounter_complete of whether each location represents
    ## the endpoint of a segment where a crossing was indicated and where a true
    ## crossing occurred. By default this will indicate no (0) for all and then
    ## update those for which a crossing is indicated.
    encounter_complete$cross_ind <- 0
    encounter_complete$cross_true <- 0
    
    ## Identify the encounters with multiple nearest barriers and assign new
    ## columns that can be used to indicate unique sequences of nearest barrier
    ## locations.
    bursts_updated_barrier <-
      encounter_complete %>% 
      dplyr::left_join(
        encounter_complete %>% 
          dplyr::group_by(burstID) %>% 
          dplyr::summarize(nbar = length(unique(bar_min))) %>% 
          sf::st_drop_geometry(),
        by = 'burstID') %>% 
      dplyr::filter(nbar > 1) %>% 
      dplyr::group_by(burstID) %>% 
      dplyr::mutate(
        bar_change = ifelse(dplyr::lag(bar_min) != bar_min, 1, 0),
        bar_change = ifelse(is.na(bar_change), 0, bar_change),
        bar_dif = cumsum(bar_change)) %>% 
      dplyr::ungroup()
    
    ## Run through those encounters for which crossings were indicated and see
    ## if they actually occurred
    bursts_updated_cross <- NULL
    encounter_complete$cross_bar <- NA
    encounter_complete$cross_x <- NA
    encounter_complete$cross_y <- NA
    for(k in unique(cross_coords$burstID)){
      ## Pull data for the individual of interest
      encounter_i <- 
        encounter_complete %>% 
        dplyr::filter(burstID == k)
      
      ## Iterate through each pair of encounter points to check for indicated
      ## and true crossings
      for(j in 1:(nrow(encounter_i)-1)){
        pts_tmp <- encounter_i[j:(j+1),]
        line_tmp <- pts_tmp %>% 
          dplyr::summarize(do_union = FALSE) %>%
          sf::st_cast(to = 'LINESTRING')
        
        ## Check for barrier intersection
        int_check <- sf::st_intersection(line_tmp, barrier)
        
        ## If there is an intersection, check for points on different sides of
        ## the barrier to indicate a true crossing
        if(nrow(int_check) > 0){
          ## Indicate that a crossing was indicated. Using j + 1 sets this to be
          ## the endpoint of the crossing (i.e., after the purported crossing
          ## has occurred), which is important when using the crossing
          ## information to split encounters into sub-bursts below.
          encounter_i$cross_ind[j+1] <- 1
          
          ## Pull the road side polygons for the intersected barrier
          poly.tmp <-
            bar_list %>% 
            purrr::pluck(int_check$Name)
          
          ## Crop the points with each polygon, setting st_arg() to avoid a
          ## warning
          sf::st_agr(pts_tmp) <- 'constant'
          enc_poly1 <- sf::st_intersection(pts_tmp, poly.tmp[1])
          enc_poly2 <- sf::st_intersection(pts_tmp, poly.tmp[2])
          
          ## If either sample size is 0, indicate a false crossing, otherwise
          ## indicate a true crossing
          encounter_i$cross_true[j+1] <- ifelse(nrow(enc_poly1) == 0 | nrow(enc_poly2) == 0, 0, 1)
          
          ## If a true crossing occurred, record the barrier that was crossed
          ## and location. If multiple barriers were crossed, combine names, if
          ## multiple locations along a barrier were crossed in the single step,
          ## return the location of the first.
          if(encounter_i$cross_true[j+1] == 1){
            encounter_i$cross_bar[j+1] <- paste(int_check$Name, collapse = '_')
            encounter_i$cross_x[j+1] <- dplyr::first(sf::st_coordinates(int_check)[,1])
            encounter_i$cross_y[j+1] <- dplyr::first(sf::st_coordinates(int_check)[,2])
          }
        }
      }
      
      ## Add column to encounter_i that indicates parts before/after crossing,
      ## which I can use to identify different sub-bursts based on crossing events.
      encounter_i$cumcross <- cumsum(encounter_i$cross_true)
      
      ## Output the results
      bursts_updated_cross <- rbind(bursts_updated_cross, encounter_i)
    }
    
    ## Incorporate the extra columns from the barrier and crossing objects into
    ## encounter_complete and use them to make a joint burstID column
    encounter_complete <-
      encounter_complete %>% 
      ## Make a unique column that is Animal.ID_ptsID and use that to merge the
      ## data
      dplyr::mutate(join_col = paste(Animal.ID, ptsID, sep = '_')) %>% 
      dplyr::left_join(
        bursts_updated_barrier %>% 
          dplyr::mutate(join_col = paste(Animal.ID, ptsID, sep = '_')) %>% 
          sf::st_drop_geometry() %>% 
          dplyr::select(join_col, bar_dif),
        by = 'join_col') %>% 
      dplyr::left_join(
        bursts_updated_cross %>% 
          dplyr::mutate(join_col = paste(Animal.ID, ptsID, sep = '_')) %>% 
          sf::st_drop_geometry() %>% 
          dplyr::select(join_col, cross_ind, cross_true, cross_bar, cross_x, cross_y, cumcross),
        by = 'join_col') %>% 
      ## Create updated columns with the desired information
      dplyr::mutate(
        ## Deal with NAs in bar_dif and cumcross, to allow summing
        bar_dif2 = ifelse(is.na(bar_dif), 0, bar_dif),
        cumcross2 = ifelse(is.na(cumcross), 0, cumcross),
        ## Sum barrier differences and cumulative crossing to get unique indicators
        ## (per burstID) of barrier/crossing events
        barcross = bar_dif2 + cumcross2) %>% 
      ## Add in the by-group sum of cumcross to each record
      dplyr::group_by(burstID) %>% 
      dplyr::mutate(sum_cumcross = sum(cumcross, na.rm = TRUE)) %>% 
      dplyr::ungroup() %>% 
      ## Add additional needed information
      dplyr::mutate(
        ## Create a new burstID as the combination of burstID and barcross, if
        ## those exist
        burstID2 = ifelse(is.na(bar_dif) & (is.na(cumcross) | sum_cumcross == 0),
                          burstID,
                          paste(burstID, barcross, sep = '_')),
        ## Create updated columns for indicated and true crossings
        cross_ind = ifelse(is.na(cross_ind.y), cross_ind.x, cross_ind.y),
        cross_true = ifelse(is.na(cross_true.y), cross_true.x, cross_true.y)) %>% 
      ## Retain only the needed columns
      dplyr::select(Animal.ID:geometry, burstID = burstID2, cross_ind, cross_true,
                    cross_bar = cross_bar.y, cross_x = cross_x.y, cross_y = cross_y.y)
    
    
    ## Rename as encounter (encounter_complete may be bigger as it includes extra
    ## points that are within tolerance)
    encounter <- encounter_complete
    
    
    
    # Classification step 2: classify events ----------------------------------

    print("classifying behaviors...") 
    ## Open progress bar
    pb <- utils::txtProgressBar(style = 3)
    
    ## Create empty object that will hold results
    event_df <- NULL
    ## Keep the encounter_i objects for use in plotting specific encounters
    encounter_i_out <- vector('list', length(unique(encounter$burstID)))
    names(encounter_i_out) <- unique(encounter$burstID)
    
    ## Run classification procedure for each encounter
    for(i in unique(encounter$burstID)) {
      ## Update progressbar
      utils::setTxtProgressBar(pb, which(unique(encounter$burstID) == i)/length(unique(encounter$burstID)))
      
      ## Subset down to the specific encounter and animal data
      encounter_i <- encounter[encounter$burstID == i, ]
      animal_i <- animal[animal$Animal.ID == encounter_i$Animal.ID[1],]
      
      ## Expand the encounter by one point on either side to include the step
      ## that brings the animal within the barrier buffer. This will also
      ## include the step that takes the animal from one nearest barrier to
      ## another within the same buffer, if applicable.
      encounter_i <-
        animal_i %>% 
        dplyr::filter(ptsID == encounter_i$ptsID[1] - 1 |
                        ptsID == encounter_i$ptsID[nrow(encounter_i)] + 1) %>% 
        ## Reorder to match encounter_i
        dplyr::select(Animal.ID:y, ptsID, season, bar_min, bar_dist_km, geometry) %>%
        ## Add the burstID
        dplyr::mutate(burstID = encounter_i$burstID[1]) %>% 
        ## Combine with existing encounter_i data
        dplyr::bind_rows(encounter_i) %>% 
        ## Rearrange to be in ptsID order
        dplyr::arrange(ptsID)
      
      ## Calculate duration and straightness
      duration <-  difftime (encounter_i$date[nrow(encounter_i)], encounter_i$date[1], units = units)
      if(round_fixes) duration <- round(duration)
      straightness_i <- strtns(encounter_i)
      
      ## Calculate turning angles of the encounter and add them to encounter_i
      encounter_i$angle <- 
        calc_angle(x = encounter_i) %>% 
        circular::circular(zero = 0) %>%
        circular::conversion.circular(units = 'degrees') %% 360
      
      ## Check whether the first step of the expanded burst crosses a road to
      ## pick up situations where the animal took a step wider than d leading to
      ## an unidentified crossing. This is only needed if a crossing is not
      ## already indicated for the second location of the burst (i.e., the
      ## endpoint of the first step). It also is only evaluated if the previous
      ## location is within 3 days of the current location, so that excessive
      ## missing data does not lead to an indication of a crossing where one is
      ## unclear.
      if(encounter_i$cross_ind[2] == 0 &
         difftime(encounter_i$date[2], encounter_i$date[1], units = units) < 72){
        pts_tmp <- encounter_i[1:2,]
        line_tmp <- pts_tmp %>% 
          dplyr::summarize(do_union = FALSE) %>%
          sf::st_cast(to = 'LINESTRING')
        sf::st_agr(pts_tmp) <- 'constant'
        int_check <- sf::st_intersection(line_tmp, barrier)
        
        ## If there's an intersection, check for points on different sides of the
        ## barrier to indicate a true crossing
        if(nrow(int_check) > 0){
          ## Indicate that a crossing was indicated. As is done above, indicate this
          ## for the endpoint of the crossing (i.e., after the purported crossing
          ## has occurred).
          encounter_i$cross_ind[2] <- 1
          
          ## Pull the road side polygons for the intersected barrier
          poly.tmp <-
            bar_list %>% 
            purrr::pluck(int_check$Name)
          
          ## Crop the points with each polygon
          enc_poly1 <- sf::st_intersection(pts_tmp, poly.tmp[1])
          enc_poly2 <- sf::st_intersection(pts_tmp, poly.tmp[2])
          
          ## If either sample size is 0, indicate a false crossing, otherwise
          ## indicate a true crossing
          encounter_i$cross_true[2] <- ifelse(nrow(enc_poly1) == 0 | nrow(enc_poly2) == 0, 0, 1)
          
          ## If a true crossing occured, record the barrier that was crossed and location
          if(encounter_i$cross_true[2] == 1){
            encounter_i$cross_bar[2] <- paste(int_check$Name, collapse = '_')
            encounter_i$cross_x[2] <- dplyr::first(sf::st_coordinates(int_check)[,1])
            encounter_i$cross_y[2] <- dplyr::first(sf::st_coordinates(int_check)[,2])
          }
        }
      }
      
      ## For the simplicity of the code, make a cross_true object that
      ## summarizes encounter_i$cross_true
      cross_true <- sum(encounter_i$cross_true, na.rm = TRUE)
      
      ## Identify the closest barrier to most points, for use below
      bar_closest <- names(which.max(table(encounter_i$bar_min)))
      
      ## Identify the season of the encounter for pulling seasonal data below
      ## and to facilitate subsequent analysis by season. If the encounter spans
      ## seasons, isolate the modal season so that we compare avg movement
      ## during the most represented season to that during the encounter.
      season_i <- names(which.max(table(encounter_i$season)))
      
      ## Identify all seasonal points outside of buffers for calculation of
      ## avg movement statistics
      animal_season_exclude <-
        animal_i %>% 
        ## Subset down to the most common season in the encounter
        dplyr::filter(season == season_i) %>% 
        ## Remove all points inside buffers
        dplyr::filter(!(ptsID %in% encounter$ptsID[encounter$Animal.ID == animal_i$Animal.ID[1]])) %>% 
        ## Add indicator of point difference. Subtract one from each so that
        ## subsequent ptsIDs get a value of zero, which will be useful for
        ## calculating the burstIDs below.
        dplyr::mutate(ptdiff = ptsID - dplyr::lag(ptsID) - 1)
      ## Set the first ptdiff value to zero
      animal_season_exclude$ptdiff[1] <- 0
      ## Reset any point differences within tolerance to zero to not split up
      ## continuous segments by missing data if it is within the specified
      ## tolerance level
      animal_season_exclude$ptdiff[animal_season_exclude$ptdiff <= tolerance] <- 0
      ## Identify continuous sections in the remaining movement data
      animal_season_exclude$continuousID <- cumsum(animal_season_exclude$ptdiff)
      
      ## For each continuous section, calculate straightness of all movements
      ## lasting the duration of the encounter (moving window of the size of the
      ## encounter) or 1000 hrs, whichever is shorter.
      duration_comp <- min(as.numeric(duration), 1000)
      straightnesses_seasonal <- NULL
      for(ii in unique(animal_season_exclude$continuousID)) {
        animal_ii <- animal_season_exclude[animal_season_exclude$continuousID == ii, ]
        duration_ii <- difftime(animal_ii$date[nrow(animal_ii)], animal_ii$date[1], units = units)
        ## Calculate straightness only if at least as long as encounter event
        if(duration_ii >= duration_comp) {
          animal_ii$indur <- !(animal_ii$date > (animal_ii$date[nrow(animal_ii)] - as.difftime(duration_comp, units = units)))
          for(iii in 1:nrow(animal_ii)){
            if(animal_ii$indur[iii]){
              straightnesses_seasonal <- 
                c(straightnesses_seasonal, 
                  animal_ii %>% 
                    dplyr::slice(iii:nrow(animal_ii)) %>% 
                    dplyr::mutate(difftime = difftime(date, date[1], units = units)) %>% 
                    dplyr::filter(difftime <= duration_comp) %>%
                    strtns(.))
            }
          }
        }
      }
      
      ## Determine the confidence interval for "average movement" as mean +/-
      ## sd*sd_multiplier
      str_mean <- mean(straightnesses_seasonal, na.rm = TRUE)
      str_sd <- stats::sd(straightnesses_seasonal, na.rm = TRUE)
      upper <- str_mean + sd_multiplier * str_sd
      lower <- str_mean - sd_multiplier * str_sd
      
      

      #### Local angle analysis of trace behavior

      ## Identify the nearest barrier point to each movement point
      bar_dist_all <- sf::st_distance(x = encounter_i, y = barrier_pts_all)
      encounter_i$bar_pt_min <- barrier_pts_all$barID[apply(bar_dist_all, 1, which.min)]
      
      ## Run through each movement step and evaluate whether the movement angle
      ## is within the CI of the nearest barrier segment
      encounter_i$bar_lcl <- NA
      for(j in 1:(nrow(encounter_i)-1)){
        ## Check if the nearest barrier points are the same point, if so, leave NA
        if(encounter_i$bar_pt_min[j] != encounter_i$bar_pt_min[j+1]){
          ## Check if the nearest barrier points are from the same barrier. If
          ## not then leave the NA.
          if(barrier_pts_all$Name[barrier_pts_all$barID == encounter_i$bar_pt_min[j]] ==
             barrier_pts_all$Name[barrier_pts_all$barID == encounter_i$bar_pt_min[j+1]]){
            ## Otherwise, pull all the barrier points between the nearest points
            ## and check for angles within the CI
            
            ## Calculate angles of the barrier segment and determine the mean
            ## and sd
            lcl_angles <-
              ## Pull all the barrier points between the nearest points
              barrier_pts_all %>% 
              dplyr::filter(barID >= min(encounter_i$bar_pt_min[j], encounter_i$bar_pt_min[j+1]) &
                       barID <= max(encounter_i$bar_pt_min[j], encounter_i$bar_pt_min[j+1])) %>% 
              ## Calculate angles of the barrier segment and determine its mean
              ## and sd
              calc_angle()
            lcl_mean <- 
              lcl_angles %>% 
              circular::circular(zero = 0) %>% 
              circular::mean.circular(na.rm = TRUE) %>% 
              circular::conversion.circular(units = 'degrees') %% 360  ## Convert this to degrees in a 360 deg range
            lcl_sd <- 
              lcl_angles %>% 
              circular::circular(zero = 0) %>% 
              circular::sd.circular(na.rm = TRUE) %>% 
              circular::deg() %% 360
            ## Calculate the upper and lower bounds for comparison. This is
            ## needed in both directions so that movement paths heading in
            ## either direction along a barrier will count as being within the
            ## barrier buffer. Keep these in the 360 degree range.
            lcl_upper <- (lcl_mean + sd_multiplier * lcl_sd) %% 360
            lcl_lower <- (lcl_mean - sd_multiplier * lcl_sd) %% 360
            lcl_upper2 <- (lcl_mean + sd_multiplier * lcl_sd + 180) %% 360
            lcl_lower2 <- (lcl_mean - sd_multiplier * lcl_sd + 180) %% 360
            ## Check whether the step-specific encounter angle lies within the
            ## calculated limits and indicate accordingly
            encounter_i$bar_lcl[j] <- 
              ifelse(angle_in_range(encounter_i$angle[j], lcl_lower, lcl_upper) |
                       angle_in_range(encounter_i$angle[j], lcl_lower2, lcl_upper2),
                     1, 0)
          }
        }
      }
      
      ## Calculate the proportion of parallel steps (1s), taking NAs into account
      lcl_prop <- sum(encounter_i$bar_lcl, na.rm = TRUE)/nrow(encounter_i)
      
      
      
      ### Classify the event ------------------------------------------------------

      ## Earlier testing showed a need to prioritize trace behavior, otherwise
      ## it tends to get classified as avg mvmt when a response to the barrier
      ## is clear.
      
      ## Check whether all bar_lcl are NA. This represents a situation in which
      ## all locations point to the same segment of the barrier and will get
      ## classified as avg movement.
      lcl_NA_check <-
        encounter_i %>% 
        summarize(burstID = unique(burstID),
                  n_tot = dplyr::n(),
                  lcl_NA = sum(is.na(bar_lcl))) %>% 
        sf::st_drop_geometry()
      
      ## Classify the encounter, assigning avg if all bar_lcl are NA
      if(lcl_NA_check$lcl_NA == lcl_NA_check$n_tot){
        classification <- 'Avg_lclNA'
        
        ## If not, check whether the trace and duration thresholds were met and
        ## assign behavior accordingly. Trace behavior involves some persistence,
        ## so to qualify as trace behavior the duration of the encounter must be
        ## at least 3 days (72 hrs), based on preliminary testing.
      } else if(lcl_prop >= 0.2 & duration >= 72){
        classification <- 'Trace'
        
        ## Classify average movement
      } else if(straightness_i >= lower & straightness_i <= upper){
        classification <- "Average_Movement"
        
        ## Classify straight encounters (trace, quick cross) 
      } else if(straightness_i > upper){
        
        ## Check whether movement is more parallel or perpendicular to the
        ## barrier by calculating the mean and sd of angles of the barrier
        ## closest to the movement path, and comparing against the mean angle of
        ## the movement path. To ensure barrier metrics reflect areas nearest
        ## the animal, only calculate the mean and sd for the area within d of
        ## the closest point to the movement path. In some cases this will
        ## encompass the entire barrier, in others it will not.
        
        ## Create a point version of the nearest barrier line that can be used
        ## below
        barrier_pts <-
          barrier_pts_all %>% 
          ## Isolate the feature nearest the movement path
          dplyr::filter(Name == bar_closest)
        
        ## Create a buffer of size d around the barrier point nearest the
        ## movement path that can be used to clip down the barrier, using the
        ## correct d value for the nearest barrier
        d_tmp <- ifelse(length(d) > 1, d[which(barrier$Name == bar_closest)], d)
        barrier_i_buf <-
          barrier_pts %>% 
          ## Calculate distance to the nearest movement point
          dplyr::mutate(
            dist_mvmt = sf::st_distance(
              x = ., 
              y = encounter_i %>%
                dplyr::filter(bar_min == bar_closest) %>% 
                dplyr::slice(which.min(.$bar_dist_km)))) %>% 
          ## Isolate the barrier point nearest to the movement path
          dplyr::slice(which.min(dist_mvmt)) %>% 
          ## Buffer this point by d
          sf::st_buffer(dist = d_tmp, nQuadSegs = 5) %>%   ## Note that nQuadSegs is set to 5 as this was the default value for rgeos::gBuffer in previous versions of BaBA and matches what was done above for the barrier
          sf::st_union()
        
        ## Extract barrier points that fall inside the buffer. Explicitly suppress constant geometry assumption warning by confirming attribute is constant throughout the geometry. See https://github.com/r-spatial/sf/issues/406 for details.
        sf::st_agr(barrier_pts) <- 'constant'   
        barrier_i <- 
          sf::st_intersection(barrier_pts, barrier_i_buf) %>% 
          ## Add in the xy coords as columns, as expected by calc_angle()
          dplyr::mutate(x = sf::st_coordinates(.)[,1],
                        y = sf::st_coordinates(.)[,2])
        
        ## Calculate angles of each barrier segment and determine their mean and sd
        barrier_i_angles <- calc_angle(barrier_i)
        barrier_i_mean <- 
          barrier_i_angles %>% 
          circular::circular(zero = 0) %>% 
          circular::mean.circular(na.rm = TRUE) %>% 
          circular::conversion.circular(units = 'degrees') %% 360  ## Convert this to degrees in a 360 deg range
        barrier_i_sd <- 
          barrier_i_angles %>% 
          circular::circular(zero = 0) %>% 
          circular::sd.circular(na.rm = TRUE) %>% 
          circular::deg() %% 360
        
        ## Calculate the upper and lower bounds for comparison. This is needed
        ## in both directions so that movement paths heading in either direction
        ## along a barrier will count as being within the barrier buffer. Keep
        ## these in the 360 degree range.
        barrier_upper <- (barrier_i_mean + sd_multiplier * barrier_i_sd) %% 360
        barrier_lower <- (barrier_i_mean - sd_multiplier * barrier_i_sd) %% 360
        barrier_upper2 <- (barrier_i_mean + sd_multiplier * barrier_i_sd + 180) %% 360
        barrier_lower2 <- (barrier_i_mean - sd_multiplier * barrier_i_sd + 180) %% 360
        
        ## Calculate the mean encounter angle
        encounter_i_mean <-
          encounter_i$angle %>% 
          circular::circular(zero = 0, units = 'degrees') %>% 
          circular::mean.circular(na.rm = TRUE) %% 360
        
        ## Check whether the mean encounter angle lies within the calculated
        ## limits and assign behavior accordingly
        if(angle_in_range(encounter_i_mean, barrier_lower, barrier_upper) |
           angle_in_range(encounter_i_mean, barrier_lower2, barrier_upper2)){
          classification <- 'Trace'
        } else {
          classification <- ifelse(cross_true > 0, 'Quick_Cross', 'Unknown')
        }
        
        ## Classify non-straight encounters (bounce, back and forth)
      } else{
        
        ## Did the animal change direction of movement after approaching the
        ## barrier?
        
        ## Split the encounter locations into two groups before and after the
        ## point closest to the barrier, including the closest point in each.
        close_pt <-
          encounter_i %>% 
          dplyr::filter(bar_min == bar_closest) %>% 
          dplyr::slice(which.min(.$bar_dist_km)) %>% 
          dplyr::pull(ptsID)
        encounter_i1 <-
          encounter_i %>% 
          dplyr::filter(ptsID <= close_pt)
        encounter_i2 <-
          encounter_i %>% 
          dplyr::filter(ptsID >= close_pt)
        
        ## Calculate the mean of the angles for each group and
        ## sd of the first group
        encounter_i1_mean <- 
          encounter_i1$angle[-nrow(encounter_i1)] %>%  ## The last angle is to the next group, so remove it from here
          circular::mean.circular(na.rm = TRUE) %% 360
        encounter_i1_sd <- 
          encounter_i1$angle[-nrow(encounter_i1)] %>% 
          circular::sd.circular(na.rm = TRUE) %>% 
          circular::deg() %% 360
        encounter_i2_mean <- 
          encounter_i2$angle %>% 
          circular::mean.circular(na.rm = TRUE) %% 360
        
        ## Calculate upper and lower limits
        upper_i1 <- (encounter_i1_mean + sd_multiplier * encounter_i1_sd) %% 360
        lower_i1 <- (encounter_i1_mean - sd_multiplier * encounter_i1_sd) %% 360
        
        ## Calculate the mean resultant length of each group
        encounter_i1_mrl <-
          encounter_i1$angle[-nrow(encounter_i1)] %>% 
          circular::rho.circular(na.rm = TRUE)
        encounter_i2_mrl <-
          encounter_i2$angle[-nrow(encounter_i2)] %>% 
          circular::rho.circular(na.rm = TRUE)
        
        ## Calculate average mean resultant length
        mrl_avg <- mean(c(encounter_i1_mrl, encounter_i2_mrl))
        
        ## Classify the result
        if(is.na(encounter_i1_mean) | is.na(encounter_i1_sd) | is.na(encounter_i2_mean)){
          ## If there is insufficient sample size to be able to run the
          ## comparison classify as unknown
          classification <- 'Unknown_insufficient_n'
          
          ## If the sd is 0 make this avg movment because there is only a single
          ## step
        } else if(encounter_i1_sd == 0){
          classification <- 'Average_sd0'
          
          ## Is the mean of the second group's angles outside the limits of the
          ## first group's angles? If so, is the avg mean resultant length
          ## greater than 0.6?
        } else if(!angle_in_range(encounter_i2_mean, lower_i1, upper_i1)){
          if(mrl_avg >= 0.6){
            ## Classify as bounce, unless it crossed a road, then make unknown
            classification <- ifelse(cross_true > 0, 'Unknown_cross_bounce', 'Bounce')
            
            ## If mrl_avg < 0.6 make back and forth
          } else{
            classification <- 'Back_n_forth'
          }
          
          ## If the second group's angles lie within the range of the first,
          ## classify as back and forth
        } else{
          classification <- 'Back_n_forth'
        }
      }
      
      

      ### Consolidate outputs -----------------------------------------------------

      ## Plot the encounters to check later, if desired
      if (export_images) {
        grDevices::png(paste0(img_path, "/", img_prefix, "_", i, "_", classification, "_", img_suffix, ".png"), width = 12, height = 12, units = "in", res = 300)
        plot(sf::st_geometry(encounter_i), main = paste0(paste(encounter_i$burstID[1], classification, sep = ' - '), '\n', paste(paste(sort(unique(encounter_i$bar_min)), collapse = '-'), season_i, sep = ' - ')),
             sub = paste0("cross = ", cross_true, ", dur =", duration, ", stri =", round(straightness_i, 2), ", str_mn = ",  round(str_mean, 2), ", str_sd = ",  round(str_sd, 2)))
        if(!is.null(img_background)) lapply(img_background, function(x) plot(sf::st_geometry(x), border = 'lightgrey', add = TRUE))
        plot(animal_i %>% 
               dplyr::filter(lubridate::year(date) == lubridate::year(encounter_i$date[1])) %>% 
               dplyr::summarize(do_union = FALSE) %>% 
               sf::st_cast(to = 'LINESTRING') %>% 
               sf::st_geometry(),
             add = TRUE)
        plot(sf::st_geometry(barrier_buffer), border = scales::alpha("red", 0.5), lty = "dashed", add = TRUE)
        plot(sf::st_geometry(barrier), col = "red", lwd = 2, add = TRUE)
        plot(sf::st_geometry(encounter_i), pch = 20, col = "cyan3", type = "o", lwd = 2, add = TRUE)
        plot(sf::st_geometry(encounter_i %>% dplyr::slice(1)), pch = 16, col = "blue", add = TRUE)
        grDevices::dev.off()
      }
      
      ## Prepare data for output
      if(classification %in% c('Trace', 'Quick_Cross', 'Unknown')){
        ang_i <- ifelse(exists('encounter_i_mean'), encounter_i_mean, NA)
        ang_mean <- ifelse(exists('barrier_i_mean'), barrier_i_mean, NA)
        ang_sd <- ifelse(exists('barrier_i_sd'), barrier_i_sd, NA)
        mrl_1 <- NA
        mrl_2 <- NA
        mrl_mean <- NA
      } else if(classification %in% c('Bounce', 'Back_n_forth', 'Unknown_cross_bounce', 'Unknown_insufficient_n', 'Average_sd0')){
        ang_i <- encounter_i2_mean
        ang_mean <- encounter_i1_mean
        ang_sd <- encounter_i1_sd
        mrl_1 <- encounter_i1_mrl
        mrl_2 <- encounter_i2_mrl
        mrl_mean <- mrl_avg
      } else{
        ang_i <- NA
        ang_mean <- NA
        ang_sd <- NA
        mrl_1 <- NA
        mrl_2 <- NA
        mrl_mean <- NA
      }
      ## Summarize barrier crossing information into a single value
      cross_bar_tmp <- 
        encounter_i$cross_bar %>%
        stats::na.omit() %>%
        as.character()
      cross_bar_tmp <- ifelse(length(cross_bar_tmp) == 0, NA, cross_bar_tmp)
      cross_x_tmp <- 
        encounter_i$cross_x %>%
        stats::na.omit() %>%
        as.character()
      cross_x_tmp <- ifelse(length(cross_x_tmp) == 0, NA, cross_x_tmp)
      cross_y_tmp <- 
        encounter_i$cross_y %>%
        stats::na.omit() %>%
        as.character()
      cross_y_tmp <- ifelse(length(cross_y_tmp) == 0, NA, cross_y_tmp)
      
      ## Combine output
      event_tmp <- 
        tibble::tibble(AnimalID = encounter_i$Animal.ID[1],
               encounter = gsub('(\\d+\\w?_\\d+)(_?\\d?)', '\\1', i),
               burstID = i,
               season = season_i,
               barrier = paste(sort(unique(encounter_i$bar_min)), collapse = '-'),
               barrier_n = paste(names(table(encounter_i$bar_min)), table(encounter_i$bar_min), sep = '-', collapse = '-'),
               barrier_min_dist = min(encounter_i$bar_dist_km),
               closest_bar = bar_closest,
               closest_dist = encounter_i %>% 
                 dplyr::filter(bar_min == bar_closest) %>% 
                 dplyr::pull(bar_dist_km) %>% 
                 min(),
               start_time = encounter_i$date[1],
               end_time = encounter_i$date[nrow(encounter_i)],
               duration,
               cross_any = sum(encounter_i$cross_ind, na.rm = TRUE),
               cross_true = cross_true,
               cross_bar = cross_bar_tmp,
               cross_x = as.numeric(cross_x_tmp),
               cross_y = as.numeric(cross_y_tmp),
               class = classification,
               str_i = straightness_i,
               str_mean,
               str_sd,
               ang_i,
               ang_mean,
               ang_sd,
               mrl_1,
               mrl_2,
               mrl_mean,
               lcl_prop,
               easting = sf::st_coordinates(encounter_i)[1, 1],
               northing = sf::st_coordinates(encounter_i)[1, 2])
      event_df <- rbind(event_df, event_tmp)
      
      ## Keep the encounter_i data
      encounter_i_out[[i]] <- encounter_i
    }
    
    ## Close progress bar
    close(pb)
    
    
    
    # Finalize data -----------------------------------------------------------

    print("creating dataframe...")
    
    ## Clean the encounter data
    encounter_final <- 
      encounter %>% 
      dplyr::filter(!duplicated(burstID)) %>% 
      dplyr::left_join(event_df %>% 
                         dplyr::select(burstID, class),
                       by = 'burstID') %>% 
      dplyr::select(Animal.ID, burstID, date, class)
    
    ## Return output as a named list
    return(list(encounters = encounter_final,
                encounter_is = encounter_i_out,
                classification = event_df))
  }
