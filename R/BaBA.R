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

#'  eturn \code{BaBA} returns a list with two components:
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
#'@references Xu W, Dejid N, Herrmann V, Sawyer H, Middleton AD. Barrier
#'  Behaviour Analysis (BaBA) reveals extensive effects of fencing on
#'  wide-ranging ungulates. J Appl Ecol.
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

