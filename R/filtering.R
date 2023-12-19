#### MigrObj filtering functions ####

#' #### Filtering object with TRACK_IDs or logical vector ####
#'
#' @examples
#' SampledTracks = sample(tracks.raw(MigrDatS1)[["LABEL"]], 15)
#' filter.w.tracks(MigrObj = MigrDatS1, track.slot = "raw", tracks = SampledTracks)
#' SampledTF = sample(c(T,F), 328, replace = T)
#' filter.w.tracks(MigrObj = MigrDatS1, track.slot = "raw", tracks = SampledTF)
#'
#' @export
#'
filter.w.tracks <- function(MigrObj, track.slot = "raw", tracks = NULL, verbose = FALSE) {
  
  # slot check!
  track.slot = match.arg(track.slot, slots.in.use(MigrObj)$in.tracks)
  # Check if target slots is available
  if (!track.slot %in% slots.in.use(MigrObj)$in.tracks) {
    print("track.slot not filled in MigrObj, please check and provide correct track.slot for filtering.")
    return(NULL)
  }
  
  
  if (is.logical(tracks)) {
    dt.tracks <- data.table::copy(slot(object = MigrObj, name = "tracks")[[track.slot]] )
    strtDim = dim(dt.tracks)
    if (strtDim[1] == length(tracks)) {
      dt.tracks <- dt.tracks[tracks, ]
      tracks.left = dt.tracks[["LABEL"]]
      fltrDim = dim(dt.tracks)
      cat("\nFiltering removed ", strtDim[1]-fltrDim[1], "Tracks.\n\n")
    } else cat("Length of the logical tracks vector did not match number of tracks!\n")
    
  } else if (!missing(tracks)) {
    # Do positive filtering of tracks by leaving the provided tracks!
    dt.tracks <- data.table::copy(slot(object = MigrObj, name = "tracks")[[track.slot]] )
    strtDim = dim(dt.tracks)
    dt.tracks <- dt.tracks[LABEL %in% tracks, ]
    tracks.left = dt.tracks[["LABEL"]]
    fltrDim = dim(dt.tracks)
    cat("\nFiltering removed ", strtDim[1]-fltrDim[1], "Tracks.\n\n")
    
  } else {
    dt.tracks <- data.table::copy(slot(object = MigrObj, name = "tracks")[[track.slot]] )
    cat("\n Argument problem, no filtering done")
    return(MigrObj)
  }
  
  # Extra info debugging
  if (verbose == TRUE) {
    cat("Tracks left:\n")
    print(tracks.left)
    print(dt.tracks)
  }
  
  # remove spots and edges with same names as well
  dt.spots <- data.table::copy(slot(object = MigrObj, name = "spots")[[track.slot]] )
  dt.edges <- data.table::copy(slot(object = MigrObj, name = "edges")[[track.slot]] )
  
  dt.spots <- dt.spots[TRACK_ID %in% tracks.left, ]
  dt.edges <- dt.edges[TRACK_ID %in% tracks.left, ]
  
  # Extra info debugging
  if (verbose == TRUE) {
    print(dt.spots)
    print(dt.edges)
  }
  
  # replace in the original object
  slot(object = MigrObj, name = "tracks")[[track.slot]] <- dt.tracks
  slot(object = MigrObj, name = "spots")[[track.slot]] <- dt.spots
  slot(object = MigrObj, name = "edges")[[track.slot]] <- dt.edges
  
  # Filter metadata
  
  #MigrObj <- upd.metadata(MigrObj)
  #MigrObj <- upd.clustering(MigrObj)
  MigrObj <- sync.MigrObj(MigrObj)
  
  return(MigrObj)
  
}



#' #### Filtering object with tracks ####
#' @name filter.tracks
#'
#' @export
#'
filter.tracks <- function(MigrObj, filter = NULL, track.slot = "raw", tracks = NULL, verbose = FALSE) {
  
  # filtering expression to string
  if (!missing(x = filter)) {
    filter <- rlang::enquo(arg = filter)
    fltrExpr <- ParseFltr(filter = filter)
    if (verbose == TRUE) {print(fltrExpr) }
  }
  
  # Build a filter variable check later
  
  # slot check!
  track.slot = match.arg(track.slot, slots.in.use(MigrObj)$in.tracks)
  # Check if target slots is available
  if (!track.slot %in% slots.in.use(MigrObj)$in.tracks) {
    print("track.slot not filled in MigrObj, please check and provide correct track.slot for filtering.")
    return(NULL)
  }
  
  # For alternative direct spot/track/edge filtering use later
  if (!missing(tracks)) {
    cat("\nFor direct Track filtering use filter.w.tracks function.\n")
  }
  
  # Filtering test
  dt.tracks <- data.table::copy(slot(object = MigrObj, name = "tracks")[[track.slot]] )
  strtDim = dim(dt.tracks)
  dt.tracks <- dt.tracks[eval(fltrExpr), ]
  tracks.left = dt.tracks[["LABEL"]]
  
  fltrDim = dim(dt.tracks)
  cat("\nFiltering removed", strtDim[1]-fltrDim[1], "Tracks and the spots and edges associated with them.\n\n")
  
  # Extra info debugging
  if (verbose == TRUE) {
    cat("Tracks left:\n")
    print(tracks.left)
    print(dt.tracks)
  }
  
  # remove spots and edges with same names as well
  dt.spots <- data.table::copy(slot(object = MigrObj, name = "spots")[[track.slot]] )
  dt.edges <- data.table::copy(slot(object = MigrObj, name = "edges")[[track.slot]] )
  
  dt.spots <- dt.spots[TRACK_ID %in% tracks.left, ]
  dt.edges <- dt.edges[TRACK_ID %in% tracks.left, ]
  
  
  # Extra info debugging
  if (verbose == TRUE) {
    print(dt.spots)
    print(dt.edges)
  }
  
  # replace in the original object
  slot(object = MigrObj, name = "tracks")[[track.slot]] <- dt.tracks
  slot(object = MigrObj, name = "spots")[[track.slot]] <- dt.spots
  slot(object = MigrObj, name = "edges")[[track.slot]] <- dt.edges
  
  MigrObj@metadata$filtering = "propagated"
  
  # Filter metadata
  #MigrObj <- upd.metadata(MigrObj)
  #MigrObj <- upd.clustering(MigrObj)
  MigrObj <- sync.MigrObj(MigrObj)
  
  return(MigrObj)
}

#' #### Filtering object with spots variables ####
#' Note: removes all tracks, spots, and edges linked to filtered out edges
#' @param MigrObj
#' @param filter expression for filtering the spots e.g. POSITION_X > 120 & POSITION_Y > 220
#' @param spot.slot raw, normalized, scale, or other used slot in spots data
#' @param spots spot names to select spots by their SPOT_ID or logical vector of same length as number of spots
#' @param propagate filtering removes also linked edges and tracks logical
#'
#' @export
#'
filter.spots <- function(MigrObj, filter, spot.slot = "raw", spots = NULL, propagate = TRUE, verbose = FALSE) {
  
  # filtering expression to string
  if (!missing(x = filter)) {
    filter <- rlang::enquo(arg = filter)
    fltrExpr <- ParseFltr(filter = filter)
    if (verbose == TRUE) {print(fltrExpr) }
  }
  
  # Build a filter variable check later!
  
  # slot check!
  spot.slot = match.arg(spot.slot, slots.in.use(MigrObj)$in.spots)
  
  # Check if target slots is available
  if (!spot.slot %in% slots.in.use(MigrObj)$in.spots) {
    print("spot.slot not filled in MigrObj, please check and provide correct spot.slot for filtering.")
    return(NULL)
  }
  
  # data.table::copy spots data.table for filtering
  dt.spots <- data.table::copy(slot(object = MigrObj, name = "spots")[[spot.slot]])
  spots.in = dim(dt.spots)[1]
  tracks.in = length(unique(MigrObj@tracks[[spot.slot]][["TRACK_ID"]]))
  
  # For alternative direct spot/track/edge filtering use later
  if (!missing(spots)) {
    
    if (is.character(spots)) {
      # Filtering spots with SPOT_ID (not LABEL)
      dt.spots.fltr <- dt.spots[SPOT_ID %in% spots, ]
      spots.lft <- dt.spots.fltr[["SPOT_ID"]]
      spots.rmv <- dt.spots[["SPOT_ID"]][!dt.spots[["SPOT_ID"]] %in% spots.lft]
      tracks.rmv <- unique(dt.spots[["TRACK_ID"]][!dt.spots[["SPOT_ID"]] %in% spots.lft])
      tracks.lft <- unique(dt.spots[["TRACK_ID"]][!dt.spots[["TRACK_ID"]] %in% tracks.rmv])
    } else if (is.logical(spots) & length(spots) == dim(dt.spots)[1]) {
      # Filtering spots
      dt.spots.fltr <- dt.spots[spots, ]
      spots.lft <- dt.spots.fltr[["SPOT_ID"]]
      spots.rmv <- dt.spots[["SPOT_ID"]][!dt.spots[["SPOT_ID"]] %in% spots.lft]
      tracks.rmv <- unique(dt.spots[["TRACK_ID"]][!dt.spots[["SPOT_ID"]] %in% spots.lft])
      tracks.lft <- unique(dt.spots[["TRACK_ID"]][!dt.spots[["TRACK_ID"]] %in% tracks.rmv])
    }
    cat("\nFiltering with spots vector!\n")
  } else if (is.expression(fltrExpr)) {
    # Filtering spots
    dt.spots.fltr <- dt.spots[eval(fltrExpr), ]
    spots.lft <- dt.spots.fltr[["SPOT_ID"]]
    spots.rmv <- dt.spots[["SPOT_ID"]][!dt.spots[["SPOT_ID"]] %in% spots.lft]
    tracks.rmv <- unique(dt.spots[["TRACK_ID"]][!dt.spots[["SPOT_ID"]] %in% spots.lft])
    tracks.lft <- unique(dt.spots[["TRACK_ID"]][!dt.spots[["TRACK_ID"]] %in% tracks.rmv])
    
  }
  
  if (propagate) {
    MigrObj = filter.w.tracks(MigrObj = MigrObj, track.slot = spot.slot, tracks = tracks.lft, verbose = FALSE)
    cat("\nFiltering removed", spots.in-length(MigrObj@spots[[spot.slot]][["SPOT_ID"]]), "spots and linked edges. \n")
    # cat("\nFiltering removed", tracks.in-length(unique(MigrObj@tracks[[spot.slot]][["TRACK_ID"]])), "tracks.\n\n")
    MigrObj@metadata$filtering = "propagated"
    
    # Filter metadata
    #MigrObj <- upd.metadata(MigrObj)
    MigrObj <- sync.MigrObj(MigrObj)
    
    return(MigrObj)
  } else {
    slot(object = MigrObj, name = "spots")[[spot.slot]] <- dt.spots.fltr
    cat("\nFiltering removed", spots.in-length(MigrObj@spots[[spot.slot]][["SPOT_ID"]]),
        "spots. \nTracks and edges were not filtered as propagation was false. \n")
    MigrObj@metadata$filtering = "non-propagated"
    
    # Filter metadata
    #MigrObj <- upd.metadata(MigrObj)
    #MigrObj <- upd.clustering(MigrObj)
    MigrObj <- sync.MigrObj(MigrObj)
    
    return(MigrObj)
  }
  
}



#' #### Filtering object with edge-data ####
#' Note: removes all tracks, spots, and edges linked to filtered out edges
#'
#' @export
#'
filter.edges <- function(MigrObj, filter, edge.slot = "raw", edges = NULL, verbose = FALSE) {
  
  # filtering expression to string
  if (!missing(x = filter)) {
    filter <- rlang::enquo(arg = filter)
    fltrExpr <- ParseFltr(filter = filter)
    if (verbose == TRUE) {print(fltrExpr) }
  }
  
  # Build a filter variable check later!
  
  # slot check!
  edge.slot = match.arg(edge.slot, slots.in.use(MigrObj)$in.edges)
  # Check if target slots is available
  if (!edge.slot %in% slots.in.use(MigrObj)$in.spots) {
    print("edge.slot not filled in MigrObj, please check and provide correct edge.slot for filtering.")
    return(NULL)
  }
  
  # For alternative direct spot/track/edge filtering use later
  if (!missing(edges)) {
    cat("\nNo direct STE filtering implemented yet!\n")
  }
  
  # Filtering edges
  dt.edges <- data.table::copy(slot(object = MigrObj, name = "edges")[[edge.slot]])
  edges.lft <- dt.edges[eval(fltrExpr), ][["EDGE_ID"]]
  edges.rmv <- dt.edges[["EDGE_ID"]][!dt.edges[["EDGE_ID"]] %in% edges.lft]
  tracks.rmv <- unique(dt.edges[["TRACK_ID"]][!dt.edges[["EDGE_ID"]] %in% edges.lft])
  tracks.lft <- unique(dt.edges[["TRACK_ID"]][!dt.edges[["TRACK_ID"]] %in% tracks.rmv])
  
  MigrObj = filter.w.tracks(MigrObj = MigrObj, track.slot = edge.slot, tracks = tracks.lft, verbose = FALSE)
  cat("Edges filtering always removes linked tracks and spots associated with them!\n\n")
  
  # Filter metadata
  #MigrObj <- upd.metadata(MigrObj)
  #MigrObj <- upd.clustering(MigrObj)
  MigrObj <- sync.MigrObj(MigrObj)
  
  return(MigrObj)
  
}


#' ParseFltr
#' fltrExpr <- ParseFltr(filter = X > 3 | Y < 14)
#' MigrDatXML_Z0_T00_C1@tracks$raw[eval(fltrExpr), ]
#'
#' 
#'
ParseFltr <- function(filter) {
  
  # filter <- rlang::enquo(arg = filter)
  
  fltr <- if (tryCatch(expr = rlang::is_quosure(x = filter), error = function(...) FALSE)) {
    filter
  } else if (is.call(x = rlang::enquo(arg = filter))) {
    rlang::enquo(arg = filter)
  } else {
    parse(text = filter)
  }
  # fltr <- as_label(fltr) # shortens expression with ... and hence does not work as intended
  # as_label() will replace as.character in the future!
  # suppressWarnings()
  fltr = as.character(fltr)[2]
  fltr <- (parse(text=fltr))
  return(fltr)
}



#' @name Apply_TMate_Filter 
#' 
#' 
Apply_TMate_Filter <- function(MigrObj) {
  # Create True False filter for spots 
  tm_filter <- MigrObj@metadata$track_filter 
  # Apply filter 
  MigrObj <- filter.tracks(MigrObj, track.slot = "raw", filter = TRACK_ID %in% tm_filter) 
  return(MigrObj) 
} 

#' @name rm_unmated_spots
#'
#'
rm_unmated_spots <- function(MigrObj) {
  # remove "ninTrack" spots, which are non-track assigned spots in the TrackMate xml data. 
  MigrObj <- filter.spots(MigrObj, filter = TRACK_ID != "ninTracks", spot.slot = "raw")
  return(MigrObj)
}




#' @name DRs.TFfilterL
#'
#'
DRs.TFfilterL <- function(DRs, TF_filter) {
  lapply(DRs, DR.TFfilter, TF_filter)
}

#' @name DR.TFfilter
#'
#'
DR.TFfilter <- function(DR, TF_filter) {
  if (class(DR)[1] == "matrix") {
    if ( dim(DR)[1] == length(TF_filter) ) {
      DR <- DR[TF_filter,]
    } else cat("Warning filtering vector length error!") 
  } 
  return(DR)
}

#' @name DRs.TFfilter
#'
#'
DRs.TFfilter <- function(DRs, TF_filter) {
  for (n in names(DRs)) {
    if (class(DRs[[n]])[1] == "matrix") {
      if ( dim(DRs[[n]])[1] == length(TF_filter) ) {
        DRs[[n]] <- DRs[[n]][TF_filter,]
      } else cat("Warning filtering vector length error!") 
    } 
  }
  return(DRs)
}



#' @name DRs.syncL
#'
#'
DRs.syncL <- function(MigrObj, labels, type = c("S","T","E")) {
  
  DRs.all = dimreducts(MigrObj)
  
  if (type == "S") {
    DRs <- DRs.all[[type]] 
    DRs <- DRs.TFfilterL(DRs = DRs, TF_filter = spots.raw(MigrObj)[["LABEL"]] %in% labels)
    DRs.all[[type]] <- DRs 
  }
  
  if (type == "T") {
    DRs <- DRs.all[[type]] 
    DRs <- DRs.TFfilterL(DRs = DRs, TF_filter = tracks.raw(MigrObj)[["LABEL"]] %in% labels)
    DRs.all[[type]] <- DRs 
  }
  
  if (type == "E") {
    DRs <- DRs.all[[type]] 
    DRs <- DRs.TFfilterL(DRs = DRs, TF_filter = edges.raw(MigrObj)[["LABEL"]] %in% labels)
    DRs.all[[type]] <- DRs 
  }
  
  dimreducts(MigrObj) <- DRs.all 
  return(MigrObj)
}


#' @name sync.MigrObj
#'
#'
#' @export
sync.MigrObj <- function(MigrObj) {
  
  # get spots labels present in all spots tables.
  spot.labels <- Reduce(intersect, lapply((spotsLS(MigrObj)[!sapply(spotsLS(MigrObj), is.null)]), `[[`, "LABEL") )
  # Check for functionality
  if (is.null(spot.labels)) {stop("No spot labels retrieved!")}
  spots.raw(MigrObj) <- spots.raw(MigrObj)[LABEL %in% spot.labels]
  tryCatch(spots.normalized(MigrObj) <- spots.normalized(MigrObj)[LABEL %in% spot.labels], error = function(e) {print("no spot norm")} )
  tryCatch(spots.scaled(MigrObj) <- spots.scaled(MigrObj)[LABEL %in% spot.labels], error = function(e) {print("no spot scaled")} )
  tryCatch(spots.meta(MigrObj) <- spots.meta(MigrObj)[LABEL %in% spot.labels], error = function(e) {print("no spot metadata")} )
  tryCatch(roi_points.raw(MigrObj) <- roi_points.raw(MigrObj)[spot.labels], error = function(e) {print("roi_points data err")} )
  tryCatch(spots.clusters(MigrObj) <- spots.clusters(MigrObj)[LABEL %in% spot.labels], error = function(e) {print("clusters data err")} )
  MigrObj <- DRs.syncL(MigrObj, labels = spot.labels, type = "S")
  #tryCatch(  (MigrObj) <-    (MigrObj)[spot.labels], error = function(e) {print("roi_points data err")} )
  
  track.labels <- Reduce(intersect, lapply((tracksLS(MigrObj)[!sapply(tracksLS(MigrObj), is.null)]), `[[`, "LABEL") )
  if (is.null(track.labels)) {
    stop("No track labels retrieved!")
  }
  tracks.raw(MigrObj) <- tracks.raw(MigrObj)[LABEL %in% track.labels]
  tryCatch(tracks.normalized(MigrObj) <- tracks.normalized(MigrObj)[LABEL %in% track.labels], error = function(e) {print("no track norm")} )
  tryCatch(tracks.scaled(MigrObj) <- tracks.scaled(MigrObj)[LABEL %in% track.labels], error = function(e) {print("no track scaled")} )
  tryCatch(tracks.meta(MigrObj) <- tracks.meta(MigrObj)[LABEL %in% track.labels], error = function(e) {print("no track metadata")} )
  tryCatch(tracks.clusters(MigrObj) <- tracks.clusters(MigrObj)[LABEL %in% track.labels], error = function(e) {print("no track clusters")} )
  MigrObj <- DRs.syncL(MigrObj, labels = tracks.labels, type = "T")
  
  edge.labels <- Reduce(intersect, lapply((edgesLS(MigrObj)[!sapply(edgesLS(MigrObj), is.null)]), `[[`, "LABEL") )
  if (is.null(edge.labels)) {
    stop("No edge labels retrieved!")
  }
  edges.raw(MigrObj) <- edges.raw(MigrObj)[LABEL %in% edge.labels]
  tryCatch(edges.normalized(MigrObj) <- edges.normalized(MigrObj)[LABEL %in% edge.labels], error = function(e) {print("no edge norm")} )
  tryCatch(edges.scaled(MigrObj) <- edges.scaled(MigrObj)[LABEL %in% edge.labels], error = function(e) {print("no edge scaled")} )
  tryCatch(edges.meta(MigrObj) <- edges.meta(MigrObj)[LABEL %in% edge.labels], error = function(e) {print("no edge meta")} )
  tryCatch(edges.clusters(MigrObj) <- edges.clusters(MigrObj)[LABEL %in% edge.labels], error = function(e) {print("no edge clusters")} )
  MigrObj <- DRs.syncL(MigrObj, labels = edge.labels, type = "E")
  
  return(MigrObj)
}
