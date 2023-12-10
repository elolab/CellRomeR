
#### data wrangling functions ####

#' Need conversion table/system for variables so that all objects will have same basic features
#'
#' Function checks the variable names in data.table according to given list and uses list names as defaults to convert names into.
#' It also reorders columns to same order as given in the lists and orders rows in TRACK_ID, FRAME/SPOT_SOURCE_ID order.
#'
#'
#'
#' @param dtIn a data.table
#' @param type type of the table
#' @param stdsS optional list of standardization names for spots
#' @param stdsT optional list of standardization names for tracks
#' @param stdsE optional list of standardization names for edges
#'
#'@examples
#' standardize_dt(readTrckMatCSVs("../data/Tcell_tracks/Exp 1 20180126 ICAM 1_spots.csv"), type = "S")
#' standardize_dt(readTrckMatCSVs("../data/Tcell_tracks/Exp 1 20180126 ICAM 1_tracks.csv"), type = "T")
#' standardize_dt(readTrckMatCSVs("../data/Tcell_tracks/Exp 1 20180126 ICAM 1_edges.csv"), type = "E")
#'
#' @export
standardize_dt <- function(dtIn, type = c("S","T","E"), stdsS = NULL, stdsT = NULL, stdsE = NULL) {


  StdP = "^LABEL|^label|ID$|id$|INDEX|^TIME|QUALITY|LOCATION|START|STOP|GAP$|^NUMBER|DURATION|^POSITION|^FRAME|VISIBILITY|LINK_COST|^EDGE_TIME$"
  # Standards is also stored in the object and could be retrieved from there, except function does not use MigrDat.

  sort.spots <- function(dt.spots) {
    dt.spots[, SPOT_ID := as.numeric(SPOT_ID)]
    setorder(dt.spots, TRACK_ID, FRAME) # TRACK_ID should be split into char and numeric part for correct ordering before sorting
    dt.spots[, SPOT_ID := as.character(SPOT_ID)]
  }

  sort.edges <- function(dt.edges) {
    dt.edges[, SPOT_SOURCE_ID := as.numeric(SPOT_SOURCE_ID)]
    setorder(dt.edges, TRACK_ID, SPOT_SOURCE_ID) # TRACK_ID should be split into char and numeric part for correct ordering before sorting
    dt.edges[, SPOT_SOURCE_ID := as.character(SPOT_SOURCE_ID)]
  }

  sort.tracks <- function(dt.tracks) {
    dt.tracks[, TRACK_ID := as.numeric(TRACK_ID)]
    setorder(dt.tracks, TRACK_ID)
    dt.tracks[, TRACK_ID := as.character(TRACK_ID)]
  }
  # sort.tracks <- function(dt.tracks) {
  #   setorder(dt.tracks, LABEL)
  # }

  # Example lists provided.
  # Standard for spots
  if (type == "S") {
    if (length(stdsS) <1 | is.null(stdsS)) {
      stdsX = list(LABEL = c("LABEL","name"),
                   SPOT_ID = c("SPOT_ID","spot_id","^ID$"),
                   TRACK_ID  = c("^TRACK_ID$","^track_id"),
                   TRACK_ID2  = c("^TRACK_ID2","^track_id2"),
                   POSITION_X = c("POSITION_X","POSIT.*X","posit.*x","Posit.*x","^X$","^x$"),
                   POSITION_Y = c("POSITION_Y","POSIT.*Y","posit.*y","Posit.*y","^Y$","^y$"),
                   POSITION_Z = c("POSITION_Z","POSIT.*Z","posit.*z","Posit.*z","^Z$","^z$"),
                   POSITION_T = c("POSITION_T","POSIT.*T$","posit.*t$","Posit.*t$"),
                   FRAME = c("FRAME", "frame","trace"),
                   TIME = c("time","TIME")
      )
    } else {stdsX = stdsS}
  }

  # Standard for tracks
  if (type == "T") {
    if (length(stdsT) <1 | is.null(stdsT)) {
      stdsX = list(LABEL = c("LABEL","name"),
                   TRACK_ID = c("TRACK_ID$"),
                   TRACK_INDEX = c("TRACK_INDEX"),
                   TRACK_X_LOCATION = c("TRACK_X_LOCATION","^X$","^x$"),
                   TRACK_Y_LOCATION = c("TRACK_Y_LOCATION","^Y$","^y$"),
                   TRACK_Z_LOCATION = c("TRACK_Z_LOCATION","^Z$","^z$"),
                   TRACK_START = c("TRACK_START"),
                   TRACK_STOP = c("TRACK_STOP"),
                   TRACK_MEAN_QUALITY = c("TRACK_MEAN_QUALITY"),
                   LONGEST_GAP = c("LONGEST_GAP"),
                   NUMBER_COMPLEX = c("NUMBER_COMPLEX"),
                   NUMBER_GAPS = c("NUMBER_GAPS"),
                   NUMBER_MERGES = c("NUMBER_MERGES"),
                   NUMBER_SPLITS = c("NUMBER_SPLITS"),
                   NUMBER_SPOTS = c("NUMBER_SPOTS"),
                   TRACK_DURATION = c("TRACK_DURATION"),
                   TRACK_DISPLACEMENT = c("TRACK_DISPLACEMENT"),
                   TRACK_MIN_SPEED = c("TRACK_MIN_SPEED"),
                   TRACK_MAX_SPEED = c("TRACK_MAX_SPEED"),
                   TRACK_MEAN_SPEED = c("TRACK_MEAN_SPEED"),
                   TRACK_MEDIAN_SPEED = c("TRACK_MEDIAN_SPEED"),
                   TRACK_STD_SPEED = c("TRACK_STD_SPEED"),
                   TOTAL_DISTANCE_TRAVELED = c("TOTAL_DISTANCE_TRAVELED"),
                   MAX_DISTANCE_TRAVELED = c("MAX_DISTANCE_TRAVELED"),
                   CONFINEMENT_RATIO = c("CONFINEMENT_RATIO"),
                   MEAN_STRAIGHT_LINE_SPEED = c("MEAN_STRAIGHT_LINE_SPEED"),
                   LINEARITY_OF_FORWARD_PROGRESSION = c("LINEARITY_OF_FORWARD_PROGRESSION"),
                   MEAN_DIRECTIONAL_CHANGE_RATE = c("MEAN_DIRECTIONAL_CHANGE_RATE")
      )
    } else {stdsX = stdsT}
  }

  # Standard for edges
  if (type == "E") {
    if (length(stdsE) <1 | is.null(stdsE)) {
      stdsX = list(LABEL = c("LABEL","name"),
                   EDGE_ID = c("EDGE_ID","Edge_ID","edge_id"),
                   TRACK_ID = c("TRACK_ID$"),
                   EDGE_INDEX = c("TRACK_INDEX"),
                   SPOT_SOURCE_ID = c("SPOT_SOURCE_ID","spot_target_id"),
                   SPOT_TARGET_ID = c("SPOT_TARGET_ID","spot_source_id"),
                   EDGE_X_LOCATION = c("EDGE_X_LOCATION","^X$","^x$"),
                   EDGE_Y_LOCATION = c("EDGE_Y_LOCATION","^Y$","^y$"),
                   EDGE_Z_LOCATION = c("EDGE_Z_LOCATION","^Z$","^z$"),
                   LINK_COST = c("LINK_COST")
      )
    } else {stdsX = stdsE}

  }

  # Column names in and out
  colmnIn = colnames(dtIn)
  colmnOut = colnames(dtIn)

  # Check and fix names
  for (var in names(stdsX)) {
    colmnOut[grep(combine_patterns(stdsX[[var]]), colmnOut)] = as.character(var)
  }

  colnames(dtIn) = colmnOut
  # Standard ordering
  standard_cols = names(stdsX)[names(stdsX) %in% colmnOut]
  nonstandard_cols = colmnOut[colmnOut %nin% names(stdsX)]
  colorder = c(standard_cols, nonstandard_cols)

  setcolorder(dtIn, colorder)

  # Convert types to what they look like.
  dtIn <- type.convert(dtIn, as.is = TRUE)

  # set some ID-numbers to characters for matching and not confusing to indices
  ToChar = c("TRACK_ID","TRACK_ID2","SPOT_SOURCE_ID","SPOT_TARGET_ID","ID","SPOT_ID")
  Cols2Char =  ToChar[ToChar %in% names(dtIn) ]
  for (col in Cols2Char) set(dtIn, j = col, value = as.character(dtIn[[col]]))
  print( data.frame(OldCols = colmnIn, NewCols = colmnOut, Order = colorder) )

  if (type == "S") {sort.spots(dtIn)}
  if (type == "T") {sort.tracks(dtIn)}
  if (type == "E") {sort.edges(dtIn)}

  return(dtIn)

}



#' #### Fetch track_ids for spots from edge table ####
#' track_ids4spotsM
#' fastest of three implementations tried
#'
#'@param dtExml is a TrackMate edges datatable imported from TMxml
#'@param dtSxml is a TrackMate spots datatable imported from TMxml
#'
#'
#'@examples
#'
#'@export
track_ids4spotsM <- function(dtExml,dtSxml) {
  # Using matching
  colsE = c("TRACK_ID","SPOT_SOURCE_ID","SPOT_TARGET_ID")
  es_ss_ts = dtExml[, .SD , .SDcols=colsE]
  for (col in colsE[2:3]) set(es_ss_ts, j = col, value = as.character(es_ss_ts[[col]]))
  es_ss_tsA <- rbindlist(list(es_ss_ts, list("ninTracks", 1e5, 1e5)  ))
  es_ss_tsA[.N]
  nomatchIdx = (dim(es_ss_tsA)[1])
  mch = dtSxml[["ID"]]
  pos1 = match(mch, es_ss_tsA[["SPOT_SOURCE_ID"]], nomatch = nomatchIdx)
  pos2 = match(mch, es_ss_tsA[["SPOT_TARGET_ID"]], nomatch = nomatchIdx)
  idx = data.table(pos1 = pos1, pos2 = pos2)[, min:= do.call(pmin, .SD)][,min]
  track_ids = es_ss_tsA[["TRACK_ID"]][idx]

  if (any(is.na(track_ids))) {
    print("Warning some spots were not matched!")
  }
  return(track_ids)
}


#' #### read data from TrackMate XML file
#'
#' fetches lists of ROI points from TrackMate XML file
#'
#' @param TMxml parsed TrackMate XML-file
#'
#' @examples
#' ROI_POINTS <- getROIpoints_TMxml(TMxml)
#'
#'
#' @export
getROIpoints_TMxml <- function(TMxml, IDs = NULL) {

  #funs
  strsplt.Space <- function(x) {
    strsplit(x,split = " ") }

  strsplt.coord <- function(x) {
    coorX = strsplit(x,split = " ")
    coorX = as.numeric(coorX)
    return(coorX) }
  coordNum <- function(x) {as.numeric(x) }

  ## Ids for the coords
  if (is.null(IDs)) {
    IDs = xpathSApply(TMxml, "//Spot", xmlGetAttr, "name")
  }


  # get and process coordinates
  coords = XML::xmlValue(x = XML::getNodeSet(TMxml, "//Spot"))
  names(coords) = IDs

  CoordSet = sapply(coords, FUN = strsplt.Space)
  CoordSet = lapply(CoordSet, FUN = coordNum)

  return(CoordSet)

}

















#### Data transfer functions ####

#' ##### speed etc. to spots #####
#'
#'
#' @export
speed2spots.raw <- function(MigrDat, NAs0 = TRUE) {

  spots <- copy(spots.raw(MigrDat))
  edges <- copy(edges.raw(MigrDat))

  spots <- spots[edges, SPEED := i.SPEED, on = c("SPOT_ID" = "SPOT_SOURCE_ID")]
  spots <- spots[edges, SPEED := i.SPEED, on = c("SPOT_ID" = "SPOT_TARGET_ID")]

  # The SPEED and DISPLACEMENT are identical! >> Speed is displacement/delta_time
  #spots <- spots[edges, DISPLACEMENT := i.DISPLACEMENT, on = c("SPOT_ID" = "SPOT_SOURCE_ID")]
  #spots <- spots[edges, DISPLACEMENT := i.DISPLACEMENT, on = c("SPOT_ID" = "SPOT_TARGET_ID")]

  spots <- spots[edges, DIRECTIONAL_CHANGE_RATE := i.DIRECTIONAL_CHANGE_RATE, on = c("SPOT_ID" = "SPOT_TARGET_ID")]
  spots <- spots[edges, DIRECTIONAL_CHANGE_RATE := i.DIRECTIONAL_CHANGE_RATE, on = c("SPOT_ID" = "SPOT_SOURCE_ID")]

  if (NAs0 == TRUE) {
    spots[is.na(SPEED), SPEED := 0]
    #spots[is.na(DISPLACEMENT), DISPLACEMENT := 0]
    spots[is.na(DIRECTIONAL_CHANGE_RATE), DIRECTIONAL_CHANGE_RATE := 0]
  }

  MigrDat@spots$raw <- spots

  return(MigrDat)
}


#'@name init.metadata
#'
#'
init.metadata <- function(MigrObj) {
  
  e = grep("QUALITY",colnames(spots.raw(MigrObj)))-1
  MigrObj@metadata[["S"]] <- spots.raw(MigrObj)[,1:e]
  MigrObj@metadata[["T"]] <- tracks.raw(MigrObj)[,1:2]
  MigrObj@metadata[["E"]] <- edges.raw(MigrObj)[,1:2]
  cat("\nMetadata data.tables for spots, tracks, and edges was added to metadata slot.")
  
  return(MigrObj)
  
}


init.clustering <- function(MigrObj) {
  
  default.clust.names <- c("kmeans","hclusts","ILoRegclusters","Seurclusters")
  
  MigrObj@clustering[["S"]] <- spots.raw(MigrObj)[,1:2]
  MigrObj@clustering[["S"]][["kmeans"]] <- 0
  MigrObj@clustering[["S"]][["hclusts"]] <- 0
  MigrObj@clustering[["S"]][["ILoRegclusters"]] <- 0
  MigrObj@clustering[["S"]][["Seurclusters"]] <- 0
  
  MigrObj@clustering[["T"]] <- tracks.raw(MigrObj)[,1:2]
  MigrObj@clustering[["T"]][["kmeans"]] <- 0
  MigrObj@clustering[["T"]][["hclusts"]] <- 0
  MigrObj@clustering[["T"]][["ILoRegclusters"]] <- 0
  MigrObj@clustering[["T"]][["Seurclusters"]] <- 0
  
  MigrObj@clustering[["E"]] <- edges.raw(MigrObj)[,1:2]
  MigrObj@clustering[["E"]][["kmeans"]] <- 0
  MigrObj@clustering[["E"]][["hclusts"]] <- 0
  MigrObj@clustering[["E"]][["ILoRegclusters"]] <- 0
  MigrObj@clustering[["E"]][["Seurclusters"]] <- 0
  
  cat("\nClustering data.tables for spots, tracks, and edges was added to clustering slot.\n")
  return(MigrObj)
}