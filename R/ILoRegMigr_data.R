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
    data.table::setorder(dt.spots, TRACK_ID, FRAME) # TRACK_ID should be split into char and numeric part for correct ordering before sorting
    dt.spots[, SPOT_ID := as.character(SPOT_ID)]
  }
  
  sort.edges <- function(dt.edges) {
    dt.edges[, SPOT_SOURCE_ID := as.numeric(SPOT_SOURCE_ID)]
    data.table::setorder(dt.edges, TRACK_ID, SPOT_SOURCE_ID) # TRACK_ID should be split into char and numeric part for correct ordering before sorting
    dt.edges[, SPOT_SOURCE_ID := as.character(SPOT_SOURCE_ID)]
  }
  
  sort.tracks <- function(dt.tracks) {
    dt.tracks[, TRACK_ID := as.numeric(TRACK_ID)]
    data.table::setorder(dt.tracks, TRACK_ID)
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
  nonstandard_cols = colmnOut[! (colmnOut %in%  names(stdsX))]
  colorder = c(standard_cols, nonstandard_cols)
  
  data.table::setcolorder(dtIn, colorder)
  
  # Convert types to what they look like.
  dtIn <- type.convert(dtIn, as.is = TRUE)
  
  # set some ID-numbers to characters for matching and not confusing to indices
  ToChar = c("TRACK_ID","TRACK_ID2","SPOT_SOURCE_ID","SPOT_TARGET_ID","ID","SPOT_ID")
  Cols2Char =  ToChar[ToChar %in% names(dtIn) ]
  for (col in Cols2Char) data.table::set(dtIn, j = col, value = as.character(dtIn[[col]]))
  
  if (type == "S") {sort.spots(dtIn)}
  if (type == "T") {sort.tracks(dtIn)}
  if (type == "E") {sort.edges(dtIn)}
  
  return(dtIn)
  
}



#' #### Fetch track_ids for spots from edge table ####
#' @name track_ids4spotsM
#'
#' @param dtExml is a TrackMate edges datatable imported from TMxml
#' @param dtSxml is a TrackMate spots datatable imported from TMxml
#' @description
#' To fetch track_ids for spots to add them to spots data.table
#' fastest of the three tested implementations 
#'
#'
#' @examples
#'
#' @export
track_ids4spotsM <- function(dtExml,dtSxml) {
  # Using matching
  colsE = c("TRACK_ID","SPOT_SOURCE_ID","SPOT_TARGET_ID")
  es_ss_ts = dtExml[, .SD , .SDcols=colsE]
  for (col in colsE[2:3]) data.table::set(es_ss_ts, j = col, value = as.character(es_ss_ts[[col]]))
  es_ss_tsA <- data.table::rbindlist(list(es_ss_ts, list("ninTracks", 1e5, 1e5)  ))
  es_ss_tsA[.N]
  nomatchIdx = (dim(es_ss_tsA)[1])
  mch = dtSxml[["ID"]]
  pos1 = match(mch, es_ss_tsA[["SPOT_SOURCE_ID"]], nomatch = nomatchIdx)
  pos2 = match(mch, es_ss_tsA[["SPOT_TARGET_ID"]], nomatch = nomatchIdx)
  idx = data.table::data.table(pos1 = pos1, pos2 = pos2)[, min:= do.call(pmin, .SD)][,min]
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
    IDs = XML::xpathSApply(TMxml, "//Spot", XML::xmlGetAttr, "name")
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
  
  spots <- data.table::copy(spots.raw(MigrDat))
  edges <- data.table::copy(edges.raw(MigrDat))
  
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
  
  # e is column index before QUALITY column. 
  e = grep("QUALITY",colnames(spots.raw(MigrObj)))-1
  MigrObj@metadata[["S"]] <- spots.raw(MigrObj)[,1:e]
  MigrObj@metadata[["T"]] <- tracks.raw(MigrObj)[,1:2]
  MigrObj@metadata[["E"]] <- edges.raw(MigrObj)[,1:2]
  
  return(MigrObj)
}

#' @name init.clustering
#' @description
#' Prepares `clustering` slot with correct number of rows for clustering data. 
#' Can be used also to wipe all or reset `clustering` slot. 
#' Should actually be a data.table for consistency. 
#' There may be a problem as fetched data does not allow data.table type of slicing. 
#'
init.clustering <- function(MigrObj) {
  
  default.clust.names <- c("kmeans","hclusts","ILoRegclusters")
  
  MigrObj@clustering[["S"]] <- spots.raw(MigrObj)[,1:2]
  MigrObj@clustering[["S"]][["kmeans"]] <- 0
  MigrObj@clustering[["S"]][["hclusts"]] <- 0
  MigrObj@clustering[["S"]][["ILoRegclusters"]] <- 0
  
  MigrObj@clustering[["T"]] <- tracks.raw(MigrObj)[,1:2]
  MigrObj@clustering[["T"]][["kmeans"]] <- 0
  MigrObj@clustering[["T"]][["hclusts"]] <- 0
  MigrObj@clustering[["T"]][["ILoRegclusters"]] <- 0
  
  MigrObj@clustering[["E"]] <- edges.raw(MigrObj)[,1:2]
  MigrObj@clustering[["E"]][["kmeans"]] <- 0
  MigrObj@clustering[["E"]][["hclusts"]] <- 0
  MigrObj@clustering[["E"]][["ILoRegclusters"]] <- 0
  
  #cat("\nClustering data.tables for spots, tracks, and edges was added to clustering slot.\n")
  return(MigrObj)
}

#' ### Get data ####
#' @name getdtcols
#' 
getdtcols <- function(MigrObj, dat.slot = "scaled", type = c("S","T","E"),
                      predef =  "none", vars = NULL, incl.pattern = NULL,
                      excld.pattern = NULL, numerics = TRUE, rnames = FALSE) {
  
  predef = match.arg(predef, c("none", "technical", "morphological", "morphplus", "clust", "coord"), several.ok = FALSE)
  
  # Standard variables like coordinates IDs etc. which we are not interested in as data 
  StdP = Std.pattern(MigrObj)
  
  # Fetch data.table
  dt <- getdt(MigrObj, dat.slot = dat.slot, type = type)
  
  # Add rownames if needed for something. Each data type should have the LABEL column for IDs 
  if (rnames) {
    rownames(dt) <- dt[["LABEL"]]
  }
  
  
  # get selected columns
  if (is.null(incl.pattern) & is.null(vars) & predef=="none" ) {
    # standard_varnames exclusion only
    data = dt.colSpt(dt, excld.pattern = excld.pattern, predef = "none", numerics = numerics)
  } else if (is.null(excld.pattern) & is.null(vars) & predef=="none" ) {
    # incl.pattern
    data = dt.colSpt(dt, excld.pattern = NULL, incl.pattern = incl.pattern, predef = "none", numerics = numerics)
  } else if (!is.null(excld.pattern) & !is.null(incl.pattern) & is.null(vars) & predef=="none" ) {
    # incl.pattern and exclusion
    data = dt.colSpt(dt, excld.pattern = excld.pattern, incl.pattern = incl.pattern, predef = "none", numerics = numerics)
  } else if (is.null(incl.pattern) & is.null(excld.pattern) & predef=="none") {
    # vars
    colz <- vars[vars %in% colnames(dt)]
    data <-  dt[,..colz]
  } else if (predef %in% c("technical", "morphological", "morphplus", "clust", "coord") ) {
    # predefined
    data = dt.colSpt(dt, excld.pattern = NULL, predef = predef, numerics = numerics)
  } else {print("Incompatible variable choosing at the moment. Check if dat.slot exists!")
    return(NULL)
  }
  return(data)
}

getlabels <- function(MigrObj,  dat.slot = dat.slot, type = type) {
  # arguments parsing
  type = match.arg(type, c("S", "T", "E"))
  # data.slot check!
  # in use could be changed to return list without in. pre...
  in.STE <- ifelse(type == "S", "in.spots", ifelse(type == "T", "in.tracks",  "in.edges"))
  lbl.STE <- ifelse(type == "S", "LABEL", ifelse(type == "T", "LABEL",  "LABEL"))
  STE <- gsub("in.", "", in.STE)
  dat.slot = match.arg(dat.slot, slots.in.use(MigrObj)[[in.STE]])
  # if error here, it means, that the data is not available, please check the data.slot and type arguments!
  # fetch data.table
  labels <- (slot(object = MigrObj, name = STE)[[dat.slot]])[[lbl.STE]]
  return(labels)
}

#### Data fetching functions ####
#'
#'
getdt <- function(MigrObj, dat.slot = "raw", type = c("S","T","E")) {
  # arguments parsing
  type = match.arg(type, c("S", "T", "E"))
  # data.slot check!
  # in use could be changed to return list without in. pre...
  in.STE <- ifelse(type == "S", "in.spots", ifelse(type == "T", "in.tracks",  "in.edges"))
  STE <- gsub("in.", "", in.STE)
  dat.slot = match.arg(dat.slot, slots.in.use(MigrObj)[[in.STE]])
  # if error here, it means, that the data is not available, please check the data.slot and type arguments!
  # fetch data.table
  dt <- data.table::copy(slot(object = MigrObj, name = STE)[[dat.slot]] )
  return(dt)
}

#' #### slots in MigrObj in use function ####
#'
#' @export
slots.in.use <- function(MigrDatObj) {
  SlotZ = slot.usage(MigrDatObj = MigrDatObj)
  in.spots = names(MigrDatObj@tracks)[SlotZ$spots]
  in.tracks = names(MigrDatObj@tracks)[SlotZ$tracks]
  in.edges = names(MigrDatObj@tracks)[SlotZ$edges]
  in.roi_points = names(MigrDatObj@roi_points)[SlotZ$roi_points]
  return(list(in.spots=in.spots, in.tracks=in.tracks, in.edges=in.edges, in.roi_points= in.roi_points))
}

#### MigrObj data auxiliary Functions ####
#' #### slots in MigrObj not NULL function ####
#'
#'
#'
#' @export
slot.usage <- function(MigrDatObj) {
  sapply(slotNames(MigrDatObj), function(s) {
    (!sapply(slot(MigrDatObj, s), FUN = is.null) )
  } )
}

#### Function to fetch dt columns ####
#' Works with some arguments only!
#' @name dt.colSpt
#' @description
#' A complex general data fetcher to be used with various functions using data. 
#' Tries to handle data fetching for different purposes including clustering and plotting. 
#' Should be split into smaller and more specific functions. 
#' Still think the patterns usage in getting is good working idea. Implementation needs re-thinking. 
#'
#' @export
dt.colSpt <- function(dt, excld.pattern = NULL, incl.pattern = NULL, predef = "none", vars = NULL, numerics = F,
                      StdP = NULL, TechP = NULL, MorphP = NULL, MorphPpos = NULL) {
  `%nin%` <- Negate(`%in%`)
  
  # If vars defined, return just dt with vars!  
  if (length(vars)>0) {
    colz <- vars[vars %in% colnames(dt)]
    dt <-  dt[,..colz]
    return(dt)
  }
  
  # Parse exclusion patterns with  "standard_varnames" option
  if (is.null(StdP)) {
    StdP = "^LABEL|^label|ID$|ID[1-9]$|id$|INDEX|^TIME|QUALITY|LOCATION|START|STOP|GAP$|^NUMBER|DURATION|^POSITION|^FRAME|VISIBILITY|LINK_COST|^EDGE_TIME$"
  }
  # predefined column group patterns
  if (is.null(TechP)) {
    TechP = "CH[1-9]$"
  }
  if (is.null(MorphP)) {
    # create check for dt type STE.
    MorphP = "^RADIUS$|^ELLIPSE|^AREA$|^PERIMETER$|^CIRCULARITY$|^SOLIDITY$|^SPEED$|^DIRECTIONAL_CHANGE_RATE$"
  }
  if (is.null(MorphPpos)) {
    # create check for dt type STE.
    MorphPpos = "^RADIUS$|^ELLIPSE|^AREA$|^PERIMETER$|^CIRCULARITY$|^SOLIDITY$"
  }
  
  # Process exclusion patterns to contain standard excluded and given exclusions.
  excld.pattern <- excld.pattern.process(excld.pattern, StdP = StdP)
  
  predef = match.arg(predef, c("coord","technical", "morphological", "morphplus","clust", "nontechnical", "none"), several.ok = TRUE)
  
  # To initiate data coordinates (dtC), get just XY-coordinates per data type using specific colnmae as check 
  if (any(predef %in% "coord")) {
    dtC <- data.table::copy(dt)
    #dtC[, grep("_ID$|^FRAME$|POSITION.*[XYZT]", colnames(dtC), invert = T):=NULL]
    if ("SPOT_ID" %in% colnames(dt)) {
      dtC <- dtC[,c("TRACK_ID", "FRAME", "POSITION_X", 'POSITION_Y', 'POSITION_T', 'SPOT_ID')]
    }
    if ("TRACK_X_LOCATION" %in% colnames(dt)) {
      dtC <- dtC[,c("LABEL", "TRACK_X_LOCATION", "TRACK_Y_LOCATION", "TRACK_Z_LOCATION")]
    }
    if ("EDGE_ID" %in% colnames(dt)) {
      dtC <- dtC[,c("TRACK_ID", "EDGE_TIME", "EDGE_X_LOCATION", "EDGE_Y_LOCATION","EDGE_Z_LOCATION","EDGE_ID")]
    }
  } else dtC = NULL
  
  # Remove variables in dt based on keep the predefined sets 
  if (any(predef %in% "technical")) {
    dtTh <- data.table::copy(dt)
    dtTh[, grep(TechP, colnames(dtTh), invert = T):=NULL]
  } else dtTh = NULL
  if (any(predef %in% "morphological")) {
    dtMh <- data.table::copy(dt)
    dtMh[, grep(MorphPpos, colnames(dtMh), invert = T):=NULL]
  } else dtMh = NULL
  if (any(predef %in% "morphplus")) {
    dtMhp <- data.table::copy(dt)
    dtMhp[, grep(MorphP, colnames(dtMhp), invert = T):=NULL]
  } else dtMhp = NULL
  
  # Any but technical and in exclusion 
  if (any(predef %in% "clust")) {
    dtClst <- data.table::copy(dt)
    #negpattern = combine_patterns(c(StdP,TechP), logic = "OR")
    negpattern = combine_patterns(c(StdP, TechP, excld.pattern), logic = "OR")
    dtClst[, grep(negpattern, colnames(dtClst), invert = F):=NULL]
    # keep just numeric columns
    num_cols <- colnames(dtClst)[sapply(dtClst, is.numeric)]
    if (length(num_cols) < dim(dtClst)[2]) {
      dtClst[,!num_cols:=NULL]
    }
  } else dtClst = NULL
  
  # Any but technical
  if (any(predef %in% "nontechnical")) {
    dtNtch <- data.table::copy(dt)
    negpattern = combine_patterns(c(StdP,TechP), logic = "OR")
    dtNtch[, grep(negpattern, colnames(dtNtch), invert = F):=NULL]
  } else dtNtch = NULL
  
  # Use inclusion patterns 
  # There was breakage with use of is.null(excld.pattern) Hence match usage. 
  # "NO_EXCLUSION" is from excld.pattern.process parsing. 
  if (!is.null(incl.pattern)  & excld.pattern=="NO_EXCLUSION" ) {
    dtF <- data.table::copy(dt)
    inclCols = grep(incl.pattern, colnames(dt), value = T, invert = F)
    # remove non numeric columns
    if (numerics) {
      # get column names  of numerical columns
      num_cols <- colnames(dt)[sapply(dt, is.numeric)]
      inclCols = inclCols[inclCols %in% num_cols]
    }
    dtF[,(colnames(dtF)[colnames(dtF) %nin% inclCols]):=NULL]
    
  } else if (is.null(incl.pattern) & !excld.pattern=="NO_EXCLUSION") {
    dtF <- data.table::copy(dt)
    dtF[, grep(excld.pattern, colnames(dtF), invert = F):=NULL]
    if (numerics) {
      # get column names  of numerical columns
      num_cols <- colnames(dtF)[sapply(dtF, is.numeric)]
      num_cols = combine_patterns(c(num_cols), logic = "OR")
      dtF[, grep(num_cols, colnames(dtF), invert = T):=NULL]
    }
    
  } else if (!is.null(incl.pattern) & !excld.pattern=="NO_EXCLUSION") {
    dtF <- data.table::copy(dt)
    exclCols = grep(excld.pattern, colnames(dt), value = T, invert = F)
    inclCols = grep(incl.pattern, colnames(dt), value = T, invert = F)
    inclCols = inclCols[inclCols %nin% exclCols]
    # remove non-numeric columns
    if (numerics) {
      # get column names  of numerical columns
      num_cols <- colnames(dt)[sapply(dt, is.numeric)]
      inclCols = inclCols[inclCols %in% num_cols]
      dtF <- dtF[,..inclCols]
    }
    #inclCols = combine_patterns(c(inclCols), logic = "OR", exacts = T)
    #dtF[, grep(inclCols, colnames(dtF), invert = T):=NULL]
  } else dtF = NULL
  
  lst = list(dtF = dtF, dtC = dtC, dtTh = dtTh, dtMh = dtMh, dtMhp = dtMhp, dtClst = dtClst, dtNtch = dtNtch)
  lst = lst[!sapply(lst, is.null)]
  
  if (length(lst) == 1) {
    return(lst[[1]])
  }
  return(lst)
}
