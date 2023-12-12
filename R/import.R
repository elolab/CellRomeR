import_XML <- function(tmXML, dataset = "MigrDatTestXML",
                         experiment = "experiment", sample = "sample",
                         condition = "normal", replicate = NULL, MigrDatObj = NULL,
                         ...){
  if(!(file.exists(tmXML)))
    stop(sprintf("Couldn't find file %s", tmXML))
  ### Timing
  strt = Sys.time()
  
  ####
  #library("foreach")
  #library("doParallel")
  
  # if (deving==TRUE) {
  # Get default names and paths
  # x = initi_readTMx(type = "tmXML")
  # }
  
  #### auxiliary functions ####
  # Get feature names from retrieved XML feature list
  getFeatNames <- function(xmllist = XML::featsXML) {
    Featslist <- c("name",unlist(xmllist))
    names(Featslist) = NULL
    return(Featslist)
  }
  
  # Fetch data
  getTMvar <- function(TMxml,var) {
    sapply(TMxml, XML::xmlGetAttr, var)
  }
  
  
  
  ##### Function #####
  # Create "TRACK_ID" vector for edges table by using tracks data.table and spots number per track "NUMBER_SPOTS"
  track_ids4edges <- function(dtTxml, SE = c("E","SE","S")) {
    
    TRACK_IDs_E = c()
    #TRACK_IDs_SE = list()
    #TRACK_IDs_SE$Tracs4Edges = TRACK_IDs_E
    #TRACK_IDs_SE$Tracs4Spots = TRACK_IDs_S
    #TRACK_IDs_S = c()
    
    xvecE = c()
    # xvecS = c()
    
    SE = match.arg(SE, c("E","SE","S"))
    
    # Create TRACK_ID vetors for edges
    if (any(SE=="SE",SE=="E")) {
      for (i in seq_along(dtTxml[["LABEL"]])) {
        n = (dtTxml[["NUMBER_SPOTS"]][i]-1)
        TRACK_IDs_Ex = rep(dtTxml[["LABEL"]][i], n)
        TRACK_IDs_E = c(TRACK_IDs_E,TRACK_IDs_Ex)
        xvecE = c(xvecE,n)
      }
      return(TRACK_IDs_E)
      TRACK_IDs_SE$Tracs4Edges = TRACK_IDs_E
    }
    
    # Create TRACK_ID vetors for spots
    if (any(SE=="SE",SE=="S")) {
      print("Track_IDs to spots matches cannot be determined based on spot IDs!")}
    # if (any(SE=="SE",SE=="S")) {
    #   for (i in seq_along(dtTxml[["LABEL"]])) {
    #     n = (dtTxml[["NUMBER_SPOTS"]][i])
    #     TRACK_IDs_Sx = rep(x = dtTxml[["LABEL"]][i], n)
    #     TRACK_IDs_S = c(TRACK_IDs_S,TRACK_IDs_Sx)
    #     xvecS = c(xvecS,n)
    #   }
    #   TRACK_IDs_SE$Tracs4Spots = "Wrong order"# TRACK_IDs_S
    # }
    
    # Spot numbers used
    print("Number of edges per track:")
    print(xvecE)
    #print(xvecS)
    
    return(TRACK_IDs_SE)
  }
  
  if (is.null(MigrDatObj)) {
    MigrDatXML = MigrDat_newobj(dataset = dataset,
                                experiment = experiment, sample = sample,
                                data_type = "TMxml",
                                condition = condition, replicate = replicate,...)
    # MigrDatXML = MigrDat_newobj(dataset = dataset, Data_type = "TMxml",...)
  } else MigrDatXML = MigrDatObj
  
  MigrDatXML@metadata$condition = condition
  
  # implement file check tryCatch() here for non-XML files
  TMxml = XML::xmlParse(tmXML)
  
  if(length(TMxml) == 0) {
    # This probably will not work, as error will break anyway!
    stop("File is not a TrackMate XML file")
  }
  
  # Get variable names in TrackMate-XML
  sublistS <- XML::getNodeSet(TMxml,"//FeatureDeclarations//SpotFeatures//Feature/@feature")
  varsS = c("ID",getFeatNames(sublistS))
  varsS = varsS[-match(c("MANUAL_SPOT_COLOR"),varsS, nomatch = 0)] # remove "MANUAL_SPOT_COLOR"
  sublistT <- XML::getNodeSet(TMxml,"//FeatureDeclarations//TrackFeatures//Feature/@feature")
  varsT = getFeatNames(sublistT)
  sublistE <- XML::getNodeSet(TMxml,"//FeatureDeclarations//EdgeFeatures//Feature/@feature")
  varsE = getFeatNames(sublistE)
  varsE = varsE[-match(c("name","MANUAL_EDGE_COLOR"), varsE, nomatch = 0)] # remove "name" and "MANUAL_EDGE_COLOR"
  varsXML = list(varsSin = varsS, varsTin = varsT, varsEin = varsE)
  
  # could implement varsdf as with the tables
  varsdf = list(warning = "varsdf not yet implemented for TMxml import.",
                varsS = "varsdf not yet implemented for TMxml import.",
                varsT = "varsdf not yet implemented for TMxml import.",
                varsE = "varsdf not yet implemented for TMxml import.")
  
  ##### Retrieve units #####
  # units
  attrList <- c("spatialunits","timeunits")
  unitVec <- sapply(attrList, function(x) XML::xpathSApply(TMxml, "//Model", XML::xmlGetAttr, x))
  unitVec <- c(unitVec,"img width pixels","img height pixels","voxeldepth","timeinterval")
  attrList <- c("pixelwidth","timeinterval","width","height","voxeldepth","timeinterval")
  valueVec <- sapply(attrList, function(x) XML::xpathSApply(TMxml, "//ImageData", XML::xmlGetAttr, x))
  calibrationDF <- data.frame(value = as.numeric(valueVec),unit = unitVec)
  # <Model spatialunits="µm" timeunits="frame">
  # voxeldepth="1.0" timeinterval="1.0"
  MigrDatXML@metadata$units = calibrationDF
  
  # Store variable names into metadata
  # Could be better to add only after retrieving the data
  # If filtered etc. needs updating
  MigrDatXML@metadata$vars = varsXML
  MigrDatXML@metadata$varsdf = varsdf
  MigrDatXML@metadata$TMxml_file = as.character(tmXML)[1] # To prevent accidentally storing raw XML data
  
  # nspots could be used for checking import!
  # nspotsTMxml = xpathSApply(TMxmlSid, "//AllSpots", xmlGetAttr, "nspots")
  
  
  ## TrackMate track filtering in metadata could be used to filer tracks.
  ## However, as this is based SPOT_NUMBER < n, TRACK_DISPLACEMENT < 10 etc.
  ## see [link](https://imagej.net/plugins/trackmate/scripting)
  trackFltXML = XML::getNodeSet(TMxml, "//FilteredTracks//TrackID")
  TRACK_IDs <- getTMvar(trackFltXML, var = "TRACK_ID")
  MigrDatXML@metadata$track_filter = TRACK_IDs
  
  ## TrackMate Log (process later)
  TMxml2 = xml2::read_xml(tmXML)
  MigrDatXML@metadata$TMlog = xml2::xml_text(xml2::xml_find_all(TMxml2, ".//Log"))
  
  # These save xmlnodesets/nodes as strings can be catted
  #MigrDatXML@metadata$imageSettingsXML = as.character(xml_find_all(TMxml2, ".//Settings"))
  #MigrDatXML@metadata$imageMeta = as.character(xml_find_all(TMxml2, ".//Settings//ImageData"))
  #MigrDatXML@metadata$TrackerSettings = as.character(xml_find_all(TMxml2, ".//Settings//TrackerSettings"))
  
  
  
  
  ## TrackMate Image metadata (process later)
  imageSettingsXML = XML::getNodeSet(TMxml, "//Settings")
  MigrDatXML@metadata$imageSettingsXML = imageSettingsXML
  
  ## TrackMate Image metadata (process later)
  #imageMetaXML = getNodeSet(TMxml, "//Settings//ImageData")
  #MigrDatXML@metadata$imageMeta = imageMetaXML
  
  ## Tracker metadata (process later)
  #imageMetaXML = getNodeSet(TMxml, "//Settings//TrackerSettings")
  #MigrDatXML@metadata$TrackerSettings = imageMetaXML
  
  # Some differences vs. imported from TrackMate-export tables
  # TMtableS and E has Track_ID to indicate in which track the spot/edge is used! Important for plotting functions!
  # Need to implement
  # - add "Track_ID" for spots and edges
  # - Change name to LABEL
  
  
  ##### Add tracks data  #####
  # Add tracks data
  #
  #
  trackXML = XML::getNodeSet(TMxml, "//Track")
  nms <- getTMvar(trackXML, var = "name")
  dtTxml <- data.table::data.table(matrix(nrow = length(nms), ncol = length(varsXML$varsTin)))
  names(dtTxml) = varsXML$varsTin
  dtTxml[,name:=nms]
  names(dtTxml)[1] = "LABEL"
  
  for (var in varsXML$varsTin[-1]) {
    dat = getTMvar(trackXML,var = var)
    dtTxml[,(var):= dat ]
  }
  
  # Note: TRACK_ID is same as TRACK_INDEX, while it would likely be better to correspond to number in LABEL (e.g. Track_1207)
  # It remains to be seen what the convention is.
  # The ids and indices are problematic to match, as can be somtimes be interpreted as indices.
  dtTxml <- standardize_dt(dtTxml, type = "T")
  
  
  MigrDatXML@tracks$raw = dtTxml
  MigrDatXML@IDs$tracks = list(LABELs = dtTxml[["LABEL"]], TRACK_IDs = dtTxml[["TRACK_ID"]]) # Need to think this still
  MigrDatXML@metadata$vars$varsT <- colnames(MigrDatXML@tracks$raw)
  
  
  ##### Add edges data  #####
  # Add edges data
  #
  #
  edgeXML = XML::getNodeSet(TMxml, "//Track//Edge")
  dat1 <- getTMvar(edgeXML, var = varsXML$varsEin[1])
  dat2 <- getTMvar(edgeXML, var = varsXML$varsEin[2])
  dtExml <- data.table::data.table(matrix(nrow = length(dat1), ncol = (length(varsXML[[3]])+2) ))
  names(dtExml) = c("LABEL","EDGE_ID", varsXML[[3]])
  Labels = paste0("ID",dat1,"_","ID",dat2)
  dtExml[,LABEL:=Labels]
  Eids = paste0(dat1,"_",dat2)
  dtExml[,EDGE_ID:=Eids]
  
  for (var in varsXML$varsEin) {
    dat = getTMvar(edgeXML,var = var)
    dtExml[,(var):= dat ]
  }
  
  # Create TRACK_ID vectors for spots and edges tables
  # This is actually the LABEL vector how is it in imported tables?
  Tracs4Edges = track_ids4edges(dtTxml, SE = "E")
  
  # Add Track_ID data
  dtExml[,TRACK_ID:=Tracs4Edges]
  
  # The ids and indices are problematic to match, as can be sometimes be interpreted as indices.
  # Should they be kept as characters? instead of numbers?
  # Hence use of LABEL is more secure way!
  # The LABEL is the actual TRACK_ID!
  # dtExml <- type.convert(dtExml, as.is = TRUE) now in standardize_dt() function
  # Also the tables are now sorted in standardize function
  dtExml <- standardize_dt(dtExml, type = "E")
  
  # Add to MigrDat object
  MigrDatXML@edges$raw = dtExml
  MigrDatXML@IDs$edges = list(LABELs = dtExml[["LABEL"]], EDGE_IDs = dtExml[["EDGE_ID"]]) # Need to think this still
  # MigrDatXML@IDs$edges = paste(dataset, dtTxml[["LABEL"]], dtTxml[["TRACK_ID"]], sep = "_") # Need to think this still
  MigrDatXML@metadata$vars$varsE <- colnames(MigrDatXML@edges$raw)
  
  
  
  ##### Add spots data  #####
  # Add spots data
  #
  #
  spotsXML <- XML::getNodeSet(TMxml,"//AllSpots//SpotsInFrame//Spot")
  IDs <- getTMvar(spotsXML,var = "ID")
  dtSxml <- data.table::data.table(matrix(nrow = length(IDs), ncol = length(varsXML$varsSin)))
  names(dtSxml) = varsXML$varsSin
  dtSxml[,ID:=IDs]
  
  for (var in varsXML$varsSin[-1]) {
    dat = getTMvar(spotsXML,var = var)
    dtSxml[,(var):= dat ]
  }
  
  # Get TRACK_IDs for the spots
  # Spots not being used in Tracks are labeled as ninTracks
  track_ids = track_ids4spotsM(dtExml,dtSxml)
  spots_nin_tracks = dtSxml[["ID"]][(track_ids=="ninTracks") ]
  
  # add to spots data.table
  dtSxml[,TRACK_ID:=track_ids]
  
  # Conversions now in standardize_dt() function
  dtSxml <- standardize_dt(dtSxml, type = "S")
  
  MigrDatXML@spots$raw = dtSxml
  MigrDatXML@IDs$spots = list(LABELs = dtSxml[["LABEL"]], SPOT_IDs = dtSxml[["SPOT_ID"]])
  MigrDatXML@metadata$vars$varsS <- colnames(MigrDatXML@spots$raw)
  
  ##### ROIs #####
  # Add ROI_points data
  #
  #
  #
  # List of ROI_points named with spot LABELs
  MigrDatXML@roi_points$raw <- getROIpoints_TMxml(TMxml)
  MigrDatXML@IDs$roi_points = names(MigrDatXML@roi_points$raw)
  
  # Add roi_point point number to spots data
  MigrDatXML@spots$raw[["roi_length"]] <- sapply(MigrDatXML@roi_points$raw, length)
  
  #update IDs list (not yet used anywhere)
  MigrDatXML@metadata$vars$varsS <- colnames(MigrDatXML@spots$raw)
  
  # Do checks
  MigrDatXML@IDs$spots$LABELs
  MigrDatXML@IDs$tracks$LABELs
  MigrDatXML@IDs$edges$LABELs
  
  MigrDatXML@IDs$spots$SPOT_IDs
  MigrDatXML@IDs$tracks$TRACK_IDs
  MigrDatXML@IDs$edges$EDGE_IDs
  
  # Add SPEED, DISPLACEMENT, and DIRECTIONAL_CHANGE_RATE to spots
  # The values from second spots are used for the first. NaN and NAs are replaced with 0.
  MigrDatXML <- speed2spots.raw(MigrDatXML, NAs0 = TRUE)
  
  # Filtering
  # Timing
  timdiff = (Sys.time() - strt)
  MigrDatXML@metadata$importduration = timdiff
  # initialize metadata
  MigrDatXML <- init.metadata(MigrDatXML)
  MigrDatXML <- init.clustering(MigrDatXML)
  
  
  #### Return MigrDat object ####
  return(MigrDatXML)
  
}
