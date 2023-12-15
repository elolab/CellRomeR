#' @include zzzz.R
#' @include ILoRegMigr_utils.R
#' @include auxiliary.R
#' @include ILoRegMigr_utils.R
#' @include ILoRegMigr_main.R
#'
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Class definitions ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# define classes here

#### Cell migration data class ####
#' #### Cell migration data class ####
#'
#'
#'
#' @export
MigrDat <- setClass(
  "MigrDat",
  slots = c(
    spots = "list",
    tracks = "list",
    edges = "list",
    roi_points  = "list",
    metadata = "list",
    dimreductions = "list",
    clustering = "list",
    IDs = "list",
    misc = "ANY"
  ),

  prototype = list(
    spots = list(raw = NULL, normalized = NULL, scaled = NULL ),
    tracks = list(raw = NULL, normalized = NULL, scaled = NULL ),
    edges = list(raw = NULL, normalized = NULL, scaled = NULL ),
    roi_points  = list(raw = NULL, coord = NULL, scaled = NULL),
    metadata = list(experiment = NULL, sample = NULL, replicate = NULL, condition = NULL),
    dimreductions = list("S" = list(PCA = NULL, UMAP = NULL, UMAP3D = NULL),
                         "T" = list(PCA = NULL, UMAP = NULL, UMAP3D = NULL),
                          "E" = list(PCA = NULL, UMAP = NULL, UMAP3D = NULL)),
    clustering = list("S" = NULL,
                      "T" = NULL,
                      "E" = NULL),
    IDs = list(spots = NULL, tracks = NULL, edges = NULL),
    misc = list()
  )
)

#### New MigrDat object ####
#'
#'
#'
#' Create new MigrDat object
#' @export

MigrDat_newobj <- function(dataset = "data", experiment = "Exp_1", sample = "Sample1",
                           data_type = c("TMxml","TMtable","MigrTable"),
                           condition = NULL, replicate = NULL
                           ) {
  # Check arguments
  stopifnot(!missing(dataset), !missing(data_type))
  data_type <- match.arg(data_type, c("TMxml","TMtable","MigrTable"))

  # Create object
  dobj = new("MigrDat")
  #dobj@spots = list(raw = NULL, filtered = NULL, normalized = NULL )
  #dobj@tracks = list(raw = NULL, filtered = NULL, normalized = NULL )
  #dobj@edges = list(raw = NULL, filtered = NULL, normalized = NULL )
  #dobj@roi_points = list(raw = NULL, coord = NULL)
  dobj@metadata[["dataset"]] = list(dataset = dataset, experiment = experiment,
                                    sample = sample, type = data_type,
                                    condition = condition, replicate = replicate)
  if (!is.null(condition)) {
    dobj@metadata$condition = condition
  }
  if (!is.null(replicate)) {
    dobj@metadata$replicate = replicate
  }

  StdP = "^LABEL|^label|ID$|ID[1-9]$|id$|INDEX|^TIME|QUALITY|LOCATION|START|STOP|GAP$|^NUMBER|DURATION|^POSITION|^FRAME|VISIBILITY|LINK_COST|^EDGE_TIME$"
  TechP = "CH[1-9]$"
  MorphP = "^RADIUS$|^ELLIPSE|^AREA$|^PERIMETER$|^CIRCULARITY$|^SOLIDITY$|^SPEED$|^DIRECTIONAL_CHANGE_RATE$"
  patterns = list(standard = StdP, technical = TechP, morphological = MorphP)

  # Store standard pattern to metadata
  dobj@metadata$Std.pattern <- StdP
  dobj@metadata$patterns <- patterns

  dobj@metadata$data.in.use = c("raw") #
  #dobj@IDs = list(spots = "NULL", tracks = "NULL", edges = "NULL")
  return(dobj)
}

#### methods for MigrDat object  ####
#' ## Limit output of the MigrDat object ####
#' Add more meaningful info later
#'
#' @export
setMethod("show",
          "MigrDat",
          function(object){
            cat("An object of class", class(x = object), "with:\n")
            cat("data of cell experiment", (object@metadata$name), "of condition ", (object@metadata$condition), "\n")
            cat("imported from", (object@metadata$TMxml_file), "\n\n")
            cat("Has following slots available: \n", slotNames(object), "\n\n")

            slotz = slots.in.use(object)

            cat("Spots available as:", slotz$in.spots,"\n"  )
            cat("Tracks available as:", slotz$in.tracks,"\n"  )
            cat("Edges available as:", slotz$in.edges,"\n"  )
            cat("Roi_points available as:", slotz$in.roi_points,"\n"  )


            cat("\nDefault PCAs are:", default.PCAs(object)[[1]], default.PCAs(object)[[2]], default.PCAs(object)[[3]] )
            cat("\nDefault UMAPs are:", default.UMAPs(object)[[1]],default.UMAPs(object)[[2]],default.UMAPs(object)[[3]] )
            cat("\nDefault UMAP3Ds are:", default.UMAP3Ds(object)[[1]], default.UMAP3Ds(object)[[2]], default.UMAP3Ds(object)[[3]],"\n" )

            cat("\nFollowing",length(names(dimreducts(object))), "dimensional reductions are available:\n", names(dimreducts(object)),"\n" )
            cat("\n",length( names(clusterings(object)) ), "clustering sets with following unique identifiers are available:\n", names(clusterings(object)),"\n\n" )

            cat("\nDefault kmeans clustering is:", default.kmeans(object)[[1]],"with", default.kmeans(object)[[2]])
            cat("\nDefault hclust clustering is:", default.hclusts(object) [[1]],"with", default.hclusts(object)[[2]])
            cat("\nDefault ILoReg clustering is:", default.ILoRegclusters(object)[[1]],"with", default.ILoRegclusters(object)[[2]] )
            #cat("\nDefault Seurat clustering is:", default.Seurclusters(object)[[1]],"with", default.Seurclusters(object)[[2]], "\n\n" )

            cat("The MigrDat object contains:\n")
            cat(nrow(object@tracks$raw), "tracks with ", ncol(object@tracks$raw),  "variables.\n")
            cat(nrow(object@edges$raw), "edges with ", ncol(object@edges$raw),  "variables.\n")
            cat(nrow(object@spots$raw), "spots with ", ncol(object@spots$raw),  "variables.\n\n")


          })

#### Assessors and getters ####


#' ## Std.pattern in MigrDat object ####
#' Fetch pattern that can be used with column names to grep standard columns
#'
#' @export
setGeneric("Std.pattern", function(x) standardGeneric("Std.pattern"))
setGeneric("Std.pattern<-", function(x, value) standardGeneric("Std.pattern<-"))

setMethod("Std.pattern", "MigrDat", function(x) {
  return(x@metadata$Std.pattern)
})
setMethod("Std.pattern<-", "MigrDat", function(x, value) {
  x@metadata$Std.pattern <- value
  x
})




#' #### Get variable names ####
#' @examples
#'
#' MigrObj.vars.raw(MigrDat)
#'
#' @export
setGeneric("MigrObj.vars.raw", function(object) {standardGeneric("MigrObj.vars.raw") })
setMethod("MigrObj.vars.raw", "MigrDat", function(object) {

  vars <- sapply(list(spots.raw(object),
                      tracks.raw(object),
                      edges.raw(object)), colnames)
  names(vars) <- c("spot.vars", "track.vars", "edge.vars")
  return(vars)

})

#' #### Get variable names ####
#' @examples
#'
#' MigrObj.vars.normalized(MigrDat)
#'
#' @export
setGeneric("MigrObj.vars.normalized", function(object) {standardGeneric("MigrObj.vars.normalized") })
setMethod("MigrObj.vars.normalized", "MigrDat", function(object) {

  vars <- sapply(list(spots.normalized(object),
                      tracks.normalized(object),
                      edges.normalized(object)), colnames)
  names(vars) <- c("spot.vars", "track.vars", "edge.vars")
  print(vars)

})


#' #### Get variable names ####
#' @examples
#'
#' MigrObj.vars.scaled(MigrDat)
#'
#' @export
setGeneric("MigrObj.vars.scaled", function(object) {standardGeneric("MigrObj.vars.scaled") })
setMethod("MigrObj.vars.scaled", "MigrDat", function(object) {

  vars <- sapply(list(spots.scaled(object),
                      tracks.scaled(object),
                      edges.scaled(object)), colnames)
  names(vars) <- c("spot.vars", "track.vars", "edge.vars")
  return(vars)

})

##### Data getters #####

#' ## Get and store whole raw data from MigrDat object ####
#' @examples
#'
#' MigrObj.raw(MigrDat)
#'
#' @export
setGeneric("MigrObj.raw", function(object) {standardGeneric("MigrObj.raw") })
setMethod("MigrObj.raw", "MigrDat", function(object) {
  print(object@spots$raw)
  print(object@tracks$raw)
  print(object@edges$raw)
  })


#' ## Get and store whole normalized data from MigrDat object ####
#' @examples
#'
#' MigrObj.normalized(MigrDat)
#'
#' @export
setGeneric("MigrObj.normalized", function(object) {standardGeneric("MigrObj.normalized") })
setMethod("MigrObj.normalized", "MigrDat", function(object) {
  print(object@spots$normalized)
  print(object@tracks$normalized)
  print(object@edges$normalized)

} )


#' #### Get and store whole spots set from MigrDat object ####
#' @examples
#'
#' spotsLS(MigrDat_XMLtest)
#'
#' @export
setGeneric("spotsLS", function(x) standardGeneric("spotsLS"))
setGeneric("spotsLS<-", function(x, value) standardGeneric("spotsLS<-"))

setMethod("spotsLS", "MigrDat", function(x) x@spots)
setMethod("spotsLS<-", "MigrDat", function(x, value) {
  x@spots <- value
  x
})


#' #### Get and store raw spots data from MigrDat object ####
#' @examples
#'
#' spots.raw(MigrDat_XMLtest)
#' spots.raw(MigrDat_XMLtest) <- new.data
#'
#' @export
setGeneric("spots.raw", function(x) standardGeneric("spots.raw"))
setGeneric("spots.raw<-", function(x, value) standardGeneric("spots.raw<-"))

setMethod("spots.raw", "MigrDat", function(x) x@spots$raw)
setMethod("spots.raw<-", "MigrDat", function(x, value) {
  x@spots$raw <- value
  x
})


#' #### Get and store normalized spots data from MigrDat object ####
#'
#'
#' @examples
#' spots.normalized(MigrDat_XMLtest)
#' spots.normalized(MigrDat_XMLtest) <- normalize(spots.raw(MigrDat_XMLtest))
#'
#' @export
setGeneric("spots.normalized", function(x) standardGeneric("spots.normalized"))
setGeneric("spots.normalized<-", function(x, value) standardGeneric("spots.normalized<-"))

setMethod("spots.normalized", "MigrDat", function(x) x@spots$normalized)
setMethod("spots.normalized<-", "MigrDat", function(x, value) {
  x@spots$normalized <- value
  x
})



#' #### Get and store scaled spots data from MigrDat object ####
#' @examples
#'
#' spots.scaled(MigrDat_XMLtest)
#'
#' @export
setGeneric("spots.scaled", function(x) standardGeneric("spots.scaled"))
setGeneric("spots.scaled<-", function(x, value) standardGeneric("spots.scaled<-"))

setMethod("spots.scaled", "MigrDat", function(x) x@spots$scaled)
setMethod("spots.scaled<-", "MigrDat", function(x, value) {
  x@spots$scaled <- value
  x
})


#' #### Get and store whole tracks set from MigrDat object ####
#' @examples
#'
#' tracksLS(MigrDat_XMLtest)
#'
#' @export
setGeneric("tracksLS", function(x) standardGeneric("tracksLS"))
setGeneric("tracksLS<-", function(x, value) standardGeneric("tracksLS<-"))

setMethod("tracksLS", "MigrDat", function(x) x@tracks)
setMethod("tracksLS<-", "MigrDat", function(x, value) {
  x@tracks <- value
  x
})


#' #### Get and store raw tracks data from MigrDat object ####
#' @examples
#'
#' tracks.raw(MigrDat_XMLtest)
#'
#' @export
setGeneric("tracks.raw", function(x) standardGeneric("tracks.raw"))
setGeneric("tracks.raw<-", function(x, value) standardGeneric("tracks.raw<-"))

setMethod("tracks.raw", "MigrDat", function(x) x@tracks$raw)
setMethod("tracks.raw<-", "MigrDat", function(x, value) {
  x@tracks$raw <- value
  x
})


#' #### Get and store normalized tracks data from MigrDat object ####
#' @examples
#'
#' tracks.normalized(MigrDat_XMLtest)
#' tracks.normalized(MigrDat_XMLtest) <- normalize(tracks.raw(MigrDat_XMLtest))
#'
#' @export
setGeneric("tracks.normalized", function(x) standardGeneric("tracks.normalized"))
setGeneric("tracks.normalized<-", function(x, value) standardGeneric("tracks.normalized<-"))

setMethod("tracks.normalized", "MigrDat", function(x) x@tracks$normalized)
setMethod("tracks.normalized<-", "MigrDat", function(x, value) {
  x@tracks$normalized <- value
  x
})


#' #### Get and store scaled tracks data from MigrDat object ####
#' @examples
#'
#' tracks.scaled(MigrDat_XMLtest)
#' tracks.scaled(MigrDat_XMLtest) <- scale(tracks.normalized(MigrDat_XMLtest))
#'
#' @export
setGeneric("tracks.scaled", function(x) standardGeneric("tracks.scaled"))
setGeneric("tracks.scaled<-", function(x, value) standardGeneric("tracks.scaled<-"))

setMethod("tracks.scaled", "MigrDat", function(x) x@tracks$scaled)
setMethod("tracks.scaled<-", "MigrDat", function(x, value) {
  x@tracks$scaled <- value
  x
})

#' #### Get and store whole edges set from MigrDat object ####
#' @examples
#'
#' edgesLS(MigrDat_XMLtest)
#'
#' @export
setGeneric("edgesLS", function(x) standardGeneric("edgesLS"))
setGeneric("edgesLS<-", function(x, value) standardGeneric("edgesLS<-"))

setMethod("edgesLS", "MigrDat", function(x) x@edges)
setMethod("edgesLS<-", "MigrDat", function(x, value) {
  x@edges <- value
  x
})



#' #### Get and store raw edges data from MigrDat object ####
#' @examples
#'
#' edges.raw(MigrDat_XMLtest)
#'
#' @export
setGeneric("edges.raw", function(x) standardGeneric("edges.raw"))
setGeneric("edges.raw<-", function(x, value) standardGeneric("edges.raw<-"))

setMethod("edges.raw", "MigrDat", function(x) x@edges$raw)
setMethod("edges.raw<-", "MigrDat", function(x, value) {
  x@edges$raw <- value
  x
})

#' #### Get and store normalized edges data from MigrDat object ####
#' @examples
#'
#' edges.normalized(MigrDat_XMLtest)
#'
#' @export
setGeneric("edges.normalized", function(x) standardGeneric("edges.normalized"))
setGeneric("edges.normalized<-", function(x, value) standardGeneric("edges.normalized<-"))

setMethod("edges.normalized", "MigrDat", function(x) x@edges$normalized)
setMethod("edges.normalized<-", "MigrDat", function(x, value) {
  x@edges$normalized <- value
  x
})

#' #### Get and store scaled edges data from MigrDat object ####
#' @examples
#'
#' edges.scaled(MigrDat_XMLtest)
#'
#' @export
setGeneric("edges.scaled", function(x) standardGeneric("edges.scaled"))
setGeneric("edges.scaled<-", function(x, value) standardGeneric("edges.scaled<-"))

setMethod("edges.scaled", "MigrDat", function(x) x@edges$scaled)
setMethod("edges.scaled<-", "MigrDat", function(x, value) {
  x@edges$scaled <- value
  x
})



#' #### Get and store whole roi_points set from MigrDat object ####
#' @examples
#'
#' roi_pointsLS(MigrDat_XMLtest)
#'
#' @export
setGeneric("roi_pointsLS", function(x) standardGeneric("roi_pointsLS"))
setGeneric("roi_pointsLS<-", function(x, value) standardGeneric("roi_pointsLS<-"))

setMethod("roi_pointsLS", "MigrDat", function(x) x@roi_points)
setMethod("roi_pointsLS<-", "MigrDat", function(x, value) {
  x@roi_points <- value
  x
})



#' #### Get and store raw roi_points data from MigrDat object ####
#' @examples
#'
#' roi_points.raw(MigrDat_XMLtest)[1:10]
#'
#' @export
setGeneric("roi_points.raw", function(x) standardGeneric("roi_points.raw"))
setGeneric("roi_points.raw<-", function(x, value) standardGeneric("roi_points.raw<-"))

setMethod("roi_points.raw", "MigrDat", function(x) x@roi_points$raw)
setMethod("roi_points.raw<-", "MigrDat", function(x, value) {
  x@roi_points$raw <- value
  x
})



#' #### Get and store normalized roi_points data from MigrDat object ####
#' @examples
#'
#' roi_points.normalized(MigrDat_XMLtest)[1:10]
#'
#' @export
setGeneric("roi_points.normalized", function(x) standardGeneric("roi_points.normalized"))
setGeneric("roi_points.normalized<-", function(x, value) standardGeneric("roi_points.normalized<-"))

setMethod("roi_points.normalized", "MigrDat", function(x) x@roi_points$normalized)
setMethod("roi_points.normalized<-", "MigrDat", function(x, value) {
  x@roi_points$normalized <- value
  x
})


#' #### Get and store scaled roi_points data from MigrDat object ####
#' @examples
#'
#' roi_points.scaled(MigrDat_XMLtest)[1:10]
#'
#' @export
setGeneric("roi_points.scaled", function(x) standardGeneric("roi_points.scaled"))
setGeneric("roi_points.scaled<-", function(x, value) standardGeneric("roi_points.scaled<-"))

setMethod("roi_points.scaled", "MigrDat", function(x) x@roi_points$scaled)
setMethod("roi_points.scaled<-", "MigrDat", function(x, value) {
  x@roi_points$scaled <- value
  x
})


#### Metadata classes ###

#' ## Get and store raw spots metadata in MigrDat object ####
#' @examples
#'
#' spots.meta(MigrDat_XMLtest)
#' spots.meta(MigrDat_XMLtest) <- new.dt
#'
#' @export
setGeneric("spots.meta", function(x) standardGeneric("spots.meta"))
setGeneric("spots.meta<-", function(x, value) standardGeneric("spots.meta<-"))

setMethod("spots.meta", "MigrDat", function(x) x@metadata[["S"]])
setMethod("spots.meta<-", "MigrDat", function(x, value) {
  x@metadata[["S"]] <- value
  x
})

#' ## Get and store raw edges metadata in MigrDat object ####
#' @examples
#'
#' edges.meta(MigrDat_XMLtest)
#' edges.meta(MigrDat_XMLtest) <- new.dt
#'
#' @export
setGeneric("edges.meta", function(x) standardGeneric("edges.meta"))
setGeneric("edges.meta<-", function(x, value) standardGeneric("edges.meta<-"))

setMethod("edges.meta", "MigrDat", function(x) x@metadata[["T"]])
setMethod("edges.meta<-", "MigrDat", function(x, value) {
  x@metadata[["T"]] <- value
  x
})

#' ## Get and store raw tracks metadata in MigrDat object ####
#' @examples
#'
#' tracks.meta(MigrDat)
#' tracks.meta(MigrDat) <- new.dt
#'
#' @export
setGeneric("tracks.meta", function(x) standardGeneric("tracks.meta"))
setGeneric("tracks.meta<-", function(x, value) standardGeneric("tracks.meta<-"))

setMethod("tracks.meta", "MigrDat", function(x) x@metadata[["E"]])
setMethod("tracks.meta<-", "MigrDat", function(x, value) {
  x@metadata[["E"]] <- value
  x
})



#### Methods dimensional reductions ####

#' ## Get and store whole DR sets from MigrDat object ####
#' @examples
#'
#' dimreducts(MigrDat_XMLtest)
#'
#' @export
setGeneric("dimreducts", function(x) standardGeneric("dimreducts"))
setGeneric("dimreducts<-", function(x, value) standardGeneric("dimreducts<-"))

setMethod("dimreducts", "MigrDat", function(x) x@dimreductions)
setMethod("dimreducts<-", "MigrDat", function(x, value) {
  x@dimreductions <- value
  x
})

#' ## Get list of names of the DR in the MigrDat object ####
#' @examples
#'
#' dimreducts.names(MigrDat)
#'
#' @export
setGeneric("dimreducts.names", function(x) standardGeneric("dimreducts.names"))
setMethod("dimreducts.names", "MigrDat", function(x) sapply(x@dimreductions, names))


#' ## Get and store whole DRs data sets from MigrDat object ####
#' @examples
#'
#' default.DRs(MigrDat)
#'
#' @export
setGeneric("default.DRs", function(x) standardGeneric("default.DRs"))
setMethod("default.DRs", "MigrDat", function(x)
  list("S" = x@dimreductions[["S"]][c("PCA","UMAP","UMAP3D")],
       "T" = x@dimreductions[["T"]][c("PCA","UMAP","UMAP3D")],
       "E" = x@dimreductions[["E"]][c("PCA","UMAP","UMAP3D")]) )

# Does not work
# setGeneric("default.DRs<-", function(x, value) standardGeneric("default.DRs<-"))
# setMethod("default.DRs<-", "MigrDat", function(x, value) {
#   x@dimreductions[["S"]] <- value[[1]]
#   x@dimreductions[["T"]] <- value[[2]]
#   x@dimreductions[["E"]] <- value[[3]] })



#' ## Get and store default DRs data set names  ####
#' @examples
#'
#' default.DRs.spots(MigrDat)
#'
#' @export
setGeneric("default.DRs.spots", function(x) standardGeneric("default.DRs.spots"))
setGeneric("default.DRs.spots<-", function(x, value) standardGeneric("default.DRs.spots<-"))

setMethod("default.DRs.spots", "MigrDat", function(x) x@dimreductions[["S"]][c("PCA","UMAP","UMAP3D")])
setMethod("default.DRs.spots<-", "MigrDat", function(x, value) {
  x@dimreductions[["S"]][c("PCA","UMAP","UMAP3D")] <- value
  x
})

#' ## Get and store default DRs data set names  ####
#' @examples
#'
#' default.DRs.tracks(MigrDat)
#'
#' @export
setGeneric("default.DRs.tracks", function(x) standardGeneric("default.DRs.tracks"))
setGeneric("default.DRs.tracks<-", function(x, value) standardGeneric("default.DRs.tracks<-"))

setMethod("default.DRs.tracks", "MigrDat", function(x) x@dimreductions[["T"]][c("PCA","UMAP","UMAP3D")])
setMethod("default.DRs.tracks<-", "MigrDat", function(x, value) {
  x@dimreductions[["T"]][c("PCA","UMAP","UMAP3D")] <- value
  x
})


#' ## Get and store default DRs data set names  ####
#' @examples
#'
#' default.DRs.edges(MigrDat)
#'
#' @export
setGeneric("default.DRs.edges", function(x) standardGeneric("default.DRs.edges"))
setGeneric("default.DRs.edges<-", function(x, value) standardGeneric("default.DRs.edges<-"))

setMethod("default.DRs.edges", "MigrDat", function(x) x@dimreductions[["E"]][c("PCA","UMAP","UMAP3D")])
setMethod("default.DRs.edges<-", "MigrDat", function(x, value) {
  x@dimreductions[["E"]][c("PCA","UMAP","UMAP3D")] <- value
  x
})

#' ## Name of the default PCA in MigrDat object ####
#' @examples
#'
#' default.PCAs(MigrDat)
#'
#' @export
setGeneric("default.PCAs", function(x) standardGeneric("default.PCAs"))
setMethod("default.PCAs", "MigrDat", function(x)
          list(spots = default.PCA.spots(x),
               tracks = default.PCA.tracks(x),
               edges = default.PCA.edges(x)) )

# setGeneric("default.PCA<-", function(x, value) standardGeneric("default.PCA<-"))
# setMethod("default.PCA<-", "MigrDat", function(x, value) {
#   x@dimreductions$PCA <- value
#   x
# })



#' ## Name of the default UMAP in MigrDat object ####
#' @examples
#'
#' default.UMAP(MigrDat)
#'
#' @export
setGeneric("default.UMAPs", function(x) standardGeneric("default.UMAPs"))
setMethod("default.UMAPs", "MigrDat", function(x)
  list(spots = default.UMAP.spots(x),
       tracks = default.UMAP.tracks(x),
       edges = default.UMAP.edges(x)) )

# setGeneric("default.UMAP<-", function(x, value) standardGeneric("default.UMAP<-"))
# setMethod("default.UMAP<-", "MigrDat", function(x, value) {
#   x@dimreductions$UMAP <- value
#   x
# })



#' ## Name of the default UMAP3D in MigrDat object ####
#' @examples
#'
#' default.UMAP3D(MigrDat)
#'
#' @export
setGeneric("default.UMAP3Ds", function(x) standardGeneric("default.UMAP3Ds"))
setMethod("default.UMAP3Ds", "MigrDat", function(x)
  list(spots = default.UMAP3D.spots(x),
       tracks = default.UMAP3D.tracks(x),
       edges = default.UMAP3D.edges(x)))

# setGeneric("default.UMAP3D<-", function(x, value) standardGeneric("default.UMAP3D<-"))
# setMethod("default.UMAP3D<-", "MigrDat", function(x, value) {
#   x@dimreductions$UMAP3D <- value
#   x
# })


#' ## Name of the default PCA in MigrDat object ####
#' @examples
#'
#' default.PCA(MigrDat)
#'
#' @export
setGeneric("default.PCA.spots", function(x) standardGeneric("default.PCA.spots"))
setGeneric("default.PCA.spots<-", function(x, value) standardGeneric("default.PCA.spots<-"))

setMethod("default.PCA.spots", "MigrDat", function(x) x@dimreductions[["S"]]$PCA)
setMethod("default.PCA.spots<-", "MigrDat", function(x, value) {
  x@dimreductions[["S"]]$PCA <- value
  x
})



#' ## Name of the default UMAP in MigrDat object ####
#' @examples
#'
#' default.UMAP(MigrDat)
#'
#' @export
setGeneric("default.UMAP.spots", function(x) standardGeneric("default.UMAP.spots"))
setGeneric("default.UMAP.spots<-", function(x, value) standardGeneric("default.UMAP.spots<-"))

setMethod("default.UMAP.spots", "MigrDat", function(x) x@dimreductions[["S"]]$UMAP)
setMethod("default.UMAP.spots<-", "MigrDat", function(x, value) {
  x@dimreductions[["S"]]$UMAP <- value
  x
})



#' ## Name of the default UMAP3D in MigrDat object ####
#' @examples
#'
#' default.UMAP3D(MigrDat)
#'
#' @export
setGeneric("default.UMAP3D.spots", function(x) standardGeneric("default.UMAP3D.spots"))
setGeneric("default.UMAP3D.spots<-", function(x, value) standardGeneric("default.UMAP3D.spots<-"))

setMethod("default.UMAP3D.spots", "MigrDat", function(x) x@dimreductions[["S"]]$UMAP3D)
setMethod("default.UMAP3D.spots<-", "MigrDat", function(x, value) {
  x@dimreductions[["S"]]$UMAP3D <- value
  x
})


#' ## Name of the default PCA in MigrDat object ####
#' @examples
#'
#' default.PCA.tracks(MigrDat)
#'
#' @export
setGeneric("default.PCA.tracks", function(x) standardGeneric("default.PCA.tracks"))
setGeneric("default.PCA.tracks<-", function(x, value) standardGeneric("default.PCA.tracks<-"))

setMethod("default.PCA.tracks", "MigrDat", function(x) x@dimreductions[["T"]]$PCA)
setMethod("default.PCA.tracks<-", "MigrDat", function(x, value) {
  x@dimreductions[["T"]]$PCA <- value
  x
})



#' ## Name of the default UMAP in MigrDat object ####
#' @examples
#'
#' default.UMAP.tracks(MigrDat)
#'
#' @export
setGeneric("default.UMAP.tracks", function(x) standardGeneric("default.UMAP.tracks"))
setGeneric("default.UMAP.tracks<-", function(x, value) standardGeneric("default.UMAP.tracks<-"))

setMethod("default.UMAP.tracks", "MigrDat", function(x) x@dimreductions[["T"]]$UMAP)
setMethod("default.UMAP.tracks<-", "MigrDat", function(x, value) {
  x@dimreductions[["T"]]$UMAP <- value
  x
})



#' ## Name of the default UMAP3D in MigrDat object ####
#' @examples
#'
#' default.UMAP3D.tracks(MigrDat)
#'
#' @export
setGeneric("default.UMAP3D.tracks", function(x) standardGeneric("default.UMAP3D.tracks"))
setGeneric("default.UMAP3D.tracks<-", function(x, value) standardGeneric("default.UMAP3D.tracks<-"))

setMethod("default.UMAP3D.tracks", "MigrDat", function(x) x@dimreductions[["T"]]$UMAP3D)
setMethod("default.UMAP3D.tracks<-", "MigrDat", function(x, value) {
  x@dimreductions[["T"]]$UMAP3D <- value
  x
})


#' ## Name of the default PCA in MigrDat object ####
#' @examples
#'
#' default.PCA.edges(MigrDat)
#'
#' @export
setGeneric("default.PCA.edges", function(x) standardGeneric("default.PCA.edges"))
setGeneric("default.PCA.edges<-", function(x, value) standardGeneric("default.PCA.edges<-"))

setMethod("default.PCA.edges", "MigrDat", function(x) x@dimreductions[["E"]]$PCA)
setMethod("default.PCA.edges<-", "MigrDat", function(x, value) {
  x@dimreductions[["E"]]$PCA <- value
  x
})



#' ## Name of the default UMAP in MigrDat object ####
#' @examples
#'
#' default.UMAP.edges(MigrDat)
#'
#' @export
setGeneric("default.UMAP.edges", function(x) standardGeneric("default.UMAP.edges"))
setGeneric("default.UMAP.edges<-", function(x, value) standardGeneric("default.UMAP.edges<-"))

setMethod("default.UMAP.edges", "MigrDat", function(x) x@dimreductions[["E"]]$UMAP)
setMethod("default.UMAP.edges<-", "MigrDat", function(x, value) {
  x@dimreductions[["E"]]$UMAP <- value
  x
})



#' ## Name of the default UMAP3D in MigrDat object ####
#' @examples
#'
#' default.UMAP3D.edges(MigrDat)
#'
#' @export
setGeneric("default.UMAP3D.edges", function(x) standardGeneric("default.UMAP3D.edges"))
setGeneric("default.UMAP3D.edges<-", function(x, value) standardGeneric("default.UMAP3D.edges<-"))

setMethod("default.UMAP3D.edges", "MigrDat", function(x) x@dimreductions[["E"]]$UMAP3D)
setMethod("default.UMAP3D.edges<-", "MigrDat", function(x, value) {
  x@dimreductions[["E"]]$UMAP3D <- value
  x
})










#### Clustering methods ####


#' ## Get and store whole clustering data sets from MigrDat object ####
#' @examples
#'
#' clusterings(MigrDat_XMLtest)
#'
#' @export
setGeneric("clusterings", function(x) standardGeneric("clusterings"))
setGeneric("clusterings<-", function(x, value) standardGeneric("clusterings<-"))

setMethod("clusterings", "MigrDat", function(x) x@clustering)
setMethod("clusterings<-", "MigrDat", function(x, value) {
  x@clustering <- value
  x
})


#' ## Get list of names of the clusterings in the MigrDat object ####
#' @examples
#'
#' clusterings.names(MigrDat)
#'
#' @export
setGeneric("clusterings.names", function(x) standardGeneric("clusterings.names"))
setMethod("clusterings.names", "MigrDat", function(x) sapply(x@clustering, names))



#' ## Get list of names of the clusterings in the MigrDat object ####
#' @examples
#'
#' clusterings.spots(MigrDat)
#'
#' @export
setGeneric("clusterings.spots", function(x) standardGeneric("clusterings.spots"))
setMethod("clusterings.spots", "MigrDat", function(x) clusterings.names(x)[["S"]] )


#' ## Get list of names of the clusterings in the MigrDat object ####
#' @examples
#'
#' clusterings.tracks(MigrDat)
#'
#' @export
setGeneric("clusterings.tracks", function(x) standardGeneric("clusterings.tracks"))
setMethod("clusterings.tracks", "MigrDat", function(x) clusterings.names(x)[["T"]])

#' ## Get list of names of the clusterings in the MigrDat object ####
#' @examples
#'
#' clusterings.edges(MigrDat)
#'
#' @export
setGeneric("clusterings.edges", function(x) standardGeneric("clusterings.edges"))
setMethod("clusterings.edges", "MigrDat", function(x) clusterings.names(x)[["E"]])



#' ## Get and store raw spots clustering in MigrDat object ####
#' @examples
#'
#' spots.clusters(MigrDat)
#' spots.clusters(MigrDat) <- new.dt
#'
#' @export
setGeneric("spots.clusters", function(x) standardGeneric("spots.clusters"))
setGeneric("spots.clusters<-", function(x, value) standardGeneric("spots.clusters<-"))

setMethod("spots.clusters", "MigrDat", function(x) x@clustering[["S"]])
setMethod("spots.clusters<-", "MigrDat", function(x, value) {
  x@clustering[["S"]] <- value
  x
})

#' ## Get and store raw edges clustering in MigrDat object ####
#' @examples
#'
#' edges.clusters(MigrDat)
#' edges.clusters(MigrDat) <- new.dt
#'
#' @export
setGeneric("edges.clusters", function(x) standardGeneric("edges.clusters"))
setGeneric("edges.clusters<-", function(x, value) standardGeneric("edges.clusters<-"))

setMethod("edges.clusters", "MigrDat", function(x) x@clustering[["T"]])
setMethod("edges.clusters<-", "MigrDat", function(x, value) {
  x@clustering[["T"]] <- value
  x
})

#' ## Get and store raw tracks clustering in MigrDat object ####
#' @examples
#'
#' tracks.clusters(MigrDat)
#' tracks.clusters(MigrDat) <- new.dt
#'
#' @export
setGeneric("tracks.clusters", function(x) standardGeneric("tracks.clusters"))
setGeneric("tracks.clusters<-", function(x, value) standardGeneric("tracks.clusters<-"))

setMethod("tracks.clusters", "MigrDat", function(x) x@clustering[["E"]])
setMethod("tracks.clusters<-", "MigrDat", function(x, value) {
  x@clustering[["E"]] <- value
  x
})



#' ## Get default clusterings ####
#' @name default.clusts.spots
#' @examples
#'
#' default.clusts.spots(MigrDat)
#'
#' @export
setGeneric("default.clusts.spots", function(x) standardGeneric("default.clusts.spots"))
setGeneric("default.clusts.spots<-", function(x, value) standardGeneric("default.clusts.spots<-"))

setMethod("default.clusts.spots", "MigrDat", function(x) x@clustering[["S"]][,c("kmeans","hclusts","ILoRegclusters","Seurclusters")])
setMethod("default.clusts.spots<-", "MigrDat", function(x, value) {
  x@clustering[["S"]][,c("kmeans","hclusts","ILoRegclusters","Seurclusters")] <- value
  x
})


#' ## Get ####
#' @name default.clusts.tracks
#' @examples
#'
#' default.clusts.tracks(MigrDat)
#'
#' @export
setGeneric("default.clusts.tracks", function(x) standardGeneric("default.clusts.tracks"))
setGeneric("default.clusts.tracks<-", function(x, value) standardGeneric("default.clusts.tracks<-"))

setMethod("default.clusts.tracks", "MigrDat", function(x) x@clustering[["S"]][,c("kmeans","hclusts","ILoRegclusters","Seurclusters")])
setMethod("default.clusts.tracks<-", "MigrDat", function(x, value) {
  x@clustering[["S"]][,c("kmeans","hclusts","ILoRegclusters","Seurclusters")] <- value
  x
})

#' ## Get ####
#' @name default.clusts.edges
#' @examples
#'
#' default.clusts.edges(MigrDat)
#'
#' @export
setGeneric("default.clusts.edges", function(x) standardGeneric("default.clusts.edges"))
setGeneric("default.clusts.edges<-", function(x, value) standardGeneric("default.clusts.edges<-"))

setMethod("default.clusts.edges", "MigrDat", function(x) x@clustering[["S"]][,c("kmeans","hclusts","ILoRegclusters","Seurclusters")])
setMethod("default.clusts.edges<-", "MigrDat", function(x, value) {
  x@clustering[["S"]][,c("kmeans","hclusts","ILoRegclusters","Seurclusters")] <- value
  x
})






#' ## Name of the default kmeans clustering in the MigrDat object ####
#' @examples
#'
#' default.kmeans(MigrDat)
#'
#' @export
setGeneric("default.kmeans", function(x) standardGeneric("default.kmeans"))
setGeneric("default.kmeans<-", function(x, value) standardGeneric("default.kmeans<-"))

setMethod("default.kmeans", "MigrDat", function(x) x@metadata$kmeans)
setMethod("default.kmeans<-", "MigrDat", function(x, value) {
  x@metadata$kmeans <- value
  x
})







#' ## Name of the default hclusts clustering in the MigrDat object ####
#' @examples
#'
#' default.hclusts(MigrDat)
#'
#' @export
setGeneric("default.hclusts", function(x) standardGeneric("default.hclusts"))
setGeneric("default.hclusts<-", function(x, value) standardGeneric("default.hclusts<-"))

setMethod("default.hclusts", "MigrDat", function(x) x@metadata$hclusts)
setMethod("default.hclusts<-", "MigrDat", function(x, value) {
  x@metadata$hclusts <- value
  x
})

#' ## Name of the default ILoRegclusters clustering in the MigrDat object ####
#' @examples
#'
#' default.ILoRegclusters(MigrDat)
#'
#' @export
setGeneric("default.ILoRegclusters", function(x) standardGeneric("default.ILoRegclusters"))
setGeneric("default.ILoRegclusters<-", function(x, value) standardGeneric("default.ILoRegclusters<-"))

setMethod("default.ILoRegclusters", "MigrDat", function(x) x@metadata$ILoRegclusters)
setMethod("default.ILoRegclusters<-", "MigrDat", function(x, value) {
  x@metadata$ILoRegclusters <- value
  x
})

#' ## Name of the default Seurclusters clustering in the MigrDat object ####
#' @examples
#'
#' default.Seurclusters(MigrDat)
#'
#' @export
setGeneric("default.Seurclusters", function(x) standardGeneric("default.Seurclusters"))
setGeneric("default.Seurclusters<-", function(x, value) standardGeneric("default.Seurclusters<-"))

setMethod("default.Seurclusters", "MigrDat", function(x) x@metadata$Seurclusters)
setMethod("default.Seurclusters<-", "MigrDat", function(x, value) {
  x@metadata$Seurclusters <- value
  x
})






#### IDs sets ####


#' ## Get and store whole IDs set from MigrDat object ####
#' @examples
#'
#' IDs.LS(MigrDat_XMLtest)
#'
#' @export
setGeneric("IDs.LS", function(x) standardGeneric("IDs.LS"))
setGeneric("IDs.LS<-", function(x, value) standardGeneric("IDs.LS<-"))

setMethod("IDs.LS", "MigrDat", function(x) x@IDs)
setMethod("IDs.LS<-", "MigrDat", function(x, value) {
  x@IDs <- value
  x
})

#' ## Get and store spots IDs and LABEls from MigrDat object ####
#' @examples
#'
#' IDs.spots(MigrDat_XMLtest)
#'
#' @export
setGeneric("IDs.spots", function(x) standardGeneric("IDs.spots"))
setGeneric("IDs.spots<-", function(x, value) standardGeneric("IDs.spots<-"))

setMethod("IDs.spots", "MigrDat", function(x) x@IDs$spots)
setMethod("IDs.spots<-", "MigrDat", function(x, value) {
  x@IDs$spots <- value
  x
})


#' ## Get and store tracks IDs and LABEls from MigrDat object ####
#' @examples
#'
#' IDs.tracks(MigrDat_XMLtest)
#'
#' @export
setGeneric("IDs.tracks", function(x) standardGeneric("IDs.tracks"))
setGeneric("IDs.tracks<-", function(x, value) standardGeneric("IDs.tracks<-"))

setMethod("IDs.tracks", "MigrDat", function(x) x@IDs$tracks)
setMethod("IDs.tracks<-", "MigrDat", function(x, value) {
  x@IDs$tracks <- value
  x
})


#' ## Get and store edges IDs and LABEls from MigrDat object ####
#' @examples
#'
#' IDs.edges(MigrDat_XMLtest)
#'
#' @export
setGeneric("IDs.edges", function(x) standardGeneric("IDs.edges"))
setGeneric("IDs.edges<-", function(x, value) standardGeneric("IDs.edges<-"))

setMethod("IDs.edges", "MigrDat", function(x) x@IDs$edges)
setMethod("IDs.edges<-", "MigrDat", function(x, value) {
  x@IDs$edges <- value
  x
})


#' ## Get and store roi_points IDs and LABEls from MigrDat object ####
#' @examples
#'
#' IDs.roi_points(MigrDat_XMLtest)
#'
#' @export
setGeneric("IDs.roi_points", function(x) standardGeneric("IDs.roi_points"))
setGeneric("IDs.roi_points<-", function(x, value) standardGeneric("IDs.roi_points<-"))

setMethod("IDs.roi_points", "MigrDat", function(x) x@IDs$roi_points)
setMethod("IDs.roi_points<-", "MigrDat", function(x, value) {
  x@IDs$roi_points <- value
  x
})



#' ## Get and store metadata info from MigrDat object ####
#' @examples
#'
#' metadata.LS(MigrDat_XMLtest)
#'
#' @export
setGeneric("metadata.LS", function(x) standardGeneric("metadata.LS"))
setGeneric("metadata.LS<-", function(x, value) standardGeneric("metadata.LS<-"))

setMethod("metadata.LS", "MigrDat", function(x) x@metadata)
setMethod("metadata.LS<-", "MigrDat", function(x, value) {
  x@metadata <- value
  x
})


#' ## Show stored TrackMate log in MigrDat object ####
#'
#'
#' @export
setGeneric("ShowTMimageSettings", function(x) standardGeneric("ShowTMimageSettings"))
setGeneric("addTMimageSettings<-", function(x, value) standardGeneric("addTMimageSettings<-"))

setMethod("ShowTMimageSettings", "MigrDat", function(x) cat(XML::toString.XMLNode(x@metadata$imageSettingsXML)))
setMethod("addTMimageSettings<-", "MigrDat", function(x, value) {
  x@metadata$imageSettingsXML <- value
  x
})


#' ## Show stored TrackMate log in MigrDat object ####
#'
#'
#' @export
setGeneric("ShowTMlog", function(x) standardGeneric("ShowTMlog"))
setGeneric("AddTMlog<-", function(x, value) standardGeneric("AddTMlog<-"))

setMethod("ShowTMlog", "MigrDat", function(x) cat(x@metadata$TMlog))
setMethod("AddTMlog<-", "MigrDat", function(x, value) {
  x@metadata$TMlog <- value
  x
})


#' ## Get and store whole spots set from MigrDat object ####
#' @examples
#'
#' spotsLS(MigrDat_XMLtest)
#'

#' @export
setGeneric("tm_track_filter", function(x) standardGeneric("tm_track_filter"))
setGeneric("tm_track_filter<-", function(x, value) standardGeneric("tm_track_filter<-"))

setMethod("tm_track_filter", "MigrDat", function(x) x@metadata$track_filter)
setMethod("tm_track_filter<-", "MigrDat", function(x, value) {
  x@metadata$track_filter <- value
  x
})




#### ILoRegMigr object ####
##### ILoRegMigr object planning ##########################################
#' Planned to contain:
#' in_slot sample_metadata sample_metadata raw data tables
#' inside metadata list(spot_metadata, track_metadata, edge_metadata)
#' inside metadata list(edge_metadata = data.frame())
#' inside IDs list(sample_IDs, spot_IDs, track_IDs, edge_IDs)
#' inside joint_tables list(joint_spots, joint_tracks, joint_edges)
#### items inside slots
#' projects.included = "list",
#' projects.active = "list",
#' active.experiments = "ANY",
#' active.slot = "ANY",
#' active.ident = "factor",
#' active.project = "character",
#' sample_metadata = "ANY",
#' spot_metadata = "ANY",
#' track_metadata = "ANY",
#' edge_metadata = "ANY",
#' sample_ID = "ANY",
#' spot_ID = "ANY",
#' track_ID = "ANY",
#' edge_ID = "ANY",
#' joint_spot = "ANY",
#' joint_track = "ANY",
#' joint_edge = 'ANY',
#' reductions.proj = "list",
#'     clusters = "list",
#' reductions.all = "list",
#' config.project = "character",
#' slot clusters dimension reduction coordinates
#### Considered
#' parameters slot to save used parameters
#' images reserved for images


#### ILoRegMigr data class ####

#' @title S4 ILoRegMigr Class
#' @description Object to store and analyze cell migration data
#' @concept ILoRegMigr object
#' @slot datasets raw data and metadata
#' @slot projects projects derived from raw data
#' @slot metadata general metadata slot containing clusterings, cell types, etc...
#' @slot active metadata slot
#' @slot IDs various IDs
#' @slot joint_tables reserved for storage of joint tables
#' @slot spatial_analysis spatial data
#' @slot reductions dimension reduction coordinates
#' @slot config Operating System info
#' @slot instructions slot for instructions
#' @slot history slot for record history to repeat analysis
#' @slot misc miscellaneous
#' @details
#'
#'
#' @export
ILoRegMigr <- setClass(
  "ILoRegMigr",
  slots = c(
    datasets = "list",
    projects = "list",
    metadata = "list",
    active = "list",
    IDs = "list",
    joint_tables = "list",
    spatial_analysis = "ANY",
    reductions = "list",
    config = "list",
    instructions = "ANY",
    history = "ANY",
    misc = "ANY"
  ),

  prototype = list(
    datasets = NULL,
    projects = NULL,
    metadata = NULL,
    active = NULL,
    IDs = NULL,
    joint_tables = NULL,
    spatial_analysis = NULL,
    reductions = NULL,
    config = NULL,
    instructions = NULL,
    history = NULL,
    misc = NULL
  )
)


#'
#'
#'
#'
#'
#'@export
ILoRegMigr.newobj <- function(project_name = "Analysis1", conditions = c("con1","cond2")) {
  # Create object
  dobj = new("ILoRegMigr")
  dobj@metadata[["project"]]  = project_name
  dobj@metadata[["conditions"]]  = conditions
  dobj@datasets = list()
  dobj@datasets[["metadata"]] = list()
  dobj@datasets[["metadata"]][["MigrDats"]] = length(dobj@datasets) -1
  dobj@IDs[["WellID_hi"]] = 0
  for (cond in conditions) {
    dobj@IDs[["CondID_Nums"]] = list(cond = 0)
  }
  dobj@IDs[["spot_ID"]] = list()
  dobj@IDs[["track_ID"]] = list()
  dobj@IDs[["edge_ID"]] = list()
  dobj@IDs[["sample_ID"]] = list()
  dobj@joint_tables[["joint_spots"]]  = list()
  dobj@joint_tables[["joint_tracks"]]  = list()
  dobj@joint_tables[["joint_edges"]]  = list()

  return(dobj)
}


#### ILoRegMigr methods ####
#' ## Set the output of the MigrDat object when called ####
#' Add more meaningful info later
#'
#' @export
setMethod("show",
          "ILoRegMigr",
          function(object){
            cat("An object of class", class(x = object), "has following slots available: \n\n")
            cat(slotNames(object), "with:\n")
            cat( (length(object@datasets) -1), "datasets of cell migration experiments using conditions:", (object@metadata[["conditions"]]) ,"\n")
          })




#' show method for ILoRegMigr class
#' @param object ILoRegMigr object
#
# setMethod(
#   f = "show",
#   signature = "ILoRegMigr",
#   definition = function(object) {
#
#     cat("An object of class",  class(object), "\n")
#
#     for(spat_unit in names(object@expression_feat)) {
#       cat('spatial units = ', spat_unit, '\n')
#       for(feat_type in unique(object@expression_feat)) {
#         cat("features = ", feat_type, "\n")
#
#         cat(
#           nrow(x = object@expression[[spat_unit]][[feat_type]][['raw']]),
#           "features across",
#           ncol(x = object@expression[[spat_unit]][[feat_type]][['raw']]),
#           "samples.\n \n"
#         )
#       }
#     }
#
#
#     cat('Steps and parameters used: \n \n')
#     print(object@parameters)
#     invisible(x = NULL)
#   }
# )

#' ILoRegMigr.newobjS3
#'
#'
#'
#'
#'
ILoRegMigr.newobjS3 <- function(dataset_name = "data1",
                                experiment_name = "Exp1",
                                condition_names = c("C1","C2")
) {
  # Create object
  dobj = list()
  dobj[["datasets"]] = list()
  dobj[["datasets"]][["metadata"]] = list()
  dobj[["datasets"]][[dataset_name]] = list()
  dobj[["datasets"]][[dataset_name]][["Metadata"]] = list()
  dobj[["metadata"]] = list()
  dobj[["active"]] = list()
  dobj[["IDs"]] = list()
  dobj[["IDs"]][["WellID_hi"]] = 0
  for (cond in condition_names) {
    dobj[["IDs"]][["CondID_Nums"]] = list(cond = 0)
  }
  dobj[["IDs"]][["spot_ID"]] = list()
  dobj[["IDs"]][["track_ID"]] = list()
  dobj[["IDs"]][["edge_ID"]] = list()
  dobj[["IDs"]][["sample_ID"]] = list()
  dobj[["joint_tables"]] = list()
  dobj[["joint_tables"]][["joint_spots"]]  = list()
  dobj[["joint_tables"]][["joint_tracks"]]  = list()
  dobj[["joint_tables"]][["joint_edges"]]  = list()
  dobj[["spatial_analysis"]] = list()
  dobj[["reductions"]] = list()
  dobj[["clusters"]] = list()
  dobj[["parameters"]] = list()
  dobj[["config"]] = list()
  dobj[["instructions"]] = list()
  dobj[["history"]] = list()
  dobj[["misc"]] = list()
  return(dobj)
}
