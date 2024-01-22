plot_rag <- function(MigrObj, clusterType = "ILoRegclusters", feature = NULL) {
  if (length(unique(MigrObj@clustering$S))==1) {
    stop("No clustering data available. Please run the clustering first.")
  }
  
  if (!is.null(feature) && !is.numeric(MigrObj@tracks[["raw"]][[feature]])) {
    stop("Please select a numeric column")
  }
  
  type = "S"
  track_dat <- MigrObj@tracks[["raw"]] 
  dat <- MigrObj@spots[["raw"]] 
  if (!is.null(feature)){
    dat$features <- track_dat[[feature]][match(dat$TRACK_ID, track_dat$LABEL)]
    dat$combinedID <- paste0(dat$features, " | ", dat$TRACK_ID)
    xlabel = paste0("CombinedID of ", feature, " and TRACK_ID")
  } else {
    dat$features = dat$TRACK_ID
    dat$combinedID = dat$TRACK_ID
    xlabel = "TRACK_ID"
  }
  dat$clusters <- clusterings(MigrObj)[["S"]][[clusterType]]
  ragplot_data <- dat[order(features, TRACK_ID),]
  
  clst.cols <- colorRampPalette(c("#3C95AC", "#E756A7", "#FFCB8A"))(length(unique(dat$clusters)))
  
  
  
  
  ggplot2::ggplot(data = dat, mapping = ggplot2::aes_string(x = "combinedID", y = "FRAME", color = "clusters")) +
    ggplot2::scale_color_manual(values = clst.cols) +
    ggplot2::geom_point(size = 1, shape = 15) + ggplot2::labs(x=xlabel) +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = 'black', colour = 'grey'),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank())
  
  
}

