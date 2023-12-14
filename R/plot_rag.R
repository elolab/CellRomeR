plot_rag <- function(MigrObj) {
  
  if (length(unique(MigrObj@clustering$S$ILoRegclusters))==1) {
    stop("No clustering data available. Please run the clustering first.")
  }
  
  type = "S"
  dat = spots.raw(MigrObj) 
  dat$clusts = clusterings(MigrObj)[[type]][["ILoRegclusters"]]
  clst.cols = pals::brewer.accent(length(unique(dat$clusts)))
  
  
  
  ggplot2::ggplot(data = dat, mapping = ggplot2::aes_string(x = "FRAME", y = "SPEED", color = "clusts")) +
    ggplot2::scale_color_manual(values = clst.cols) +
    ggplot2::geom_point(size = 1, shape = 15) +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = 'black', colour = 'grey'),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank())

  
}

