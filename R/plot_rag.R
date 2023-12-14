plot_rag <- function(MigrObj) {
  
  if (length(unique(MigrObj@clustering$S$ILoRegclusters))==1) {
    stop("No clustering data available. Please run the clustering first.")
  }
  
  type = "S"
  dat = spots.raw(MigrObj) 
  dat$clusts = clusterings(MigrObj)[[type]][["ILoRegclusters"]]
  clst.cols = pals::brewer.accent(length(unique(dat$clusts)))
  
  
  
  ggplot(data = dat, mapping = aes_string(x = "FRAME", y = "SPEED", color = "clusts")) +
    scale_color_manual(values = clst.cols) +
    geom_point(size = 1, shape = 15) +
    theme(panel.background = element_rect(fill = 'black', colour = 'grey'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

  
}

