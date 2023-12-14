plot_rag <- function(MigrObj, clusterType = "ILoRegclusters") {
  
  type = "S"
  dat = spots.raw(MigrObj) 
  dat$clusts = clusterings(MigrObj)[[type]][[clusterType]]
  clst.cols = pals::brewer.accent(length(unique(dat$clusts)))
  
  
  
  ggplot(data = dat, mapping = aes_string(x = "FRAME", y = "SPEED", color = "clusts")) +
    scale_color_manual(values = clst.cols) +
    geom_point(size = 1, shape = 15) +
    theme(panel.background = element_rect(fill = 'black', colour = 'grey'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

  
}

