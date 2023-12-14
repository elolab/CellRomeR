# Plots violin plot of features grouped by clusters

plot_violin <- function(MigrObj, feature = NULL, data.slot = "raw", type = "S") {
  
  if (length(unique(MigrObj@clustering$S$ILoRegclusters))==1) {
    stop("No clustering data available. Please run the clustering first.")
  }
  
  if (is.null(feature)) {
    stop("Please specify a feature to plot")
  }
  
  violin_data <- data.frame(
    cluster = clusterings(MigrObj)[[type]][["ILoRegclusters"]],
    feature = getdt(MigrObj, data.slot, type)[[feature]]
  )
  
  colors <- colorRampPalette(c("#3C95AC", "#E756A7", "#FFCB8A"))(length(unique(violin_data$cluster)))
  
  
  ggplot2::ggplot(violin_data, ggplot2::aes(x=cluster, y=feature, fill=cluster)) +
    ggplot2::scale_fill_manual(values=colors) + 
    ggplot2::geom_violin() + ggplot2::theme_classic() + ggplot2::labs(y=feature)
    
  
}

