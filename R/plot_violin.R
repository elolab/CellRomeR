# Plots violin plot of features grouped by clusters

plot_violin <- function(MigrObj, feature = NULL, data.slot="raw", type="S") {
  
  if (length(unique(MigrObj@clustering$S$ILoRegclusters))==1) {
    stop("No clustering data available. Please run the clustering first.")
  }
  
  if (is.null(feature)) {
    stop("Please specify a feature to plot")
  }
  
  if (type=="ES") {
    cluster <- clusterings(MigrObj)[["S"]][["ILoRegclusters"]]
    cluster.id <- getdt(MigrObj, data.slot, "S")[["SPOT_ID"]]
    values <- getdt(MigrObj, data.slot, "E")[[feature]]
    values.id <- getdt(MigrObj, data.slot, "E")[["SPOT_SOURCE_ID"]]
    cluster <- cluster[match(values.id, cluster.id)]
    violin_data <- data.frame(cluster, values)
  } else {
    cluster <- clusterings(MigrObj)[[type]][["ILoRegclusters"]]
    values <- getdt(MigrObj, data.slot, type)[[feature]]
    violin_data <- data.frame(cluster, values)
  }
  
  colors <- colorRampPalette(c("#3C95AC", "#E756A7", "#FFCB8A"))(length(unique(violin_data$cluster)))
  
  ggplot2::ggplot(violin_data, ggplot2::aes(x=cluster, y=values, fill=cluster)) +
    ggplot2::scale_fill_manual(values=colors) + 
    ggplot2::geom_violin() + ggplot2::theme_classic() + ggplot2::labs(y=feature)
    
}

