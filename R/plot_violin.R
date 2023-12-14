# Plots violin plot of features grouped by clusters

plot_violin <- function(MigrObj, feature = NULL, data.slot = "raw", type = "S", clusterType = "ILoRegclusters") {
  if (is.null(feature)) {
    stop("Please specify a feature to plot")
  }
  
  violin_data <- data.frame(
    cluster = clusterings(MigrObj)[[type]][[clusterType]],
    feature = getdt(MigrObj, data.slot, type)[[feature]]
  )
  
  ggplot2::ggplot(violin_data, ggplot2::aes(x=cluster, y=feature, fill=cluster)) +
    ggplot2::geom_violin() + ggplot2::theme_classic() + ggplot2::labs(y=feature)
  
}

