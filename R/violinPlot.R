# Plots violin plot of features grouped by clusters

plot_violin <- function(MigrObj, dat.slot = NULL, feature = NULL, type = NULL, clusterType = "ILoRegclusters") {
  violin_data <- data.frame(
    cluster = clusterings(MigrDatXML)[[type]][[clusterType]],
    feature = getdt(MigrObj, dat.slot, type)[[feature]]

  )
  
  
  ggplot(violin_data, aes(x=cluster, y=feature, fill=cluster)) +
    geom_violin() + theme_classic() + labs(y=feature)
  
  
}

# Plots violin plot of features grouped by clusters horizontally
plot_hviolin <- function(MigrObj, dat.slot = NULL, feature = NULL, type = NULL, clusterType = "ILoRegclusters") {
  violin_data <- data.frame(
    cluster = clusterings(MigrDatXML)[[type]][[clusterType]],
    feature = getdt(MigrObj, dat.slot, type)[[feature]]
    
  )
  
  
  ggplot(violin_data, aes(x=cluster, y=feature, fill=cluster)) +
    geom_violin() + theme_classic() + coord_flip() + labs(y=feature)
  
  
}
