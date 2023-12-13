# Plots violin plot of features grouped by clusters

plotViolin <- function(clusters, features) {
  violin_data <- data.frame(
    cluster = clusters,
    feature = features
  )
  
  
  ggplot(violin_data, aes(x=cluster, y=feature, fill=cluster)) +
    geom_violin() + theme_classic()
  
  
}

# Plots violin plot of features grouped by clusters horizontally
plotHviolin <- function(clusters, features) {
  violin_data <- data.frame(
    cluster = clusters,
    feature = features
  )
  
  
  ggplot(violin_data, aes(x=cluster, y=feature, fill=cluster)) +
    geom_violin() + theme_classic() + coord_flip()
  
  
}