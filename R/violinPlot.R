# Plots violin plot of features grouped by clusters

CellRomeR.violin <- function(clusters, features) {
  violin_data <- data.frame(
    cluster = clusters,
    feature = features
  )
  
  
  p <- ggplot(violin_data, aes(x=cluster, y=feature, fill=cluster)) +
    geom_violin() + theme_classic()
  
  return(p)
  
}

# Plots violin plot of features grouped by clusters horizontally
CellRomeR.hviolin <- function(clusters, features) {
  violin_data <- data.frame(
    cluster = clusters,
    feature = features
  )
  
  
  p <- ggplot(violin_data, aes(x=cluster, y=feature, fill=cluster)) +
    geom_violin() + theme_classic() + coord_flip()
  
  
  return(p)
  
}