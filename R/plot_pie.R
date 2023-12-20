# Plots number of spots per cluster

plot_pie <- function(MigrObj) {
  
  if (length(unique(MigrObj@clustering$S$ILoRegclusters))==1) {
    stop("No clustering data available. Please run the clustering first.")
  }
  
  clusters <- MigrObj@clustering$S$ILoRegclusters
  names <- unique(MigrObj@clustering$S$ILoRegclusters)
  colors <- colorRampPalette(c("#3C95AC","#E756A7","#FFCB8A"))(length(names))
  names(colors) <- names
  
  tbl <- table(clusters)
  labels <- paste(names(tbl), " (", round((tbl/sum(tbl))*100,1), "%)", sep="")
  pie(tbl, labels=labels, col=colors[names(tbl)], radius=0.5, cex=0.5)
  
}

