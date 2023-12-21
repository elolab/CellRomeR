plot_umap <- function(MigrObj, uniq=NULL) {
  
  if (length(unique(MigrObj@clustering$S$ILoRegclusters))==1) {
    stop("No clustering data available. Please run the clustering first.")
  }
  
  if (is.null(uniq)) {
    stop("Please specify a named clustering")
  }
  
  clusters <- clusterings(MigrObj)[["S"]]$ILoRegclusters
  tsne <- MigrObj@dimreductions[["S"]][[paste0(uniq, "_UMAP")]]
  
  colors <- colorRampPalette(c("#3C95AC","#E756A7","#FFCB8A"))(length(unique(clusters)))
  plot(tsne[,1], tsne[,2], bty="n", pch=19, cex=0.1, col=colors[clusters], axes=FALSE, xlab="", ylab="")
  
  legend("bottom", legend=unique(clusters), pch=21, ncol=ceiling(length(unique(clusters))/2), pt.bg=colors, bty="n")
  
}