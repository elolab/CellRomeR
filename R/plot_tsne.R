plot_tsne <- function(MigrObj, uniq=NULL) {
  
  if (is.null(uniq)) {
    stop("Please specify a named clustering")
  }
  
  clusters <- clusterings(MigrObj)[["S"]]$ILoRegclusters
  tsne <- MigrObj@dimreductions[["S"]][[paste("CH1_15", "_TSNE")]]
  
  colors <- colorRampPalette(c("gold","royalblue","forestgreen","darkmagenta"))(length(unique(clusters)))
  plot(tsne[,1], tsne[,2], bty="n", pch=19, cex=0.1, col=colors[clusters], axes=FALSE, xlab="", ylab="")
  
  legend("bottom", legend=unique(clusters), pch=21, ncol=ceiling(length(unique(clusters))/2), pt.bg=colors, bty="n")
  
}
