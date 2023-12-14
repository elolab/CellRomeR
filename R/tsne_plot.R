tSNE_plot <- function(MigrObj, type, uniq, colors) {
  clusters <- clusterings(MigrDatXML)[["S"]]$ILoRegclusters
  tsne <- MigrDatXML@dimreductions[["S"]][[paste("CH1_15", "_TSNE")]] 
  tsne.x = tsne[, 1]
  tsne.y = tsne[, 2]
  plot(tsne.x, tsne.y, bty="n", pch=19, cex=0.1, col=color_clusters[clusters], axes=FALSE, xlab="", ylab="")
}
