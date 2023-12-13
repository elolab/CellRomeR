
plot_heatmap <- function(x) {
  
  if (is.null(x@clustering$S$ILoRegclusters)) {
    stop("No clustering data available. Please run the clustering first.")
  }
  
  clusters <- unique(x@clustering$S$ILoRegclusters)
  markers <- which(sapply(x@spots$raw, typeof)=="double" | sapply(x@spots$raw, typeof)=="integer")
  
  data.heatmap <- data.frame(x@spots$raw)
  data.heatmap <- t(apply(data.heatmap[,markers], 2, function(y) tapply(y, migrdata@clustering$S$ILoRegclusters, function(z)  median(z, na.rm=TRUE))))
  data.heatmap <- data.heatmap[which(apply(data.heatmap,1,sd)>0),]
  
  breaks <- seq(-2, 2, by=0.1)
  col <- colorRampPalette(c("darkmagenta", "white", "royalblue"))(length(breaks)-1)
  heatmap(data.heatmap, scale="row", col=col, breaks=breaks)
  
}





