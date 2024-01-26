
plot_heatmap <- function(MigrObj, feature_vector=NULL) {
  
  if (length(unique(MigrObj@clustering$S$ILoRegclusters))==1) {
    stop("No clustering data available. Please run the clustering first.")
  }
  if (!is.null(feature_vector) && all(feature_vector %in% which(sapply(MigrObj@spots$raw, typeof)=="double" | sapply(MigrObj@spots$raw, typeof)=="integer"))) {
    stop("Please select a numeric column")
  }
  
  clusters <- unique(MigrObj@clustering$S$ILoRegclusters)
  if(length(feature_vector) == 0){
    markers <- which(sapply(MigrObj@spots$raw, typeof)=="double" | sapply(MigrObj@spots$raw, typeof)=="integer")
  } else {
    markers = feature_vector  
  }
  
  
  data.heatmap <- data.frame(MigrObj@spots$raw)
  data.heatmap <- t(apply(data.heatmap[,markers], 2, function(y) tapply(y, migrdata@clustering$S$ILoRegclusters, function(z)  median(z, na.rm=TRUE))))
  data.heatmap <- data.heatmap[which(apply(data.heatmap,1,sd)>0),]
  
  breaks <- seq(-2, 2, by=0.1)
  col <- colorRampPalette(c("darkmagenta", "white", "royalblue"))(length(breaks)-1)
  
  pheatmap::pheatmap(data.heatmap, scale="row", col=col, breaks=breaks)


  
}





