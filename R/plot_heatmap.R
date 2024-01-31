
plot_heatmap <- function(MigrObj, features=NULL) {
  
  if (length(unique(MigrObj@clustering$S$ILoRegclusters))==1) {
    stop("No clustering data available. Please run the clustering first.")
  }
  if (!is.null(features) && all(features %in% which(sapply(MigrObj@spots$raw, typeof)=="double" | sapply(MigrObj@spots$raw, typeof)=="integer"))) {
    stop("Please select a numeric column")
  }
  
  clusters <- unique(MigrObj@clustering$S$ILoRegclusters)
  if(length(features) == 0){
    markers <- which(sapply(MigrObj@spots$raw, typeof)=="double" | sapply(MigrObj@spots$raw, typeof)=="integer")
  } else {
    markers = features  
  }
  
  
  data.heatmap <- data.frame(MigrObj@spots$raw)
  data.heatmap <- t(apply(data.heatmap[,markers], 2, function(y) tapply(y, MigrObj@clustering$S$ILoRegclusters, function(z)  median(z, na.rm=TRUE))))
  data.heatmap <- data.heatmap[which(apply(data.heatmap,1,sd)>0),]
  
  breaks <- seq(-2, 2, by=0.1)
  col <- colorRampPalette(c("#FCB317", "white", "#AB1F91"))(length(breaks)-1)
  
  pheatmap::pheatmap(data.heatmap, scale="row", col=col, breaks=breaks)


  
}





