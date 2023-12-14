plot.MigrDat <- function(x, ...) {
  
  n <- length(unique(x@clustering$S$ILoRegclusters))
  colors <- colorRampPalette(c("#3C95AC","#E756A7","#FFCB8A"))(n)[x@clustering$S$ILoRegclusters]
  
  x.coords <- x@spots$raw$POSITION_X
  y.coords <- x@spots$raw$POSITION_Y
  plot(x.coords, y.coords, pch=NA, ...)
  
  ord <- order(x@spots$raw$POSITION_T, decreasing=FALSE)
  for(track in unique(x@spots$raw$TRACK_ID)) {
    sel <- which(x@spots$raw$TRACK_ID[ord]==track)
    lines(x.coords[ord][sel], y.coords[ord][sel], col="lightgrey")
  }
  
  points(x.coords, y.coords, pch=21, bg=colors[ord], cex=0.5)
  
}
