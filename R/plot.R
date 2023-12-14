plot.MigrDat <- function(MigrObj, ...) {
  
  n <- length(unique(MigrObj@clustering$S$ILoRegclusters))
  colors <- colorRampPalette(c("#3C95AC","#E756A7","#FFCB8A"))(n)[MigrObj@clustering$S$ILoRegclusters]
  
  x <- MigrObj@spots$raw$POSITION_X
  y <- MigrObj@spots$raw$POSITION_Y
  plot(x, y, pch=NA, ...)
  
  ord <- order(MigrObj@spots$raw$POSITION_T, decreasing=FALSE)
  for(track in unique(MigrObj@spots$raw$TRACK_ID)) {
    sel <- which(MigrObj@spots$raw$TRACK_ID[ord]==track)
    lines(x[ord][sel], y[ord][sel], col="lightgrey")
  }
  
  points(x, y, pch=21, bg=colors[ord], cex=0.5)
  
}
