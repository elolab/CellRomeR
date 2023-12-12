plot.MigrDat <- function(x) {
  
  n <- length(unique(migrdata@clustering$S$ILoRegclusters))
  colors <- rainbow(n)[migrdata@clustering$S$ILoRegclusters]
  
  x <- migrdata@spots$raw$POSITION_X
  y <- migrdata@spots$raw$POSITION_Y
  plot(x, y, pch=NA)
  
  ord <- order(migrdata@spots$raw$POSITION_T, decreasing=FALSE)
  for(track in unique(migrdata@spots$raw$TRACK_ID)) {
    sel <- which(migrdata@spots$raw$TRACK_ID[ord]==track)
    lines(x[ord][sel], y[ord][sel], col="lightgrey")
  }
  
  points(x, y, pch=21, bg=colors[ord], cex=0.5)
  
}
