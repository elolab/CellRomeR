plot.MigrDat <- function(x, cex=1, ...) {
  
  n <- length(unique(x@clustering$S$ILoRegclusters))
  colors <- colorRampPalette(c("#3C95AC","#E756A7","#FFCB8A"))(n)[x@clustering$S$ILoRegclusters]
  
  x.coords <- x@spots$raw$POSITION_X
  y.coords <- x@spots$raw$POSITION_Y
  plot(x.coords, y.coords, pch=NA, bty="n", xaxt="n", yaxt="n", ann=FALSE, ...)
  
  ord <- order(x@spots$raw$POSITION_T, decreasing=FALSE)
  for(track in unique(x@spots$raw$TRACK_ID)) {
    sel <- which(x@spots$raw$TRACK_ID[ord]==track)
    lines(x.coords[ord][sel], y.coords[ord][sel], col="lightgrey")
  }
  for(track in unique(x@spots$raw$TRACK_ID)) {
    sel <- which(x@spots$raw$TRACK_ID[ord]==track)
    sel.a <- sel[-length(sel)]
    sel.b <- sel[length(sel)]
    points(x.coords[ord][sel.a], y.coords[ord][sel.a], pch=21, col=adjustcolor(colors[ord][sel.a], offset=c(-0.3, -0.3, -0.3, 0)), bg=colors[ord][sel.a], cex=0.5*cex)
    points(x.coords[ord][sel.b], y.coords[ord][sel.b], pch=21, col=adjustcolor(colors[ord][sel.b], offset=c(-0.3, -0.3, -0.3, 0)), bg=colors[ord][sel.b], cex=cex)
  }
  
}
