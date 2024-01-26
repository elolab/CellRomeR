plot_tsne <- function(MigrObj, uniq=NULL, color.by="ILoRegclusters", 
                      spots.slot="raw", point.size=1) {
  
  if (length(unique(MigrObj@clustering$S$ILoRegclusters))==1) {
    stop("No clustering data available. Please run the clustering first.")
  }
  
  if (is.null(uniq)) {
    stop("Please specify a named clustering.")
  }
  
  if (color.by=="ILoRegclusters") {
    var.plot <- clusterings(MigrObj)[["S"]]$ILoRegclusters
  } else if (color.by %in% colnames(MigrObj@spots[[spots.slot]])) {
    var.plot <- MigrObj@spots[[spots.slot]][[color.by]]
  } else {
    stop("'color.by' variable not found. Please check the variable provided.")
  }
  
  tsne <- MigrObj@dimreductions[["S"]][[paste0(uniq, "_tsne")]]
  
  if (is.factor(var.plot) | is.character(var.plot)) {
    var.plot <- as.factor(var.plot)
    colors <- colorRampPalette(c("#3C95AC","#E756A7","#FFCB8A"))(length(unique(var.plot)))
    names(colors) <- levels(var.plot)
    color.legend <- levels(var.plot)
    par(mar = c(10, 4.5, 4.5, 4.5)) # b, l, t, r
    plot(tsne[,1], tsne[,2], bty="n", pch=21, cex=point.size, col="white", bg=colors[var.plot],  
         xlab="tsne 1", ylab="tsne 2")
    box(bty = "l")
    legend("bottom", title = color.by, legend=color.legend, pch=21,
           ncol=ceiling(length(unique(var.plot))/2),
           pt.bg=colors, bty="n", inset=c(0,-0.45),
           xpd=TRUE)
  }
  if (is.numeric(var.plot)) {
    var.order <- findInterval(var.plot, sort(var.plot))
    colors <- colorRampPalette(c("blue", "white", "red"))(length(var.plot))
    color.legend <- range(var.plot)
    color.legend[3] <- ((color.legend[2] - color.legend[1])/2) + color.legend[1]
    color.legend <- round(sort(color.legend), 1)
    par(mar = c(10, 4.5, 4.5, 4.5)) # b, l, t, r
    plot(tsne[,1], tsne[,2], bty="n", pch=19, cex=point.size, col=colors[var.order],  
         xlab="tsne 1", ylab="tsne 2")
    box(bty = "l")
    legend("bottom", title = color.by, legend=color.legend, pch=21,  
           col="black", pt.bg=c("blue", "white", "red"), bty="n", 
           inset=c(0,-0.45), xpd=TRUE, ncol=3, cex=1)
  }
}