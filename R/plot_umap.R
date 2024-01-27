plot_umap <- function(MigrObj, uniq="ILoReg", color.by="ILoRegclusters", spots.slot="raw", cex=0.1) {
  
  if (length(unique(MigrObj@clustering$S$ILoRegclusters))==1) {
    stop("No clustering data available. Please run the clustering first.")
  }
  
  if (is.null(uniq)) {
    stop("Please specify a named clustering")
  }
  
  umap <- MigrObj@dimreductions[["S"]][[paste0(uniq, "_UMAP")]]
  
  if (color.by=="ILoRegclusters") {
    clusters <- clusterings(MigrObj)[["S"]]$ILoRegclusters
    colors <- colorRampPalette(c("#3C95AC","#E756A7","#FFCB8A"))(length(unique(clusters)))
    plot(umap[,1], umap[,2], bty="n", pch=19, cex=cex, col=colors[clusters], axes=FALSE, xlab="", ylab="")
    legend("bottom", legend=unique(clusters), pch=21, ncol=ceiling(length(unique(clusters))/2), pt.bg=colors, bty="n")
  } else if (color.by %in% colnames(MigrObj@spots[[spots.slot]])) {
    values <- MigrObj@spots[[spots.slot]][[color.by]]
    var.order <- findInterval(values, sort(values))
    colors <- colorRampPalette(c("#FCB317", "white", "#AB1F91"))(length(values))
    color.legend <- range(values)
    color.legend[3] <- ((color.legend[2] - color.legend[1])/2) + color.legend[1]
    color.legend <- round(sort(color.legend), 1)
    col=colors[var.order]
    plot(umap[,1], umap[,2], bty="n", pch=19, cex=cex, col=colors[var.order], axes=FALSE, xlab="", ylab="")
  } else {
    stop("'color.by' variable not found. Please check the variable provided.")
  }

}
