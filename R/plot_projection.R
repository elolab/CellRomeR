plot_projection_factory <- function(projection_string, slotname_suffix = paste0("_",projection_string,))
{
  function(MigrObj, uniq, color.by="ILoRegclusters", 
           spots.slot="raw", point.size=1) {
    
    if (length(unique(MigrObj@clustering$S$ILoRegclusters))==1) {
      stop("No clustering data available. Please run the clustering first.")
    }

    create_dimname <-
      function(uniq)
        paste0(uniq,slutname_suffix)
    if(!create_dimname(uniq) %in% names(MigrObj@dimreductions[["S"]]))
             {
               stop(sprintf("Couldn't find %s in MigrObj,(derived from  %s). Try a different value for `uniq'",
                            
                            create_dimname(uniq),
                            uniq))
               
             }
    
    if (color.by=="ILoRegclusters") {
      var.plot <- clusterings(MigrObj)[["S"]]$ILoRegclusters
    } else if (color.by %in% colnames(MigrObj@spots[[spots.slot]])) {
      var.plot <- MigrObj@spots[[spots.slot]][[color.by]]
    } else {
      stop("'color.by' variable not found. Please check the variable provided.")
    }
    
    projection <- MigrObj@dimreductions[["S"]][[create_dimname(uniq)]]
    
    if (is.factor(var.plot) || is.character(var.plot)) {
      var.plot <- as.factor(var.plot)
      colors <- colorRampPalette(c("#3C95AC","#E756A7","#FFCB8A"))(length(unique(var.plot)))
      names(colors) <- levels(var.plot)
      color.legend <- levels(var.plot)
      par(mar = c(10, 4.5, 4.5, 4.5)) # b, l, t, r
      plot(projection[,1], projection[,2], bty="n", pch=21, cex=point.size, col="white", bg=colors[var.plot],  
           xlab=sprintf("%s 1",projection_string), ylab=sprintf("%s 2",projection_string))
      box(bty = "l")
      legend("bottom", title = color.by, legend=color.legend, pch=21,
             ncol=ceiling(length(unique(var.plot))/2),
             pt.bg=colors, bty="n", inset=c(0,-0.45),
             xpd=TRUE)
    }
    if (is.numeric(var.plot)) {
      var.order <- findInterval(var.plot, sort(var.plot))
      colors <- colorRampPalette(c("#FCB317", "white", "#AB1F91"))(length(var.plot))
      color.legend <- range(var.plot)
      color.legend[3] <- ((color.legend[2] - color.legend[1])/2) + color.legend[1]
      color.legend <- round(sort(color.legend), 1)
      par(mar = c(10, 4.5, 4.5, 4.5)) # b, l, t, r
      plot(projection[,1], projection[,2], bty="n", pch=19, cex=point.size, col=colors[var.order],
           xlab=sprintf("%s 1",projection_string), ylab=sprintf("%s 2",projection_string))
      box(bty = "l")
      legend("bottom", title = color.by, legend=color.legend, pch=21,  
             col="black", pt.bg=c("#FCB317", "white", "#AB1F91"), bty="n", 
             inset=c(0,-0.45), xpd=TRUE, ncol=3, cex=1)
    }
  }}

plot_umap <- plot_projection_factory("UMAP")
plot_tsne <- plot_projection_factory("tsne")
