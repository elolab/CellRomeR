plot_umap <- function(MigrObj, uniq="ILoReg", color.by="ILoRegclusters", spots.slot="raw", cex=0.1,use.opinionated.layout = TRUE) {
  matetype <- "S" # one of "S","T","E"
  is_clustering_name  <- function(x)
    x %in% names(MigrObj@clustering[[matetype]])
  select_clustering <- function(x)
    clusterings(MigrObj)[[matetype]][[x]]
  
  
  is_matetype_column <- function(x){
    x %in% colnames(
      switch(EXPR = matetype,
             S = 
               MigrObj@spots,
             E =
               MigrObj@edges,
             T =
               MigrObj@tracks)[[spots.slot]])

  }
  
  select_matetype_column <- function(x){
      switch(EXPR = matetype,
             S = 
               MigrObj@spots,
             E =
               MigrObj@edges,
             T =
               MigrObj@tracks)[[spots.slot]][[x]]
  }
  
  find_color_sequence_discrete <- function(categories)
  {
        colorRampPalette(c("#3C95AC","#E756A7","#FFCB8A"))(length(unique(categories)))
  }
  
  
  map_colors_discrete <- function(categories){
    find_color_sequence_discrete(categories)[categories]
  }

  find_color_sequence_continuous <- function(values)
    colorRampPalette(c("#FCB317", "white", "#AB1F91"))(length(values))
  
  map_colors_continuous <- function(values)
  {
    var.order <- findInterval(values, sort(values))
    colors <- find_color_sequence_continuous(values)
    colors[var.order]  
  }
  # default plotting function
  # where colormapf is a fucntion in the form of map_colors_{discrete,continuous}
  fplot <- function(umap,values,colormapf,...)
    plot(umap[,1], umap[,2], bty="n", pch=19, cex=cex, col=colormapf(values), axes=FALSE, xlab="", ylab="",...)
  
  plot_continuous_minimal <- function(umap,values,...)
    fplot(umap,values,map_colors_continuous,...)
    
  plot_discrete_minimal <- function(umap,values,...)
  {
    fplot(umap,values,map_colors_discrete,...)
    legend("bottom", legend=unique(values), pch=21, ncol=ceiling(length(unique(values))/2), pt.bg=find_color_sequence_discrete(values), bty="n")
  }
  
  plot_discrete_opinionated <- function(umap, values)
  {
    layout(matrix(c(1,2), ncol=2), widths=c(4,1))
    oldpar <-     par(mar=c(5, 4, 4, 2))
    on.exit(par(oldpar))
    fplot(umap, values,map_colors_discrete, main = color.by,)
    par(mar=c(5,0,4,1))
    plot.new()  # Empty plot space for legend 
    legend("center", legend=unique(values), pch=21, pt.bg=find_color_sequence_discrete(values), bty="n", cex=1.5)  # Increase legend dot size
  }
  # find minimum, maximum and the average of the two
  find_min_midp_max <- function(values)
  {
    summ <- range(values)
    summ[3] <- ((summ[2] - summ[1])/2) + summ[1]
    summ <- round(sort(summ), 1)
    summ
  }
    
  plot_continuous_opinionated <- function(umap,values)
  {
    layout(matrix(c(1,2), ncol=2), widths=c(4.5,0.5))
    oldpar <- par(mar=c(5,1,4,1))
    on.exit(par(oldpar))
    #par(mfrow=c(1, 2), mar=c(5, 4, 4, 2) + 0.1)
    plot_continuous_minimal(umap, values, main = color.by)
    par(mar=c(5, 1, 4, 1) + 0.1)
    colors <- find_color_sequence_continuous(values)
    min_midp_max <- find_min_midp_max(values)
    image(1, seq_along(colors), t(matrix(seq_along(colors), ncol=1)), col=colors, axes=FALSE)
    axis(4, at=seq(1, length(colors), length.out=length(min_midp_max)), labels=min_midp_max)
  }
  

  if (is_clustering_name(color.by) && length(unique(select_clustering(color.by))) == 1) {
    stop("No clustering data available. Please run the clustering first.")
  }
  
  if (is.null(uniq)) {
    stop("Please specify a named clustering")
  }
  
  umap <- MigrObj@dimreductions[[matetype]][[paste0(uniq, "_UMAP")]]
  
  if (is_clustering_name(color.by)) {
    (if(use.opinionated.layout) plot_discrete_opinionated else plot_discrete_minimal)(umap, select_clustering(color.by))
  } else if(is_matetype_column(color.by)){
    (if(use.opinionated.layout) plot_continuous_opinionated else plot_continuous_minimal)(umap, select_matetype_column(color.by))
  }
  else {
    stop("'color.by' is not one of: ",
         paste0(c(
           colnames(
             switch(matetype,
                    S = 
                      MigrObj@spots,
                    E =
                      MigrObj@edges,
                    T =
                      MigrObj@tracks)[[spots.slot]]),
           names(MigrObj@clustering[[matetype]])),
           collapse=","
           ))}}
