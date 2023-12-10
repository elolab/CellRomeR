#### Plotting ####

##### Palettes #####
cbPalette8 <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
                "#D55E00", "#CC79A7")
colrz.clust = as.vector(pals::alphabet())
colrz.clust2= c("#E21900B2","#C840C8B2",'#8DD3C7',"#1900FFB2",'#FFFFB3','#BEBADA',
                '#FB8072','#80B1D3','#FDB462','#B3DE69','#FCCDE5',
                '#D9D9D9','#BC80BD','#CCEBC5','#FFED6F','#e31a1c',"#009E73",'#8e0152',
                '#276419',"#00F1F1B2","#F0BC95")
colrz.clust3 = as.vector(c(pals::alphabet(), pals::alphabet2(), pals::brewer.accent(32))) # up to XX clusters
colrz.clust4 = c(as.vector(pals::polychrome(36)), "#C8C8C8B2","#1900FFB2",
                 "#E21900B2", "#C840C8B2", "#00F1F1B2", "#010101B2" )

cbPalette <- c('#a6cee3','#fb9a99','#ffff99','#1f78b4','#b2df8a','#fdbf6f',
               '#33a02c','#ff7f00','#e31a1c','#cab2d6','#6a3d9a','#b15928')
colz.hto8 = c(c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
                "#D55E00", "#CC79A7")) # safe 8 colors
colz.hto9 = c("#C8C8C8B2","#1900FFB2", "#E21900B2", "#C840C8B2", "#00F1F1B2",
              "#010101B2",'#FFFF99','#BC80BD','#CCEBC5')
colz.hto12 = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
               '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928')
colz.hto12b = c('#8DD3C7','#FFFFB3','#BEBADA','#FB8072','#80B1D3','#FDB462','#B3DE69',
                '#FCCDE5','#D9D9D9','#BC80BD','#CCEBC5','#FFED6F')


#### Themes ####

## DimPlot themes
themeDPgr = ggplot2::theme(panel.background = ggplot2::element_rect(fill = "darkgrey",
                                                                    colour = "darkgrey", size = 0.25, linetype = "solid"),
                           panel.grid.major = ggplot2::element_line(size = 0.25, linetype = 'solid', colour = "white"),
                           panel.grid.minor = ggplot2::element_line(size = 0.15, linetype = 'solid', colour = "white") )
themeDP = ggplot2::theme(panel.background = ggplot2::element_rect(fill = "darkgrey", colour = "darkgrey", size = 0.25, linetype = "solid"))

catcolors <- function(mobj, group.by) {
  cats <- unique(sobj[[group.by]][[1]])

  lencats <- length(cats)
  if (is.null(group.by)) {
    colors <- as.vector(pals::alphabet2())
  } else if (lencats > 26) {
    colors <- colrz.clust4
  } else if (lencats > 21 & lencats < 27) {
    colors <- colrz.clust
  } else if (lencats > 12 & lencats < 22) {
    colors <- colrz.clust2
  } else if (lencats < 13) {
    colors <- cbPalette
  } else colors <- as.vector(pals::alphabet2())
  return(colors)
}

colrz <- function(colrzX = NA, names = NA) {

  # under dev. only up to 22 now
  if (is.na(colrzX)) {
    colrz.clust = c('#7f3b08','#b35806','#e08214','#fdb863','#fee0b6','#f7f7f7',
                    '#d8daeb','#b2abd2','#8073ac','#542788','#2d004b', '#8e0152',
                    '#c51b7d','#de77ae','#f1b6da','#fde0ef','#f7f7f7','#e6f5d0',
                    '#b8e186','#7fbc41','#4d9221','#276419')
  }  else { colrz.clust = colrz }
  # need to add matching names for colors Need update after annotation!
  if (!is.na(names)) {
    colrz.clust = colrz.clust[1:length(names)]
    names(colrz.clust) = names
  } else names(colrz.clust) = 1:length(colrz.clust)
  return(colrz.clust)
}

spectrl <- function(x) {
  ggplot2::scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))
}

##### Data boxplots #####
IRmigr.boxplot <- function(MigrObj, dat.slot = "raw", type = "STE", log2 = TRUE,
                           vars = NULL, incl.pattern = NULL, legend = TRUE,
                           predef =  "none", excld.pattern = NULL) {

  predef = match.arg(predef, c("none", "technical", "morphological", "clust"), several.ok = FALSE)

  data <- getdtcols(MigrObj, dat.slot = dat.slot, type = type,
                    predef = predef,
                    vars = vars, incl.pattern = incl.pattern,
                    excld.pattern = excld.pattern,
                    numerics = TRUE, rnames = TRUE)
  # data
  Bxplot = suppressWarnings(
    ggplot2::ggplot(reshape2::melt(data), ggplot2::aes(x=variable, y=value, group=variable)) +
      ggplot2::geom_boxplot(aes(fill=variable)) +
      ggplot2::theme(axis.text.x=ggplot2::element_text(angle=60,hjust=1)) +
      if (log2 == TRUE) {ggplot2::scale_y_continuous(trans='log2')} )
  if (!legend) {
    Bxplot <- Bxplot + ggplot2::theme(legend.position="none")
  }
  return(Bxplot)
}

IRmigr.corrplots <- function(MigrObj, dat.slot = "raw", type = "STE", log2 = TRUE,
                             vars = NULL, incl.pattern = NULL,
                             predef =  "none", excld.pattern = NULL, numerics = TRUE) {

  predef = match.arg(predef, c("none", "technical", "morphological", "clust"), several.ok = FALSE)

  data <- getdtcols(MigrObj, dat.slot = dat.slot, type = type,
                    predef = predef,
                    vars = vars, incl.pattern = incl.pattern,
                    excld.pattern = excld.pattern,
                    numerics = TRUE, rnames = TRUE)

  data[,STE:= type]
  corrplots = suppressWarnings(GGally::ggpairs(data,ggplot2::aes(alpha = 0.5)))

  return(corrplots)
}


#### Plots spots as tracks ####

##### Spots plotting #####

##### For categorical coloring #####
SpotPlotCut <- function(MigrObj, dat.slot = "scaled",
                        feature = NULL, cut = 6,
                        clusters = NULL,
                        facetting = NULL, fcut = 9,
                        grp = "TRACK_ID2",
                        point.size = 0.3,
                        legend = FALSE, saves = FALSE,
                        savename = "Tracks_Plot.png") {
  type = "S"
  dats = getdtcols(MigrObj, dat.slot = dat.slot, type = type, predef =  "coord")

  if (!is.null(grp)) {
    dats$grp = as.factor(getdtcols(MigrObj, dat.slot = dat.slot, type = type, vars = grp)[[1]] )
  } else if  ("TRACK_ID2" %in% getvars(MigrObj, dat.slot, type) ) {
    dats$grp = as.factor(getdtcols(MigrObj, dat.slot = dat.slot, type = type, vars = "TRACK_ID2")[[1]] )
  } else dats$grp = as.factor(dats$TRACK_ID)

  # Add feature to plotting table
  if (!is.null(feature)) {
    if (feature %in% getvars(MigrObj, dat.slot, type)) {
      dats <- cut2dt(MigrObj, dt = dats, feature = feature, dat.slot = dat.slot, cut = cut, type = type)
      colnames(dats)[grep(feature,colnames(dats))] <- "colors"
    } else if (feature %in% colnames(MigrObj@metadata[[type]]) ) {
      dats$colors = as.factor(MigrObj@metadata[[type]][[feature]])
    }
  } else if (!is.null(grp) & is.null(feature) & is.null(clusters)) {
    dats$colors = as.factor(getdtcols(MigrObj, dat.slot = dat.slot, type = type, vars = grp)[[1]] )
  } else if (!is.null(clusters) & is.null(grp) & is.null(feature)) {
    dats$colors = factor(getclst(MigrObj, clst = "", incl.patterns = clusters, type = "S", defaults = F))
  } else dats$colors <- rep("uknwn", dim(dats)[1])

  # Facetting
  if (!is.null(facetting)) {
    if (facetting %in% getvars(MigrObj, dat.slot, type) ) {
      facet = getdtcols(MigrObj, dat.slot = dat.slot, type = type, vars = facetting)[[1]]
      if (length(unique(facet))<43) {
        facet = as.factor(facet)
      } else if (is.numeric(facet)) {
        facet = as.factor(ggplot2::cut_interval(facet,n = fcut) )
      } else stop("Incompatible variable for facetting")
    } else if (facetting %in% colnames(MigrObj@metadata[["S"]])) {
      facet = MigrObj@metadata[[type]][[facetting]]
      if (length(unique(facet))<43) {
        facet = as.factor(facet)
      } else if (is.numeric(facet)) {
        facet = as.factor(ggplot2::cut_interval(facet,n = fcut) )
      } else stop("Incompatible variable for facetting")
    }
    if (nlevels(facet) < 43) {
      dats$facet = facet
    } else stop("Too many levels for facetting")
  }

  xrng = c( floor(min(dats$POSITION_X)) -1, ceiling(max(dats$POSITION_X))+1 )
  yrng = c( floor(min(dats$POSITION_Y)) -1, ceiling(max(dats$POSITION_Y))+1 )


  ncolorz = nlevels(factor(dats$colors))

  if (ncolorz < 20) {
    colorz = pals::brewer.accent(ncolorz)
  } else colorz = pals::viridis(ncolorz)

  plt1 <- ggplot2::ggplot(dats, ggplot2::aes(x=POSITION_X, y=POSITION_Y, colour = colors, group = grp)) +
    ggplot2::geom_point(size = point.size) + ggplot2::geom_path() +
    ggplot2::scale_color_manual(values = colorz) +
    ggplot2::scale_x_continuous(limits = xrng) + ggplot2::scale_y_continuous(limits = yrng)  +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = 'black', colour = 'grey'),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank() )

  if (!legend) {
    plt1 <- plt1 + ggplot2::theme(legend.position = "none")
  }

  if (!is.null(facetting)) {
    plt1 <- plt1 + ggplot2::facet_wrap(~ facet, ncol = 3)
  }

  # Name of the plot...
  customtitle = 1
  if (customtitle == 1) {
    title = paste("Spots colored with", feature)
    plt1 <- plt1 + ggplot2::ggtitle(title)
  }

  if (saves) {
    width = ((3*(xrng[2]-xrng[1])/100)+0.2)
    height = ((3*(yrng[2]-yrng[1])/100)+0.2)

    ggplot2::ggsave(path = "../plots/", filename = "Tracks_Plot.png", plot = plt1, device = NULL,
                    width = width, height = height, units = "in", dpi = 300)
  }
  return(plt1)
}

# ##### Spots plotting simple #####
SpotzPlot <- function(dats = spotz_new, grp = "TRACK_ID2",
                      feature = "POSITION_T", clustering = NULL,
                      facetting = NULL, point.size = 0.3) {

  dats$colors = as.factor(dats[[feature]])
  dats$grp = dats[[grp]]
  xrng = c( floor(min(dats$POSITION_X)) -1, ceiling(max(dats$POSITION_X))+1 )
  yrng = c( floor(min(dats$POSITION_Y)) -1, ceiling(max(dats$POSITION_Y))+1 )

  if (!is.null(facetting)) {
    facet = as.factor(dats[[facetting]])
    if (nlevels(facet) < 30) {
      dats$facet = facet
    } else stop("Too many levels for facetting")
  }
  ncolorz = nlevels(dats$colors)
  colorz = pals::brewer.accent(ncolorz)
  plt1 <- ggplot2::ggplot(dats, ggplot2::aes(x=POSITION_X, y=POSITION_Y, colour = colors, group = grp)) +
    ggplot2::geom_point(size = point.size) + ggplot2::geom_path() +
    ggplot2::scale_color_manual(values = colorz) +
    ggplot2::scale_x_continuous(limits = xrng) +
    ggplot2::scale_y_continuous(limits = yrng)  +
    ggplot2::ggtitle(title)  +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = 'black', colour = 'grey'),
                   panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank())
  if (!legend) {
    plt1 <- plt1 + ggplot2::theme(legend.position = "none")
  }

  if (!is.null(facetting)) {
    plt1 <- plt1 + ggplot2::facet_wrap(~ facetting, ncol = 3)
  }
  # Name of the plot...
  customtitle = 1
  if (customtitle == 1) {
    title = paste("Spots colored with", feature)
    plt1 <- plt1 + ggplot2::ggtitle(title)
  }
  print(plt1)
}


# ##### Spots clusters plotting v1 #####
SpotzClstPlot <- function(dats = spotz_new, grp = "TRACK_ID2",
                          feature = "POSITION_T", clustering = NULL,
                          legend = FALSE,
                          facetting = NULL, point.size = 0.3) {

  available_feats <- colnames(dats)

  if (length(clustering) == dim(dats)[1]) {
    dats[["clusts"]] = clustering
  } else stop("Dimensions did not match.")

  if (is.null(clustering) & !is.null(feature)) {
    dats$colors = as.factor(dats[[feature]])
  } else if (!is.null(clustering) & is.null(feature)) {
    dats$colors = as.factor(dats[["clusts"]])
  }

  dats$grp = dats[[grp]]
  xrng = c( floor(min(dats$POSITION_X)) -1, ceiling(max(dats$POSITION_X))+1 )
  yrng = c( floor(min(dats$POSITION_Y)) -1, ceiling(max(dats$POSITION_Y))+1 )

  if (!is.null(facetting)) {
    if (facetting %in% available_feats) {
      facet = as.factor(dats[[facetting]])
      if (nlevels(facet) < 30) {
        dats$facet = facet
      } else stop("Too many levels for facetting")
    } else if (!facetting %in% available_feats & !is.null(clustering)) {
      facet = as.factor(dats[["clusts"]])
      dats$facet = facet
    } else dats$facet = "0"
  }
  ncolorz = nlevels(dats$colors)
  colorz = pals::brewer.accent(ncolorz)
  plt1 <- ggplot2::ggplot(dats, ggplot2::aes(x=POSITION_X, y=POSITION_Y, colour = colors, group = grp)) +
    ggplot2::geom_point(size = point.size) + ggplot2::geom_path() +
    ggplot2::scale_color_manual(values = colorz) +
    ggplot2::scale_x_continuous(limits = xrng) +
    ggplot2::scale_y_continuous(limits = yrng)  +
    ggplot2:: ggtitle(title)  +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = 'black', colour = 'grey'),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank())
  if (!legend) {
    plt1 <- plt1 + ggplot2::theme(legend.position = "none")
  }
  if (!is.null(facetting)) {
    plt1 <- plt1 + ggplot2::facet_wrap(~ facet, ncol = 3)
  }
  # Name of the plot...
  customtitle = 1
  if (customtitle == 1) {
    title = paste("Spots colored with", feature, "or", facetting)
    plt1 <- plt1 + ggplot2::ggtitle(title)
  }
  return(plt1)
}

#### UMAPs ####

# ##### 2D umap plot function with save option #####
UMAP2DplotF <- function(dats = fact.spots, x = "x", palette = colrz.clust, filename = NULL) {

  if (length(palette) < length(levels(dats[[x]])) ) {
    palette = pals::brewer.accent(length(levels(dats[[x]])))
  }

  if (!is.null(filename)) {
    png(filename = filename, width = 1600, height = 900, units = "px", pointsize = 6)
  }
  print(
    ggplot2::ggplot(data = dats, mapping = ggplot2::aes_string(x = "UMAP1", y = "UMAP2", color = x)) +
      ggplot2::geom_point(size=2.5) +
      #scale_color_gradient2(low = "white", mid = "darkorange", high = "black") +
      ggplot2::scale_color_manual(values = palette) +
      ggplot2::theme(legend.position = "bottom",
            legend.text=element_text(size=16) )
  )
  if (!is.null(filename)) {
    dev.off()
  }
}

##### 2D dim plot function with save option #####
DimplotF <- function(dats = fact.spots, x = "x", palette = colrz.clust,
                     dimreduction = "UMAP", dims = c(1,2),
                     point.size = 3, legend.text.size = 12,
                     filename = NULL) {


  if (length(palette) < length(levels(dats[[x]])) ) {
    palette = pals::brewer.accent(length(levels(dats[[x]])))
  }

  if (!is.null(filename)) {
    png(filename = filename, width = 1600, height = 900, units = "px", pointsize = point.size)
  }

  if (grepl("UMAP", dimreduction) ) {
    dims = paste0("UMAP", dims)
    print(
      ggplot2::ggplot(data = dats, mapping = ggplot2::aes_string(x = dims[1], y = dims[2], color = x)) +
        ggplot2::geom_point(size=point.size) +
        ggplot2::scale_color_manual(values = palette) +
        ggplot2::theme(legend.position = "right",
                       legend.text=ggplot2::element_text(size=legend.text.size) )
    )
  }
  if (grepl("PCA",dimreduction) ) {
    dims = paste0("PCA", dims)
    print(
      ggplot2::ggplot(data = dats, mapping = ggplot2::aes_string(x = dims[1], y = dims[2], color = x)) +
        ggplot2::geom_point(size=point.size) +
        ggplot2::scale_color_manual(values = palette) +
        ggplot2::theme(legend.position = "right",
                       legend.text=ggplot2::element_text(size=legend.text.size) )
    )
  }
  if (!is.null(filename)) {
    dev.off()
  }
}

##### 2D umap plot saving function #####
UMAP2DplotC <- function(filename = "filename", dats = fact.spots, x = "x", palette = palette) {

  png(filename = filename, width = 1600, height = 900, units = "px", pointsize = 6)
  print(
    ggplot2::ggplot(data = dats, mapping = ggplot2::aes_string(x = "UMAP_1", y = "UMAP_2", color = x)) +
      ggplot2::geom_point(size=2.5) +
      ggplot2::scale_color_gradient2(low = "white", mid = "darkorange", high = "black") +
      ggplot2::theme(legend.position = "bottom",
                     ggplot2::legend.text=element_text(size=16) )
  )
  dev.off()
}

#### Utils ####
# Factorize set of columns
factorizeD <- function(datnx = datn1, dat.cols = NULL,
                       vaxlst = c("SAMPLE","Exp","^spdc","xcam","kmeans","hclust","Ilo"),
                       speeds = c(2:6), sort = FALSE, ...) {

  if (is.null(dat.cols)) {
    dat.cols = clstcols(datn = datnx, vaxlst = vaxlst)[[2]]
  }
  # variables to factorize
  vars = Filter(Negate(is.na), names(datnx)[dat.cols])
  if (sort) {
    vars = sort(vars)
  }
  for (var in vars) {
    datnx[[var]] = as.factor(datnx[[var]])
  }
  # add speeds
  if (is.numeric(speeds) ) {
    datnx = SpeedCuts(datn1 = datnx, datxy2 = NULL, speeds = speeds)
  }
  return(datnx)
}

# Select data for clusters and speed cats
clstcols <- function(datn = dataX, vaxlst = c("^spdc","^kmeans","^hclust","^Ilo"), verbose = FALSE) {
  # colsz = NULL
  varsz = list(vars = NULL, varidxs = NULL)
  for (var in vaxlst) {
    vars = grep(var, names(datn), value = T)
    varidxs = grep(var, names(datn), value = F)

    varsz$vars = unique(c(varsz$vars, vars))
    varsz$varidxs = unique(c(varsz$varidxs, varidxs))

    # print(varidxs)
    # colsz = c(colsz,varidxs)
  }
  if (verbose) {
    print(varsz)
  }
  return(varsz)
}

#### Add speed classes to data ####
SpeedCuts <- function(datn1 = tbl3, datxy2 = NULL, speeds = c(2:6)) {

  vars = names(datn1)

  if (!is.numeric(speeds) ) {
    print("Speed categories need to be numeric!")
    break()
  }

  if (is.null(datxy2)) {
    for (n in speeds) {
      if (exists("trckspeed", datn1)) {
        spdcatT = paste0("spdcatT",n)
        datn1[[spdcatT]] = ggplot2::cut_interval(datn1[["trckspeed"]], n = n)
      }

      if (exists("edgspeed", datn1)) {
        spdcatS = paste0("spdcatE",n)
        datn1[[spdcatS]] = ggplot2::cut_interval(datn1[["edgspeed"]], n = n)
      }

      if (exists("edgspdRaw", datn1)) {
        spdcatS = paste0("spdcatS",n)
        datn1[[spdcatS]] = ggplot2::cut_interval(datn1[["edgspdRaw"]], n = n)
      }

      if (exists("TRACK_MEAN_SPEED", datn1)) {
        spdcatS = paste0("spdcatTrkMn",n)
        datn1[[spdcatS]] = ggplot2::cut_interval(datn1[["TRACK_MEAN_SPEED"]], n = n)
      }
      if (exists("trkspdRaw", datn1)) {
        spdcatS = paste0("spdcatTrkMn",n)
        datn1[[spdcatS]] = ggplot2::cut_interval(datn1[["trkspdRaw"]], n = n)
      }
    }
    return(datn1)
  }

  if (!is.null(datxy2)) {

    for (n in speeds) {
      if (exists("trckspeed", datn1)) {
        spdcatT = paste0("spdcatT",n)
        datxy2[[spdcatT]] = ggplot2::cut_interval(datn1[["trckspeed"]], n = n)
      }

      if (exists("edgspeed", datn1)) {
        spdcatS = paste0("spdcatE",n)
        datxy2[[spdcatS]] = ggplot2::cut_interval(datn1[["edgspeed"]], n = n)
      }

      if (exists("edgspdRaw", datn1)) {
        spdcatS = paste0("spdcatS",n)
        datxy2[[spdcatS]] = ggplot2::cut_interval(datn1[["edgspdRaw"]], n = n)
      }

      if (exists("TRACK_MEAN_SPEED", datn1)) {
        spdcatS = paste0("spdcatTrkMn",n)
        datxy2[[spdcatS]] = ggplot2::cut_interval(datn1[["TRACK_MEAN_SPEED"]], n = n)
      }

      if (exists("trkspdRaw", datn1)) {
        spdcatS = paste0("spdcatTrkMn",n)
        datxy2[[spdcatS]] = ggplot2::cut_interval(datn1[["trkspdRaw"]], n = n)
      }

    }
    return(datxy2)
  }
}

