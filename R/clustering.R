#### Main clustering function ####

clustering <- function(MigrObj, dat.slot = "raw", type = c("S","T","E"),
                       clust.label = "base",
                       kILoReg = 5,
                       predef =  c("none", "technical", "morphological", "morphplus", "clust"),
                       vars = NULL, incl.pattern = NULL,
                       excld.pattern = NULL,
                       scale = TRUE, 
                       set.default = TRUE, threads = 0,...) {
  stopifnot(kILoReg > 1)
  type <- match.arg(type)
  predef <- match.arg(predef)
  
  # set thread count
  if (threads==0) {
    threads=parallel::detectCores()
  }
  
  # standard variable names we are not interested in
  StdP = CellRomeR::Std.pattern(MigrObj)
  
  data <- CellRomeR::getdtcols(MigrObj, dat.slot = dat.slot, type = type,
                                predef = predef,
                                vars = vars, incl.pattern = incl.pattern,
                                excld.pattern = excld.pattern, numerics = TRUE)
  
  # ILoReg
  mtx = as.matrix(data)
  rownames(mtx) <- CellRomeR::getlabels(MigrObj, dat.slot = dat.slot, type = type)
  scemg <- CellRomeR::cRomeR.ILoReg(mtx, kILoReg = kILoReg, type = type, threads = threads, scale = scale, ...)
  
  ILoRegClsts = list()
  for (k in kILoReg) {
    scemg <- ILoReg::SelectKClusters(scemg, K=k)
    ilonm = paste0("IloRegK",k)
    ILoRegClsts[[ilonm]] = scemg@metadata$iloreg$clustering.manual
  }
  
  l <- SingleCellExperiment::reducedDims(scemg)
  MigrObj@dimreductions[[type]][[paste(clust.label, "_UMAP")]] <- l$UMAP
  MigrObj@dimreductions[[type]][[paste(clust.label, "_TSNE")]] <- l$TSNE
  
  for (nm in names(ILoRegClsts)) {
    ILoRegname = paste0(nm, "_", type,"_",clust.label)
    MigrObj@clustering[[type]][[ILoRegname]] <- ILoRegClsts[[nm]]
  }
  if (set.default) {
    MigrObj@clustering[[type]][["ILoRegclusters"]] <- ILoRegClsts[[nm]]
    default.ILoRegclusters(MigrObj)[[type]] <- ILoRegname
  }
  
  return(MigrObj)
  
}

#### Individual clustering functions ####

cRomeR.ILoReg <- function(mtx, kILoReg, scale=TRUE, seed = 1917, L = 50,
                          K = 15, d = 0.3, r = 500, C = 0.3,
                          icp.batch.size = 4000,
                          reg.type = "L1", threads = 0, type = c("S","T","E")) {
  # Checks
  type <- match.arg(type)
  
  # Data  
  if (scale) {
    mtx = scale(mtx)
  }
  
  vars = colnames(mtx)
  
  # Get rid of sce-object later
  scemg = SingleCellExperiment::SingleCellExperiment(assays = list(logcounts = t(mtx)))
  scemg <- ILoReg::PrepareILoReg(scemg)
  scemg@metadata$iloreg$vars = vars
  
  if (is.numeric(seed)) {
    set.seed(seed)
  }
  
  if (type %in% c("S","E")) {
    scemg <- ILoReg::RunParallelICP(object = scemg, k = K,
                                    d = d, L = L,
                                    r = r, C = C,
                                    icp.batch.size = icp.batch.size,
                                    reg.type = reg.type, threads = threads)
    
  } else if (type == "T") {
    scemg <- ILoReg::RunParallelICP(object = scemg, k = K,
                                    d = d, L = L,
                                    r = r, C = C,
                                    icp.batch.size = Inf,
                                    reg.type = reg.type, threads = threads)
    
  } else {
    scemg <- ILoReg::RunParallelICP(object = scemg, k = K,
                                    d = d, L = L,
                                    r = r, C = C,
                                    icp.batch.size = icp.batch.size,
                                    reg.type = reg.type, threads = threads)
  }
  
  scemg <- ILoReg::RunPCA(scemg, p=50, scale = FALSE)
  scemg <- ILoReg::HierarchicalClustering(scemg)
  scemg <- ILoReg::RunUMAP(scemg)
  scemg <- ILoReg::RunTSNE(scemg)
  
  return (scemg)
  
}

cRomeR.kmeansX <- function(mtx, kmeans, scale = TRUE) {
  
  if (scale) {
    mtx = scale(mtx)
  }
  kmeansclst = list()
  for (k in kmeans) {
    if (k>0) kmeansdat = stats::kmeans(mtx, centers = k, nstart = 15, iter.max = 2500, algorithm="MacQueen")
    knms = paste0("kmeans",k)
    kmeansclst[[knms]] = kmeansdat$cluster
  }
  return(kmeansclst)
}

cRomeR.kmeans <- function(mtx, k, scale=TRUE) {
  
  if (scale) {
    mtx = scale(mtx)
  }
  if (k>0 & !scale) {
    kmeansdat = stats::kmeans(mtx, centers = k, nstart = 15, iter.max = 2500, algorithm="MacQueen")  #
    knms = paste0("kmeans",k)
    return(kmeansdat$cluster)
  } else return(NULL)
}

cRomeR.hclust <- function(mtx, khclust, scale=TRUE) {
  
  if (scale) {
    mtx = scale(mtx)
  }
  
  d = dist(mtx, method = "euclidean")
  # clustering
  h_euclidean = fastcluster::hclust(d, method = "ward.D2")
  hclusts = list()
  for (k in khclust) {
    if (k>0) {
      # euclidean distances
      hclustk = cutree(h_euclidean, k = k)
      hcnms = paste0("hclust",k)
      hclusts[[hcnms]] = hclustk
    }
  }
  return(hclusts)
}
