#### Main clustering function ####

clustering <- function(MigrObj, dat.slot = "raw", type = c("S", "T", "E"),
                              uniq = "ILoReg",
                              kILoReg = 5,
                              predef =  c("none","morphological","morphplus","technical"),
                              vars = NULL, incl.pattern = NULL,
                              excld.pattern = NULL,
                              scale = TRUE, 
                              set.default = TRUE, threads = 0,...) {
  stopifnot(kILoReg > 0)
  

  # set thread count
  if (threads==0) {
    threads <- parallel::detectCores()
  }

  predef <- match.arg(predef)
  type <- match.arg(type)

  # standard variable names we are not interested in
  StdP <- Std.pattern(MigrObj)
  

  data <- getdtcols(MigrObj, dat.slot = dat.slot, type = type,
                    predef = predef,
                    vars = vars, incl.pattern = incl.pattern,
                    excld.pattern = excld.pattern, numerics = TRUE)
  if(any(dim(data) == 0))
    stop(sprintf("Did not select columns from the %s table,consider changing one of `predef`, `incl.pattern`,`vars`,excld.pattern`", switch(type, S  = "Spots", E = "Edges", T = "Tracks")))
  
    

  # ILoReg
  mtx <- as.matrix(data)
  rownames(mtx) <- getlabels(MigrObj, dat.slot = dat.slot, type = type)
  scemg <- IRMigr.ILoReg(mtx, K = max(kILoReg), type = type, threads = threads, 
                         scale = scale, icp.batch.size = ceiling(nrow(mtx)/4), ...)
  ILoRegClsts <-
    Map(function(k)
      ILoReg::SelectKClusters(scemg, K=k)@metadata$iloreg$clustering.manual,
      kILoReg)
  names(ILoRegClsts) <- sprintf("IloRegK%d",kILoReg)
      
      
  # Simplification of DR-storage would be like this 
  MigrObj@dimreductions[[type]] <- SingleCellExperiment::reducedDims(scemg)
  
  l <- SingleCellExperiment::reducedDims(scemg)
  MigrObj@dimreductions[[type]][[paste0(uniq, "_UMAP")]] <- l$UMAP
  #MigrObj@dimreductions[[type]][[paste0(uniq, "_TSNE")]] <- l$TSNE

  for (nm in names(ILoRegClsts)) {
    ILoRegname = paste0(nm, "_", type,"_",uniq)
    MigrObj@clustering[[type]][[ILoRegname]] <- ILoRegClsts[[nm]]
  }
  if (set.default) {
    if(length(ILoRegClsts) > 1)
      warning(sprintf("More than one clustering, setting default to = %s", ILoRegname))
    
    MigrObj@clustering[[type]][["ILoRegclusters"]] <- ILoRegClsts[[nm]]
    default.kmeans(MigrObj)[[type]] <- ILoRegname
  }
    
  return(MigrObj)

}

#### Individual clustering functions ####

IRMigr.ILoReg <- function(mtx, scale=TRUE, seed = 1917, L = 50,
                          K = 5, d = 0.3, r = 10, C = 0.3,
                          icp.batch.size = 1000,
                          reg.type = "L1", threads = 0, 
                          type = c("S", "T", "E"),
                          # the `p` parameter for ILOreg::runPCA
                          # for testing purpopses
                          p = 25) {
  
  if (scale) {
    mtx = scale(mtx)
  }
  type <- match.arg(type)
  
  #ILoRegClsts = data.frame()
  ILoRegClsts = list()
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
  
  scemg <- ILoReg::RunPCA(scemg, p=p, scale = FALSE)
  scemg <- ILoReg::HierarchicalClustering(scemg)
  # Faster implementation of the UMAP 
  SingleCellExperiment::reducedDim(scemg, "UMAP") <- uwot::umap(X = SingleCellExperiment::reducedDim(scemg, "PCA"))
  # scemg <- ILoReg::RunUMAP(scemg)
  # scemg <- ILoReg::RunTSNE(scemg)
  
  return (scemg)
  
}


IRMigr.kmeansX <- function(mtx, kmeans, scale = FALSE) {

  if (scale) {
    mtx = scale(mtx)
  }
  kmeansclst = list()
  for (k in kmeans) {
    if (k>0) kmeansdat = kmeans(mtx, centers = k, nstart = 15, iter.max = 2500, algorithm="MacQueen")
    knms = paste0("kmeans",k)
    kmeansclst[[knms]] = kmeansdat$cluster
  }
  return(kmeansclst)
}

IRMigr.kmeans <- function(mtx, k, scale=FALSE) {

  if (scale) {
    mtx = scale(mtx)
  }
  if (k>0 & !scale) {
    kmeansdat = kmeans(mtx, centers = k, nstart = 15, iter.max = 2500, algorithm="MacQueen")  #
    knms = paste0("kmeans",k)
    return(kmeansdat$cluster)
  } else return(NULL)
}

IRMigr.hclust <- function(mtx, khclust, scale=FALSE) {

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
