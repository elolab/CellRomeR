#### Main clustering function ####

IRmigr.clustering <- function(MigrObj, dat.slot = "scaled", type = "STE",
                              uniq = "base",
                              kmeans = 0, khclust = 0, kILoReg = 0,
                              predef =  "none",
                              vars = NULL, incl.pattern = NULL,
                              excld.pattern = NULL,
                              scale = FALSE, center = FALSE,
                              set.default = TRUE, threads = 0,...) {
  startTime <- Sys.time()


  # set thread count
  if (threads==0) {
    threads=parallel::detectCores()
  }

  predef = match.arg(predef, c("none", "technical", "morphological", "clust"), several.ok = FALSE)

  # standard variable names we are not interested in
  StdP = Std.pattern(MigrObj)

  data <- getdtcols(MigrObj, dat.slot = dat.slot, type = type,
                    predef = predef,
                    vars = vars, incl.pattern = incl.pattern,
                    excld.pattern = excld.pattern, numerics = TRUE)
  #data[,STE:= type]

  # dimreducts <- getdtcols(MigrObj, dat.slot = dat.slot, type = type,
  #                         predef = "coord", numerics = FALSE)

  # if (is.null(MigrObj@clustering[[uniq]])) {
  #   MigrObj@clustering[[uniq]] <- list()
  # }

  # kmean
  if (kmeans[1]>1) {
    kmeansclsts <- IRMigr.kmeansX(data, kmeans = kmeans, scale = scale)

    for (nm in names(kmeansclsts)) {
      kmeansname = paste0(nm, "_", type,"_",uniq)
      MigrObj@clustering[[type]][[kmeansname]] <- kmeansclsts[[nm]]
    }
    if (set.default) {
      MigrObj@clustering[[type]][["kmeans"]] <- kmeansclsts[[nm]]
      default.kmeans(MigrObj)[[type]] <- kmeansname
    }

  } else kmeansclsts = list()

  # Hierarchical hclust
  if (khclust[1]>1) {
    hclusts <- IRMigr.hclust(data, khclust = khclust, scale = scale)
    for (nm in names(hclusts)) {
      hclustname = paste0(nm, "_", type,"_",uniq)
      MigrObj@clustering[[type]][[hclustname]] <- hclusts[[nm]]
    }
    if (set.default) {
      MigrObj@clustering[[type]][["hclusts"]] <- hclusts[[nm]]
      default.kmeans(MigrObj)[[type]] <- hclustname
      #
    }
  } else hclusts = list()


  # ILoReg
  if (kILoReg[1]>1) {
    mtx = as.matrix(data)
    rownames(mtx) <- getlabels(MigrObj, dat.slot = dat.slot, type = type)
    ILoRegclsts <- IRMigr.ILoReg(mtx, kILoReg = kILoReg, type = type, threads = threads, scale = scale, ...)

    for (nm in names(ILoRegclsts)) {
      ILoRegname = paste0(nm, "_", type,"_",uniq)
      MigrObj@clustering[[type]][[ILoRegname]] <- ILoRegclsts[[nm]]
    }
    if (set.default) {
      MigrObj@clustering[[type]][["ILoRegclusters"]] <- ILoRegclsts[[nm]]
      default.kmeans(MigrObj)[[type]] <- ILoRegname
    }
  } else ILoRegclsts = list()

  cat("\nClustering results were stored with ", uniq," identifier in clusterings slot.")

  timing=Sys.time() - startTime
  # cat("\nClustering took: ",timing,"\n\n" )
  cat("\nClustering took: \n" )
  cat(Sys.time() - startTime,"\n\n")
  #return(MigrObj)
  return(MigrObj)

}

#### Individual clustering functions ####

IRMigr.ILoReg <- function(mtx, kILoReg, scale=TRUE, seed = 1917, L = 50,
                             K = 15, d = 0.3, r = 500, C = 0.3,
                             icp.batch.size = 4000,
                             reg.type = "L1", threads = 0, type = "STE") {

  if (scale) {
    mtx = scale(mtx)
  }

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
    scemg <- ILoReg::RunParallelICP(object = scemg, k = 10,
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


  for (k in kILoReg) {
    scemg <- ILoReg::SelectKClusters(scemg, K=k)
    ilonm = paste0("IloRegK",k)

    ILoRegClsts[[ilonm]] = scemg@metadata$iloreg$clustering.manual
  }
  return(ILoRegClsts)
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
  # h_euclidean = hclust(d, method = "ward.D2")
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
