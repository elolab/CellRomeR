\name{clustering}
\title{Clusteriing of Trackmate Data}

\description{Invoke ILoReg clustering on TrackMate data}
\usage{
clustering(MigrObj, dat.slot = "raw", type = "STE",
                       uniq = "base",
                       # a numeric vector of k's
                       # iloreg will be invoked for each of its elements.
                       kILoReg = c(8),
                       predef =  c("none", "technical", "morphological", "clust"),
                       vars = NULL, incl.pattern = NULL,
                       excld.pattern = NULL,
                       scale = FALSE, 
                       set.default = TRUE, threads = 0,\dots)
}
		       




