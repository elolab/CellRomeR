\name{clustering}
\alias{clustering}
\title{Clustering of Trackmate Data}

\description{Invoke ILoReg clustering on TrackMate data}

\usage{
clustering(MigrObj, dat.slot = "raw", type = c("S","T","E"),
            uniq = "ILoReg", kILoReg = 5, predef =  c("none", "morphological", "morphplus", "technical"),
            vars = NULL, incl.pattern = NULL, excld.pattern = NULL,
            scale = TRUE, set.default = TRUE, threads = 0,\dots)
}


\arguments{
  \item{MigrObj}{A MigrDat object}
  \item{dat.slot}{Future compatability parameter, currently ignored.}
  \item{type}{One of S,T,E for \emph{S}pots, \emph{T}racks, or
  \emph{E}dges, which to cluster one.}
  \item{uniq}{An Identifier for this clustering for later functions.}
  \item{kILoReg}{Initial number of cluster for ILoReg, can be length \gt 1}
  \item{predef}{Set to one of \code{morphological},\code{morphplus},
  \code{technical} to select these clustering, see see \dQuote{TrackMate Column Selection}}
  \item{vars}{A character vector of Trackmate data columnnames to select
  for clustering, see \dQuote{TrackMate Column Selection}}
  \item{incl.pattern}{a character vector of regular expressions, determining what columns to use for clustering }
  \item{excld.pattern}{character vector of regular expression of which
  columns to \emph{not} use for clustering, see see \dQuote{TrackMate Column Selection}}
  \item{scale}{Whether or not to \code{\link{scale}} the data before clustering.}
  \item{set.default}{Whether or not to set the default clustering to this result. (default TRUE)}
  \item{threads}{Number of threads (default 0, i.e. all available)}
  \item{\dots}{Additional arguments passed to \code{\link[ILoReg]{RunParallelICP}}
\details{
  \subsection{TrackMate Column Selection}
  The arguments \code{predef},\code{inc.pattern},
  \code{excl.pattern}, \code{vars} all determine the columns of
  trackate data to consider for clustering.
}
\value{
  The input MigrObj with its \code{clustering} fields set to the result
  of ILoreg clustering, based on the selected columns.
}  

\seealso{
  \code{\link{scale}}
  
  \code{\link[ILoReg]{RunParallelICP}}
}
      





