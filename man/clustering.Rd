\name{clustering}
\alias{clustering}
\title{Clustering of Trackmate Data}

\description{Invoke ILoReg clustering on TrackMate data}

\usage{
clustering(MigrObj, dat.slot = "raw", type = c("S","T","E"),
            uniq = "IloReg", kILoReg = 0, predef =  c("none", "morphological", "morphplus", "technical"),
            vars = NULL, incl.pattern = NULL, excld.pattern = NULL,
            scale = TRUE, set.default = TRUE, threads = 0,\dots)
}

\arguments{
	\item{MigrObj}{A MigrDat object}
	\item{dat.slot}{Data slot}
	\item{type}{Type}
	\item{uniq}{Identifier for clustering}
	\item{kILoReg}{Initial number of cluster for ILoReg}
	\item{predef}{predef}
	\item{vars}{vars}
	\item{incl.pattern}{Inclusion pattern}
	\item{excld.pattern}{Exclusion pattern}
	\item{scale}{A logical (deafult TRUE) indicating scaling}
	\item{set.default}{A logical (deafult FALSE)}
	\item{threads}{Number of threads (deafult 0, i.e. all available)}
	\item{\dots}{Additional items passed to ILoReg}
}





