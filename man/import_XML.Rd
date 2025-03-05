\name{import_XML}
\alias{import_XML}
\title{Load TrackMate XML exported files}
\description{This function reads in the XML exported files from TrackMate, and converts it to various `data.table` which are split by the types i.e "Spots", "Tracks" or "Edges".
}

\usage{
import_XML(tmXML, dataset = "MigrDatTestXML",
          experiment = "experiment", sample = "sample",
          condition = "normal", replicate = NULL, MigrDatObj = NULL)
}

\arguments{
\item{tmXML}{Input XML file from TrackMate export, or an
  \code{http:}/\code{ftp:}/\code{file:} url pointing to such a file.
  \code{https:} is
  not supported by R's \code{XML}, so just leave out the \code{s}}
\item{dataset}{The name to give the dataset within the `MigrDat` object.
One `MigrDat` can hold multiple datasets}
\item{MigrDatObj}{Optional `MigrDat` object to add the dataset to.}
\item{sample}{Optional name for sample}
\item{condition}{Optional name for condition}
\item{experiment}{Optional name for experiment}
\item{replicate}{Optional name for replicate}
}

\details{
\subsection{Optional Metadata Parameters}{
The parameters \code{sample},\code{condition},\code{experiment} and \code{replicate} are
labels you can give to describe the data. It doesn't do much useful at the moment, but it's honest work.
}
}

\examples{
## Import the "file_path.xml" file
\dontrun{MigrObject <- import_XML(tmXML = "file_path.xml")}
}

\keyword{file}









