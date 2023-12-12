\name{import_XML}
\alias{import_XML}
\title{Load TrackMate XML exported files}
\description{This function reads in the XML exported files from TrackMate, and converts it to various `data.table` which are split by the types i.e "Spots", "Tracks" or "Edges".
}
\newcommand{\optional}{\item{#1}{(Optional) Metadata, see Details.}}
\usage{
import_XML(tmXML,
           dataset = "MigrDatTestXML", 
           experiment = "experiment", sample = "S-test",
           condition = "Test", MigrDatObj = NULL)
}
\arguments{
  \item{tmXML}{Input XML file from TrackMate export}
  \item{dataset}{The name to give the dataset within the `MigrDat`
        object. One `MigrDat` can hold multiple datasets}
  \item{MigrDatObj}{WHEN IS THIS NEEDED?}
  \optional{sample}
  \optional{condition}
  \optional{experiment}
  \optional{replicate}
}
\details{
\subsection{Optional Metadata Parameters}{
The parameters \code{sample},\code{condition},\code{experiment} and \code{replicate} are
labels you can give to describe the data. It doesn't do much useful at the moment, but it's honest work.
}
}
\seealso{
  \code{\link{MigrDat}}.
}
\examples{
## Import the "file_path.xml" file
MigrObject <- CellRomeR::import_XML(tmXML = "file_path.xml", 
                           dataset = "Dataset1", 
                           experiment = "ExpT", sample = "S-test",
                           condition = "Test", MigrDatObj = NULL)
}
\keyword{file}









