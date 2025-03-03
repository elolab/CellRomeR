\name{plot_umap}
\alias{plot_umap}
\title{Plot a UMAP plot}
\description{This function plots a UMAP from clustered migration data object.
}

\usage{
plot_umap(MigrObj, uniq="ILoReg", color.by="ILoRegclusters",
spots.slot="raw", cex=0.1,use.opinionated.layout = TRUE)
}

\arguments{
\item{MigrObj}{MigrDat object, that has been through \code{\link{clustering}}}
\item{uniq}{A label given to it by the \code{uniq} parameter in the \code{\link{clustering}} method}
\item{color.by}{Indicates the field to color the dots by, either a
  column name of the trackmate spots, edges, tracks data, or one of
   \code{ILoReg<K>_<uniq>_<type>}
  where \code{<K>} is a kILoReg parameter given to
  \code{\link{clustering}}, and \code{<uniq>} and \code{<type>} are as
  in that function.
  Alterinatively, leave to the default to use the value of clustering
  that was set as the default.
}
\item{spots.slot}{Spots slot data,as in dat.slot in \code{\link{clustering}}, currently should be left unchanged}
\item{cex}{Point size of dots}
\item{use.opinionated.layout}{Whether or not to use the more spaceconsuming layout}
}

\examples{
 \dontrun{
  migob <- import_XML("file.XML") |>
           Apply_TMate_Filter() |>
           clustering(predef = "morphological")
  plot_umap(migob)}
}









