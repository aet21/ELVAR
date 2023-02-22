\docType{data}
\name{purmodInfo}
\alias{purmodInfo}
\title{Purity, modularity and cluster number values of 100 EVA runs for
  different choices of alpha purity parameter.}
\format{Contains three matrices, one for purity, one for modularity and
  one for cluster number, each defined over 100 EVA-runs (rows) and for nine
  different choices of alpha parameter (columns). In addition, it also
  provides the estimated number of clusters using deterministic Louvain.
}

\usage{
data(purmodInfo)
}
\description{
 The purity, modularity and cluster number values were estimated by
 running the EVA  clustering algorithm on the lung-tissue Cd4+ T-cell
 scRNA-Seq dataset,  using age as the attribute. EVA was run a 100 times
 for each of nine  different choices of alpha (0.1,0.2,...0.9).
 The higher alpha is the more purity contributes to the optimization function.
}

\details{
Upon loading the data object, we find the following:
\itemize{
  \item \code{pur.m}: the purity matrix
  \item \code{mod.m}: the modularity matrix
  \item \code{nC.m}: the cluster number matrix
  \item \code{nc0}: the number of inferred clusters with Louvain
}
}
\references{

}
\keyword{datasets}
