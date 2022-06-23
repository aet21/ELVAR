\docType{data}
\name{purmodInfo}
\alias{purmodInfo}
\title{Purity and modularity values of 100 EVA runs for different
  choices of alpha parameter}
\format{Contains two matrices, one for purity and another for
  modularity, each defined over 100 EVA-runs (rows) and for nine
  different choices of alpha parameter (columns).
}

\usage{
data(purmodInfo)
}
\description{
 The purity and modularity values were estimated by running the EVA
 clustering algorithm on the lung-tissue Cd4+ T-cell scRNA-Seq dataset,
 using age as the attribute. EVA was run a 100 times for each of nine
 different choices of alpha (0.1,0.2,...0.9). The higher alpha is the
 more purity contributes to the optimization function.
}

\details{
Upon loading the data object, we find the following:
\itemize{
  \item \code{pur.m}: the purity matrix
  \item \code{mod.m}: the modularity matrix
}
}
\references{

}
\keyword{datasets}
