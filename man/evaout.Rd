\docType{data}
\name{evaout}
\alias{evaout}
\title{The output of one EVA run at alpha=0.8}
\format{For convenience, we provide the output of one EVA algorithm run
  on the Cd4+ T-cell scRNA-Seq dataset at alpha=0.8. The format is as
  specified by the output of the \code{Eva_partitions} function.
}

\usage{
data(evaout);
}
\description{
  A list object containing the community memberships of all cells, the
  estimated purity and modularity values.
}

\details{
Upon loading the data object, we find the following:
\itemize{
  \item \code{eva.o}: the output from EVA-clustering algorithm
}
}
\references{

}
\keyword{datasets}
