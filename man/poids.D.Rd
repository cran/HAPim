\name{poids.D}

\alias{poids.D}

\title{poids.D}

\description{
The function can be viewed as an internal function. 
The user does not have to call it by himself.
}

\usage{
poids.D(dist.test, pos.QTL, res.structure)
}


\arguments{
  \item{dist.test}{results provided by \code{distance.test()} function,a list of n elements (n=number of markers-1) which contains the distance of each marker and ranged in two consecutive markers interval.}
  \item{pos.QTL}{numeric value.}
  \item{res.structure}{results provided by \code{structure()} function, list of numeric objects.}
}


\value{
the value returned is a list of numeric objects.
}

\references{publication to be submitted:C.Cierco-Ayrolles,S.Dejean,A.Legarra,H.Gilbert,
T.Druet,F.Ytournel,D.Estivals,N.Oumouhou and B.Mangin,
Combining linkage analysis and linkage disequilibrium for QTL fine mapping in animal pedigrees.}

\author{S. Dejean, N. Oumouhou, D. Estivals, B. Mangin }

\seealso{\code{\link{structure}},\code{\link{distance.test}} }

\keyword{models}
