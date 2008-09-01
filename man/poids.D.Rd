\name{poids.D}

\alias{poids.D}

\title{poids.D}

\description{
This function calculates the probability of no recombinaison between loci for each Bennett's desequilibrium.
The function can be viewed as an internal function. 
The user does not have to call it by himself.
}

\usage{
poids.D(dist.test, pos.QTL, res.structure)
}


\arguments{
  \item{dist.test}{results provided by \code{distance.test()} function, list of numeric objects. }
  \item{pos.QTL}{numeric value, interval of a test position}
  \item{res.structure}{results provided by \code{structure.hap()} function, list of objects.}
}


\value{
The returned  value is a list of numeric objects.
}


\references{publication to be submitted: C. Cierco-Ayrolles, S. Dejean, A. Legarra, H. Gilbert,
T. Druet, F. Ytournel, D. Estivals, N. Oumouhou and B. Mangin.
Combining linkage analysis and linkage disequilibrium for QTL fine mapping in animal pedigrees.}




\author{S. Dejean, N. Oumouhou, D. Estivals, B. Mangin }

\seealso{\code{\link{structure.hap}}, \code{\link{distance.test}} }

\keyword{models}
