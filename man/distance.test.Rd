\name{distance.test}

\alias{distance.test}

\title{test positions}

\description{
The function gives test positions compared to the first marker and ranged in two consecutive markers interval.
It can be viewed as an internal function. 
The user does not have to call it by himself.
}

\usage{
distance.test(position, dist.marq)
}


\arguments{
  \item{position}{numeric vector of test positions.
   Positions of "position" vector have to be the same as positions of "PLA" matrix.}
  \item{dist.marq}{results provided by \code{distance.marqueurs()} function, numeric vector 
of length=number of markers which gives the distance of each marker compared to the first marker.}
}


\value{
the value returned is a list of n elements (n=number of markers-1) which contains the distance 
of each marker and ranged in two consecutive markers interval.
}

\references{publication to be submitted:C.Cierco-Ayrolles,S.Dejean,A.Legarra,H.Gilbert,
T.Druet,F.Ytournel,D.Estivals,N.Oumouhou and B.Mangin,
Combining linkage analysis and linkage disequilibrium for QTL fine mapping in animal pedigrees.}

\author{S. Dejean, N. Oumouhou, D. Estivals, B. Mangin }

\seealso{\code{\link{distance.marqueurs}} }

\keyword{models}

 
