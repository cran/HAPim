\name{distance.test}

\alias{distance.test}

\title{localisation of test positions}

\description{
The function gives the distance of test positions from  the first marker and sorts them out by marker interval.
It can be viewed as an internal function. 
The user does not have to call it by himself.
}

\usage{
distance.test(position, dist.marq)
}


\arguments{
  \item{position}{numeric vector of test positions.}
  \item{dist.marq}{results provided by \code{distance.marqueurs()} function, numeric vector 
of length=number of markers which gives the distance of each marker from the bigining of the chromosome.}
}


\value{
the returned value  is a list of n numeric vectors (n=number of markers-1). The $j$th vector contains the distance (from the first marker) of test positions belonging to the $j$th marker interval.
}


\references{publication to be submitted: C. Cierco-Ayrolles, S. Dejean, A. Legarra, H. Gilbert,
T. Druet, F. Ytournel, D. Estivals, N. Oumouhou and B. Mangin.
Combining linkage analysis and linkage disequilibrium for QTL fine mapping in animal pedigrees.}


\author{S. Dejean, N. Oumouhou, D. Estivals, B. Mangin }

\seealso{\code{\link{distance.marqueurs}} }

\keyword{models}

 
