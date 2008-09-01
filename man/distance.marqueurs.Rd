\name{distance.marqueurs}

\alias{distance.marqueurs}

\title{distance of markers}

\description{
The function gives the distance of each marker to the first marker.
It can be viewed as an internal function. 
The user does not have to call it by himself.
}

\usage{
distance.marqueurs(map)
}

\arguments{
  \item{map}{numeric vector of length=(number of markers-1) giving the distance 
between two consecutive markers on all the chromosome.}
}

\value{
the returned   valueis a numeric vector of length=number of markers containing the positions of each marker from the beginning of the chromosome.
}



\references{publication to be submitted: C. Cierco-Ayrolles, S. Dejean, A. Legarra, H. Gilbert,
T. Druet, F. Ytournel, D. Estivals, N. Oumouhou and B. Mangin.
Combining linkage analysis and linkage disequilibrium for QTL fine mapping in animal pedigrees.}




\author{S. Dejean, N. Oumouhou, D. Estivals, B. Mangin }

\keyword{models}
