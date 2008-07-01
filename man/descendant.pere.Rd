\name{descendant.pere}

\alias{descendant.pere}

\title{The first and last indexes of individuals per sire}

\description{
The function defines, for each sire, the first and last indexes of its progeny.
It can be viewed as an internal function. 
The user does not have to call it by himself.
}

\usage{
descendant.pere(genea)
}

\arguments{
  \item{genea}{numeric matrix (number of individuals*2) which contains 
individual index and corresponding sire index of each individual.}
}

\value{
 The value returned is a numeric matrix (number of sires*2) which gives for each sire, the first and last indexes of its progeny.
}

\references{publication to be submitted:C.Cierco-Ayrolles,S.Dejean,A.Legarra,H.Gilbert,
T.Druet,F.Ytournel,D.Estivals,N.Oumouhou and B.Mangin,
Combining linkage analysis and linkage disequilibrium for QTL fine mapping in animal pedigrees.}

\author{S. Dejean, N. Oumouhou, D. Estivals, B. Mangin }

\keyword{models}

