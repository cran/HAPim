\name{depart.LD}

\alias{depart.LD}

\title{starting values for the optimization of HAPimLD method}

\description{
the function calculates the starting values of the performance mean and the error variance for the optimization of HAPimLD method.
It can be viewed as an internal function. 
The user does not have to call it by himself.
}

\usage{
depart.LD(perf, CD)
}

\arguments{
  \item{perf}{numeric vector of length=number of individuals which contains the performances of individuals.}
  \item{CD}{numeric vector of length=number of individuals which contains the CD of individuals.}
}

\value{
The value returned is a numeric vector of length=2 which contains estimates of the performance mean and the error variance.
}

\references{publication to be submitted: C. Cierco-Ayrolles, S. Dejean, A. Legarra, H. Gilbert,
T. Druet, F. Ytournel, D. Estivals, N. Oumouhou and B. Mangin.
Combining linkage analysis and linkage disequilibrium for QTL fine mapping in animal pedigrees.}

\author{S. Dejean, N. Oumouhou, D. Estivals, B. Mangin }

\keyword{models}


