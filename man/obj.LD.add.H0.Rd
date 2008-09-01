\name{obj.LD.add.H0}

\alias{obj.LD.add.H0}

\title{log-likelihood value of HAPimLD method under H0}

\description{
The function calculates the log-likelihood value of HAPimLD method under hypothesis HO.
It can be viewed as an internal function. 
The user does not have to call it by himself.
}

\usage{
obj.LD.add.H0(perf, CD)
}

\arguments{
  \item{perf}{numeric vector of length = number of individuals which contains the performances of individuals.}
  \item{CD}{numeric vector of length=number of individuals which contains the CD of individuals.  var(perf$_i$)=s/CD$^2_i$}
}

\value{
The  returned value is the log-likelihood value of HAPimLD method under hypothesis HO.
}

\references{publication to be submitted: C. Cierco-Ayrolles, S. Dejean, A. Legarra, H. Gilbert,
T. Druet, F. Ytournel, D. Estivals, N. Oumouhou and B. Mangin.
Combining linkage analysis and linkage disequilibrium for QTL fine mapping in animal pedigrees.}



\author{S. Dejean, N. Oumouhou, D. Estivals, B. Mangin }

\keyword{models}


