\name{obj.LD.add}

\alias{obj.LD.add}

\title{log-likelihood value of HAPimLD method under H1}

\description{
The function calculates the log-likelihood value of HAPimLD method under alternative hypothesis H1.
It can be viewed as an internal function. 
The user does not have to call it by himself.
}

\usage{
obj.LD.add(param,don)
}

\arguments{
  \item{param}{numeric vector of length=5 containing parameters to maximize:

   -error variance

   -Q allele effect

   -evolution time

   -frequency of Q allele at time t=0

   -performances mean}


  \item{don}{list of 8 objects. Some objects are results provided by \code{corresp()}, \code{pi.hap()}, \code{structure.hap()}, \code{poids.D()} functions.}
}


\value{
The  returned value is the log-likelihood value of HAPimLD method under hypothesis H1.
}

\references{publication to be submitted: C. Cierco-Ayrolles, S. Dejean, A. Legarra, H. Gilbert,
T. Druet, F. Ytournel, D. Estivals, N. Oumouhou and B. Mangin.
Combining linkage analysis and linkage disequilibrium for QTL fine mapping in animal pedigrees.}

\author{S. Dejean, N. Oumouhou, D. Estivals, B. Mangin }

\seealso{\code{\link{depart.LD}}, \code{\link{corresp}}, \code{\link{pi.hap}}, \code{\link{structure.hap}}, \code{\link{poids.D}} }

\keyword{models}
