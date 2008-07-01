\name{esp.freq.hap}

\alias{esp.freq.hap}

\title{esp.freq.hap}

\description{
The function can be viewed as an internal function. 
The user does not have to call it by himself.
}

\usage{
esp.freq.hap(hap.assoc, piQ.t0, timeT, pi.hap, res.structure, poids.D)
}

\arguments{
  \item{hap.assoc}{numeric value.}
  \item{piQ.t0}{frequency of Q allele at time t=0.}
  \item{timeT}{time.}
  \item{pi.hap}{provided by \code{pi.hap()} function, list of numeric objects.}
  \item{res.structure}{provided by \code{structure()} function, list of numeric objects.}
  \item{poids.D}{provided by \code{poids.D()} function, list of numeric objects.}
}


\value{
The value returned is a numeric vector.
}

\references{publication to be submitted:C.Cierco-Ayrolles,S.Dejean,A.Legarra,H.Gilbert,
T.Druet,F.Ytournel,D.Estivals,N.Oumouhou and B.Mangin,
Combining linkage analysis and linkage disequilibrium for QTL fine mapping in animal pedigrees.}

\author{S. Dejean, N. Oumouhou, D. Estivals, B. Mangin }

\seealso{\code{\link{pi.hap}},\code{\link{structure}},\code{\link{poids.D}} }

\keyword{models}
