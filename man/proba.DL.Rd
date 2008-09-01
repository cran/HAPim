\name{proba.DL}

\alias{proba.DL}

\title{conditional probabilities}

\description{
The function calculates conditional probabilities of having Q allele given marker information due to linkage disequilibrium.
}

\usage{
proba.DL(piQ.t0, esp.freq.hap, res.structure, pi.hap, res.corresp)
}

\arguments{
  \item{piQ.t0}{frequency of Q allele at time t=0.}
  \item{esp.freq.hap}{results provided by \code{esp.freq.hap()} function, numeric vector.}
  \item{res.structure}{results provided by \code{structure.hap()} function, list of numeric objects.}
  \item{pi.hap}{provided by \code{pi.hap()} function, list of numeric objects.}
  \item{res.corresp}{results provided by \code{corresp()} function, list of 2 numeric objects.}
}


\value{
The returned value  is a numeric vector of length=number of individuals 
which contains conditional probabilities of having Q allele given marker information due to linkage disequilibrium.
}
\references{publication to be submitted: C. Cierco-Ayrolles, S. Dejean, A. Legarra, H. Gilbert,
T. Druet, F. Ytournel, D. Estivals, N. Oumouhou and B. Mangin.
Combining linkage analysis and linkage disequilibrium for QTL fine mapping in animal pedigrees.}



\author{S. Dejean, N. Oumouhou, D. Estivals, B. Mangin }

\seealso{\code{\link{esp.freq.hap}}, \code{\link{structure.hap}}, \code{\link{pi.hap}}, \code{\link{corresp}} }

\keyword{models}
