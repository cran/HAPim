\name{obj.haplomax.add}

\alias{obj.haplomax.add}

\title{Analysis of variance}

\description{
The function gives an analysis of variance with a sire effect and a haplotype effect. 
It can be viewed as an internal function. 
The user does not have to call it by himself.
}


\usage{
obj.haplomax.add(perf, CD, assoc, res.structure, pi.hap, cor.pere, cor.mere)
}


\arguments{
  \item{perf}{numeric vector of length=number of individuals which contains the performances of individuals.}
  \item{CD}{numeric vector of length=number of individuals which contains the CD of individuals.}
  \item{assoc}{numeric value.}
  \item{res.structure}{results provided by \code{structure()} function, list of numeric objects.}
  \item{pi.hap}{provided by \code{pi.hap()} function, list of numeric objects.}
  \item{cor.pere}{provided by \code{corresp()} function, list of numeric objects.}
  \item{cor.mere}{provided by \code{corresp()} function, list of numeric objects.}
}


\value{
The value returned is a F-test value.
}

\references{publication to be submitted:C.Cierco-Ayrolles,S.Dejean,A.Legarra,H.Gilbert,
T.Druet,F.Ytournel,D.Estivals,N.Oumouhou and B.Mangin,
Combining linkage analysis and linkage disequilibrium for QTL fine mapping in animal pedigrees.}

\author{S. Dejean, N. Oumouhou, D. Estivals, B. Mangin }

\seealso{\code{\link{corresp}},\code{\link{pi.hap}},\code{\link{structure}},\code{\link{proba.DL}} }

\keyword{models}
