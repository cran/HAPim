\name{obj.LDL.add.H0}

\alias{obj.LDL.add.H0}

\title{log-likelihood value of the HAPIM LDL model under H0}

\description{
The function calculates the log-likelihood value of the HAPIM LDL model under hypothesis H0.
It can be viewed as an internal function. 
The user does not have to call it by himself.
}

\usage{
obj.LDL.add.H0(moyenne.pere, perf, CD, desc.pere)
}

\arguments{
  \item{moyenne.pere}{results provided by \code{moyenne.pere()} function, numeric vector of length=number of sires which contains the mean of performances of individuals per sire.}
  \item{perf}{numeric vector of length = number of individuals which contains the performances of individuals.}
  \item{CD}{numeric vector of length=number of individuals which contains the CD of individuals.}
  \item{desc.pere}{results provided by \code{descendant.pere()} function,numeric matrix (number of sires*2) which gives for each sire, the first and last indexes of its progeny.}
}

\value{
The value returned is the log-likelihood value of HAPIML LDL model under hypothesis HO.
}

\references{publication to be submitted:C.Cierco-Ayrolles,S.Dejean,A.Legarra,H.Gilbert,
T.Druet,F.Ytournel,D.Estivals,N.Oumouhou and B.Mangin,
Combining linkage analysis and linkage disequilibrium for QTL fine mapping in animal pedigrees.}

\author{S. Dejean, N. Oumouhou, D. Estivals, B. Mangin }

\seealso{\code{\link{moyenne.pere}},\code{\link{descendant.pere}} }

\keyword{models}
