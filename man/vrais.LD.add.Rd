\name{vrais.LD.add}

\alias{vrais.LD.add}

\title{log likelihood of HAPIM LD model under H1}

\description{
The function calculates the value of log likelihood of HAPIM LD model under alternative hypothesis H1.
It can be viewed as an internal function. 
The user does not have to call it by himself.
}

\usage{
vrais.LD.add(mu, alpha.Q, s, CD, perf, DL.d)
}

\arguments{
  \item{mu}{mean of performances of individuals.}
  \item{alpha.Q}{QTL effect.}
  \item{s}{variance of the error.}
  \item{CD}{numeric vector of length=number of individuals which contains the CD of individuals.}
  \item{perf}{numeric vector of length=number of individuals which contains the performances of individuals.}
  \item{DL.d}{results provided by \code{proba.DL.diplotype()} function, numeric matrix (number of individuals*4) containing probabilities due to linkage disequilibrium for a diplotype.}
}

\value{
The value returned is the value of log likelihood of HAPIM LD model under alternative hypothesis H1.
}

\references{publication to be submitted:C.Cierco-Ayrolles,S.Dejean,A.Legarra,H.Gilbert,
T.Druet,F.Ytournel,D.Estivals,N.Oumouhou and B.Mangin,
Combining linkage analysis and linkage disequilibrium for QTL fine mapping in animal pedigrees.}

\author{S. Dejean, N. Oumouhou, D. Estivals, B. Mangin }

\seealso{\code{\link{proba.DL.diplotype}} }

\keyword{models}
