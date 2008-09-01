\name{vrais.LDL.add.pere}

\alias{vrais.LDL.add.pere}

\title{intra-sire log likelihood of HAPimLDL method under H1.}

\description{
The function calculates the value of intra-sire log likelihood of HAPimLDL method under alternative hypothesis H1.
It can be viewed as an internal function. 
The user does not have to call it by himself.
}

\usage{
vrais.LDL.add.pere(moyenne.pere, alpha.Q, s, CD, perf, PLA, LD.m, LD.chrom1,

LD.chrom2, mean.gene)
}


\arguments{
  \item{moyenne.pere}{results provided by \code{moyenne.pere()} function, mean of half-sib family performances.}
  \item{alpha.Q}{parameter of QTL effect.}
  \item{s}{parameter of the error variance.}
  \item{CD}{numeric vector of length=number of individuals which contains the CD of individuals.  var(perf$_i$)=s/CD$^2_i$}
  \item{perf}{numeric vector of length=number of individuals which contains the performances of individuals.}
  \item{PLA}{numeric vector of length=number of individuals 
which contains transmission probabilities at a single test position}
  \item{LD.m}{results provided by \code{proba.DL()} function, numeric vector of length=number of individuals 
which contains probabilities due to linkage disequilibrium on dam.}
  \item{LD.chrom1}{results provided by \code{proba.DL()} function, numeric vector of length=number of individuals 
which contains probabilities due to linkage disequilibrium on the first chromosome of sire.}
  \item{LD.chrom2}{results provided by \code{proba.DL()} function, numeric vector of length=number of individuals 
which contains probabilities due to linkage disequilibrium on the second chromosome of sire.}
  \item{mean.gene}{parameter of performance mean.}
}


\value{
The returned value  is the value of intra-sire log likelihood of HAPimLDL method under alternative hypothesis H1.
}
\references{publication to be submitted: C. Cierco-Ayrolles, S. Dejean, A. Legarra, H. Gilbert,
T. Druet, F. Ytournel, D. Estivals, N. Oumouhou and B. Mangin.
Combining linkage analysis and linkage disequilibrium for QTL fine mapping in animal pedigrees.}


\author{S. Dejean, N. Oumouhou, D. Estivals, B. Mangin }


\seealso{\code{\link{moyenne.pere}}, \code{\link{proba.DL}} }


\keyword{models}


