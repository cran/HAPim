\name{vrais.LDL.add.pere}

\alias{vrais.LDL.add.pere}

\title{intra-sire log likelihood of HAPIM LDL model under H1.}

\description{
The function calculates the value of intra-sire log likelihood of HAPIM LDL model under alternative hypothesis H1.
It can be viewed as an internal function. 
The user does not have to call it by himself.
}

\usage{
vrais.LDL.add.pere(moyenne.pere, alpha.Q, s, CD, perf, PLA, LD.m, LD.chrom1, LD.chrom2, mean.gene)
}


\arguments{
  \item{moyenne.pere}{results provided by \code{moyenne.pere()} function,
mean of performances of half-brothers.}
  \item{alpha.Q}{QTL effect.}
  \item{s}{variance of the error.}
  \item{CD}{numeric vector of length=number of individuals which contains the CD of individuals.}
  \item{perf}{numeric vector of length=number of individuals which contains the performances of individuals.}
  \item{PLA}{numeric matrix (number of individuals*number of test positions)
 which contains transmission probabilities i.e probability that the progeny receives the first chromosone from its sire given a marker,
see Jean-Michel Elsen, Brigitte Mangin, Bruno Goffinet, Didier Boichard, Pascale Le Roy, Alternative models for QTL detection, Genetics (1999) 217.}
  \item{LD.m}{results provided by \code{proba.DL()} function, numeric vector of length=number of individuals 
which contains probabilities due to linkage disequilibrium concerning dam.}
  \item{LD.chrom1}{results provided by \code{proba.DL()} function, numeric vector of length=number of individuals 
which contains probabilities due to linkage disequilibrium on the first chromosome of sire.}
  \item{LD.chrom2}{results provided by \code{proba.DL()} function, numeric vector of length=number of individuals 
which contains probabilities due to linkage disequilibrium on the second chromosome of sire.}
  \item{mean.gene}{general mean of performances of individuals.}
}


\value{
The value returned is the value of intra-sire log likelihood of HAPIM LDL model under alternative hypothesis H1.
}

\references{publication to be submitted:C.Cierco-Ayrolles,S.Dejean,A.Legarra,H.Gilbert,
T.Druet,F.Ytournel,D.Estivals,N.Oumouhou and B.Mangin,
Combining linkage analysis and linkage disequilibrium for QTL fine mapping in animal pedigrees.}

\author{S. Dejean, N. Oumouhou, D. Estivals, B. Mangin }


\seealso{\code{\link{moyenne.pere}},\code{\link{proba.DL}} }


\keyword{models}


