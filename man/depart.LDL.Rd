\name{depart.LDL}

\alias{depart.LDL}

\title{starting values for the optimization of HAPIM LDL model}

\description{
The function calculates the starting value of the variance of the error and the starting value of the QTL effect for the optimization of HAPIM LDL model.
It can be viewed as an internal function. 
The user does not have to call it by himself.
}

\usage{
depart.LDL(moyenne.pere, perf, CD, PLA, desc.pere)
}

\arguments{
  \item{moyenne.pere}{results provided by \code{moyenne.pere()} function, numeric vector of length=number of sires which contains the average performance of individuals per sire.}
  \item{perf}{numeric vector of length=number of individuals which contains the performances of individuals.}
  \item{CD}{numeric vector of length=number of individuals which contains the CD of individuals.}
  \item{PLA}{numeric matrix (number of individuals*number of test positions) which contains transmission probabilities
i.e probability that the progeny receives the first chromosone from its sire given a marker,see Jean-Michel Elsen, 
Brigitte Mangin, Bruno Goffinet, Didier Boichard, Pascale Le Roy, Alternative models for QTL detection, Genetics (1999) 217.}
  \item{desc.pere}{results provided by \code{descendant.pere()} function, numeric matrix (number of sires*2)
which gives for each sire, the first and last indexes of its progeny.}
}

\value{
The value returned is a numeric vector of length=2 which contains the variance of the error and the QTL effect.
}

\references{publication to be submitted:C.Cierco-Ayrolles,S.Dejean,A.Legarra,H.Gilbert,
T.Druet,F.Ytournel,D.Estivals,N.Oumouhou and B.Mangin,
Combining linkage analysis and linkage disequilibrium for QTL fine mapping in animal pedigrees.}


\author{S. Dejean, N. Oumouhou, D. Estivals, B. Mangin }

\seealso{\code{\link{moyenne.pere}},\code{\link{descendant.pere}} }

\keyword{models}


