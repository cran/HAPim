\name{depart.LDL}

\alias{depart.LDL}

\title{starting values for the optimization of HAPimLDL method}

\description{
The function calculates the starting value of the error variance  and the starting value of the QTL effect for the optimization of HAPimLDL method.
It can be viewed as an internal function. 
The user does not have to call it by himself.
}

\usage{
depart.LDL(moyenne.pere, perf, CD, PLA, desc.pere)
}

\arguments{
  \item{moyenne.pere}{results provided by \code{moyenne.pere()} function, mean of half-sib family performances.}
  \item{perf}{numeric vector of length=number of individuals which contains the performances of individuals.}
  \item{CD}{numeric vector of length=number of individuals which contains the CD of individuals. var(perf$_i$)=s/CD$^2_i$}
  \item{PLA}{numeric vector (number of individuals) 
which contains transmission probabilities at a single test position.}
  \item{desc.pere}{results provided by \code{descendant.pere()} function, numeric matrix (number of sires x 2)
which gives for each sire, the first and last indexes of its progeny.}
}

\value{
The returned value  is a numeric vector of length=2 which contains estimates of the error variance and the Q allele effect.
}

\references{publication to be submitted: C. Cierco-Ayrolles, S. Dejean, A. Legarra, H. Gilbert,
T. Druet, F. Ytournel, D. Estivals, N. Oumouhou and B. Mangin.
Combining linkage analysis and linkage disequilibrium for QTL fine mapping in animal pedigrees.}


\author{S. Dejean, N. Oumouhou, D. Estivals, B. Mangin }

\seealso{\code{\link{moyenne.pere}}, \code{\link{descendant.pere}} }

\keyword{models}


