\name{moyenne.pere}

\alias{moyenne.pere}

\title{average of  performances per sire}

\description{
The function computes the empirical performance mean per sire.
It can be viewed as an internal function. 
The user does not have to call it by himself.
}

\usage{
moyenne.pere(perf, CD, desc.pere)
}

\arguments{
  \item{perf}{numeric vector of length=number of individuals 
which contains the performances of individuals.}
  \item{CD}{numeric vector of length=number of individuals which
 contains the CD of individuals.  var(perf$_i$)=error variance/CD$^2_i$}
  \item{desc.pere}{results provided by \code{descendant.pere()} function, 
numeric matrix (number of sires x 2) which gives for each sire, the first and last indexes of its progeny.}
}


\value{
The returned value  is a vector of length=number of sires which contains the empirical performance mean per sire.
}

\references{

publication to be submitted: C. Cierco-Ayrolles, S. Dejean, A. Legarra, H. Gilbert,
T. Druet, F. Ytournel, D. Estivals, N. Oumouhou and B. Mangin.
Combining linkage analysis and linkage disequilibrium for QTL fine mapping in animal pedigrees.

}

\author{S. Dejean, N. Oumouhou, D. Estivals, B. Mangin }

\seealso{\code{\link{descendant.pere}} }

\keyword{models}

