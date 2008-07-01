\name{proba.DL.diplotype}

\alias{proba.DL.diplotype}

\title{probabilities due to linkage disequilibrium}

\description{
The function calculates probabilities due to linkage disequilibrium for a diplotype.
It can be viewed as an internal function. 
The user does not have to call it by himself.
}

\usage{
proba.DL.diplotype(DL.chrom1, DL.chrom2)
}

\arguments{
  \item{DL.chrom1}{results provided by \code{proba.DL()} function, 
numeric vector of length=number of individuals which contains probabilities due to linkage disequilibrium on the first chromosome given a genotype.}
  \item{DL.chrom2}{results provided by \code{proba.DL()} function,
numeric vector of length=number of individuals which contains probabilities due to linkage disequilibrium on the second chromosome given a genotype.}
   }


\value{
The value returned is a numeric matrix (number of individuals*4)
containing probabilities due to linkage disequilibrium for a diplotype.
}

\references{publication to be submitted:C.Cierco-Ayrolles,S.Dejean,A.Legarra,H.Gilbert,
T.Druet,F.Ytournel,D.Estivals,N.Oumouhou and B.Mangin,
Combining linkage analysis and linkage disequilibrium for QTL fine mapping in animal pedigrees.}

\author{S. Dejean, N. Oumouhou, D. Estivals, B. Mangin }

\seealso{\code{\link{proba.DL}} }

\keyword{models}
