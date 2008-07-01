\name{freq.all}

\alias{freq.all}

\title{allelic frequencies}

\description{
The function calculates allelic frequencies by marker given a set of haplotypes.
It can be viewed as an internal function. 
The user does not have to call it by himself.
}

\usage{
freq.all(hap)
}

\arguments{
  \item{hap}{numeric matrix (number of individuals*number of markers) which contains the haplotype of individuals.}
}

\value{
The value returned is a list of n elements (n=number of markers) which 
gives the allelic frequencies by marker and ranged in ascending order of allele number.
}

\references{publication to be submitted:C.Cierco-Ayrolles,S.Dejean,A.Legarra,H.Gilbert,
T.Druet,F.Ytournel,D.Estivals,N.Oumouhou and B.Mangin,
Combining linkage analysis and linkage disequilibrium for QTL fine mapping in animal pedigrees.}

\author{S. Dejean, N. Oumouhou, D. Estivals, B. Mangin }

\keyword{models}



