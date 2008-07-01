\name{recode.hap}

\alias{recode.hap}

\title{recoding of haplotypes}

\description{
The function recodes alleles of haplotypes of individuals in consecutive numeric values.
It can be viewed as an internal function. 
The user does not have to call it by himself.
}

\usage{
recode.hap(hap, all.marq)
}


\arguments{
  \item{hap}{character matrix (number of individuals*number of markers) 
which contains the haplotype of individuals.}
  \item{all.marq}{results provided by \code{allele.marq()} function,
list of alleles (coded with character values) for each marker given a haplotype.}
}


\value{
the value returned is a numeric matrix (number of individuals*number of markers)
 which contains the recoded haplotype of individuals.
}

\references{publication to be submitted:C.Cierco-Ayrolles,S.Dejean,A.Legarra,H.Gilbert,
T.Druet,F.Ytournel,D.Estivals,N.Oumouhou and B.Mangin,
Combining linkage analysis and linkage disequilibrium for QTL fine mapping in animal pedigrees.}

\author{S. Dejean, N. Oumouhou, D. Estivals, B. Mangin }

\seealso{\code{\link{allele.marq}} }

\keyword{models}
