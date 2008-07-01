\name{retrouve.all}

\alias{retrouve.all}

\title{Alleles of reconstruted haplotype}

\description{
The function finds alleles for a reconstructed haplotype.
It can be viewed as an internal function. 
The user does not have to call it by himself.
}

\usage{
retrouve.all(assoc, res.structure, all.marq)
}

\arguments{
  \item{assoc}{numeric value.}
  \item{res.structure}{results provided by \code{structure()} function, list of numeric objects.}
  \item{all.marq}{results provided by \code{allele.marq()} function,list of alleles for each marker given a haplotype.}
}

\details{
It aims at giving for each test position the corresponding sequence of alleles of reconstruted haplotype.
}

\value{
The value returned is a concatened character sequence of alleles given a reconstructed haplotype.
}

\references{publication to be submitted:C.Cierco-Ayrolles,S.Dejean,A.Legarra,H.Gilbert,
T.Druet,F.Ytournel,D.Estivals,N.Oumouhou and B.Mangin,
Combining linkage analysis and linkage disequilibrium for QTL fine mapping in animal pedigrees.}

\author{S. Dejean, N. Oumouhou, D. Estivals, B. Mangin }

\seealso{\code{\link{structure}},\code{\link{allele.marq}} } 

\keyword{models}
