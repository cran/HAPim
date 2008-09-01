\name{retrouve.all}

\alias{retrouve.all}

\title{ haplotype reconstruction}

\description{
The function finds alleles of the associated haplotype.
It can be viewed as an internal function. 
The user does not have to call it by himself.
}

\usage{
retrouve.all(assoc, res.structure, all.marq)
}

\arguments{
  \item{assoc}{numeric value, associated haplotype.}
  \item{res.structure}{results provided by \code{structure.hap()} function, list of objects.}
  \item{all.marq}{results provided by \code{allele.marq()} function, list of alleles for each marker.}
}


\value{
The returned value  is a character which corresponds to the concatenation of  alleles of the associated haplotype.
}



\references{publication to be submitted: C. Cierco-Ayrolles, S. Dejean, A. Legarra, H. Gilbert,
T. Druet, F. Ytournel, D. Estivals, N. Oumouhou and B. Mangin.
Combining linkage analysis and linkage disequilibrium for QTL fine mapping in animal pedigrees.}




\author{S. Dejean, N. Oumouhou, D. Estivals, B. Mangin }

\seealso{\code{\link{structure.hap}}, \code{\link{allele.marq}} } 

\keyword{models}
