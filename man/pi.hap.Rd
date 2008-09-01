\name{pi.hap}

\alias{pi.hap}

\title{haplotype frequencies}

\description{
This function computes the frequencies of each haplotype.
The function can be viewed as an internal function. 
The user does not have to call it by himself.
}

\usage{
pi.hap(freq.marq, res.structure)
}


\arguments{
  \item{freq.marq}{results provided by \code{freq.all()} function, list of numeric objects which contains the allele frequencies.}
  \item{res.structure}{results provided by \code{structure.hap()} function, list of objects.}
}


\value{
The  returned  value is a list of numeric vectors containing haplotypes freqencies estimated under marker independency.
}



\references{publication to be submitted: C. Cierco-Ayrolles, S. Dejean, A. Legarra, H. Gilbert,
T. Druet, F. Ytournel, D. Estivals, N. Oumouhou and B. Mangin.
Combining linkage analysis and linkage disequilibrium for QTL fine mapping in animal pedigrees.}



\author{S. Dejean, N. Oumouhou, D. Estivals, B. Mangin }

\seealso{\code{\link{freq.all}}, {\code{\link{structure.hap}}} }

\keyword{models}
