\name{pi.hap.NI}

\alias{pi.hap.NI}

\title{haplotype frequencies}

\description{
This function computes the observed frequencies of each haplotype.
The function can be viewed as an internal function. 
The user does not have to call it by himself.
}

\usage{
pi.hap.NI(res.structure, cor.hap)
}

\arguments{
  \item{res.structure}{results provided by \code{structure.hap()} function, list of objects.}
  \item{cor.hap}{results provided by \code{corresp()} function, list of numeric objects.}
}


\value{
The returned  value is a list of numeric vectors containing haplotypes freqencies estimated on haplotypes with non missing information.
}
\references{publication to be submitted: C. Cierco-Ayrolles, S. Dejean, A. Legarra, H. Gilbert,
T. Druet, F. Ytournel, D. Estivals, N. Oumouhou and B. Mangin.
Combining linkage analysis and linkage disequilibrium for QTL fine mapping in animal pedigrees.}



\author{S. Dejean, N. Oumouhou, D. Estivals, B. Mangin }

\seealso{\code{\link{structure.hap}} }

\keyword{models}


