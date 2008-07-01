\name{corresp}

\alias{corresp}

\title{correspondance}

\description{
The function can be viewed as an internal function. 
The user does not have to call it by himself.
}

\usage{
corresp(hap.o, res.structure)
}

\arguments{
  \item{hap.o}{numeric matrix (n*p).}
  \item{res.structure}{results provided by \code{structure()} function, list of numeric objects.}
}


\value{
  A list containing the following components:
  \item{corresp}{numeric matrix (n*p).}
  \item{assoc}{numeric vector.}
}

\references{publication to be submitted:C.Cierco-Ayrolles,S.Dejean,A.Legarra,H.Gilbert,
T.Druet,F.Ytournel,D.Estivals,N.Oumouhou and B.Mangin,
Combining linkage analysis and linkage disequilibrium for QTL fine mapping in animal pedigrees.}

\author{S. Dejean, N. Oumouhou, D. Estivals, B. Mangin }

\seealso{\code{\link{structure}} }

\keyword{models}
