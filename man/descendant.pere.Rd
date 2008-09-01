\name{descendant.pere}

\alias{descendant.pere}

\title{first and last indexes of each sire progeny}

\description{
The function defines, for each sire, the first and last indexes of its progeny.
It can be viewed as an internal function. 
The user does not have to call it by himself.
}

\usage{
descendant.pere(genea)
}

\arguments{
   \item{genea}{numeric matrix (number of individuals x 2) 
which contains the progeny index and its father index.}
}

\details{
Progeny data are prealably ordered by family}

\value{
 The returned value  is a numeric matrix (number of sires x 2) which gives for each sire, the first and last indexes of its progeny.
}

\references{publication to be submitted: C. Cierco-Ayrolles, S. Dejean, A. Legarra, H. Gilbert,
T. Druet, F. Ytournel, D. Estivals, N. Oumouhou and B. Mangin.
Combining linkage analysis and linkage disequilibrium for QTL fine mapping in animal pedigrees.}



\author{S. Dejean, N. Oumouhou, D. Estivals, B. Mangin }

\keyword{models}

