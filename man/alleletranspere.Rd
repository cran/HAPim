\name{alleletranspere}

\alias{alleletranspere}

\title{Allele transmitted by the sire to its progeny}

\description{
The function provides the value of the allele transmitted by the sire to its progeny given a marker and a transmission probability.
It can be viewed as an internal function. 
The user does not have to call it by himself.
}

\usage{
alleletranspere(all1.pere,all2.pere,all1.fils,all2.fils,PLA,seuil)
}

\arguments{
  \item{all1.pere}{the value of allele on the first chromosome of the sire given a marker.}
  \item{all2.pere}{the value of allele on the second chromosome of the sire given a marker.}
  \item{all1.fils}{the value of allele on the first chromosome of the progeny given a marker.}
  \item{all2.fils}{the value of allele on the second chromosome of the progeny given a marker.}
  \item{PLA}{transmission probability i.e probability that the progeny receives the first chromosone 
from its sire given a marker,see Jean-Michel Elsen, Brigitte Mangin, Bruno Goffinet, Didier Boichard, 
Pascale Le Roy, Alternative models for QTL detection, Genetics (1999) 217.}
  \item{seuil}{threshold / transmission probability>=seuil.}
}

\value{
The value returned is the value of the allele transmitted by the sire given a marker.
}

\references{publication to be submitted:C.Cierco-Ayrolles,S.Dejean,A.Legarra,H.Gilbert,
T.Druet,F.Ytournel,D.Estivals,N.Oumouhou and B.Mangin,
Combining linkage analysis and linkage disequilibrium for QTL fine mapping in animal pedigrees.}


\author{S. Dejean, N. Oumouhou, D. Estivals, B. Mangin }


\keyword{models}


