\name{obj.LDL.add}

\alias{obj.LDL.add}

\title{log-likelihood value of the HAPIM LDL model under H1}

\description{
The function calculates the log-likelihood value of the HAPIM LDL model under alternative hypothesis H1.
It can be viewed as an internal function. 
The user does not have to call it by himself.
}

\usage{
obj.LDL.add(param,don)
}

\arguments{
  \item{param}{numeric vector of length=5 containing containing parameters to maximize:

   -variance of the error

   -QTL effect

   -frequency of Q allele at time t=0

   -mean of performances of individuals

   The values of the variance and the QTL effect are provided by \code{depart.LDL()} function.}
  \item{don}{list of 12 objects.Some objects of don are results provided by \code{descendant.pere()},
  \code{moyenne.pere()},\code{corresp()},\code{pi.hap()},\code{structure()},\code{poids.D()} functions.}
}

\value{
The value returned is the log-likelihood value of the HAPIM LDL model under alternative hypothesis H1.
}

\references{publication to be submitted:C.Cierco-Ayrolles,S.Dejean,A.Legarra,H.Gilbert,
T.Druet,F.Ytournel,D.Estivals,N.Oumouhou and B.Mangin,
Combining linkage analysis and linkage disequilibrium for QTL fine mapping in animal pedigrees.}

\author{S. Dejean, N. Oumouhou, D. Estivals, B. Mangin }

\seealso{\code{\link{depart.LDL}},\code{\link{descendant.pere}},\code{\link{moyenne.pere}},\code{\link{corresp}},
\code{\link{pi.hap}},\code{\link{structure}},\code{\link{poids.D}} }

\keyword{models}

 
