\name{vrais.LD.add}

\alias{vrais.LD.add}

\title{log likelihood of HAPIimLD method under H1}

\description{
The function calculates the value of log likelihood of HAPimLD method under alternative hypothesis H1.
It can be viewed as an internal function. 
The user does not have to call it by himself.
}

\usage{
vrais.LD.add(mu, alpha.Q, s, CD, perf, DL.d)
}

\arguments{
  \item{mu}{parameter of the  mean  performances.}
  \item{alpha.Q}{parameter of QTL effect.}
  \item{s}{parameter of the error variance.}
  \item{CD}{numeric vector of length=number of individuals which contains the CD of individuals. var(perf$_i$)=s/CD$^2_i$}
  \item{perf}{numeric vector of length=number of individuals which contains the performances of individuals.}
  \item{DL.d}{results provided by \code{proba.DL.diplotype()} function, numeric matrix (number of individuals x 4) containing probabilities due to linkage disequilibrium for a diplotype.}
}

\value{
The returned value  is the value of log likelihood of HAPimLD method under alternative hypothesis H1.
}


\references{publication to be submitted: C. Cierco-Ayrolles, S. Dejean, A. Legarra, H. Gilbert,
T. Druet, F. Ytournel, D. Estivals, N. Oumouhou and B. Mangin.
Combining linkage analysis and linkage disequilibrium for QTL fine mapping in animal pedigrees.}


\author{S. Dejean, N. Oumouhou, D. Estivals, B. Mangin }

\seealso{\code{\link{proba.DL.diplotype}} }

\keyword{models}
