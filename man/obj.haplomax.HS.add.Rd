\name{obj.haplomax.HS.add}

\alias{obj.haplomax.HS.add}

\title{Analysis of variance of the HaploMmax method in half-sib family design}

\description{
The function computes the  regression analysis with a dose haplotype effect and a sire effect in a design of half-sib families.
It can be viewed as an internal function. 
The user does not have to call it by himself.
}

\usage{
obj.haplomax.HS.add(genea, perf, CD, assoc, res.structure, pi.hap, cor.pere, 

cor.mere)
}


\arguments{
  \item{genea}{numeric matrix (number of individuals x 2) which contains individual index and corresponding sire index of each individual.}
  \item{perf}{numeric vector of length=number of individuals which contains the performances of individuals.}
  \item{CD}{numeric vector of length=number of individuals which contains the CD of individuals. var(perf$_i$)=error variance/CD$^2_i$}
  \item{assoc}{numeric value, associated haplotype}
  \item{res.structure}{provided by \code{structure.hap()} function, list of  objects.}
  \item{pi.hap}{provided by \code{pi.hap()} function, list of numeric objects.}
  \item{cor.pere}{provided by \code{corresp()} function, list of numeric objects.}
  \item{cor.mere}{provided by \code{corresp()} function, list of numeric objects.}
}

\value{
The returned  value is an object of aov class containing the dose haplotype + sire effect regression.
}

\references{publication to be submitted: C. Cierco-Ayrolles, S. Dejean, A. Legarra, H. Gilbert,
T. Druet, F. Ytournel, D. Estivals, N. Oumouhou and B. Mangin.
Combining linkage analysis and linkage disequilibrium for QTL fine mapping in animal pedigrees.}

\author{S. Dejean, N. Oumouhou, D. Estivals, B. Mangin }

\seealso{\code{\link{corresp}}, \code{\link{pi.hap}}, \code{\link{structure.hap}} }

\keyword{models}
