\name{haplomax.add}

\alias{haplomax.add}

\title{HaploMax method in unrelated population}

\description{
The function computes an analysis of variance with a dose haplotype effect. 
}

\usage{
haplomax.add(hap.trans.pere, hap.trans.mere, perf, CD, map, marq.hap)
}


\arguments{
  \item{hap.trans.pere}{character matrix (number of individuals x number of markers) 
which provides, for each individual, the haplotype transmitted by its father.}
  \item{hap.trans.mere}{character matrix (number of individuals x number of markers) 
which provides, for each individual, the haplotype transmitted by its mother.}
  \item{perf}{numeric vector of length=number of individuals which contains the performances of individuals.}
  \item{CD}{numeric vector of length=number of individuals which contains the CD of individuals. var(perf$_i$)=error variance/CD$^2_i$}
  \item{map}{numeric vector of length=(number of markers-1) which contains the distance in Morgan between two consecutive markers on the chromosome.}
  \item{marq.hap}{number of markers of the mutated haplotype}
}

\details{
Individual information have to be ranged in the same order in hap.trans.mere, hap.trans.pere, perf, CD.

All distances are assumed to be Haldame's distance in Morgan.

Test positions are located on the middles of marq.hap marker sliding windows.
}


\value{
The value returned is a data frame which contains 5 columns:

-Test positions

-Value of Fisher test 

-Mutated (i.e. associated to Q allele) haplotype 

-Estimate of the error variance 

-Estimate of the Q allele effect 

}

\references{
publication to be submitted: C. Cierco-Ayrolles, S. Dejean, A. Legarra, H. Gilbert,
T. Druet, F. Ytournel, D. Estivals, N. Oumouhou and B. Mangin.
Combining linkage analysis and linkage disequilibrium for QTL fine mapping in animal pedigrees.}


\author{S. Dejean, N. Oumouhou, D. Estivals, B. Mangin }

\examples{

data(data.test)
map=data.test[[1]]
hap.trans.mere=data.test[[2]]
hap.trans.pere=data.test[[3]]
perf=data.test[[6]]
CD=data.test[[7]]

# In this example,marker positions are: {0, 0.010, 0.020, 0.030, 0.040, 0.050, 0.060, 
# 0.070, 0.080, 0.090 }. 
# we use a 2 markers-associated haplotype.
marq.hap=2

haplomax=haplomax.add(hap.trans.pere,hap.trans.mere,perf,CD,map, marq.hap)

haplomax

}

\keyword{models}


