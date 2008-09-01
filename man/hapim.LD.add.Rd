\name{hapim.LD.add}

\alias{hapim.LD.add}

\title{HAPimLD method}

\description{
HAPimLD is a method of QTL(Quantitative Trait Loci) detection developed 
by Boitard et al. (2006).
It is an interval-mapping method designed for unrelated individuals with no family information. It is based on a maximum-likelihood calculation and makes use of linkage disequilibrium through a Wright-Fisher modelisation of the population evolution.
}

\usage{
hapim.LD.add(hap.trans.pere, hap.trans.mere, perf, CD, map, position, 

temps.depart, perfectLD, marq.hap.left)
}

\arguments{
  \item{hap.trans.pere}{character matrix (number of individuals x number of markers) 
which provides, for each individual, the haplotype transmitted by its father.}
  \item{hap.trans.mere}{character matrix (number of individuals x number of markers) 
which provides, for each individual, the haplotype transmitted by its mother.}
  \item{perf}{numeric vector of length=number of individuals which contains the performances of individuals.}
  \item{CD}{numeric vector of length=number of individuals which contains the CD of individuals. var(perf$_i$)=error variance/CD$^2_i$}
  \item{map}{numeric vector of length=(number of markers-1) which contains the distance in Morgan between two consecutive markers on the chromosome.}
  \item{position}{numeric vector which contains the distance in Morgan of test positions from the beginning of the chromosome (first marker).}
  \item{temps.depart}{numeric value which provides a start value for the evolution  time of the population.}
  \item{perfectLD}{need to be equal to TRUE: linkage disequilibrium is complete between mutated haplotype and Q allele at time 0.}
  \item{marq.hap.left}{(number of markers of the mutated haplotype)/2.}
}

\details{
Individual information have to be ranged in the same order in hap.trans.mere, hap.trans.pere, perf, CD.

All distances are assumed to be Haldame's distance in Morgan.
}

\value{
The returned value  is a data frame which contains 8 columns:

-Test positions

-Value of Likelihood Ratio Test (LRT) 

-Mutated (i.e. associated to Q allele) haplotype 

-Estimate of the error variance 

-Estimate of the Q allele effect 

-Estimate of the time of population  evolution 

-Estimate of the Q allele frequency  at time t=0

-Estimate of the performance mean  
}

\references{
Boitard et al. Linkage disequilibrium interval mapping of quantitative trait loci.{\it BMC Genomics} (2006) 7:54.


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


# In this example,marker positions are: 0, 0.010, 0.020, 0.030, 0.040, 0.050, 0.060, 
# 0.070, 0.080, 0.090. 
# We want to test the presence/absence of a QTL between 2 consecutive markers, so

position=c(0.005,0.015,0.025,0.035,0.045,0.055,0.065,0.075,0.085)

# we use a 2 markers-associated haplotype.
marq.hap.left=1

# We assume an evolution of 50 generations.
temps.depart=50
perfectLD=TRUE


hapim.LD.add=hapim.LD.add(hap.trans.pere,hap.trans.mere,perf,CD,map,position,

temps.depart,perfectLD,marq.hap.left)

hapim.LD.add

}

\keyword{models}


