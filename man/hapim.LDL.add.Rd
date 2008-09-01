\name{hapim.LDL.add}

\alias{hapim.LDL.add}

\title{HAPimLDL method}

\description{
HAPimLDL is a method of QTL (Quantitative Trait Loci) detection  for a design of half-sib families.
It is an interval-mapping method which uses family information 
and combines linkage analysis and linkage disequilibrium. 
It is based on a maximum-likelihood calculation and makes use of linkage disequilibrium through a Wright-Fisher modelisation of the population evolution.
}

\usage{
hapim.LDL.add(hap.chrom1.pere, hap.chrom2.pere, hap.trans.mere, perf, CD, 

genea, PLA, map, position, temps.depart, perfectLD, marq.hap.left)
}

\arguments{
  \item{hap.chrom1.pere}{character matrix (number of sires x number of markers)
which gives the haplotype of the first chromosome for each sire.}
  \item{hap.chrom2.pere}{character matrix (number of sires x number of markers)
 which gives the haplotype of the second chromosome for each sire.}
  \item{hap.trans.mere}{numeric matrix (number of individuals x number of markers) 
which provides, for each individual, the haplotype transmitted by its mother.}
  \item{perf}{numeric vector of length=number of individuals which contains the performances of individuals.}
  \item{CD}{numeric vector of length=number of individuals which contains the CD of individuals. var(perf$_i$)=error variance/CD$^2_i$}
  \item{genea}{numeric matrix (number of individuals x 2) 
which contains the progeny index and its father index.}
  \item{PLA}{numeric matrix (number of individuals x number of test positions)
 which contains transmission probabilities i.e probability that the progeny receives the first chromosome from its father at the test positions given marker information,
see J.-M. Elsen, B. Mangin, B. Goffinet, D. Boichard, P. Le Roy.
Alternative models for QTL detection in livestock. I. General introduction. {\it Genet. Sel. Evol.} 31 (1999) 213-224.}
  \item{map}{numeric vector of length=(number of markers-1) which contains the distance in Morgan between two consecutive markers on the chromosome.}
  \item{position}{numeric vector which contains the distance in Morgan of test positions from the beginning of the chromosome (first marker).}
  \item{temps.depart}{numeric value which provides a start value for the evolution  time of the population.}
  \item{perfectLD}{need to be equal to TRUE: linkage disequilibrium is complete between mutated haplotype and Q allele at time 0.}
  \item{marq.hap.left}{(number of markers of the mutated haplotype)/2.}


}

\details{
Progeny information have to be ranged in the same order in genea, hap.trans.mere, perf, CD and PLA.

Columns of PLA have to correspond to test positions.  

Sire information have to be ranged in the same order in unique(genea[,2]), hap.chrom1.pere and hap.chrom2.pere.

All distances are assumed to be Haldame's distance in Morgan.
}

\value{
The returned value is a data frame which contains 8 columns:

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

publication to be submitted: C. Cierco-Ayrolles, S. Dejean, A. Legarra, H. Gilbert,
T. Druet, F. Ytournel, D. Estivals, N. Oumouhou and B. Mangin.
Combining linkage analysis and linkage disequilibrium for QTL fine mapping in animal pedigrees.

}

\author{S. Dejean, N. Oumouhou, D. Estivals, B. Mangin }

\examples{
data(data.test)
map=data.test[[1]]
hap.trans.mere=data.test[[2]]
hap.chrom1.pere=data.test[[4]]
hap.chrom2.pere=data.test[[5]]
perf=data.test[[6]]
CD=data.test[[7]]
PLA=data.test[[8]]
genea=data.test[[9]]


# In this example, marker positions are : {0, 0.010, 0.020, 0.030, 0.040, 0.050, 0.060, 
# 0.070, 0.080, 0.090 }. 
# We want to test the presence/absence of a QTL between 2 consecutive markers, so
# transmission probabilities are given for the middle of each interval and
 
position=c(0.005,0.015,0.025,0.035,0.045,0.055,0.065,0.075,0.085)

# we use a 2 markers-associated haplotype.
marq.hap.left=1

# We assume an evolution of 50 generations.
temps.depart=50
perfectLD=TRUE



hapim.LDL.add=hapim.LDL.add(hap.chrom1.pere,hap.chrom2.pere,hap.trans.mere,perf,CD,genea,

PLA,map,position, temps.depart,perfectLD,marq.hap.left)

hapim.LDL.add

}

\keyword{models}


