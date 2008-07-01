\name{hapim.LDL.add}

\alias{hapim.LDL.add}

\title{HAPim LDL model}

\description{
HAPim LDL is a method of detection of QTLs (Quantitative Trait Loci) for animal pedigrees.
It is an interval-mapping method which uses family information of individuals 
and combines linkage and linkage disequilibrium. 
It is based on a maximum-likelihood calculation.
}

\usage{
hapim.LDL.add(hap.chrom1.pere, hap.chrom2.pere, hap.trans.mere, perf, CD, genea, PLA, map, position, temps.depart, perfectLD, marq.hap.left)
}

\arguments{
  \item{hap.chrom1.pere}{character matrix (number of sires*number of markers)
which gives the haplotype of the first chromosome for each sire.}
  \item{hap.chrom2.pere}{character matrix (number of sires*number of markers)
 which gives the haplotype of the second chromosome for each sire.}
  \item{hap.trans.mere}{numeric matrix (number of individuals*number of markers) 
which provides, for each individual, the haplotype transmitted by the dam.}
  \item{perf}{numeric vector of length=number of individuals which contains the performances of individuals.}
  \item{CD}{numeric vector of length=number of individuals which contains the CD of individuals.}
  \item{genea}{numeric matrix (number of individuals*2) 
which contains individual index and corresponding sire index of each individual.}
  \item{PLA}{numeric matrix (number of individuals*number of test positions)
 which contains transmission probabilities i.e probability that the progeny receives the first chromosone from its sire given a marker,
see Jean-Michel Elsen, Brigitte Mangin, Bruno Goffinet, Didier Boichard, Pascale Le Roy, Alternative models for QTL detection, Genetics (1999) 217.}
  \item{map}{numeric vector of length=(number of markers-1) 
giving the distance between two consecutive markers on all the chromosome.}
  \item{position}{numeric vector of test positions.
   Positions of "position" vector have to be the same as positions of "PLA".}
  \item{temps.depart}{time of evolution of the population.}
  \item{perfectLD}{TRUE:linkage disequilibrium is total between mutated haplotype and QTL at time 0.}
  \item{marq.hap.left}{(number of markers of the reconstructed haplotype)/2.}
}

\details{
Individuals have to be ranged in the same order in genea,hap.trans.mere,perf,CD,PLA.
}

\value{
The value returned is a data frame which contains 8 columns:

-Test positions

-Value of Likelihood Ratio Test (LRT) 

-Reconstructed haplotype coded with letters

-Estimation of the variance of the error 

-Estimation of the QTL effect 

-Estimation of the time of evolution of the population 

-Estimation of the frequency of Q allele at time t=0

-Estimation of the mean of performances of individuals 
}

\references{publication to be submitted:C.Cierco-Ayrolles,S.Dejean,A.Legarra,H.Gilbert,
T.Druet,F.Ytournel,D.Estivals,N.Oumouhou and B.Mangin,
Combining linkage analysis and linkage disequilibrium for QTL fine mapping in animal pedigrees.}

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

########################################################################################
# warning:positions in "position" vector have to correspond to positions in "PLA" matrix.
########################################################################################
# In this example:
#	-positions test are expressed in Morgan.
#	-markers are biallelic (SNPs).
#	-marker positions are : {0, 0.010, 0.020, 0.030, 0.040, 0.050, 0.060, 0.070, 0.080, 0.090 }. 
# We want to test the presence/absence of a QTL between 2 consecutive markers.
# Transmission probabilities are calculated for the middle of each interval between 2 consecutive markers.
# So, positions of "position" vector are as follow:
 
position=c(0.005,0.015,0.025,0.035,0.045,0.055,0.065,0.075,0.085)

# we use a 2 markers-reconstructed haplotype
marq.hap.left=1

# We assume temps.depart=50 which corresponds to time of evolution of the population.
temps.depart=50

# We consider the linkage disequilibrium is total
perfectLD=TRUE


hapim.LDL.add=hapim.LDL.add(hap.chrom1.pere,hap.chrom2.pere,hap.trans.mere,perf,CD,genea,PLA,map,position,temps.depart,perfectLD,marq.hap.left)

hapim.LDL.add

}

\keyword{models}


