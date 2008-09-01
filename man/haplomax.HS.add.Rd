\name{haplomax.HS.add}

\alias{haplomax.HS.add}

\title{HaploMax method in half-sib family design}

\description{
The function computes an analysis of variance with a sire effect and a dose haplotype effect.
}

\usage{
haplomax.HS.add(hap.chrom1.pere, hap.chrom2.pere, hap.trans.pere, 

hap.trans.mere, perf, CD, genea, map, marq.hap)
}


\arguments{

  \item{hap.chrom1.pere}{character matrix (number of sires x number of markers)
which gives the haplotype of the first chromosome for each sire.}
  \item{hap.chrom2.pere}{character matrix (number of sires x number of markers)
 which gives the haplotype of the second chromosome for each sire.}
  \item{hap.trans.pere}{numeric matrix (number of individuals x number of markers) 
which provides, for each individual, the haplotype transmitted by its father.}
  \item{hap.trans.mere}{numeric matrix (number of individuals x number of markers) 
which provides, for each individual, the haplotype transmitted by its mother.}
  \item{perf}{numeric vector of length=number of individuals which contains the performances of individuals.}
  \item{CD}{numeric vector of length=number of individuals which contains the CD of individuals. var(perf$_i$)=error variance/CD$^2_i$}
  \item{genea}{numeric matrix (number of individuals x 2) 
which contains the progeny index and its father index.}
  \item{map}{numeric vector of length=(number of markers-1) which contains the distance in Morgan between two consecutive markers on the chromosome.}
 \item{marq.hap}{number of markers of the mutated haplotype.}
}
\details{
Progeny information have to be ranged in the same order in genea, hap.trans.pere, hap.trans.mere, perf and CD.
 

Sire information have to be ranged in the same order in unique(genea[,2]), hap.chrom1.pere and hap.chrom2.pere.

All distances are assumed to be Haldame's distance in Morgan.

Test positions are located on the middles of marq.hap marker sliding windows.
}

\value{
The returned value is a data frame which contains 5 columns:

-Test positions

-Value of Fisher test 

-Mutated (i.e. associated to Q allele) haplotype 

-Estimate of the error variance 

-Estimate of the Q allele effect 

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
hap.trans.pere=data.test[[3]]
hap.chrom1.pere=data.test[[4]]
hap.chrom2.pere=data.test[[5]]
perf=data.test[[6]]
CD=data.test[[7]]
genea=data.test[[9]]

# In this example, marker positions are : 0, 0.010, 0.020, 0.030, 0.040, 0.050, 0.060,
# 0.070, 0.080, 0.090. 
# we use a 2 markers-associated haplotype

marq.hap=2

haplomax.HS=haplomax.HS.add(hap.chrom1.pere,hap.chrom2.pere,hap.trans.pere,hap.trans.mere,

perf,CD,genea,map,marq.hap)

haplomax.HS

}

\keyword{models}
