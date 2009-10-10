\name{data.test}

\alias{data.test}

\docType{data}

\title{data.test}

\description{
"data.test" is a list of 9 elements.
We consider a 10 cM chromosomal region with 10 equally spaced biallelic markers and a design of 5 half-sib families, each sire having 10 sons.
Data are prepared for a QTL detection in the middle of each marker intervalle.
}

\usage{data(data.test)}


\format{
  A list containing the following components:

 - \code{map}: vector (9) with numerical values

 - \code{hap.trans.mere}: matrix (50 x 10) with character values 

 - \code{hap.trans.pere}: matrix (50 x 10) with character values 

 - \code{hap.chrom1.pere}: matrix (5 x 10) with character values 

 - \code{hap.chrom2.pere}: matrix (5 x 10) with character values 

 - \code{perf}: vector (50) with numerical values
 
 - \code{CD}: vector (50) with numerical values

 - \code{PLA}: matrix (50 x 9) with numerical values

 - \code{genea}: matrix (50 x 2) with numerical values
}

\examples{

data(data.test)

#distance between two consecutive markers on the chromosome
 map=data.test[[1]]
 map

#haplotype transmitted by dams
#son information (lines) are ordered following genea[,1]
 hap.trans.mere=data.test[[2]]
 hap.trans.mere

#haplotype transmitted by sires
#son information (lines) are ordered following genea[,1]
 hap.trans.pere=data.test[[3]]
 hap.trans.pere

#haplotype of the first chromosome for each sire
#sire information (lines) are ordered following unique(genea[,2])
 hap.chrom1.pere=data.test[[4]]
 hap.chrom1.pere

#haplotype of the second chromosome for each sire
#sire information (lines) are ordered following unique(genea[,2])
 hap.chrom2.pere=data.test[[5]]
 hap.chrom2.pere

#performances of sons
#son information are ordered following genea[,1]
 perf=data.test[[6]]
 perf

#CD of sons
#son information are ordered following genea[,1]
 CD=data.test[[7]]
 CD

#transmission probabilities of sons for each test position
#son information (lines) are ordered following genea[,1]
 PLA=data.test[[8]]
 PLA

#son index and  index of his father
 genea=data.test[[9]]
 genea

}

\keyword{datasets}


