\name{data.test}

\alias{data.test}

\docType{data}

\title{data.test}

\description{
"data.test" is a set of data.It is a list containing 9 elements.
We consider a 10 cM chromosomal region with 10 equally spaced biallelic markers.So,map density is equal to 1 cM.
We assume a set of 5 sires, each sire having 10 sons with 10 different dams.
}

\usage{data(data.test)}


\format{
  A list containing the following components:
  \item \code{map}: vector (9) with numerical values
  \item \code{hap.trans.mere}: matrix (50*10) with character values coded with "A","B","NA"
  \item \code{hap.trans.pere}: matrix (50*10) with character values coded with "A","B","NA"
  \item \code{hap.chrom1.pere}: matrix (5*10) with character values coded with "A","B"
  \item \code{hap.chrom2.pere}: matrix (5*10) with character values coded with "A","B"
  \item \code{perf}: vector (9) with numerical values
  \item \code{CD}: vector (50) with numerical values
  \item \code{PLA}: vector (50) with numerical values
  \item \code{genea}: matrix (50*2) with numerical values
}

\examples{

data(data.test)

#distance between two consecutive markers on all the chromosome
 map=data.test[[1]]
 map

#haplotype transmitted by the dam
 hap.trans.mere=data.test[[2]]
 hap.trans.mere

#haplotype transmitted by the sire
 hap.trans.pere=data.test[[3]]
 hap.trans.pere

#haplotype of the first chromosome for each sire
 hap.chrom1.pere=data.test[[4]]
 hap.chrom1.pere

#haplotype of the second chromosome for each sire
 hap.chrom2.pere=data.test[[5]]
 hap.chrom2.pere

#performances of sons
 perf=data.test[[6]]
 perf

#CD of sons
 CD=data.test[[7]]
 CD

#transmission probabilities of sons
 PLA=data.test[[8]]
 PLA

#genealogy of sons
 genea=data.test[[9]]
 genea

}

\keyword{datasets}


