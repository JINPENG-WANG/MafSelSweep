\name{MafSelSweep-package}
\alias{MafSelSweep-package}
\alias{MafSelSweep}
\docType{package}
\title{
\packageTitle{MafSelSweep}
}
\description{
\packageDescription{MafSelSweep}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{MafSelSweep}
\packageIndices{MafSelSweep}
find.maf.block(snp.data,maf.threshold,snp.number,out.file.name)
}
\author{
\packageAuthor{MafSelSweep}
Jinming Huang, Jinpeng Wang
Maintainer: \packageMaintainer{MafSelSweep}
Jinpeng Wang wangjinpeng0225@163.com
}
\keyword{ package }

\examples{
define.selsweep.region(snp.data,size,out.file.name)
snp.data is an object in the class of "data.frame". It must contain 7 columns.
The first column contains the Chromosome number.
The second column contains the SNP identifiers which are unique for every SNP.
The third and fourth columns contain some characters.
The fifth column contains the frequency of this SNP in the population.
The sixth column contains scores of  SNPs.
The seventh column contains the position of the SNP on the chromosome.
}
