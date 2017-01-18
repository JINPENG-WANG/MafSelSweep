\name{find.maf.block}
\alias{find.maf.block}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
define selective sweep regions 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
find.maf.block(snp.data, maf.threshold, snp.number, out.file.name)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snp.data}{
%%     ~~Describe \code{snp.data} here~~
}
  \item{maf.threshold}{
%%     ~~Describe \code{maf.threshold} here~~
}
  \item{snp.number}{
%%     ~~Describe \code{snp.number} here~~
}
  \item{out.file.name}{
%%     ~~Describe \code{out.file.name} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (snp.data, maf.threshold, snp.number, out.file.name) 
{
    if (missing(snp.data)) {
        stop("An data.frame object containing SNP data must be provided!")
    }
    else {
        column_count = ncol(snp.data)
        if (column_count != 7) {
            stop("Input data must be in the required format. See the example input file!")
        }
    }
    if (missing(maf.threshold)) {
        maf.threshold <- 0.05
    }
    if (missing(snp.number)) {
        snp.number <- 5
    }
    if (missing(out.file.name)) {
        out.file.name = "MafSelSweep.out"
    }
    chroms <- snp.data[, 1]
    chroms <- chroms[!duplicated(chroms)]
    for (chr in chroms) {
        tyu <- subset(snp.data, snp.data[, 1] == chr)
        row.num <- length(tyu[, 1])
        tyu$flag <- tyu[, 5] <= maf.threshold
        tyu$flag[is.na(tyu$flag)] <- FALSE
        tyu$chromosome <- NA
        tyu$start <- NA
        tyu$end <- NA
        tyu$blocksize <- NA
        tyu$snpnumber <- NA
        count1 = 1
        tyu
        while (count1 <= row.num) {
            if (tyu$flag[count1]) {
                sum = 1
                count2 = count1 + 1
                condition1 = 1
                while (condition1) {
                  if (tyu$flag[count2]) {
                    sum = sum + 1
                    count2 = count2 + 1
                  }
                  else {
                    condition1 = 0
                    chr.name = paste("chr", chr, sep = "")
                    start = tyu[count1, 7]
                    end = tyu[count2 - 1, 7]
                    blocksize = end - start + 1
                    for (index in count1:(count2 - 1)) {
                      if (sum >= snp.number) {
                        tyu$chromosome[index] = chr.name
                        tyu$start[index] = start
                        tyu$end[index] = end
                        tyu$blocksize[index] = blocksize
                        tyu$snpnumber[index] = sum
                      }
                    }
                  }
                }
                count1 = count2
            }
            else {
                count1 = count1 + 1
            }
        }
        snp.data$Chromosome[snp.data[, 1] == chr] <- tyu$chromosome
        snp.data$Start[snp.data[, 1] == chr] <- tyu$start
        snp.data$End[snp.data[, 1] == chr] <- tyu$end
        snp.data$Block_size[snp.data[, 1] == chr] <- tyu$blocksize
        snp.data$SNP_number[snp.data[, 1] == chr] <- tyu$snpnumber
    }
    write.table(snp.data, file = out.file.name, sep = "\t", row.names = F, 
        quote = F)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
