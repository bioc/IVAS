\name{MsqtlFinder}
\alias{MsqtlFinder}
\title{
  Find SQTLs in multiple genes.
}
\description{
  This function enables one to analyze multiple genes using multi-thread version of the foreach function and joins output results from sqtlfinder function. Moreover, it calculates the FDR using P-values of the matrix result data. 
}
\usage{
MsqtlFinder(expdata = NULL, snpdata = NULL, snplocus = NULL, GTFdata = NULL,
met = NULL, Ncor = 1, bplotout = NULL)
}
\arguments{
  \item{expdata}{
  Dataframe of expression data.
}
  \item{snpdata}{
  Dataframe of genotype data.
}
  \item{snplocus}{
  Locus of SNP markers in the snpdata.
}
  \item{GTFdata}{
  The transcriptDb object in the GnomicFeatures package.
}
  \item{met}{
  The option for statistical models.("lm" : analysis using linear regression model, "glm" : analysis using generalized linear mixed model,and "both" : "lm" and "glm").
}
  \item{Ncor}{
  The number of cores for multi-threads.
}
  \item{bplotout}{
  A directory saving boxplots
}
}
\value{
  This function returns the result matrix including SNP markers ID, chromosome number, alternative exons range, intron ranges, alternative type, P value, information of differential median values of expression ratio among genotypes ("sig" if differential median > 0.1 and "not sig" otherwise), gene names, methods ("lm" or "glm").
}
\references{
  Lawrence M, Huber W, Pages H, Aboyoun P, Carlson M, Gentleman R, Morgan M, and Carey V. Software for Computing and Annotating Genomic Ranges. PLoS Computational Biology, 9, e1003118. 2013.
  
  Benjamini, Yoav, Hochberg, and Yosef. Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society, Series B 57, 289-300. 1995.
}

\author{
Seonggyun Han, Sangsoo Kim
}
\seealso{
\code{\link[foreach]{foreach}},
\code{\link[GenomicRanges]{GRanges}},
\code{\link{p.adjust}}
}

\examples{
sampleDB <- system.file("extdata", "sampleDB", package="IVAS")
sample.Txdb <- loadDb(sampleDB)
data(sampleexp)
data(samplesnp)
data(samplesnplocus)
#final.result <- MsqtlFinder(sampleexp,samplesnp,samplesnplocus,sample.Txdb,"lm",1)
}