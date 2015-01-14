\name{sqtlfinder}
\alias{sqtlfinder}
\title{
  Find SQTLs.
}
\description{
  Find significant SNPs using the calSignificant function.
}
\usage{
sqtlfinder(altInvalue = NULL, overapvalue = NULL, expdata = NULL, snpdata = NULL, method = NULL)
}
\arguments{
  \item{altInvalue}{
  A list data set from the findAlternative function.
}
  \item{overapvalue}{
  A matrix data with SNPs in the flanking introns of alternative exons and ranges of those SNPs from findOversnp function.
}
  \item{expdata}{
  Expression data of samples.
}
  \item{snpdata}{
  Genotype data of samples.
}
  \item{method}{
  The option for statistical models and boxplot.("lm" : analysis using linear regression model, "glm" : analysis using generalized linear mixed model, "both" : "lm" and "glm",and "boxplot" : for writing boxplot).
}
}
\value{
  The lm or glm method returns matrix data including SNP markers ID, chromosome number, alternative exons range, intron ranges, alternative type, P value, information of differential median values of expression ratio among genotypes ("sig" if differential median > 0.1 and "not sig" otherwise), a gene name, methods ("lm" or "glm"),and strand information of the gene.
  The boxplot method returns matrix data with relative ratio values and genotypes of samples.
}
\author{
Seonggyun Han, Sangsoo Kim
}

\examples{
sampleDB <- system.file("extdata", "sampleDB", package="IVAS")
sample.Txdb <- loadDb(sampleDB)
data(samplesnplocus)
data(sampleexp)
data(samplesnp)
filtered.txdb <- chrseparate(sample.Txdb,6)
trans.exon.range <- exonsBy(filtered.txdb,by="tx")
trans.intron.range <- intronsByTranscript(filtered.txdb)
txTable <- select(filtered.txdb, keys=names(trans.exon.range),
columns=c("TXID","TXNAME","GENEID","TXSTART","TXEND"), keytype="TXID")
ch.snp.locus <- as.matrix(samplesnplocus[samplesnplocus[,2] == 6,])
ch.snps <- matrix(ch.snp.locus[is.element(ch.snp.locus[,1],rownames(samplesnp)),],ncol=3,byrow=FALSE)
ch.snps.range <- GRanges(seqnames=Rle(6),ranges=IRanges(start=as.integer(ch.snps[,3]),
end=as.integer(ch.snps[,3])),metadata=ch.snps[,1])
Altvalue <- findAlternative("ENSG00000196735",txTable,trans.exon.range,trans.intron.range,6)
overlapsnp <- findOversnp(Altvalue,ch.snps.range)
sqtl.result <- sqtlfinder(Altvalue,overlapsnp,sampleexp,samplesnp,"lm")
}