\name{saveBplot}
\alias{saveBplot}
\title{
 Save boxplots
}
\description{
  Save boxplots
}
\usage{
saveBplot(sig.sqtl = NULL, expdata = NULL, snpdata = NULL,
snplocus = NULL, GTFdata = NULL, outdir = NULL)

}
\arguments{
  \item{sig.sqtl}{
  A matrix of significant SQTLs from the sqtlfinder function
}
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
  \item{outdir}{
  A directory saving boxplots
}
}

\value{
  This function draws the boxplot
}

\author{
Seonggyun Han, Sangsoo Kim
}

\seealso{
\code{\link{boxplot}}
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
saveBplot(sqtl.result[1,],sampleexp,samplesnp,samplesnplocus,filtered.txdb,"./result")
}