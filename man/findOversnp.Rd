\name{findOversnp}
\alias{findOversnp}
\title{
  Find SNPs which belongs to alternative exons and flanking introns of them.
}
\description{
    Find SNPs which belong to alternative exons and flanking introns of them.
}
\usage{
findOversnp(altInvalue = NULL, snprange = NULL)
}
\arguments{
  \item{altInvalue}{
  A list data set from the findAlternative function.
}
  \item{snprange}{
  A matrix of SNP ranges.
}
}

\value{
This function returns a matrix with SNPs in alternative exons and flanking introns and ranges of those SNPs.
}

\author{
Seonggyun Han, Sangsoo Kim
}
\seealso{
\code{\link[IRanges]{findOverlaps}}
}

\examples{
sampleDB <- system.file("extdata", "sampleDB", package="IVAS")
sample.Txdb <- loadDb(sampleDB)
data(samplesnplocus)
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
}