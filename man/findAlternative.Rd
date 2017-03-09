\name{findAlternative}
\alias{findAlternative}
\title{
    Find alternative exons of a gene.
}
\description{
    Search alternative exons among transcript isoforms from a single gene.
}
\usage{
    findAlternative(geneid = NULL, txTable = NULL, totalExrange = NULL,
        totalInrange = NULL, one.chr = NULL)
}
\arguments{
    \item{geneid}{
        Ensembl gene name.
    }
    \item{txTable}{
        The matrix of transcripts including transcript IDs, Ensembl gene names, Ensembl transcript names, transcript start sites, and transcript end sites.
    }
    \item{totalExrange}{
        A list of GRanges objects including total exon ranges in each transcript resulted from the \code{\link{exonsBy}} function in \pkg{GenomicFeatures}.
    }
    \item{totalInrange}{
        A list of GRanges objects including total intron ranges in each transcript resulted from the \code{\link{intronsByTranscript}} function in \pkg{GenomicFeatures}.
    }
    \item{one.chr}{
        The chromosome number that you would like to test
    }
}

\value{
    \item{alterIntron}{A GRanges object with flanking introns of alternative exons}
    \item{tableBygene}{An information table of transcripts including transcript IDs, Ensembl gene names, Ensembl transcript names, transcript start sites, and transcript end sites.}
    \item{exonRange}{All exons locus of a gene}
    \item{intronRange}{All intron locus of a gene}
}
\references{
    Lawrence M, Huber W, Pages H, Aboyoun P, Carlson M, Gentleman R, Morgan M, and Carey V. Software for Computing and Annotating Genomic Ranges. PLoS Computational Biology, 9, e1003118. 2013.
}
\author{
    Seonggyun Han, Sangsoo Kim
}
\seealso{
    \code{\link[GenomicRanges]{GRanges}},
    \code{\link[IRanges]{IRanges}}
}

\examples{
    sampleDB <- system.file("extdata", "sampleDB", package="IVAS")
    sample.Txdb <- loadDb(sampleDB)
    filtered.txdb <- chrseparate(sample.Txdb,19)
    trans.exon.range <- exonsBy(filtered.txdb,by="tx")
    trans.intron.range <- intronsByTranscript(filtered.txdb)
    txTable <- select(filtered.txdb, keys=names(trans.exon.range),
                                        columns=c("TXID","TXNAME","GENEID","TXSTART","TXEND"), keytype="TXID")
    Altvalue <- findAlternative("ENSG00000170889",txTable,trans.exon.range,trans.intron.range,19)
}
