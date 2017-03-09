\name{chrseparate}
\alias{chrseparate}
\title{
    Separate a TxDb object based on a chromosome.
}
\description{
    With the \code{\link{isActiveSeq}} method in \pkg{GenomicFeatures} package, this function filters the TxDb object in the \pkg{GenomicFeatures} package based on a single chromosome.
}
\usage{
    chrseparate(GTFdb = NULL, chrname = NULL)
}
\arguments{
    \item{GTFdb}{
        The TxDb object in the \pkg{GenomicFeatures} package.
    }
    \item{chrname}{
        The chromosome number you would like to select from TxDb
    }
}

\value{
    This function returns the TxDb limited to the chromosome number that you want.
}
\references{
    Lawrence M, Huber W, Pages H, Aboyoun P, Carlson M, Gentleman R, Morgan M, and Carey V. Software for Computing and Annotating Genomic Ranges. PLoS Computational Biology, 9, e1003118. 2013.
}
\author{
    Seonggyun Han, Sangsoo Kim
}
\seealso{
    \code{\link[GenomicFeatures]{isActiveSeq}},
    \code{\link[GenomeInfoDb]{seqinfo}}
}
\examples{
    sampleDB <- system.file("extdata", "sampleDB", package="IVAS")
    sample.Txdb <- loadDb(sampleDB)
    filtered.txdb <- chrseparate(sample.Txdb,19)
}