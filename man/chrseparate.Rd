\name{chrseparate}
\alias{chrseparate}
\title{
Separate a TranscriptDb object based on a chromosome.
}
\description{
With the isActiveSeq method in GenomicFeatures package, this function filters the transcriptDb object in the GenomicFeatures package based on a single chromosome.
}
\usage{
chrseparate(transdb = NULL, chrname = NULL)
}
\arguments{
  \item{transdb}{
  The transcriptDb object in the GnomicFeatures package.
}
  \item{chrname}{
  The chromosome number you would like to selecet from TrnascriptDb
}
}

\value{
This function returns the TrnascriptDb limited to the chromosome number that you want.
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
filtered.txdb <- chrseparate(sample.Txdb,6)
}