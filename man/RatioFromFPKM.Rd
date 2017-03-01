\name{RatioFromFPKM}
\alias{RatioFromFPKM}
\title{
  Estimate relative expression ratio.
  }
\description{
  With the FPKM expression data set of transcripts, this function estimates relative expression ratio between transcripts with and without alternatively spliced exons based on splicing models of the ASdb object
  }
\usage{
  RatioFromFPKM(GTFdb = NULL, ASdb = NULL, 
  Total.expdata = NULL, CalIndex = NULL, Ncor = 1, out.dir = NULL)
  }
\arguments{
  \item{GTFdb}{
    A TxDb object in the \pkg{GenomicFeatures} package.
  }
  \item{ASdb}{
    A ASdb object including "SplicingModel" slot from the \code{\link{Splicingfinder}} function.
  }
  \item{Total.expdata}{
    A Data frame of expression data.
  }
  \item{CalIndex}{
    An index number of a splicing model in the ASdb object which will be calculated.
  }
  \item{Ncor}{
    The number of cores for multi-threads function.
  }
  \item{out.dir}{
    A directory saving results.
  }
}

\value{
  This function returns ASdb including the Ratio slot. The Ratio slot has expression ratios between transcripts with and without alternatively spliced exons for each alternative splicing type (i.e. Exon skipping, Alternative splice site, Intron retention). The values of the Ratio slot is as follow;
  \item{ES}{a data frame of splicing models for Exon skipping, including; index number, gene name, chromosome number, strand information of the gene, alternative exons range, downstream exon range, upstream exon range, splicing type, and ratios of individuals.}
  \item{ASS}{a data frame of splicing models for Alternative splice site, including; index number, gene name, chromosome number, strand information of the gene, alternative exons range, neighbor exon range, splicing type, and ratios of individuals.}
  \item{IR}{a data frame of splicing models for Intron retention, including; index number, gene name, chromosome number, strand information of the gene, alternative exons range, reteined exon range, splicing type, and ratios of individuals.}
  
}

\author{
  Seonggyun Han, Sangsoo Kim
}

\seealso{
  \code{\link[GenomicFeatures]{isActiveSeq}},
  \code{\link[GenomeInfoDb]{seqinfo}},
  \code{\link{Splicingfinder}}
}

\examples{
  sampleDB <- system.file("extdata", "sampleDB", package="IVAS")
  sample.Txdb <- loadDb(sampleDB)
  data(sampleexp)
  ASdb <- Splicingfinder(sample.Txdb)
  ASdb <- RatioFromFPKM(sample.Txdb,ASdb,sampleexp)
}