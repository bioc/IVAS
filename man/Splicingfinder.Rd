\name{Splicingfinder}
\alias{Splicingfinder}
\title{
  Find alternatively spliced exons
  }
\description{
  Find alternatively spliced exons based on GTF reference transcript models.
  }
\usage{
  Splicingfinder(GTFdb = NULL , txTable = NULL , calGene = NULL , Ncor = 1 , out.dir = NULL)
  }
\arguments{
  \item{GTFdb}{
    A TxDb object in the \pkg{GenomicFeatures} package.
    }
  \item{txTable}{
    A matrix of transcripts including transcript IDs, gene names, transcript names, transcript start sites, and transcript end sites based on a GTF reference transcript model file.
    }
  \item{calGene}{
    An interest of a gene that will be tested. If calGene is inputted by a single gene, the splicing pattern for the only gene are tested. if not, the splicing patterns for total of genes are tested 
    }
  \item{Ncor}{
    The number of cores for multi-threads.
    }
  \item{out.dir}{
    A directory saving results.
    }
  }

\value{
  This function returns ASdb including the SplicingModel slot. The SplicingModel slot has alternative splicing model for each alternative splicing type (i.e. Exon skipping, Alternative splice site, Intron retention). The values of the Ratio slot is as follow;
  \item{ES}{a data frame of splicing models for Exon skipping, including; index number, gene name, chromosome number, strand information of the gene, alternative exons range, downstream exon range, upstream exon range, and splicing type.}
  \item{ASS}{a data frame of splicing models for Alternative splice site, including; index number, gene name, chromosome number, strand information of the gene, alternative exons range, neighbor exon range, and splicing type.}
  \item{IR}{a data frame of splicing models for Intron retention, including; index number, gene name, chromosome number, strand information of the gene, alternative exons range, reteined exon range, and splicing type.}
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
  ASdb <- Splicingfinder(sample.Txdb)
  }