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
    An interest of a gene that will be tested. If calGene is inputted by a single gene, the splicing pattern for the only gene is tested. If not, the splicing patterns for total of genes are tested 
    }
  \item{Ncor}{
    The number of cores for multi-threads.
    }
  \item{out.dir}{
    An output directory.
    }
  }

\value{
  ASdb with the slot (labeled by "SplicingModel") containing results from the the \code{\link{Splicingfinder}} function. The "Splicingfinder" slot contains a list object and each element of the list object returns the results assigned to three elements, which is of each alternative splicing type (i.e. Exon skipping, Alternative splice site, Intron retention). Three elements are as follows;
  \item{ES}{A data frame for the result of Exon skipping, consisting of the columns named as follows; Index (index number), EnsID (gene name), Nchr (chromosome name), 1stEX (alternatively spliced target exon), 2ndEX (second alternatively spliced target exon which is the other one of the mutually exclusive spliced exons), DownEX (downstream exon range), UpEX (upstream exon range), 1st_des (alternatively spliced target exons in a representative exon), 2nd_des (second alternatively spliced target exons in a representative exon), Do_des (downstream exons in a representative exon), Up_des (upstream exons in a representative exon), and Types (splicing type).}
  \item{ASS}{A data frame for the result of Alternative splice site, consisting of the columns named as follows; Index (index number), EnsID (gene name), Nchr (chromosome nam), ShortEX (shorter spliced target exon), LongEX (longer spliced target exon), NeighborEX (neighboring down or upstream exons), Short_des (shorter spliced target exons in a representative exon), Long_des (longer spliced target exons in a representative exon), Neighbor_des (neighboring down or upstream exons in a representative exon), and Types (splicing type).}
  \item{IR}{A data frame for the result of Intron retention, consisting of the columns named as follows; Index (index number), EnsID (gene name), Nchr (chromosome name), RetainEX (retained intron exon), DownEX (downstream exon range), UpEX (upstream exon range), Retain_des (retained intron exons in a representative exon), Do_des (downstream exons in a representative exon), Up_des (upstream exons in a representative exon), and Types (splicing type).}
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
