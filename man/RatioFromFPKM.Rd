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
        A data frame of expression data.
    }
    \item{CalIndex}{
        An index number in the ASdb object which will be tested in this function.
    }
    \item{Ncor}{
        The number of cores for multi-threads function.
    }
    \item{out.dir}{
        An output directory.
    }
}

\value{
        ASdb with the slot (labeled by "Ratio") containing results from the the \code{\link{RatioFromFPKM}} function. The "Ratio" slot contains a list object and each element of the list object returns the results assigned to three elements, which is of each alternative splicing type (i.e. Exon skipping, Alternative splice site, Intron retention). Three elements are as follows;
        \item{ES}{A data frame for the result of Exon skipping, consisting of the columns named as follows; Index (index number), EnsID (gene name), Nchr (chromosome name), 1stEX (alternatively spliced target exon), 2ndEX (second alternatively spliced target exon which is the other one of the mutually exclusive spliced exons), DownEX (downstream exon range), UpEX (upstream exon range), Types (splicing type), and names of individuals.}
        \item{ASS}{A data frame for the result of Alternative splice sites, consisting of the columns named as follows; Index (index number), EnsID (gene name), Nchr (chromosome name), ShortEX (shorter spliced target exon), LongEX (longer spliced target exon), NeighborEX (neighboring down or upstream exons), Types (splicing type), and names of individuals.}
        \item{IR}{A data frame for the result of Intron retention, consisting of the columns named as follows; Index (index number), EnsID (gene name), Nchr (chromosome name), RetainEX (retained intron exon), DownEX (downstream exon range), UpEX (upstream exon range), Types (splicing type), and names of individuals.}
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
