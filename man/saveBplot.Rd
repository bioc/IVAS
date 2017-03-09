\name{saveBplot}
\alias{saveBplot}
\title{
    Save boxplots
}
\description{
    Save boxplots
}
\usage{
    saveBplot(ASdb=ASdb,Total.snpdata=NULL,Total.snplocus=NULL,
                CalIndex=NULL,out.dir=NULL)
}
\arguments{
    \item{ASdb}{
        A ASdb object including "sQTLs" slot from the \code{\link{sQTLsFinder}} function.
    }
    \item{Total.snpdata}{
        A data frame of genotype data.
    }
    \item{Total.snplocus}{
        A data frame containing locus information of SNP markers in the snpdata.
    }
    \item{CalIndex}{
        n index number in the ASdb object which will be tested in this function.
    }
    \item{out.dir}{
        An output directory.
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
    data(sampleexp)
    data(samplesnp)
    data(samplesnplocus)
    ASdb <- Splicingfinder(sample.Txdb)
    ASdb <- RatioFromFPKM(sample.Txdb,ASdb,sampleexp)
    ASdb <- sQTLsFinder(ASdb,samplesnp,samplesnplocus,method="lm")
    saveBplot(ASdb=ASdb,Total.snpdata=samplesnp,Total.snplocus=samplesnplocus,CalIndex="ASS7",out.dir="./result")
}
