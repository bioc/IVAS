\name{ASdb-class}
\alias{ASdb-class}
\docType{class}
\title{
    ASdb s4 class - a container for results from functions of the \pkg{IVAS} package.
}
\description{
    This class is the main object for storing results of the present package.
}

\seealso{
    \code{\link{Splicingfinder}},
    \code{\link{RatioFromFPKM}},
    \code{\link{sQTLsFinder}}
}

\note{
    An ASdb object stores information of alternative splicing patterns, expression ratios between transcripts with and without alternative target exons, and significant sQTLs from the functions of the present package. This ASdb object can be populated further slots during the analysis using functions for the analysis. Typically, an ASdb object can be created when the function \code{\link{Splicingfinder}} completes to define alternative splicing patterns. After creation, the ASdb contains the slot labeled as "SplicingModel", and the slot includes a list object named by "ES", "ASS", and "IR" (alternative splicing exons are saved separately in each element of the list based on their splicing pattern types; "ES": Exon skipping, "ASS": Alternative splice site, and "IR": Intron retention). In the next analysis step, further result slots can be added. The function \code{\link{RatioFromFPKM}} can add the "Ratio" slot containing expression ratio for each alternative splicing pattern based on the "SplicingModel" slot of the present class and for each individual from a matrix of FPKM values. Then, the result of the \code{\link{sQTLsFinder}} function can be saved by adding the "sQTLs" slot including significance of association between the expression ratios, which is stored in the "Ratio" slot of the present class, and SNPs for each alternative splicing exon.
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
    ASdb
}
