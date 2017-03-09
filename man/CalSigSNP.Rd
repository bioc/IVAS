\name{CalSigSNP}
\alias{CalSigSNP}
\title{
    Calculate significance SNPs
}
\description{
    This function performs linear regression test to identify significance associations between expression ratio and genotypes using the \code{\link{lm}} function.
}
\usage{
    CalSigSNP(ratio.mat=NULL,snp.mat=NULL,overlapsnp=NULL,
        each.snplocus=NULL,chr,each.gene=NULL,GroupSam=NULL,method="lm")
}
\arguments{
    \item{ratio.mat}{
        A data frame consisting of expression ratio of an alternatively spliced exon.
    }
    \item{snp.mat}{
        A data frame of genotype data.
    }
    \item{overlapsnp}{
        A data frame containing SNPs which is within an alternatively spliced exon and its flanking introns.
    }
    \item{each.snplocus}{
        A data frame consisting of locus information of SNP markers in the snpdata.
    }
    \item{chr}{
        The chromosome number that you would like to test in this function.
    }
    \item{each.gene}{
        The gene name that you would like to test in this function
    }
    \item{GroupSam}{
        A list object of a group of each sample.
    }
    \item{method}{
        The option for statistical models and boxplot.("lm" : analysis using linear regression model, "glm" : analysis using generalized linear mixed model, "both" : "lm" and "glm", and "boxplot" : for writing boxplot).
    }
}

\value{
    The lm or glm method returns matrix including; SNP marker IDs, P values, information of differential median values of expression ratio among genotypes ("sig" if differential median > 0.1 and "not sig" otherwise), gene names, and methods ("lm" or "glm").
    The boxplot method returns matrix with relative ratio values and genotypes of samples.
}
\references{
    Chambers, J. M. (1992) Linear models. Chapter 4 of Statistical Models in S eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.
}

\author{
    Seonggyun Han, Sangsoo Kim
}

\seealso{
    \code{\link[stats]{lm}},
    \code{\link[lme4]{glmer}}
}

\examples{
    sampleDB <- system.file("extdata", "sampleDB", package="IVAS")
    sample.Txdb <- loadDb(sampleDB)
    data(sampleexp)
    data(samplesnp)
    data(samplesnplocus)
    ASdb <- Splicingfinder(sample.Txdb)
    ASdb <- RatioFromFPKM(sample.Txdb,ASdb,sampleexp,CalIndex="ASS7")
    ratio.mat <- slot(ASdb,"Ratio")$ASS
    ratio.mat <- rbind(ratio.mat[,grep("NA",colnames(ratio.mat))])
    each.snp <- rbind(samplesnp[rownames(samplesnp) == "rs3810232",])
    each.snplocus <- rbind(samplesnplocus[samplesnplocus[,"SNP"] == "rs3810232",])
    overlapsnp <- rbind(c(snp="rs3810232",locus="54704760"))
    CalSigSNP(ratio.mat,as.matrix(each.snp),overlapsnp,each.snplocus,"19","ENSG00000170889",method="lm")
}




