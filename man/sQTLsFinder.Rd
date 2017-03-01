\name{sQTLsFinder}
\alias{sQTLsFinder}
\title{
  Find SQTLs.
  }
\description{
  Find significant SNPs using the calSignificant function.
  }
\usage{
  sQTLsFinder(ASdb = NULL , Total.snpdata = NULL , Total.snplocus = NULL , 
  GroupSam = NULL , method = "lm" , CalIndex = NULL , Ncor = 1 , out.dir = NULL)
  }
\arguments{
  \item{ASdb}{
    A ASdb object including "SplicingModel" and "Ratio" slots from the \code{\link{Splicingfinder}} and \code{\link{RatioFromFPKM}} funtions, respectively.
    }
  \item{Total.snpdata}{
    A data frame of genotype data.
    }
  \item{Total.snplocus}{
    A data frame containing locus information of SNP markers in the snpdata.
    }
  \item{GroupSam}{
    A list object of a conditions for each individual. If GroupSam is not NULL ,the odds ratio and its confidence intervals between conditions are calculated.
    }
  \item{method}{
    An option for statistical models and boxplot.("lm" : analysis using linear regression model, "glm" : analysis using generalized linear mixed model, "both" : "lm" and "glm",and "boxplot" : for writing boxplot).
    }
  \item{CalIndex}{
    An index number of a splicing model in the ASdb object which will be calculated.
    }
  \item{Ncor}{
    The number of cores for multi-threads function.
    }
  \item{out.dir}{
    A directory saving results
    }
  }

\value{
  The lm or glm method returns ASdb including the sQTLs slot. The SplicingModel slot has alternative splicing model for each alternative splicing type (i.e. Exon skipping, Alternative splice site, Intron retention). The values of the Ratio slot is as follow;
  \item{ES}{a data frame of splicing models for Exon skipping, including; SNP markers ID, index number, gene name, strand information of the gene, chromosome number, alternative exons range, neighbor exon range, splicing type, P value, information of differential median values of expression ratio among genotypes ("Diff" if differential median > 0.1 and "NonDiff" otherwise), methods ("lm" or "glm").}
  \item{ASS}{a data frame of splicing models for Alternative splice site, including; SNP markers ID, index number, gene name, strand information of the gene, chromosome number, alternative exons range, neighbor exon range, splicing type, P value, information of differential median values of expression ratio among genotypes ("Diff" if differential median > 0.1 and "NonDiff" otherwise), methods ("lm" or "glm").}
  \item{IR}{a data frame of splicing models for Intron retention, including; SNP markers ID, index number, gene name, strand information of the gene, chromosome number, alternative exons range, neighbor exon range, splicing type, P value, information of differential median values of expression ratio among genotypes ("Diff" if differential median > 0.1 and "NonDiff" otherwise), methods ("lm" or "glm").}
  The boxplot method returns matrix data with relative ratio values and genotypes of samples.
  }

\references{
  Chambers, J. M. (1992) Linear models. Chapter 4 of Statistical Models in S eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.
  Breslow, N.E. Clayton, D.G. (1993). Approximate Inference in Generalized Linear Mixed Models. Journal of the American Statistical Association 88.
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
  ASdb <- RatioFromFPKM(sample.Txdb,ASdb,sampleexp)
  ASdb <- sQTLsFinder(ASdb,samplesnp,samplesnplocus,method="lm")
  }