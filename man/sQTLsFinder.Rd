\name{sQTLsFinder}
\alias{sQTLsFinder}
\title{
  Find SQTLs.
  }
\description{
  Find significant SNPs using the calSignificant function.
  }
\usage{
  sQTLsFinder(ASdb, Total.snpdata = NULL , Total.snplocus = NULL , 
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
    A list object of a conditions for each individual. If GroupSam is not NULL, the odds ratio and its confidence intervals between conditions are calculated.
    }
  \item{method}{
    An option for statistical models and boxplot.("lm" : analysis using linear regression model, "glm" : analysis using generalized linear mixed model, "both" : "lm" and "glm",and "boxplot" : for writing boxplot).
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
  ASdb with the slot (labeled by "sQTLs") containing results from the the \code{\link{sQTLsFinder}} function. The "Splicingfinder" slot contains a list object and each element of the list object returns the results assigned to three elements, which is of each alternative splicing type (i.e. Exon skipping, Alternative splice site, Intron retention). Three elements are as follows;
  \item{ES}{A data frame for the result of Exon skipping, consisting of the columns named as follows; Index (index number), EnsID (gene name), Nchr (chromosome name), 1stEX (alternatively spliced target exon), 2ndEX (second alternatively spliced target exon which is the other one of the mutually exclusive spliced exons), DownEX (downstream exon range), UpEX (upstream exon range), Types (splicing type), pByGeno (P-values of "lm" or "glm" test for association PSI values and genotypes), FdrByGeno (pByGeno), diff ("diff" if differential median > 0.1 and "Nondiff" otherwise), pByGroups (P-values of chi-square test for association between genotypes of two groups), fdrByGroups (FDR values for the pByGroups column), OR (odds ratio), lowCI(low confidence interval), highCI(high confidence interval), and methods ("lm" or "glm").}
  \item{ASS}{A data frame for the result of Alternative splice sites, consisting of the columns named as follows; Index (index number), EnsID (gene name), Nchr (chromosome nam), ShortEX (shorter spliced target exon), LongEX (longer spliced target exon), NeighborEX (neighboring down or upstream exons), Types (splicing type), pByGeno (P-values of "lm" or "glm" test for association PSI values and genotypes), FdrByGeno (pByGeno), diff ("diff" if differential median > 0.1 and "Nondiff" otherwise), pByGroups (P-values of chi-square test for association between genotypes of two groups), fdrByGroups (FDR values for the pByGroups column), OR (odds ratio), lowCI(low confidence interval), highCI(high confidence interval), and methods ("lm" or "glm").}
  \item{IR}{A data frame for the result of Intron retention, consisting of the columns named as follows; Index (index number), EnsID (gene name), Nchr (chromosome name), RetainEX (retained intron exon), DownEX (downstream exon range), UpEX (upstream exon range), Types (splicing type), pByGeno (P-values of "lm" or "glm" test for association PSI values and genotypes), FdrByGeno (pByGeno), diff ("diff" if differential median > 0.1 and "Nondiff" otherwise), pByGroups (P-values of chi-square test for association between genotypes of two groups), fdrByGroups (FDR values for the pByGroups column), OR (odds ratio), lowCI(low confidence interval), highCI(high confidence interval), and methods ("lm" or "glm").}
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