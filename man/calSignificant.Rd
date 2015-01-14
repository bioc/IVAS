\name{calSignificant}
\alias{calSignificant}
\title{
Calculates P.values by using two statistical models.
}
\description{
calculates P.values for association between expressions and genotypes by using the linear regression model and/or generalized linear mixed model.

}
\usage{
calSignificant(tx.gene=NULL, total.locus=NULL, exon.locus=NULL, intron.locus=NULL, 
info.strand=NULL, overapvalue=NULL, chrnum=NULL, expdata=NULL, snpdata=NULL, method=NULL)
}
\arguments{
  \item{tx.gene}{
  The matrix of transcripts including transcript IDs, Ensembl gene names, Ensembl transcript names, transcript start sites, and transcript end sites.
}
  \item{total.locus}{
  Ranges including alternative exons and flanking introns in a gene.
}
  \item{exon.locus}{
  All exon locus of a single gene.
}
  \item{intron.locus}{
  All intron locus of a single gene.
}
  \item{info.strand}{
  The strand informatrion of a single gene (forward strand = "+", reverse strand = "-").
}
  \item{overapvalue}{
  Snps located in the alternative exons and the flanking introns
}
  \item{chrnum}{
  The chromosome number of a single gene.
}
  \item{expdata}{
  Dataframe of expression data.
}
  \item{snpdata}{
  Datafame of genotype data.
}
  \item{method}{
  The option for statistical models and boxplot.("lm" : analysis using linear regression model, "glm" : analysis using generalized linear mixed model, "both" : "lm" and "glm", and "boxplot" : for writing boxplot).
}
}

\value{
  The lm or glm method returns matrix including; SNP marker IDs, Chromosome numbers, alternative exons ranges, Intron rangess, alternative types, P values, information of differential median values of expression ratio among genotypes ("sig" if differential median > 0.1 and "not sig" otherwise), gene names, and methods ("lm" or "glm").
  The boxplot method returns matrix with relative ratio values and genotypes of samples.
}
\references{
Chambers, J. M. (1992) Linear models. Chapter 4 of Statistical Models in S eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.
Breslow, N.E. Clayton, D.G. (1993). Approximate Inference in Generalized Linear Mixed Models. Journal of the American Statistical Association 88.
}

\author{
Seonggyun Han, Sangsoo Kim
}

\seealso{
\code{\link[IRanges]{findOverlaps}},
\code{\link[stats]{lm}},
\code{\link[lme4]{glmer}}
}
