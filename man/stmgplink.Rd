\name{stmgplink}
\alias{stmgplink}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Smooth-threshold multivariate genetic prediction for genome-wide association or whole-genome sequencing data in PLINK format
}
\description{
Build prediction model from training data and predict test data phenotype through smooth-threshold multivariate genetic prediction (STMGP) method.
Data must be in PLINK binary format and marginal test p-values (i.e. test for each variant) are computed by PLINK software,
which enables rapid computation even for data having very large number of variants.
An optimal p-value cutoff is selected by Cp-type criterion.
Both quantitative and binary phenotypes are acceptable, in which data must be in PLINK fam file format.
}
\usage{
stmgplink(trainbed, testbed = NULL, gamma = 1, taun = NULL,
  lambda = 1, Z = NULL, Zte = NULL, plink = "plink --noweb",
  maf = 0.01, hwe = 1e-04, geno = 0.1, fout = "stp",
  trainfout = "train", testfout = "test", ll = 50, maxal = NULL, alc = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trainbed}{
A training data file name in PLINK-binary format;
Binary phenotype must be coded as 1 or 2 in PLINK fam file.
}
  \item{testbed}{
A test data file name in PLINK-binary format; NULL (default) means unspecified.
}
  \item{gamma}{
gamma parameter; \code{gamma=1} is default as suggested in Ueki and Tamiya (2016).
}
  \item{taun}{
tau parameter divided by (sample size) (allowed to be a vector object); NULL (default) specifies \code{tau=n/log(n)} as suggested in Ueki and Tamiya (2016).
}
  \item{lambda}{
lambda parameter (default=1).
}
  \item{Z}{
A covariate file name for training data (PLINK's covariate format); NULL (default) means unspecified.
}
  \item{Zte}{
A covariate file name for test data (PLINK's covariate format); NULL (default) means unspecified.
}
  \item{plink}{
PLINK command, e.g. "plink2", or "./plink --noweb" (default is \code{plink --noweb}); PLINK must be installed.
}
  \item{maf}{
Minor allele frequency (MAF) cutoff for \code{--maf} option in PLINK.
}
  \item{hwe}{
Hardy-Weinberg equilibrium (HWE) cutoff for \code{--hwe} option in PLINK.
}
  \item{geno}{
Missing call rate cutoff for \code{--geno} option in PLINK (default=0.1).
}
  \item{fout}{
An output file name (default="stp").
}
  \item{trainfout}{
An output file name for training data (default="train").
}
  \item{testfout}{
An output file name for test data (default="test").
}
  \item{ll}{
Number of candidate p-value cutoffs for search (default=50) as determined by
  \code{10^seq( log10(maxal),log10(5e-8), length=ll)}.
}
  \item{maxal}{
Maximum p-value cutoff for search (default=NULL); If most variants are null \code{maxal}*(number of variants) gives approximate number of filtered variants, which is useful for rapid computation even for data with large number of variants.
}
  \item{alc}{
User-specified candidate p-value cutoffs for search; \code{ll} option is effective if \code{alc=NULL}.
}
}
\details{
See Ueki and Tamiya (2016).
}
\value{
  \item{Muhat}{Estimated phenotypes from linear model evaluated at each candidate tuning parameters (\code{al} and \code{tau}) whose size is of (sample size) x (length of \code{al}) x (length of \code{tau}).}
  \item{gdf}{Generalized degrees of freedom (GDF, Ye 1998) whose size is of (length of \code{al}) x (length of \code{tau}).}
  \item{sig2hat}{Error variance estimates (=1 for binary traits) whose size is of (length of \code{al}) x (length of \code{tau}).}
  \item{df}{Number of nonzero regression coefficients whose size is of (length of \code{al}) x (length of \code{tau}).}
  \item{al}{Candidate p-value cutoffs for search.}
  \item{lopt}{An optimal tuning parameter indexes for \code{al} and \code{gamma} selected by Cp-type criterion, \code{CP}}
  \item{BA}{Estimated regression coefficient matrix whose size is of (1 + number of columns of \code{Z} + number of columns of \code{X}) x (length of \code{al})) x (length of \code{tau})); the first element, the second block and third block correspond to intercept, \code{Z} and \code{X}, respectively.}
  \item{Loss}{Loss (sum of squared residuals or -2*loglikelihood) whose size is of (length of \code{al}) x (length of \code{tau}).}
  \item{sig2hato}{An error variance estimate (=1 for binary traits) in computing the variance term of Cp-type criterion.}
  \item{tau}{Candidate tau parameters for search.}
  \item{CP}{Cp-type criterion whose size is of (length of \code{al}) x (length of \code{tau}).}
  \item{PE}{PLINK fam file for training data with additional column including the predicted phenotype from linear model.}
  \item{PEte}{PLINK fam file for test data with additional column including the predicted phenotype from linear model estimated from training data.}
  \item{nonzero}{Variants with nonzero regression coefficients at the optimal parameter in PLINK file.}
}
\references{

Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira M, Bender D, Maller J, Sklar P, de Bakker P, Daly MJ, Sham PC. (2007) PLINK: A tool set for whole-genome and population-based linkage analyses. Am J Hum Genet 81:559-75.

Chang CC, Chow CC, Tellier LCAM, Vattikuti S, Purcell SM, Lee JJ. (2015) Second-generation PLINK: rising to the challenge of larger and richer datasets. GigaScience 4.

}
\examples{

\dontrun{
wd = system.file("extdata",package="stmgp")

# quantitive traits
# training data (plink format)
trainbed = paste(wd,"train",sep="/")
# test data (plink format)
testbed = paste(wd,"test",sep="/")
# number of SNPs
#n.snp = length(readLines(paste(trainbed,".bim",sep="")))
n.snp = 80000

# model building from training data without covariates
sp1 = stmgplink(trainbed=trainbed,maxal=3000/n.snp)
head(sp1$PE)

# model building from training data to predict test data without covariates
sp2 = stmgplink(trainbed=trainbed,testbed=testbed,maxal=3000/n.snp)
head(sp2$PEte)
sp2$nonzero

# covariates files
# covariates (PLINK's covariate format, training data)
Z = as.matrix(read.table(paste(wd,"train.cov",sep="/"),h=TRUE)[,-(1:2)])
# covariates (PLINK's covariate format, test data)
Zte = as.matrix(read.table(paste(wd,"test.cov",sep="/"),h=TRUE)[,-(1:2)])

# model building from training data
sp3 = stmgplink(trainbed=trainbed,maxal=3000/n.snp,Z=Z)
head(sp3$PE)

# model building from training data and predicting test data
sp4 = stmgplink(trainbed=trainbed,testbed=testbed,maxal=3000/n.snp,Z=Z,Zte=Zte)
head(sp4$PEte)
sp4$nonzero
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
