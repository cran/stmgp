\name{stmgeplink}
\alias{stmgeplink}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Smooth-threshold multivariate genetic prediction (incorporating gene-environment interactions) for genome-wide association or whole-genome sequencing data in PLINK format
}
\description{
Build prediction model from training data and predict test data phenotype through smooth-threshold multivariate genetic prediction (STMGP) method incorporating gene-environment (GxE) interactions,
in which GxE interaction effects are linearly added to the STMGP model with marginal effects. 
Data must be in PLINK binary format and marginal test p-values (i.e. test for each variant) are computed by PLINK software,
which enables rapid computation even for data having very large number of variants.
An optimal p-value cutoff is selected by Cp-type criterion.
Both quantitative and binary phenotypes are acceptable, in which data must be in PLINK fam file format or in a separate file (PLINK format, i.e. FID and IID are needed).
Environment variables need be in covariate file by specifying the column names. 
}
\usage{
stmgeplink(trainbed, Z ,Enames ,Zte = NULL, testbed = NULL, gamma = 1, taun = NULL,
   lambda = 1, plink = "plink --noweb", maf = 0.01, hwe = 1e-4, geno = 0.1,
   fout = "stp", trainfout = "train", testfout = "test", ll = 50, maxal = NULL,
   alc = NULL, tdir = NULL, Znames = NULL, trainphenofile = NULL, testphenofile = NULL, 
   phenoname = NULL, Assc = FALSE, AsscGE = FALSE, centerE = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trainbed}{
A training data file name in PLINK binary format or a vector of three file names of .bed, .bim, .fam;
Binary phenotype must be coded as 1 or 2 in PLINK fam file.
Missing values in .fam file are -9 as usual.
}
  \item{Z}{
A covariate file name for training data including environment variables (PLINK format, i.e. FID and IID are needed) or data matrix, missing values are "-9".
}
  \item{Enames}{
        Vector of (column) names of environment variables in the covariate file specified in \code{Z}.
}
  \item{Zte}{
A covariate file name for test data including environment variables (PLINK format, i.e. FID and IID are needed) or data matrix, missing values are "-9"; NULL (default) means unspecified.
}
  \item{testbed}{
A test data file name in PLINK binary format or a vector of three file names of .bed, .bim, .fam; NULL (default) means unspecified.
Missing values in .fam are -9 as usual.
}
  \item{gamma}{
gamma parameter; \code{gamma=1} is default as suggested in Ueki and Tamiya (2016).
}
  \item{taun}{
tau parameter divided by (sample size) (allowed to be a vector object; optimal parameter is chosen by Cp); NULL (default) specifies \code{tau=n/log(n)^0.5} as suggested in Ueki and Tamiya (2016).
}
  \item{lambda}{
lambda parameter (default=1).
}
  \item{plink}{
PLINK command, e.g. "plink2", "./plink --noweb", or "plink1.9 --memory 100000" (default is \code{plink --noweb}) where options can be added; PLINK must be installed.
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
Maximum p-value cutoff for search (default=NULL); If most variants are null \code{maxal}*(number of variables to be selected) gives approximate number of filtered variants, which is useful for rapid computation even for data with large number of variants.
}
  \item{alc}{
User-specified candidate p-value cutoffs for search; \code{ll} option is effective if \code{alc=NULL}.
}
  \item{tdir}{
Location of temporary files (default=\code{tempdir()}).
}
  \item{Znames}{
Name(s) of covariate used; NULL (default) means unspecified.
}
  \item{trainphenofile}{
A phenotype file name for training data (PLINK format, i.e. FID and IID are needed) with header columns (i.e. FID, IID, phenoname1, phenoname2, ...) missing values are "-9"; NULL (default) means unspecified.
}
  \item{testphenofile}{
A phenotype file name for test data (PLINK format, i.e. FID and IID are needed) with header columns (i.e. FID, IID, phenoname1, phenoname2, ...) missing values are "-9"; NULL (default) means unspecified.
}
  \item{phenoname}{
Phenotype name in \code{trainphenofile}; NULL (default) means unspecified but users should provide if \code{trainphenofile} is provided.
}
  \item{Assc}{
Whether the marginal association result is stored or not (default=FALSE).
}
  \item{AsscGE}{
Whether the gene-environment interaction result is stored or not (default=FALSE).
}
  \item{centerE}{
Whether environment variables are centered or not by subtracting their mean (default=TRUE).
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
  \item{lopt}{An optimal tuning parameter indexes for \code{al} and \code{tau} selected by Cp-type criterion, \code{CP}}
  \item{BA}{Estimated regression coefficient matrix whose size is of (1 + number of columns of \code{Z} + number of columns of \code{X}) x (length of \code{al})) x (length of \code{tau})); the first element, the second block and third block correspond to intercept, \code{Z} and \code{X}, respectively.}
  \item{Loss}{Loss (sum of squared residuals or -2*loglikelihood) whose size is of (length of \code{al}) x (length of \code{tau}).}
  \item{sig2hato}{An error variance estimate (=1 for binary traits) used in computing the variance term of Cp-type criterion.}
  \item{tau}{Candidate tau parameters for search.}
  \item{CP}{Cp-type criterion whose size is of (length of \code{al}) x (length of \code{tau}).}
  \item{PE}{PLINK .fam file for training data with additional column including the predicted phenotype from linear model.}
  \item{PEte}{PLINK .fam file for test data with additional column including the predicted phenotype from linear model estimated from training data.}
  \item{nonzero}{Variants with nonzero regression coefficients at the optimal parameter in PLINK file.}
  \item{DataTr}{Training dataset used (\code{y}, \code{X} and \code{Z})}
  \item{lapprox}{lapprox values for gene-environment interaction tests (discrepancy from 1 suggests model misspecification)}
  \item{ASSC}{Marginal association result from PLINK if \code{Assc} is TRUE}
  \item{ASSCGE}{Gene-environment interaction result (approx test) from PLINK if \code{AsscGE} is TRUE}
  \item{ASSCGEa}{Gene-environment interaction result (all tests) from PLINK if \code{AsscGE} is TRUE}
}
\references{

Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira M, Bender D, Maller J, Sklar P, de Bakker P, Daly MJ, Sham PC. (2007) PLINK: A tool set for whole-genome and population-based linkage analyses. Am J Hum Genet 81:559-75.

Chang CC, Chow CC, Tellier LCAM, Vattikuti S, Purcell SM, Lee JJ. (2015) Second-generation PLINK: rising to the challenge of larger and richer datasets. GigaScience 4.

Ueki M, Fujii M, Tamiya G. (2019) Quick assessment for systematic test statistic inflation/deflation due to null model misspecifications in genome-wide environment interaction studies. PLoS ONE 14: e0219825.

}
\examples{

\dontrun{
wd = system.file("extdata",package="stmgp")


# quantitative traits
# training data (plink format)
trainbed = paste(wd,"train",sep="/")
# test data (plink format)
testbed = paste(wd,"test",sep="/")
# number of SNPs
#n.snp = length(readLines(paste(trainbed,".bim",sep="")))
n.snp = 80000


#> head(read.table(paste0(trainbed,".fam")))
# training sample .fam file (quantitative phenotype in the 6th column)
#       V1      V2 V3 V4 V5    V6
#1 id1_100 id2_100  0  0  1 -1.23
#2 id1_101 id2_101  0  0  1  1.48
#3 id1_102 id2_102  0  0  1 -4.27
#4 id1_103 id2_103  0  0  1 -2.61
#5 id1_104 id2_104  0  0  1 -0.27
#6 id1_105 id2_105  0  0  1 -0.50


#> head(read.table(paste0(testbed,".fam")))
# test sample .fam file
# (quantitative phenotype in the 6th column but missing (i.e. "-9") allowed)
#     V1    V2 V3 V4 V5    V6
#1 id1_0 id2_0  0  0  1 -0.59
#2 id1_1 id2_1  0  0  1  1.11
#3 id1_2 id2_2  0  0  1 -2.45
#4 id1_3 id2_3  0  0  1  0.11
#5 id1_4 id2_4  0  0  1 -1.17
#6 id1_5 id2_5  0  0  1  2.08


# using covariates files

Zf = paste(wd,"train.cov",sep="/")
Ztef = paste(wd,"test.cov",sep="/")


#> head(read.table(Zf,header=TRUE))   # covariate for training sample
#      FID     IID        COV1 COV2 qphen bphen
#1 id1_340 id2_340  1.27203944    1 -2.47     2
#2 id1_430 id2_430 -0.44144482    1 -0.71     2
#3 id1_199 id2_199 -0.18200011    1 -3.42     2
#4 id1_473 id2_473  0.03965880    0  0.32     1
#5 id1_105 id2_105  0.20418279    0 -0.50     2
#6 id1_188 id2_188 -0.04838519    0  2.98     1


#> head(read.table(Ztef,header=TRUE))   # covariate for test sample
#     FID    IID       COV1 COV2
#1 id1_80 id2_80 -0.2057512    0
#2 id1_53 id2_53 -0.8627601    1
#3 id1_59 id2_59 -0.2973529    1
#4 id1_71 id2_71  1.4728727    1
#5 id1_92 id2_92  3.5614472    0
#6 id1_25 id2_25  0.5135032    1


# model building from training data
# (incorporating COV1xG interaction as well as two covariates, COV1 and COV2)
sge1p0 = stmgeplink(trainbed=trainbed,Z=Zf,Enames="COV1",
                 maxal=5000/n.snp,Znames=c("COV1","COV2"))
head(sge1p0$PE)

# model building from training data and predicting test data
# (incorporating COV1xG interaction as well as two covariates, COV1 and COV2)
sge1p = stmgeplink(trainbed=trainbed,testbed=testbed,Z=Zf,Zte=Ztef,Enames="COV1",
                 maxal=5000/n.snp,Znames=c("COV1","COV2"))
head(sge1p$PEte)
head(sge1p$nonzero)


# model building from training data
# (incorporating COV1xG and COV2xG interactions as well as two covariates, COV1 and COV2)
sge12p0 = stmgeplink(trainbed=trainbed,Z=Zf,Enames=c("COV1","COV2"),
                 maxal=5000/n.snp,Znames=c("COV1","COV2"))
head(sge12p0$PE)

# model building from training data and predicting test data
# (incorporating COV1xG and COV2xG interactions as well as two covariates, COV1 and COV2)
sge12p = stmgeplink(trainbed=trainbed,testbed=testbed,Z=Zf,Zte=Ztef,Enames=c("COV1","COV2"),
                 maxal=5000/n.snp,Znames=c("COV1","COV2"))
head(sge12p$PEte)
head(sge12p$nonzero)





#### binary traits ####
# training data (plink format)
trainbed = paste(wd,"train",sep="/")
# test data (plink format)
testbed = paste(wd,"test",sep="/")
# number of SNPs
#n.snp = length(readLines(paste(trainbed,".bim",sep="")))
n.snp = 80000


#> head(read.table(paste0(trainbed,"b.fam"))) 
# training sample .fam file (binary phenotype (1 or 2) in the 6th column)
#       V1      V2 V3 V4 V5 V6
#1 id1_100 id2_100  0  0  1  2
#2 id1_101 id2_101  0  0  1  1
#3 id1_102 id2_102  0  0  1  2
#4 id1_103 id2_103  0  0  1  2
#5 id1_104 id2_104  0  0  1  2
#6 id1_105 id2_105  0  0  1  2


#> head(read.table(paste0(testbed,"b.fam")))  
# test sample .fam file (binary phenotype (1 or 2) in the 6th column
# but missing (i.e. "-9") allowed)
#     V1    V2 V3 V4 V5 V6
#1 id1_0 id2_0  0  0  1  2
#2 id1_1 id2_1  0  0  1  1
#3 id1_2 id2_2  0  0  1  2
#4 id1_3 id2_3  0  0  1  1
#5 id1_4 id2_4  0  0  1  2
#6 id1_5 id2_5  0  0  1  1


# using covariates files

# model building from training data
# (incorporating COV1xG interaction as well as two covariates, COV1 and COV2)
sge1p0b = stmgeplink(trainbed=paste0(trainbed,c(".bed",".bim","b.fam")),
       Z=paste(wd,"train.cov",sep="/"),Enames="COV1",
       maxal=5000/n.snp,Znames=c("COV1","COV2"))
head(sge1p0b$PE)

# model building from training data and predicting test data
# (incorporating COV1xG interaction as well as two covariates, COV1 and COV2)
sge1pb = stmgeplink(trainbed=paste0(trainbed,c(".bed",".bim","b.fam")),
       testbed=paste0(testbed,c(".bed",".bim","b.fam")),
       Z=paste(wd,"train.cov",sep="/"),Zte=paste(wd,"test.cov",sep="/"),Enames="COV1",
       maxal=5000/n.snp,Znames=c("COV1","COV2"))
head(sge1pb$PEte)
head(sge1pb$nonzero)


# model building from training data
# (incorporating COV1xG and COV2xG interactions as well as two covariates, COV1 and COV2)
sge12p0b = stmgeplink(trainbed=paste0(trainbed,c(".bed",".bim","b.fam")),
       Z=paste(wd,"train.cov",sep="/"),Enames=c("COV1","COV2"),
       maxal=5000/n.snp,Znames=c("COV1","COV2"))
head(sge12p0b$PE)

# model building from training data and predicting test data
# (incorporating COV1xG and COV2xG interactions as well as two covariates, COV1 and COV2)
sge12pb = stmgeplink(trainbed=paste0(trainbed,c(".bed",".bim","b.fam")),
       testbed=paste0(testbed,c(".bed",".bim","b.fam")),
       Z=paste(wd,"train.cov",sep="/"),Zte=paste(wd,"test.cov",sep="/"),
       Enames=c("COV1","COV2"),maxal=5000/n.snp,Znames=c("COV1","COV2"))
head(sge12pb$PEte)
head(sge12pb$nonzero)






}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}% use one of  RShowDoc("KEYWORDS")
\keyword{regression}% __ONLY ONE__ keyword per line
