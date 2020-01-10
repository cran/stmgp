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
Both quantitative and binary phenotypes are acceptable, in which data must be in PLINK fam file format or in a separate file (PLINK format, i.e. FID and IID are needed).
}
\usage{
stmgplink(trainbed, testbed = NULL, gamma = 1, taun = NULL,
  lambda = 1, Z = NULL, Zte = NULL, plink = "plink --noweb",
  maf = 0.01, hwe = 1e-04, geno = 0.1, fout = "stp",
  trainfout = "train", testfout = "test", ll = 50, maxal = NULL, alc = NULL,
  tdir = NULL, Znames=NULL, trainphenofile=NULL, phenoname=NULL, pSum = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trainbed}{
A training data file name in PLINK binary format or a vector of three file names of .bed, .bim, .fam;
Binary phenotype must be coded as 1 or 2 in PLINK fam file.
Missing values in .fam file are -9 as usual.
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
  \item{Z}{
A covariate file name for training data (PLINK format, i.e. FID and IID are needed) or data matrix, missing values are "-9"; NULL (default) means unspecified.
}
  \item{Zte}{
A covariate file name for test data (PLINK format, i.e. FID and IID are needed) or data matrix, missing values are "-9"; NULL (default) means unspecified.
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
Maximum p-value cutoff for search (default=NULL); If most variants are null \code{maxal}*(number of variants) gives approximate number of filtered variants, which is useful for rapid computation even for data with large number of variants.
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
  \item{phenoname}{
Phenotype name in \code{trainphenofile}; NULL (default) means unspecified but users should provide if \code{trainphenofile} is provided.
}
  \item{pSum}{
User-specified p-values matrix from other studies that are independent of the study data (optional, default=NULL), a matrix object with rows for each variant and columns for each study (multiple studies are capable) where rownames must be specified from SNP IDs that exist in PLINK .bim file. Missing p-values must be coded as NA. Summary p-values are combined with p-values in the study data by the Fisher's method.
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



pSum = read.table(paste0(wd,"/pSum.txt"));
rownames(pSum) = pSum[,2]; pSum = pSum[,7:8,drop=F]

# summary p-values format
# ("rownames(pSum)" for SNP ID should be provided; "NA" means unavailable)
#> head(pSum,15)  # summary p-values from two studies
#                   V7        V8
#rs12255619         NA        NA
#rs11252546         NA        NA
#rs7909677          NA        NA
#rs10904494         NA        NA
#rs11591988         NA        NA
#rs4508132          NA        NA
#rs10904561         NA        NA
#rs7917054          NA        NA
#rs7906287          NA        NA
#rs12775579         NA        NA
#rs4495823  0.08731436 0.2150108
#rs11253478 0.24030258 0.8726241
#rs9419557  0.49243856 0.3823173
#rs9286070  0.31277506 0.8521706
#rs9419560          NA 0.7604134



# model building from training data without covariates
sp1 = stmgplink(trainbed=trainbed,maxal=5000/n.snp)
head(sp1$PE)

# model building from training data to predict test data without covariates
sp2 = stmgplink(trainbed=trainbed,testbed=testbed,maxal=5000/n.snp)
head(sp2$PEte)
head(sp2$nonzero)

# model building from training data to predict test data without covariates
# (using 1 pSum)
sp2p = stmgplink(trainbed=trainbed,testbed=testbed,maxal=5000/n.snp,
	pSum=pSum[,1,drop=F])
head(sp2p$PEte)
head(sp2p$nonzero)

# model building from training data to predict test data without covariates
# (using 2 pSum)
sp2pp = stmgplink(trainbed=trainbed,testbed=testbed,maxal=5000/n.snp,pSum=pSum)
head(sp2pp$PEte)
head(sp2pp$nonzero)


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


#> head(read.table(Zfte,header=TRUE))   # covariate for test sample
#     FID    IID       COV1 COV2
#1 id1_80 id2_80 -0.2057512    0
#2 id1_53 id2_53 -0.8627601    1
#3 id1_59 id2_59 -0.2973529    1
#4 id1_71 id2_71  1.4728727    1
#5 id1_92 id2_92  3.5614472    0
#6 id1_25 id2_25  0.5135032    1


# model building from training data
sp3 = stmgplink(trainbed=trainbed,maxal=5000/n.snp,Z=Zf,Znames=c("COV1","COV2"))
head(sp3$PE)

# model building from training data and predicting test data
sp4 = stmgplink(trainbed=trainbed,testbed=testbed,maxal=5000/n.snp,Z=Zf,Zte=Ztef,
                 Znames=c("COV1","COV2"))
head(sp4$PEte)
head(sp4$nonzero)

# model building from training data and predicting test data (using 1 pSum)
sp4p = stmgplink(trainbed=trainbed,testbed=testbed,maxal=5000/n.snp,
                 Z=Zf,Zte=Ztef,Znames=c("COV1","COV2"),pSum=pSum[,1,drop=F])
head(sp4p$PEte)
head(sp4p$nonzero)

# model building from training data and predicting test data (using 2 pSum)
sp4pp = stmgplink(trainbed=trainbed,testbed=testbed,maxal=5000/n.snp,
                  Z=Zf,Zte=Ztef,Znames=c("COV1","COV2"),pSum=pSum)
head(sp4pp$PEte)
head(sp4pp$nonzero)



# no summary p-values
cor(sp2$PE[,6:7],use="pair")[1,2]  # training (no covariate)
cor(sp2$PEte[,6:7],use="pair")[1,2]  # test (no covariate)

cor(sp4$PE[,6:7],use="pair")[1,2]  # training (covariates)
cor(sp4$PEte[,6:7],use="pair")[1,2]  # test (covariates)


# 1 summary p-values
cor(sp2p$PE[,6:7],use="pair")[1,2]  # training (no covariate)
cor(sp2p$PEte[,6:7],use="pair")[1,2]  # test (no covariate)

cor(sp4p$PE[,6:7],use="pair")[1,2]  # training (covariates)
cor(sp4p$PEte[,6:7],use="pair")[1,2]  # test (covariates)


# 2 summary p-values
cor(sp2pp$PE[,6:7],use="pair")[1,2]  # training (no covariate)
cor(sp2pp$PEte[,6:7],use="pair")[1,2]  # test (no covariate)

cor(sp4pp$PE[,6:7],use="pair")[1,2]  # training (covariates)
cor(sp4pp$PEte[,6:7],use="pair")[1,2]  # test (covariates)




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



# model building from training data without covariates
sp1b = stmgplink(trainbed=paste0(trainbed,c(".bed",".bim","b.fam")),maxal=5000/n.snp)
head(sp1b$PE)

# model building from training data to predict test data without covariates
sp2b = stmgplink(trainbed=paste0(trainbed,c(".bed",".bim","b.fam")),
       testbed=paste0(testbed,c(".bed",".bim","b.fam")),maxal=5000/n.snp)
head(sp2b$PEte)
head(sp2b$nonzero)

# model building from training data to predict test data without covariates
# (using 1 pSum)
sp2bp = stmgplink(trainbed=paste0(trainbed,c(".bed",".bim","b.fam")),
        testbed=paste0(testbed,c(".bed",".bim","b.fam")),maxal=5000/n.snp,
	pSum=pSum[,1,drop=F])
head(sp2bp$PEte)
head(sp2bp$nonzero)

# model building from training data to predict test data without covariates
# (using 2 pSum)
sp2bpp = stmgplink(trainbed=paste0(trainbed,c(".bed",".bim","b.fam")),
         testbed=paste0(testbed,c(".bed",".bim","b.fam")),maxal=5000/n.snp,pSum=pSum)
head(sp2bpp$PEte)
head(sp2bpp$nonzero)


# using covariates files

# model building from training data
sp3b = stmgplink(trainbed=paste0(trainbed,c(".bed",".bim","b.fam")),
       maxal=5000/n.snp,Z=paste(wd,"train.cov",sep="/"),Znames=c("COV1","COV2"))
head(sp3b$PE)

# model building from training data and predicting test data
sp4b = stmgplink(trainbed=paste0(trainbed,c(".bed",".bim","b.fam")),
       testbed=paste0(testbed,c(".bed",".bim","b.fam")),maxal=5000/n.snp,Z=Zf,Zte=Ztef,
       Znames=c("COV1","COV2"))
head(sp4b$PEte)
head(sp4b$nonzero)

# model building from training data and predicting test data (using 1 pSum)
sp4bp = stmgplink(trainbed=paste0(trainbed,c(".bed",".bim","b.fam")),
        testbed=paste0(testbed,c(".bed",".bim","b.fam")),maxal=5000/n.snp,
        Z=Zf,Zte=Ztef,Znames=c("COV1","COV2"),pSum=pSum[,1,drop=F])
head(sp4bp$PEte)
head(sp4bp$nonzero)

# model building from training data and predicting test data (using 2 pSum)
sp4bpp = stmgplink(trainbed=paste0(trainbed,c(".bed",".bim","b.fam")),
         testbed=paste0(testbed,c(".bed",".bim","b.fam")),maxal=5000/n.snp,
         Z=Zf,Zte=Ztef,Znames=c("COV1","COV2"),pSum=pSum)
head(sp4bpp$PEte)
head(sp4bpp$nonzero)



# no summary p-values
cor(sp2b$PE[,6:7],use="pair")[1,2]  # training (no covariate)
cor(sp2b$PEte[,6:7],use="pair")[1,2]  # test (no covariate)

cor(sp4b$PE[,6:7],use="pair")[1,2]  # training (covariates)
cor(sp4b$PEte[,6:7],use="pair")[1,2]  # test (covariates)


# 1 summary p-values
cor(sp2bp$PE[,6:7],use="pair")[1,2]  # trainig (no covariate)
cor(sp2bp$PEte[,6:7],use="pair")[1,2]  # test (no covariate)

cor(sp4bp$PE[,6:7],use="pair")[1,2]  # training (covariates)
cor(sp4bp$PEte[,6:7],use="pair")[1,2]  # test (covariates)


# 2 summary p-values
cor(sp2bpp$PE[,6:7],use="pair")[1,2]  # training (no covariate)
cor(sp2bpp$PEte[,6:7],use="pair")[1,2]  # test (no covariate)

cor(sp4bpp$PE[,6:7],use="pair")[1,2]  # training (covariates)
cor(sp4bpp$PEte[,6:7],use="pair")[1,2]  # test (covariates)




}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}% use one of  RShowDoc("KEYWORDS")
\keyword{regression}% __ONLY ONE__ keyword per line
