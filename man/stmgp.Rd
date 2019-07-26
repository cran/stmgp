\name{stmgp}
\alias{stmgp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Smooth-threshold multivariate genetic prediction
}
\description{
Smooth-threshold multivariate genetic prediction (STMGP) method,
which is based on the smooth-threshold estimating equations (Ueki 2009).
Variable selection is performed based on marginal association test p-
values (i.e. test of nonzero slope parameter in univariate regression for each predictor variable) with an optimal p-value cutoff selected by a Cp-type criterion.
Quantitative and binary phenotypes are modeled via linear and logistic regression, respectively.
}
\usage{
stmgp(y, X, Z = NULL, tau, qb, maxal, gamma = 1, ll = 50,
  lambda = 1, alc = NULL, pSum = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
A response variable, either quantitative or binary (coded 0 or 1); Response type is specified by \code{qb}.
}
  \item{X}{
Predictor variables subjected to variable selection.
}
  \item{Z}{
Covariates; \code{Z=NULL} means unspecified.
}
  \item{tau}{
tau parameter (allowed to be a vector object); NULL (default) specifies \code{tau=n/log(n)^0.5} as suggested in Ueki and Tamiya (2016).
}
  \item{qb}{
Type of response variable, \code{qb="q"} and \code{"b"} specify quantitative and binary traits, respectively.
}
  \item{maxal}{
Maximum p-value cutoff for search.
}
  \item{gamma}{
gamma parameter; \code{gamma=1} is default as suggested in Ueki and Tamiya (2016).
}
  \item{ll}{
Number of candidate p-value cutoffs for search (default=50) as determined by
  \code{10^seq( log10(maxal),log10(5e-8), length=ll)}.
}
  \item{lambda}{
lambda parameter (default=1).
}
  \item{alc}{
User-specified candidate p-value cutoffs for search; \code{ll} option is effective if \code{alc=NULL}.
}
  \item{pSum}{
User-specified p-values matrix from other studies that are independent of the study data (optional, default=NULL), a matrix object having rows with the same size of \code{X} and columns for each study (multiple studies are capable). Missing p-values must be coded as NA. Summary p-values are combined with p-values in the study data by the Fisher's method.
}
}
\details{
See Ueki and Tamiya (2016).
}
\references{
Ye J. (1988) On measuring and correcting the effects of data mining and
model selection. J Am Stat Assoc 93:120-31.

Ueki M. (2009) A note on automatic variable selection using smooth-threshold estimating equations. Biometrika 96:1005-11.
}
\value{
  \item{Muhat}{Estimated phenotypic values from linear model evaluated at each candidate tuning parameters (\code{al} and \code{tau}) whose size is of (sample size) x (length of \code{al}) x (length of \code{tau}).}
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
}
\examples{

\dontrun{


set.seed(22200)

wd = system.file("extdata",package="stmgp")

D = read.table(unzip(paste(wd,"snps.raw.zip",sep="/"),exdir=tempdir()),header=TRUE)

X = D[,-(1:6)]
X = (X==1) + 2*(X==2)
p = ncol(X)
n = nrow(X)
ll = 30
p0 = 50; b0 = log(rep(1.2,p0))
iA0 = sample(1:p,p0)
Z = as.matrix(cbind(rnorm(n),runif(n)))  # covariates
eta = crossprod(t(X[,iA0]),b0) - 4 + crossprod(t(Z),c(0.5,0.5))


# quantitative trait
mu = eta
sig = 1.4
y = mu + rnorm(n)*sig
STq = stmgp(y,X,Z,tau=n*c(1),qb="q",maxal=0.1,gamma=1,ll=ll)
boptq = STq$BA[,STq$lopt[1],STq$lopt[2]]  # regression coefficient in selected model
nonzeroXq = which( boptq[(1+ncol(Z))+(1:p)]!=0 )  # nonzero regression coefficient
# check consistency
cor( STq$Muhat[,STq$lopt[1],STq$lopt[2]], crossprod(t(cbind(1,Z,X)),boptq) )
cor( STq$Muhat[,STq$lopt[1],STq$lopt[2]], eta)  # correlation with true function
# proportion of correctly identified true nonzero regression coefficients
length(intersect(which(boptq[-(1:(ncol(Z)+1))]!=0),iA0))/length(iA0)


# binary trait
mu = 1/(1+exp(-eta))
Y = rbinom(n,size=1,prob=mu)
STb = stmgp(Y,X,Z,tau=n*c(1),qb="b",maxal=0.1,gamma=1,ll=ll)
boptb = STb$BA[,STb$lopt[1],STb$lopt[2]]  # regression coefficient in selected model
nonzeroXb = which( boptb[(1+ncol(Z))+(1:p)]!=0 )  # nonzero regression coefficient
# check consistency
cor( STb$Muhat[,STb$lopt[1],STb$lopt[2]], crossprod(t(cbind(1,Z,X)),boptb) )
Prob = 1/(1+exp(-STb$Muhat[,STb$lopt[1],STb$lopt[2]]))  # Pr(Y=1) (logistic regression)
cor( STb$Muhat[,STb$lopt[1],STb$lopt[2]], eta)  # correlation with true function
# proportion of correctly identified true nonzero regression coefficients 
length(intersect(which(boptb[-(1:(ncol(Z)+1))]!=0),iA0))/length(iA0)



# simulated summary p-values
pSum = cbind(runif(ncol(X)),runif(ncol(X)));
pSum[iA0,1] = pchisq(rnorm(length(iA0),5,1)^2,df=1,low=F); # study 1 summary p-values
pSum[iA0,2] = pchisq(rnorm(length(iA0),6,1)^2,df=1,low=F); # study 2 summary p-values
pSum[sample(1:length(pSum),20)] = NA
head(pSum)


# quantitative trait using summary p-values
STqs = stmgp(y,X,Z,tau=n*c(1),qb="q",maxal=0.1,gamma=1,ll=ll,pSum=pSum)
boptqs = STqs$BA[,STqs$lopt[1],STqs$lopt[2]]  # regression coefficient in selected model
nonzeroXqs = which( boptqs[(1+ncol(Z))+(1:p)]!=0 )  # nonzero regression coefficient
# check consistency
cor( STqs$Muhat[,STqs$lopt[1],STqs$lopt[2]], crossprod(t(cbind(1,Z,X)),boptqs) )
cor( STqs$Muhat[,STqs$lopt[1],STqs$lopt[2]], eta)  # correlation with true function
# proportion of correctly identified true nonzero regression coefficients 
length(intersect(which(boptqs[-(1:(ncol(Z)+1))]!=0),iA0))/length(iA0)



# binary trait using summary p-values
STbs = stmgp(Y,X,Z,tau=n*c(1),qb="b",maxal=0.1,gamma=1,ll=ll,pSum=pSum)
boptbs = STbs$BA[,STbs$lopt[1],STbs$lopt[2]]  # regression coefficient in selected model
nonzeroXbs = which( boptbs[(1+ncol(Z))+(1:p)]!=0 )  # nonzero regression coefficient
# check consistency
cor( STbs$Muhat[,STbs$lopt[1],STbs$lopt[2]], crossprod(t(cbind(1,Z,X)),boptbs) )
Prob = 1/(1+exp(-STbs$Muhat[,STbs$lopt[1],STbs$lopt[2]]))  # Pr(Y=1) (logistic regression)
cor( STbs$Muhat[,STbs$lopt[1],STbs$lopt[2]], eta)  # correlation with true function
# proportion of correctly identified true nonzero regression coefficients 
length(intersect(which(boptbs[-(1:(ncol(Z)+1))]!=0),iA0))/length(iA0)





}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Regression}% use one of  RShowDoc("KEYWORDS")
