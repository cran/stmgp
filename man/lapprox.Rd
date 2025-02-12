\name{lapprox}
\alias{lapprox}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute lapprox
}
\description{
Compute the lapprox statistic developed in Ueki et al. (2019) to assess inflation/deflation in gene-environment interaction test before genome-wide scan.
No genotype data is used in computing the lapprox.
The statistic far from 1 implies inflation/deflation in genome-wide gene-environment interaction test (Kraft joint test), suggesting that the null model can be misspecified, where the null model assumes unrelated samples without population stratification.
Quantitative and binary phenotypes are modeled via linear and logistic regression, respectively.
}
\usage{
lapprox(Z, X, y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{
Covariates to be adjusted for in regression model. Intercept term (i.e. 1 vector) should be included.
}
  \item{X}{
Covariates to be multipled by genetic variants (G), which should be a subset of X. Environment variables can be included in X. The terms to be tested are X[,1]*G, X[,2]*G, ..., X[,ncol(X)]*G. See example below.
}
  \item{y}{
A response variable, either quantitative or binary (coded 0 or 1).
}
}
\details{
See Ueki et al. (2019).
}
\references{
Kraft P, Yen YC, Stram DO, Morrison J, Gauderman WJ (2007) Exploiting gene-environment interaction to detect genetic associations. Human Heredity 63:111-9. 

Ueki M, Fujii M, Tamiya G (2019) Quick assessment for systematic test statistic inflation/deflation due to null model misspecifications in genome-wide environment interaction studies. PLoS ONE 14: e0219825.
}
\value{
Value of lapprox. Deviation from 1 implies inflation/deflation in genome-wide gene-environment interaction test.
}
\examples{

\dontrun{


set.seed(3222)

# lapprox both for quantitative and binary phenotypes
# input:
# Z: covariate matrix including the intercept
# X: matrix for joint GxE interaction effects (marginal effect must be represented by the intercept)
# y: phenotype vector (0/1 code needed for binary phenotype)


n = 400  # sample size
q = 2  # number of covariates
Z = matrix(rnorm(n*q),n,q)  # covariates, the first component is used as the environment variable
mu = -1.5 + Z[,1]
y = mu + rnorm(n)  # quantitative phenotype
Y = rbinom(n,size=1,prob=1/(1+exp(-mu*5)))  # binary phenotype


#Example usage for quantitative trait

# lapprox for joint GxE interaction test (correctly specified)
# y ~ b0*1 + b1*Z[,1] + b2*Z[,2] + b3*G*1 + b4*G*Z[,1]
# H0: b3 = b4 = 0 is tested
lapprox(Z=cbind(1,Z),X=cbind(1,Z[,1]),y=y)
# [1] 1.012104

# lapprox for joint GxE interaction test (misspecified)
# y ~ b0*1 + b1*Z[,1]*Z[,1] + b2*Z[,2]*Z[,2] + b3*G*1 + b4*G*Z[,1]*Z[,1]
# H0: b3 = b4 = 0 is tested
lapprox(Z=cbind(1,Z*Z),X=cbind(1,Z[,1]*Z[,1]),y=y)
# [1] 1.969929

# lapprox for marginal association test
# y ~ b0*1 + b1*Z[,1] + b2*Z[,2] + b3*G*1
# H0: b3 = 0 is tested is tested
lapprox(Z=cbind(1,Z),X=matrix(1,nrow(Z),1),y=y)
# [1] 1.010101


#Example usage for binary trait

# lapprox for joint GxE interaction test (correctly specified)
# y ~ b0*1 + b1*Z[,1] + b2*Z[,2] + b3*G*1 + b4*G*Z[,1]
# H0: b3 = b4 = 0 is testedd
lapprox(Z=cbind(1,Z),X=cbind(1,Z[,1]),y=Y)
# [1] 0.9388417

# lapprox for joint GxE interaction test (misspecified)
# y ~ b0*1 + b1*Z[,1]*Z[,1] + b2*Z[,2]*Z[,2] + b3*G*1 + b4*G*Z[,1]*Z[,1]
# H0: b3 = b4 = 0 is tested
lapprox(Z=cbind(1,Z*Z),X=cbind(1,Z[,1]*Z[,1]),y=Y)
# [1] 1.150531

# lapprox for marginal association test
# y ~ b0*1 + b1*Z[,1] + b2*Z[,2] + b3*G*1
# H0: b3 = 0 is tested
lapprox(Z=cbind(1,Z),X=matrix(1,nrow(Z),1),y=Y)
# [1] 1.078894


}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}% use one of  RShowDoc("KEYWORDS")
\keyword{regression}% __ONLY ONE__ keyword per line
