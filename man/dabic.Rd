\name{dabic}
\alias{dabic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute daBIC
}
\description{
Compute the deflation-adjusted BIC (daBIC) developed in Ueki (2025) for model selection of K-means clustering.
}
\usage{
dabic(x, cluster)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Variables that were used for K-means clustering. No standardization is needed.
}
  \item{cluster}{
Integer variable indicating cluster membership for each sample obtained from K-means clustering.
}
}
\details{
See Ueki (2025).
}
\references{
Ueki M (2025) A deflation-adjusted Bayesian information criterion for selecting the number of clusters in K-means clustering. Comput Stat Data Anal 209:108170.
}
\value{
\item{lik}{Value of log-likelihood.}
\item{bic}{Value of BIC.}
\item{bica}{Value of daBIC.}
}
\examples{

\dontrun{



# simulated data
set.seed(222)

nn = rep(200,3)
n = sum(nn)

x = rbind( cbind(rnorm(nn[1],-5,1),rnorm(nn[1],-5,1)), 
           cbind(rnorm(nn[2],0,1),rnorm(nn[2],0,1)), 
           cbind(rnorm(nn[3],5,1),rnorm(nn[3],5,1)) )


# maximum number of clusters to search
Kmax = 10

bic = bica = numeric(Kmax)
clusterbic = clusterbica = rep(1,n)

da1 = dabic(x,rep(1,n))  # daBIC

bic[1] = da1$bic  # K=1
bica[1] = da1$bica  # K=1


for(kk in 2:Kmax){
	km = kmeans(scale(x),centers=kk,nstart=25)  # K-means clustering
	dakk = dabic(x,km$cluster)  # daBIC
	bic[kk] = dakk$bic  # K
	bica[kk] = dakk$bica  # K
	if(min(bic[1:(kk-1)])>bic[kk]) clusterbic = km$cluster  # update cluster membership
	if(min(bica[1:(kk-1)])>bica[kk]) clusterbica = km$cluster  # update cluster membership
}


# Optimal K from BIC
which.min(bic)

# Optimal K from daBIC
which.min(bica)

# plot x colored with optimal cluster membership obtained from daBIC
plot(x,col=clusterbica)



}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}% use one of  RShowDoc("KEYWORDS")
\keyword{regression}% __ONLY ONE__ keyword per line
