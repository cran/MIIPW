% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bgnml.R
\name{bgnml}
\alias{bgnml}
\title{Bayesian analysis of generalised linear mixed model for Poisson outcome variable with two random effect}
\usage{
bgnml(m, n, n.chains, data)
}
\arguments{
\item{m}{starting column number}

\item{n}{ending column number}

\item{n.chains}{number of MCMC chains}

\item{data}{dataset with integer entries(including NA values), first row represents time points.
In this function we are using observations at four timepoints, function takes first four columns of the countdata
for the example given in the package}
}
\value{
posterior distribution result of parameters
}
\description{
provides bayesian analysis of generalised linear mixed model with log link function 
for categorical response using MCMC
}
\details{
The response variable \eqn{Y_{ij}} follows poisson distribution
,mean and variance given random effects \eqn{E(Y_{ij}|b_{1i},b_{2i})=Var(Y_{ij}|b_{1i},b_{2i})}
with link function \eqn{log(\mu_{ij})=\beta_1+b_{1i}+(\beta_2+b_{2i})(X_{ij}-\beta_3)}
where i is the ith subject and j is the timepoint.
}
\examples{
##
data(countdata)
bgnml(m=1,n=2,n.chains = 1,data=countdata)
##
}
\references{
Broemeling, Lyle D. Bayesian methods for repeated measures. CRC Press, 2015.
}
\author{
Atanu Bhattacharjee,Bhrigu Kumar Rajbongshi and Gajendra K Vishwakarma
}
