% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/brgan.R
\name{brgan}
\alias{brgan}
\title{Bayesian analysis of mean response model with autoregressive  covariance
matrix}
\usage{
brgan(m, n, n.chains, data)
}
\arguments{
\item{m}{starting column number}

\item{n}{ending column number}

\item{n.chains}{number of MCMC chains}

\item{data}{dataset with first column is age(there are missing values in age),and columns other than age are observation at
four different timepoints, where timepoints are the respective column names.}
}
\value{
posterior distribution result of the parameters
}
\description{
Bayesian analysis of mean response over time and age (quadratic trends)
using MCMC for AR1 covariance structure,where age follows normal distribution
}
\details{
The model for the response is \deqn{Y_{ij}=\beta_1+\beta_{2}t_{ij}
+\beta_{3}t_{ij}^2+\beta_{4} age_i+e_{ij}} \deqn{e_{ij}=\rho e_{ij-1}+u_{ij}} ,\eqn{u_{ij}\sim N(0,1/\tau);\rho} is the correlation
coefficient where i refers to ith individual and j is the timepoint. Missing values of covariate age is
imputed assuming age follows normal distribution.
}
\examples{
 ##
 data(agedata)
 brgan(m=1,n=3,n.chains=1,data=agedata)
 ##
}
\author{
Atanu Bhattacharjee,Bhrigu Kumar Rajbongshi and Gajendra K Vishwakarma
}
