% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CIreal.R
\name{CIreal}
\alias{CIreal}
\alias{CIreal.default}
\alias{CIreal.hierarchical}
\title{Confidence intervals for the natural (i.e., real) parameters}
\usage{
CIreal(m, alpha = 0.95, covs = NULL, parms = NULL)

\method{CIreal}{default}(m, alpha = 0.95, covs = NULL, parms = NULL)

\method{CIreal}{hierarchical}(m, alpha = 0.95, covs = NULL, parms = NULL)
}
\arguments{
\item{m}{A \code{momentuHMM}, \code{momentuHierHMM}, \code{miHMM}, or \code{miSum} object}

\item{alpha}{Significance level of the confidence intervals. Default: 0.95 (i.e. 95\% CIs).}

\item{covs}{Data frame consisting of a single row indicating the covariate values to be used in the calculations. By default, no covariates are specified.}

\item{parms}{Optional character vector indicating which groups of real parameters to calculate confidence intervals for (e.g., 'step', 'angle', 'gamma', 'delta', etc.). Default: NULL, in which case confidence intervals are calculated for all groups of parameters in the model.}
}
\value{
A list of the following objects:
\item{...}{List(s) of estimates ('est'), standard errors ('se'), and confidence intervals ('lower', 'upper') for the natural parameters of the data streams}
\item{gamma}{List of estimates ('est'), standard errors ('se'), and confidence intervals ('lower', 'upper') for the transition probabilities}
\item{delta}{List of estimates ('est'), standard errors ('se'), and confidence intervals ('lower', 'upper') for the initial state probabilities}
\item{hierGamma}{A hierarchical data structure \code{\link[data.tree]{Node}} including a list of estimates ('est'), standard errors ('se'), and confidence intervals ('lower', 'upper') for the transition probabilities for each level of the hierarchy (only applies if \code{m} is a hierarchical model object)}
\item{hierDelta}{A hierarchical data structure \code{\link[data.tree]{Node}} including a list of estimates ('est'), standard errors ('se'), and confidence intervals ('lower', 'upper') for the initial state probabilities for each level of the hierarchy (only applies if \code{m} is a hierarchical model object)}
}
\description{
Computes the standard errors and confidence intervals on the real (i.e., natural) scale of the data stream probability distribution parameters,
as well as for the transition probabilities parameters. If covariates are included in the probability distributions or TPM formula, the mean values
of non-factor covariates are used for calculating the natural parameters. For any covariate(s) of class 'factor', then the value(s) from the first observation 
in the data are used.
}
\details{
For any covariates that are not specified using \code{covs}, the means of the covariate(s) are used 
(unless the covariate is a factor, in which case the first factor in the data is used).
}
\examples{
# m is a momentuHMM object (as returned by fitHMM), automatically loaded with the package
m <- example$m

ci1<-CIreal(m)

# specify 'covs'
ci2<-CIreal(m,covs=data.frame(cov1=mean(m$data$cov1),cov2=mean(m$data$cov2)))

all.equal(ci1,ci2)

}
