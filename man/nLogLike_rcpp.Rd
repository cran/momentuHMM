% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{nLogLike_rcpp}
\alias{nLogLike_rcpp}
\title{Negative log-likelihood}
\usage{
nLogLike_rcpp(
  nbStates,
  covs,
  data,
  dataNames,
  dist,
  Par,
  aInd,
  zeroInflation,
  oneInflation,
  stationary,
  knownStates,
  betaRef,
  mixtures
)
}
\arguments{
\item{nbStates}{Number of states,}

\item{covs}{Covariates,}

\item{data}{A \code{\link{momentuHMMData}} object of the observations,}

\item{dataNames}{Character vector containing the names of the data streams,}

\item{dist}{Named list indicating the probability distributions of the data streams.}

\item{Par}{Named list containing the state-dependent parameters of the data streams, matrix of regression coefficients 
for the transition probabilities ('beta'), and initial distribution ('delta').}

\item{aInd}{Vector of indices of the rows at which the data switches to another animal}

\item{zeroInflation}{Named list of logicals indicating whether the probability distributions of the data streams are zero-inflated.}

\item{oneInflation}{Named list of logicals indicating whether the probability distributions of the data streams are one-inflated.}

\item{stationary}{\code{false} if there are time-varying covariates in \code{formula} or any covariates in \code{formulaDelta}. If \code{true}, the initial distribution is considered
equal to the stationary distribution. Default: \code{false}.}

\item{knownStates}{Vector of values of the state process which are known prior to fitting the
model (if any). Default: NULL (states are not known). This should be a vector with length the number
of rows of 'data'; each element should either be an integer (the value of the known states) or NA if
the state is not known.}

\item{betaRef}{Indices of reference elements for t.p.m. multinomial logit link.}

\item{mixtures}{Number of mixtures for the state transition probabilities}
}
\value{
Negative log-likelihood
}
\description{
Computation of the negative log-likelihood (forward algorithm - written in C++)
}
