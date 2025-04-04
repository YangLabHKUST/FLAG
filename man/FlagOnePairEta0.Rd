% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FlagOnePairEta0.R
\name{FlagOnePairEta0}
\alias{FlagOnePairEta0}
\title{FLAG for one pair of random variables fixing eta as zero, using likelihood-ratio test.}
\usage{
FlagOnePairEta0(
  Y,
  X,
  Gamma_beta,
  Gamma_e,
  eps = 1e-07,
  max.iter = 5000,
  crit.loglik = 1e-04
)
}
\arguments{
\item{Y}{Matrix, with size n*2.}

\item{X}{Matrix, with size n*(p-2).}

\item{Gamma_beta}{Matrix, with size 2*2.}

\item{Gamma_e}{Matrix, with size 2*2.}

\item{eps}{Numeric, a small term to avoid numerical problems, default to be 1e-4.}

\item{max.iter}{Integer, the maximum number of iterations, default to be 5000.}

\item{crit.loglik}{Numeric, the criteria of the change ratio of log likelihood to stop.}
}
\value{
List,
the list of log likelihood during iterations,
Numeric, the log likelihood of the last iterations,
the estimated Gamma_beta matrix with size 2\emph{2, in the random effects model,
the estimated Gamma_epsilon matrix with size 2}2, in the random effects model,
the estimated 2*2 submatrix of the precision matrix.
}
\description{
FLAG for one pair of random variables fixing eta as zero, using likelihood-ratio test.
}
\examples{
\donttest{
FlagOnePairEta0(Y, X, Gamma_beta, Gamma_e)
}
}
