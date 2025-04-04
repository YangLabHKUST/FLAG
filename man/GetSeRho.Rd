% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/InferWald.R
\name{GetSeRho}
\alias{GetSeRho}
\title{Get the standard error of rho.}
\usage{
GetSeRho(Gamma_e, Gamma_e_cov)
}
\arguments{
\item{Gamma_e}{Matrix, with size 2*2.}

\item{Gamma_e_cov}{Matrix, with size 2*2.}
}
\value{
Numeric.
}
\description{
Get the standard error of rho.
}
\examples{
\donttest{
GetSeRho(Gamma_e,Gamma_e_cov)
}
}
