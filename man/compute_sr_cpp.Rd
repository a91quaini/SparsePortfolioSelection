% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{compute_sr_cpp}
\alias{compute_sr_cpp}
\title{Compute Portfolio Sharpe Ratio}
\usage{
compute_sr_cpp(weights, mu, sigma, do_checks = FALSE)
}
\arguments{
\item{weights}{A numeric vector of portfolio weights.}

\item{mu}{A numeric vector of expected returns.}

\item{sigma}{A numeric covariance matrix.}

\item{do_checks}{Logical flag indicating whether to perform input checks (default = false).}
}
\value{
A double representing the Sharpe ratio.
}
\description{
Given a vector of portfolio weights, a vector of mean returns, and a covariance matrix,
this function computes the Sharpe ratio defined as
\deqn{\frac{w^T \mu}{\sqrt{w^T \Sigma w}}.}
}
