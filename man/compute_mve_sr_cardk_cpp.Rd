% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{compute_mve_sr_cardk_cpp}
\alias{compute_mve_sr_cardk_cpp}
\title{Compute Mean-Variance Efficient Sharpe Ratio with Cardinality K}
\usage{
compute_mve_sr_cardk_cpp(mu, sigma, max_card, max_comb = 0L, do_checks = FALSE)
}
\arguments{
\item{mu}{Mean vector.}

\item{sigma}{Coveriance matrix.}

\item{max_card}{Maximum investment cardinality (from 1 up to the number of assets).}

\item{max_comb}{Maximum number of combinations to consider. If 0 (default),
all combinations are computed.}

\item{do_checks}{Logical flag indicating whether to perform input checks (default = \code{FALSE}).}
}
\value{
A list with \code{sr} (the optimal Sharpe ratio) and \code{selection} (the asset indices of the optimal selection).
}
\description{
This function takes as inputs the mean vector \code{mu}, the covariance matrix \code{sigma},
the maximum active portfolio cardinality \code{max_card},
and the maximum number of combinations per cardinality to evaluate \code{max_comb}.
With these inputs, it searches over all combinations of assets with cardinality from 1 up to \code{max_card}
and computes the Sharpe ratio defined as
\eqn{\mu^T \, \sigma^{-1}\, \mu}.
It returns the highest Sharpe ratio found along with the associated asset selection.
}
