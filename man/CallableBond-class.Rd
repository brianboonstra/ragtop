% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/instruments.R
\docType{class}
\name{CallableBond-class}
\alias{CallableBond-class}
\alias{CallableBond}
\title{Callable (and putable) corporate or government bond.}
\description{
When a bond is \emph{callable}, the issuer may choose to pay the call
price to the bond holder and end the life of the contract.
}
\details{
When a bond is \emph{putable}, the bond holder may choose to force the
issuer pay the put price to the bond holder thus ending the life of the contract.
}
\section{Fields}{

\describe{
\item{\code{calls}}{A data.frame of details for each call.  It should have the columns \code{call_price} and \code{effective_time}.}

\item{\code{puts}}{A data.frame of details for each put.  It should have the columns \code{put_price} and \code{effective_time}.}
}}

\section{Methods}{

\describe{
\item{\code{critical_times()}}{Important times in the life of this instrument for simulation and grid solvers}
}}

\concept{bond}
\concept{callable}
\concept{putable}
