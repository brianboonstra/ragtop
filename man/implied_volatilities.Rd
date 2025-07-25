% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calibration.R
\name{implied_volatilities}
\alias{implied_volatilities}
\title{Implied volatilities of european-exercise options under Black-Scholes or a jump-process extension}
\usage{
implied_volatilities(
  option_price,
  callput,
  S0,
  K,
  r,
  time,
  const_default_intensity = 0,
  divrate = 0,
  borrow_cost = 0,
  dividends = NULL,
  relative_tolerance = 1e-06,
  max.iter = 100,
  max_vola = 4
)
}
\arguments{
\item{option_price}{Present option values (may be a vector)}

\item{callput}{\code{1} for calls, \code{-1} for puts (may be a vector)}

\item{S0}{initial underlying price (may be a vector)}

\item{K}{strike (may be a vector)}

\item{r}{risk-free interest rate (may be a vector)}

\item{time}{Time from \code{0} until expiration (may be a vector)}

\item{const_default_intensity}{hazard rate of underlying default (may be a vector)}

\item{divrate}{A continuous rate for dividends and other cashflows such as foreign interest rates (may be a vector)}

\item{borrow_cost}{A continuous rate for stock borrow costs (may be a vector)}

\item{dividends}{A \code{data.frame} with columns \code{time}, \code{fixed},
and \code{proportional}.  Dividend size at the given \code{time} is
then expected to be equal to \code{fixed + proportional * S / S0}.  Fixed
dividends will be converted to proprtional for purposes of this algorithm.}

\item{relative_tolerance}{Relative tolerance in option price to achieve before halting the search}

\item{max.iter}{Number of iterations to try before abandoning the search}

\item{max_vola}{Maximum volatility to try in the search}
}
\value{
Scalar volatilities
}
\description{
Find default-free volatilities based on known interest rates and hazard rates, using
  a given option price.
}
\seealso{
Other Implied Volatilities: 
\code{\link{american_implied_volatility}()},
\code{\link{equivalent_bs_vola_to_jump}()},
\code{\link{equivalent_jump_vola_to_bs}()},
\code{\link{fit_variance_cumulation}()},
\code{\link{implied_jump_process_volatility}()},
\code{\link{implied_volatilities_with_rates_struct}()},
\code{\link{implied_volatility}()},
\code{\link{implied_volatility_with_term_struct}()}

Other European Options: 
\code{\link{black_scholes_on_term_structures}()},
\code{\link{blackscholes}()},
\code{\link{implied_volatilities_with_rates_struct}()},
\code{\link{implied_volatility}()},
\code{\link{implied_volatility_with_term_struct}()}

Other Equity Independent Default Intensity: 
\code{\link{american}()},
\code{\link{american_implied_volatility}()},
\code{\link{black_scholes_on_term_structures}()},
\code{\link{blackscholes}()},
\code{\link{equivalent_bs_vola_to_jump}()},
\code{\link{equivalent_jump_vola_to_bs}()},
\code{\link{implied_volatilities_with_rates_struct}()},
\code{\link{implied_volatility}()},
\code{\link{implied_volatility_with_term_struct}()}
}
\concept{Equity Independent Default Intensity}
\concept{European Options}
\concept{Implied Volatilities}
