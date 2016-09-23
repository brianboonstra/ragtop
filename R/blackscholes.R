## ragtop -- convertibles pricing in R
##
## Copyright (C) 2016  Brian Boonstra <ragtop@boonstra.org>
##
## This file is part of the ragtop package for GNU R.
## It is made available under the terms of the GNU General Public
## License, version 2, or at your option, any later version,
## incorporated herein by reference.
##
## This program is distributed in the hope that it will be
## useful, but WITHOUT ANY WARRANTY; without even the implied
## warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
## PURPOSE.  See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License
## along with ragtop.  If not, see <http://www.gnu.org/licenses/>.

library(futile.logger)


#' Constant CALL for defining option contracts
#'
#' @export CALL
CALL = 1

#' Constant PUT for defining option contracts
#'
#' @export PUT
PUT = -1

#' Vectorized Black-Scholes pricing of european-exercise options
#'
#' Price options according to the famous Black-Scholes formula, with the
#' optional addition of a jump-to-default intensity and discrete dividends.
#'
#' Note that if the \code{default_intensity} is set larger than zero then
#'  put-call parity still holds.  Greeks are reduced according to cumulated default
#'  probability.
#'
#' All inputs must either be scalars or have the same nonscalar shape.
#' @param callput 1 for calls, -1 for puts
#' @param S0 initial underlying price
#' @param K strike
#' @param r risk-free interest rate
#' @param time Time from \code{0} until expiration
#' @param vola Default-free volatility of the underlying
#' @param default_intensity hazard rate of underlying default
#' @param divrate A continuous rate for dividends and other cashflows such as foreign interest rates
#' @param borrow_cost A continuous rate for stock borrow costs
#' @param dividends A \code{data.frame} with columns \code{time}, \code{fixed},
#'   and \code{proportional}.  Dividend size at the given \code{time} is
#'   then expected to be equal to \code{fixed + proportional * S / S0}.  Fixed
#'   dividends will be converted to proportional for purposes of this algorithm.
#' @return A list with elements \describe{
#'   \item{\code{Price}}{The present value(s)}
#'   \item{\code{Delta}}{Sensitivity to underlying price}
#'   \item{\code{Vega}}{Sensitivity to volatility}
#' }
#' @family European Options
#' @family Equity Independent Default Intensity
#' @examples
#' blackscholes(callput=-1, S0=100, K=90, r=0.03, time=1, # -1 is a PUT
#'              vola=0.5, default_intensity=0.07)
#' @export blackscholes
blackscholes = function(callput, S0, K, r, time, vola,
                        default_intensity=0, divrate=0, borrow_cost=0,
                        dividends=NULL)
{
  if (!is.blank(dividends)) {
    # Divs are present. Sum any relevant ones and get present value.
    # Form an initial sum of the correct shape.
    div_sum = 0 * S0 + 0 * K + 0 * r + 0 * time + 0 * vola + 0 * callput
    included_ix = (dividends$time > 0)  & (dividends$time <= time)
    relevant_divs = dividends[included_ix,c('time', 'fixed', 'proportional')]
    if (nrow(relevant_divs) > 0) {
      # Discount without consideration for default intensity.  For options
      #  the absence of dividends in case of default is handled by terminal
      #  likelihood.
      # Subtracting present value properly adjusts the terminal distibution
      #  in case of proportional dividends.  For fixed dividends, only a
      #  grid scheme can properly handle it.  We just pretend they are
      #  proportional here.
      div_amt = time_adj_dividends(relevant_divs, 0, r,
                                   h=0, S=S0, S0=S0)
      div_sum = div_sum + div_amt
      flog.info("Found %s dividends summing to %s in present value",
                nrow(relevant_divs), div_amt,
                name="ragtop.blackscholes")
      S0 = S0 - div_sum
    }
  }
  flog.debug("blackscholes():  callput %s\nS0=%s\n  K=%s\n  r=%s\n  time=%s\n  vola=%s\n  default_intensity=%s\n  divrate=%s\n  borrow_cost=%s",
             callput, S0, K, r, time, vola,
             default_intensity, divrate, borrow_cost,
             name="ragtop.blackscholes")
  sd = vola*sqrt(time)
  q = divrate + borrow_cost - default_intensity
  d1 = log(S0/K)+(r-q)*time+0.5*sd^2
  d1 = d1/sd
  d2 = d1-sd
  v = callput*(S0*exp(-q*time) * stats::pnorm(callput*d1)-K*exp(-r*time) * stats::pnorm(callput*d2))
  delta = exp(-q*time)*callput * stats::pnorm(callput*d1)
  vega = S0*exp(-q*time) * stats::dnorm(d1)*sqrt(time)*abs(callput)  # Include abs(callput) to properly vectorize
  surv_prob = exp(-default_intensity*time)
  default_value = max(0, -callput*K*exp(-r*time))
  surv_delta = 0
  surv_vega = 0
  ans = list(Price=surv_prob*v + (1-surv_prob)*default_value,
             Delta=surv_prob*delta,
             Vega=surv_prob*vega)
  ans
}

#' Black-Scholes pricing of european-exercise options with term structure arguments
#'
#' Price an option according to the famous Black-Scholes formula, with the
#' optional addition of a jump-to-default intensity and discrete dividends.  Volatility
#' and rates may be provided as constants or as 2+ parameter functions with
#' first argument \code{T} corresponding to maturity and second argument \code{t} corresponding to
#' model date.
#'
#' Any term structures will be converted to equivalent constant arguments by calling
#' them with the arguments \code{(time, 0)}.
#'
#' @inheritParams blackscholes
#' @param dividend_rate A continuous rate for dividends and other cashflows such as foreign interest rates
#' @param discount_factor_fcn A function for computing present values to
#'   time \code{t} of various cashflows occurring during this timestep, with
#'   arguments \code{T}, \code{t}
#' @param const_volatility A constant to use for volatility in case \code{variance_cumulation_fcn}
#'  is not given
#' @param const_short_rate A constant to use for the instantaneous interest rate in case \code{discount_factor_fcn}
#'  is not given
#' @param const_default_intensity A constant to use for the instantaneous default intensity in case \code{default_intensity_fcn}
#'  is not given
#' @param survival_probability_fcn A function for probability of survival, with
#'   arguments \code{T}, \code{t} and \code{T>t}.  E.g. with
#'   a constant volatility \eqn{s} this takes the form \eqn{(T-t)s^2}.
#' @param variance_cumulation_fcn A function for computing total stock variance
#'   occurring during this timestep, with arguments \code{T}, \code{t}.  E.g. with
#'   a constant volatility \eqn{s} this takes the form \eqn{(T-t)s^2}.
#' @family European Options
#' @family Equity Independent Default Intensity
#' @examples
#' black_scholes_on_term_structures(callput=-1, S0=100, K=90, time=1,
#'                                  discount_factor_fcn = function(T, t, ...) {
#'                                    exp(-0.03 * (T - t))
#'                                  },
#'                                  survival_probability_fcn = function(T, t, ...) {
#'                                    exp(-0.07 * (T - t))
#'                                  },
#'                                  variance_cumulation_fcn = function(T, t) {
#'                                    0.45 ^ 2 * (T - t)
#'                                  })
#' @import futile.logger
#' @export black_scholes_on_term_structures
black_scholes_on_term_structures = function(callput, S0, K, time,
           const_volatility=0.5, const_short_rate=0, const_default_intensity=0,
           discount_factor_fcn = function(T, t, ...){exp(-const_short_rate*(T-t))},
           survival_probability_fcn = function(T, t, ...){exp(-const_default_intensity*(T-t))},
           variance_cumulation_fcn = function(T, t){const_volatility^2*(T-t)},
           dividends=NULL,
           borrow_cost=0.0,
           dividend_rate=0.0)
{
  if (length(time)>1) {
    stop("Only a single expiration time may be provided.")
  }
  if (is.na(time) || time<=0) {
    stop("Expiration time must be strictly positive.")
  }
  flog.debug("black_scholes_on_term_structures():  callput %s\nS0=%s\n  K=%s\n  time=%s\n  divrate=%s\n  borrow_cost=%s",
             callput, S0, K, time, dividend_rate, borrow_cost,
             name="ragtop.blackscholes")
  vola = sqrt(variance_cumulation_fcn(time, 0)/time)
  std_time = time
  r = -log(discount_factor_fcn(time,0))/time
  h = -log(survival_probability_fcn(time,0))/time
  flog.debug("black_scholes_on_term_structures() finds constant equivalents r=%s, h=%s and vola=%s",
             r, h, vola,
             name="ragtop.blackscholes")
  bs = blackscholes(callput=callput, S0=S0, K=K, r=r, time=std_time, vola=vola,
                    default_intensity=h, divrate=dividend_rate, borrow_cost=borrow_cost,
                    dividends=dividends)
  flog.debug("black_scholes_on_term_structures() finds P=%s D=%s Vg=%s",
             bs$Price, bs$Delta, bs$Vega,
             name="ragtop.blackscholes")
  bs
}
