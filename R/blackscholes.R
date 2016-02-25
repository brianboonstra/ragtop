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

#' Vectorized Black-Scholes pricing of european-exercise options
#'
#' Price options according to the famous Black-Scholes formula, with the
#' optional addition of a jump-to-default intensity.
#'
#' Note that if the \code{default_intensity} is set larger than zero then
#'  put-call parity still holds.  Greeks are reduced according to cumulated default
#'  probability.
#'
#' All inputs must either be scalars or have the same nonscalar shape.
#' @param callput 1 for calls, -1 for puts
#' @param S0 initial underlying price
#' @param default_intensity hazard rate of stock default
#' @return A list with elements \describe{
#'   \item{\code{Price}}{The present value(s)}
#'   \item{\code{Delta}}{Sensitivity to underlying price}
#'   \item{\code{Vega}}{Sensitivity to volatility}
#' }
#' @export blackscholes
blackscholes = function(callput, S0, K, r, time, vola, default_intensity=0, divrate=0, borrow_cost=0)
{
  sd = vola*sqrt(time)
  q = divrate + borrow_cost - default_intensity
  d1 = log(S0/K)+(r-q)*T+0.5*sd^2
  d1 = d1/sd
  d2 = d1-sd
  v = callput*(S0*exp(-q*T)*pnorm(callput*d1)-K*exp(-r*T)*pnorm(callput*d2))
  delta = exp(-q*T)*callput*pnorm(callput*d1)
  vega = S0*exp(-q*T)*dnorm(d1)*sqrt(T)*abs(callput)  # Include abs(callput) to properly vectorize
  surv_prob = exp(-default_intensity*T)
  default_value = max(0, -callput*K*exp(-r*T))
  surv_delta = 0
  surv_vega = 0
  ans = list(Price=surv_prob*v + (1-surv_prob)*default_value,
             Delta=surv_prob*delta,
             Vega=surv_prob*vega)
  ans
}
