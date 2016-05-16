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

#' Create a variance cumulation function from a volatility term structure
#'
#' Given a volatility term structure, create a corresponding variance
#'   cumulation function.  The function assumes piecewise constant
#'   forward volatility, with the final such forward volatility
#'   extending to infinity.
#'
#' @param vols_df A data.frame with numeric columns \code{time} (in
#'   increasing order) and \code{volatility} (not decreasing so quickly
#'   as to give negative forward variance)
#' @return A function taking two time arguments, which returns the cumulated
#'   variance from the second to the first
#' @export variance_cumulation_from_vols
variance_cumulation_from_vols = function(vols_df)
{
  N = nrow(vols_df)
  cumulated_variances = c(0, vols_df$volatility^2 * vols_df$time)
  if (any(cumulated_variances<0)) {
    stop("Nonsensical negative forward variance")
  }
  augmented_t = c(0, vols_df$time)
  time_diffs = diff(augmented_t)
  max_t = max(vols_df$time)
  last_vol = vols_df$volatility[[N]]
  fwd_variances = diff(cumulated_variances)/time_diffs
  cumul_var_0 = function(x) {
    if (x==0) {
      cmvar = 0
    } else if (x>max_t) {
      cmvar = cumulated_variances[[N]] + fwd_variances[[N]] * (x-max_t)
    } else {
      k = findInterval(x,augmented_t)  # Will not be larger than N
      dt = (x-augmented_t[[k]])
      cmvar = cumulated_variances[[k]] + fwd_variances[[k]] * dt
      flog.debug("Found k=%s, using dt=%s applied to fwd variance %s and prev var %s",
                 k, dt, fwd_variances[[k]], cumulated_variances[[k]])
    }
    cmvar
  }
  cumul_var = function(T,t=0) {
    cumul_var_0(T) - cumul_var_0(t)
  }
  cumul_var
}

#' Create a discount factor function from a yield curve
#'
#' Use a piecewise constant approximation to the given spot curve to
#'  generate a function capable of returning corresponding discount factors
#'
#' @param yield_curve A data.frame with numeric columns \code{time} (in
#'   increasing order) and \code{rate} (in natural units)
#' @return A function taking two time arguments, which returns the discount factor from the second to the first
#' @export spot_to_df_fcn
spot_to_df_fcn = function(yield_curve) {
  yield_curve$dfs = exp(-yield_curve$time*yield_curve$rate)
  later = yield_curve$dfs[2:length(yield_curve$dfs)]
  sooner = yield_curve$dfs[1:(length(yield_curve$dfs)-1)]
  fwd_rates = -log(later/sooner)/diff(yield_curve$time)
  yield_curve$fwd_rate = fwd_rates[[length(fwd_rates)]]
  yield_curve$fwd_rate[1:(length(yield_curve$fwd_rate)-1)] = fwd_rates
  ycdf = function(x) {
    n = findInterval(x,yield_curve$time)
    dt = (x-yield_curve$time[[n]])
    yield_curve$dfs[[n]] * exp(-yield_curve$fwd_rate[[n]]*dt)
    }
  treasury_df_fcn = function(T,t=0,...) {ycdf(T)/ycdf(t)}
  treasury_df_fcn
}

#' Get a US Treasury curve discount factor function
#'
#' @param on_date Date for which to query Quandl for the curve
#' @return A function taking two time arguments, which returns the discount factor from the second to the first
#' @export Quandl_df_fcn_UST_raw
Quandl_df_fcn_UST_raw = function(on_date='2016-04-18') {
  yield_curve_elems = Quandl::Quandl("USTREASURY/YIELD", start_date=on_date, end_date=on_date)
  yield_curve_elems$Date = NULL
  yc_rates = as.numeric(yield_curve_elems)/100  # Values are reported as percent
  yield_curve = data.frame(time=c(0, 30/360, 90/360, 1/2, 1,2,3,5,7,10,20,30), rate=c(0,yc_rates))
  spot_to_df_fcn(yield_curve)
}

#' Get a US Treasury curve discount factor function
#'
#' This is a caching wrapper for \code{\link{Quandl_df_fcn_UST_raw}}
#'
#' @param on_date Date for which to query Quandl for the curve
#' @return A function taking two time arguments, which returns the discount factor from the second to the first
#' @export Quandl_df_fcn_UST
Quandl_df_fcn_UST = R.cache::addMemoization(Quandl_df_fcn_UST_raw)


