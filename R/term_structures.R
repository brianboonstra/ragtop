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
#' @examples
#' vc = variance_cumulation_from_vols(
#'   data.frame(time=c(0.1,2,3),
#'   volatility=c(0.2,0.5,1.2)))
#' vc(1.5, 0)
#'
#' @export variance_cumulation_from_vols
variance_cumulation_from_vols = function(vols_df)
{
  N = nrow(vols_df)
  cumulated_variances = c(0, vols_df$volatility^2 * vols_df$time)
  if (any(cumulated_variances<0)) {
    stop("Nonsensical variance")
  }
  fwd_variances = diff(cumulated_variances)
  if (any(fwd_variances<0)) {
    stop("Nonsensical negative forward variance")
  }
  augmented_t = c(0, vols_df$time)
  time_diffs = diff(augmented_t)
  max_t = max(vols_df$time)
  last_vol = vols_df$volatility[[N]]
  vols_df$fwd_vols = sqrt(fwd_variances/time_diffs)
  cumul_var_0 = function(x) {
    if (x==0) {
      cmvar = 0
    } else if (x>=max_t) {
      cmvar = cumulated_variances[[N+1]] + vols_df$fwd_vols[[N]]^2 * (x-max_t)
      flog.debug("Found %s was beyond max_t %s, N=%s, using time diff %s from last anchor time max_t=%s applied to fwd vol %s and prev var %s",
                 x, max_t, N, (x-max_t), max_t, vols_df$fwd_vols[[N]], cumulated_variances[[N]],
                 name='ragtop.term_structures.variance_cumulation_from_vols')
    } else {
      k = findInterval(x, augmented_t)  # Will not be larger than N
      dt = x - augmented_t[[k]]
      if (dt<0) {
        stop("Negative time interval after call to findInterval() in variance_cumulation_from_vols()")
      }
      cmvar = cumulated_variances[[k]] + vols_df$fwd_vols[[k]]^2 * dt
      flog.debug("Found k=%s, using time diff %s from prev anchor time %s applied to fwd vol %s and prev var %s",
                 k, dt, augmented_t[[k]], vols_df$fwd_vols[[k]], cumulated_variances[[k]],
                 name='ragtop.term_structures.variance_cumulation_from_vols')
    }
    cmvar
  }
  cumul_var = function(T,t=0) {
    cv = cumul_var_0(T) - cumul_var_0(t)
    if ((T>t) && (cv<=0)) {
      stop("Nonsensical cumulative variance ", cv, " from t=", t, " to T=", T)
    }
    cv
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
#' @examples
#' disct_fcn = ragtop::spot_to_df_fcn(
#'   data.frame(time=c(1, 5, 10, 15),
#'              rate=c(0.01, 0.02, 0.03, 0.05)))
#' print(disct_fcn(1, 0.5))
#' @export spot_to_df_fcn
spot_to_df_fcn = function(yield_curve) {
  yield_curve$dfs = exp(-yield_curve$time*yield_curve$rate)
  later = yield_curve$dfs[2:length(yield_curve$dfs)]
  sooner = yield_curve$dfs[1:(length(yield_curve$dfs)-1)]
  fwd_rates = -log(later/sooner)/diff(yield_curve$time)
  yield_curve$fwd_rate = fwd_rates[[length(fwd_rates)]]
  yield_curve$fwd_rate[1:(length(yield_curve$fwd_rate)-1)] = fwd_rates
  ycdf = function(x) {
    loc_df = NA
    n = findInterval(x, yield_curve$time)
    if (n>0) {
      dt = (x-yield_curve$time[[n]])
      loc_df = yield_curve$dfs[[n]] * exp(-yield_curve$fwd_rate[[n]]*dt)
    } else {
      loc_df = exp(-yield_curve$rate[[1]]*x)
    }
    loc_df
  }
  treasury_df_fcn = function(T,t=0,...) {ycdf(T)/ycdf(t)}
  treasury_df_fcn
}

#' Get a US Treasury curve discount factor function
#'
#' @param on_date Date for which to query Quandl for the curve
#' @return A function taking two time arguments, which returns the discount factor from the second to the first
#' @export Quandl_df_fcn_UST_raw
Quandl_df_fcn_UST_raw = function(on_date) {
  if (is.element('R.cache', utils::installed.packages()[,1])) {
    yield_curve_elems = Quandl::Quandl("USTREASURY/YIELD", start_date=on_date, end_date=on_date)
    yield_curve_elems$Date = NULL
    yc_rates = as.numeric(yield_curve_elems)/100  # Values are reported as percent
    yield_curve = data.frame(time=c(0, 30/360, 90/360, 1/2, 1,2,3,5,7,10,20,30), rate=c(0,yc_rates))
    df_frame = spot_to_df_fcn(yield_curve)
  } else {
    flog.error('Quandl package not available for treasury curve queries')
    df_frame = data.frame()
  }
  df_frame
}

#' Get a US Treasury curve discount factor function
#'
#' This is a caching wrapper for \code{\link{Quandl_df_fcn_UST_raw}}
#'
#' @param ... Arguments passed to \code{\link{Quandl_df_fcn_UST_raw}}
#' @param envir Environment passed to \code{\link{Quandl_df_fcn_UST_raw}}
#' @return A function taking two time arguments, which returns the discount factor from the second to the first
#' @export Quandl_df_fcn_UST
Quandl_df_fcn_UST = function(...,envir=parent.frame()) {
  Quandl_df_fcn_UST_raw(...)
}
if (is.element('R.cache', utils::installed.packages()[,1])) {
  Quandl_df_fcn_UST = R.cache::addMemoization(Quandl_df_fcn_UST_raw)
}

#' Convert output of BondValuation::AnnivDates to inputd for Bond
#'
#' The BondValuation package provides day count convention treatments superior
#'  to quantmod or any other R package known (as of May 2019).  This function
#'  takes output from BondValuation::AnnivDates(...) and parses it into
#'  notionals, maturity time, and coupon times and sizes.
#'
#' Note: volatilities used in `ragtop` must have compatible time units to these times.
#'
#' @param anvdates Output of BondValuation::AnnivDates(), which must have included a `Coup` argument so that the resulting list contains an entry for `PaySched`
#' @param as_of Date or time from whose perspective times should be computed
#' @param normalization_factor Factor by which raw R time differences should be multiplied.  If volatilites are going to be annualized, then this should typically be 365 or so.
#' @return A list with some of the arguments appropriate for defining a Bond as follows:
#'             maturity - maturity
#'             notional - notional amount
#'             coupons - `data.frame` with `payment_time`, `payment_size`
#' @export detail_from_AnnivDates
detail_from_AnnivDates = function(anvdates, as_of=Sys.time(), normalization_factor=365.25) {
  as_of = as.POSIXct(as_of)
  payment_time = as.numeric(as.POSIXct(anvdates$PaySched$CoupDates) - as_of)/normalization_factor
  payment_size = anvdates$PaySched$CoupPayments
  coupons = data.frame(payment_time=payment_time, payment_size=payment_size)
  maturity_time = as.numeric(as.POSIXct(anvdates$Traits$Mat) - as_of)/normalization_factor
  parity = anvdates$Traits$Par
  list(notional=parity, coupons=coupons, maturity=maturity_time)
}
