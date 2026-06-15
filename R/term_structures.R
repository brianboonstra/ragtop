## ragtop -- convertibles pricing in R
##
## Copyright (C) 2016-2025  Brian Boonstra <ragtop@boonstra.org>
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
#' @export
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
    if (any(x < 0)) {
      stop("Negative time passed to cumulative variance function in variance_cumulation_from_vols()")
    }
    # findInterval lands each x in [augmented_t[k], augmented_t[k+1]); for
    #  nonnegative x it is in 1..N+1.  The x==0 case (k=1) and the x>=max_t
    #  case (k=N+1) both fall out of the one formula once we clamp the forward
    #  vol index to N, since vols_df$fwd_vols has only N entries.
    k = findInterval(x, augmented_t)
    fwd_ix = pmin(k, N)
    cmvar = cumulated_variances[k] + vols_df$fwd_vols[fwd_ix]^2 * (x - augmented_t[k])
    flog.debug("Cumulative variance at times %s (anchor indexes %s) is %s",
               toString(x), toString(k), toString(cmvar),
               name='ragtop.term_structures.variance_cumulation_from_vols')
    cmvar
  }
  cumul_var = function(T,t=0) {
    cv = cumul_var_0(T) - cumul_var_0(t)
    if (any((T>t) & (cv<=0))) {
      stop("Nonsensical cumulative variance ", toString(cv), " from t=", toString(t), " to T=", toString(T))
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
#' @export
spot_to_df_fcn = function(yield_curve) {
  yield_curve$dfs = exp(-yield_curve$time*yield_curve$rate)
  later = yield_curve$dfs[2:length(yield_curve$dfs)]
  sooner = yield_curve$dfs[1:(length(yield_curve$dfs)-1)]
  fwd_rates = -log(later/sooner)/diff(yield_curve$time)
  yield_curve$fwd_rate = fwd_rates[[length(fwd_rates)]]
  yield_curve$fwd_rate[1:(length(yield_curve$fwd_rate)-1)] = fwd_rates
  ycdf = function(x) {
    n = findInterval(x, yield_curve$time)
    loc_df = numeric(length(x))
    # Times at or beyond the first curve knot (n>0) use the piecewise-constant
    #  forward rate from the bracketing knot; earlier times use the first spot rate.
    on_curve = n > 0
    loc_df[on_curve] = yield_curve$dfs[n[on_curve]] *
      exp(-yield_curve$fwd_rate[n[on_curve]] * (x[on_curve] - yield_curve$time[n[on_curve]]))
    loc_df[!on_curve] = exp(-yield_curve$rate[[1]] * x[!on_curve])
    loc_df
  }
  treasury_df_fcn = function(T,t=0,...) {ycdf(T)/ycdf(t)}
  treasury_df_fcn
}

#' Get a US Treasury curve discount factor function
#'
#' @param on_date Date for which to query for the curve, year-month-day format
#' @return A function taking two time arguments, which returns the discount factor from the second to the first
#' @export
treasury_df_raw = function(on_date) {
  if (requireNamespace('treasury', quietly = TRUE)) {
    on_date = lubridate::ymd(on_date)
    yield_curve_df = treasury::tr_yield_curve(date=format(on_date, "%Y%m"))
    maturities_used = c("1 month", "3 month", "6 month", "1 year", "2 year", "3 year", "5 year", "7 year", "10 year", "20 year", "30 year")
    ix_rows = ((yield_curve_df$date == format(on_date, "%Y-%m-%d") )
               & (yield_curve_df$maturity %in% maturities_used))
    order_rows = match(yield_curve_df$maturity[ix_rows], maturities_used)
    yield_curve_elems = yield_curve_df$rate[ix_rows][order(order_rows)]
    yc_rates = as.numeric(yield_curve_elems)/100  # Values are reported as percent
    yield_curve = data.frame(time=c(0, 30/360, 90/360, 1/2, 1,2,3,5,7,10,20,30), rate=c(0,yc_rates))
    df_frame = spot_to_df_fcn(yield_curve)
  } else {
    flog.error('treasury package not available for treasury curve queries')
    df_frame = data.frame()
  }
  df_frame
}

#' Get a US Treasury curve discount factor function
#'
#' This is a caching wrapper for \code{\link{treasury_df_raw}}
#'
#' @param ... Arguments passed to \code{\link{treasury_df_raw}}
#' @param envir Environment passed to \code{\link{treasury_df_raw}}
#' @return A function taking two time arguments, which returns the discount factor from the second to the first (see \code{spot_to_df_fcn})
#' @export
treasury_df = function(...,envir=parent.frame()) {
  treasury_df_raw(...)
}
if (requireNamespace('R.cache', quietly = TRUE)) {
  treasury_df = R.cache::addMemoization(treasury_df_raw)
}

#' Convert output of BondValuation::AnnivDates to inputs for Bond
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
#' @export
detail_from_AnnivDates = function(anvdates, as_of=Sys.time(), normalization_factor=365.25) {
  as_of = as.POSIXct(as_of)
  payment_time = as.numeric(as.POSIXct(anvdates$PaySched$CoupDates) - as_of)/normalization_factor
  payment_size = anvdates$PaySched$CoupPayments
  coupons = data.frame(payment_time=payment_time, payment_size=payment_size)
  maturity_time = as.numeric(as.POSIXct(anvdates$Traits$Mat) - as_of)/normalization_factor
  parity = anvdates$Traits$Par
  list(notional=parity, coupons=coupons, maturity=maturity_time)
}
