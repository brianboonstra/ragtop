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

#' Shift grid values to account for changes in stock price from dividends paid
#'
#' Derivative values on the grid are presumed to be set according to
#' ex-dividend stock prices.  We interpolate them to cum-dividend stock prices
#' using a natural spline.
#'
#' @return Grid values for the derivatives in the same shape
#'   as \code{grid_values_before_shift}, matched to cum-divided stock prices
shift_for_dividends = function(grid_values_before_shift, stock_prices, div_sum)
{
  shifted_prices = stock_prices - div_sum
  # Interpolate and extrapolate using a natural cubic spline
  # In theory a linear interpolation is fine but in practice
  # a smoother interpolator gives better results with little
  # computational penalty, so we use splines
  #TODO: We can actually let grid_values_before_shift be a matrix
  # with columns for multiple layers, but we would then want to apply()
  # the spline interpolation over each column
  spl = stats::spline(
    x = stock_prices, y = grid_values_before_shift,
    xout = shifted_prices,
    method = "natural"
  ) # natural method allows linear extrapolation
  new_grid_values = spl$y
  new_grid_values
}

#' Find the total of all dividends, according to stock price
#'
#' For each of the N elements of S/h find the sum of the
#' given M dividends, discounted to t_final by r and h
#' The dividends are expected to be represented by a data frame with
#' columns 'time', 'fixed' and 'proportional' where the proportional amount
#' is a multiplier of the stock price over S0 to form conditional dividend size
#' Returns a vector with one dividend sum per entry in S
#'
#' @return The sum of dividend value, corrected for default intensity \code{h} and
#'   short rate \code{r}, in the same shape as stock prices \code{S}
time_adj_dividends = function(relevant_divs, t_final, r, h, S, S0)
{
  if (length(h)==1 && length(S)>1) {
    # Force h and S to have the same shape
    h = h + 0*S
  }
  inv_discount_factor = exp(+(t_final - relevant_divs$time) %o% (r + h)) # positive/inverse because we are carrying forward
  fixed_amt = inv_discount_factor * relevant_divs$fixed
  proportional_contrib = relevant_divs$proportional * inv_discount_factor
  proportional_amt = t(t(proportional_contrib) * (S / S0))
  div_amt_by_price = fixed_amt + proportional_amt
  if (is.null(dim(div_amt_by_price)) || length(dim(div_amt_by_price))==1) {
    div_sum = div_amt_by_price
  } else {
    div_sum = colSums(div_amt_by_price)
  }
  div_sum
}

#' Adjust grid values according to any dividends paid during a time interval
#'
#' @param grid_values A set of derivative prices adapted to the shape of \code{S}
#' @param r The short rate of interest during the time interval
#' @param h The default intensity during the time interval
#' @param t The smallest time in the time interval
#' @param dt The distance from \code{t} to the largest time in the time
#'   interval.  That is, our time interval is \code{(t, t+dt])}
#' @param dividends A \code{data.frame} with columns \code{time}, \code{fixed},
#'   and \code{proprtional}.  Dividend size at the given \code{time} is
#'   then expected to be equal to \code{fixed + proportional * S / S0}
#' @param S A set of stock prices from which we can infer the size of
#'   proportional dividends
#' @return A set of adjusted grid values in the same shape as \code{grid_values}
adjust_for_dividends = function(grid_values, t, dt, r, h, S, S0, dividends)
{
  gs1 = size_in_dimension(grid_values,1)
  gs2 = size_in_dimension(grid_values,2)
  rS = nrow(S)
  cS = ncol(S)
  flog.info("adjust_for_dividends() grid_values %s by %s, S %s by %s",
            gs1, gs2, rS, cS)
  num_layers = size_in_dimension(grid_values,2)
  num_grid_vals_per_layer = size_in_dimension(grid_values,1)
  num_underlying_vals = size_in_dimension(S,1)
  if ( num_grid_vals_per_layer!=num_underlying_vals ) {
    stop(paste("grid values length",num_grid_vals_per_layer,
               "must match underlying prices length",num_underlying_vals))
  }
  if (!is.blank(dividends)) {
    # Divs are present. Sum any relevant ones
    div_sum = 0 * S
    included_ix = (dividends$time > t)  & (dividends$time <= t + dt)
    relevant_divs = dividends[included_ix,c('time', 'fixed', 'proportional')]
    if (nrow(relevant_divs) > 0) {
      flog.info("Found %s dividends from t=%s to t+dt=%s",
                nrow(relevant_divs), t, t + dt)
      div_sum = time_adj_dividends(relevant_divs, t + dt, r, h, S, S0)
    }
    if (any(div_sum != 0)) {
      if (is.null(dim(grid_values)) || length(dim(grid_values))==1 || dim(grid_values)[2]==1) {
        flog.info("Grid values for div adjustment have just one dimension: %s",dput(grid_values))
        grid_values = shift_for_dividends(grid_values, S, div_sum)
      } else {
        flog.info("Grid values for div adjustment have %s layers", num_layers)
        grid_values = apply(grid_values, 2, shift_for_dividends, S, div_sum) # MARGIN=2 means to apply column-wise
      }
    }
  }
  grid_values
}

#' Dividend rate equivalent to discrete dividends
#' @return A scalar consistent with a continuous dividend model having the same
#'  present values as the discrete ones found in \code{dividends}
effective_dividend_rate = function(Tmax, S0, discount_factor_fcn, dividends, t=0)
{
  included_ix = (dividends$time > t)  & (dividends$time <= Tmax)
  relevant_divs = dividends[included_ix,c('time', 'fixed', 'proportional')]
  relevant_divs['discount_factor'] = discount_factor_fcn(relevant_divs$time, t=t)
  # We are not varying stock price so proportional dividends are taken
  # at face value
  relevant_divs['total_size'] = relevant_divs['fixed'] + relevant_divs['proportional']
  relevant_divs['present_value'] = relevant_divs['total_size'] * relevant_divs['discount_factor']
  total_present_value = sum(relevant_divs['fwd_value'])
  if (total_present_value >= S) {
    stop('Specified dividends are so large they imply a negative stock price')
  }
  q = -log(1-total_present_value/S0)
  q
}

#' Rolled up value of coupons to a given time
#'
#' @param t Time to which coupon discounting should be applied, and after which coupons will not be included
#' @param coupons_df A data.frame of details for each coupon.  It should have the
#' fields \code{payment_time} and \code{payment_size}.
#' @param discount_factor_fcn A function used to form values of coupons as of \code{t}.  It should take
#'   two arguments, a target time and a source time.  This function will always be
#'   given \code{t} as its second argument and coupon payment times in the first argument.
#' @return The total present value as of time \code{t} for coupons paid since the \code{model_t}
value_from_prior_coupons = function(t, coupons_df, discount_factor_fcn, model_t=0)
{
  coups = coupons_df[(coupons_df$payment_time<=t) & (coupons_df$payment_time>model_t),]
  disc_factors = discount_factor_fcn(coups$payment_time, t)
  pvs = coups$payment_size / disc_factors
  sum(pvs)
}

#' Coupon present value according to an acceleration schedule from terms and conditions
#'
#' Compute "present" value as of time \code{t} for coupons that
#' would otherwise have been paid up to time \code{acceleration_t}, in the
#' case of accelerated coupon provisions for forced conversions (or
#' sometimes even unforced ones).
#'
#' @param t Time to which coupon discounting should be applied, and before which coupons will not be included
#' @param coupons_df A data.frame of details for each coupon.  It should have the
#' fields \code{payment_time} and \code{payment_size}.
#' @param discount_factor_fcn A function used to form values of coupons as of \code{t}.  It should take
#'   two arguments, a target time and a source time.  This function will always be
#'   given \code{t} as its first argument and coupon payment times in the second argument.
#' @param acceleration_t Coupons with payment time after \code{acceleration_t} will not be included.
#' @return Sum of values as of \code{t}
accelerated_coupon_value = function(t, coupons_df, discount_factor_fcn, acceleration_t=Inf)
{
  coups = coupons_df[(coupons_df$payment_time<=acceleration_t) & (coupons_df$payment_time>t),]
  disc_factors = discount_factor_fcn(t, coups$payment_time)
  pvs = coups$payment_size / disc_factors
  sum(pvs)
}

#' Total coupon value in the case of exercising a conversion option or call
#'
#' Compute "present" value as of time \code{t} of past coupons paid, plus
#' the sum of coupons that would otherwise
#' have been paid up to time \code{acceleration_t}, in the
#' case of calls or of accelerated coupon provisions for
#' conversions (which usually apply only to "forced" conversions).
#'
#' @inheritParams accelerated_coupon_value
#' @param discount_factor_fcn A function used to form values of past coupons as of \code{t}.  It should take
#'   two arguments, a target time and a source time.  This function will always be
#'   given \code{t} as its first argument and coupon payment times in the second argument.
#' @param acceleration_discount_factor_fcn A function used to form values of accelerated future coupons as of \code{t}.  It should take
#'   two arguments, a target time and a source time.  This function will always be
#'   given \code{t} as its first argument and coupon payment times in the second argument.
coupon_value_at_exercise = function(t, coupons_df, discount_factor_fcn,
                                    accelerate_future_coupons=FALSE,
                                    acceleration_discount_factor_fcn=discount_factor_fcn,
                                    acceleration_t=Inf)
{
  if (accelerate_future_coupons) {
    accel_value = accelerated_coupon_value(t, coupons_df,
                                           discount_factor_fcn=acceleration_discount_factor_fcn,
                                           acceleration_t=acceleration_t)
  } else {
    accel_value = 0
  }
  past_value = value_from_prior_coupons(t, coupons_df,
                                        discount_factor_fcn=discount_factor_fcn,
                                        model_t=model_t)
  past_value + accel_value
}
