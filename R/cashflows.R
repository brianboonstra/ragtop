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

#' Shift a set of grid values for dividends paid, using spline interpolation
#'
#' @param grid_values_before_shift Values on grid before accounting for
#'   expected dividends
#' @param stock_prices Stock prices for which to shift the grid
#' @param div_sum Sum of dividend values at each grid point
#' @return An object like \code{grid_values_before_shift} with entries shifted
#'   according to the dividend sums
#' @family Dividends
shift_for_dividends = function(grid_values_before_shift, stock_prices, div_sum)
{
  ## In theory a linear interpolation is fine but in practice
  ## a smoother interpolator gives better results with little
  ## computational penalty, so we use splines
  shifted_prices = stock_prices - div_sum
  # Interpolate and extrapolate using a natural cubic spline
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

#' Find the sum of time-adjusted dividend values
#'
#' For each of the N elements of \code{S/h} find the sum of the
#' given M dividends, discounted to \code{t_final} by \code{r} and \code{h}
#'
#' @param relevant_divs A \code{data.frame} with columns \code{time}, \code{fixed},
#'   and \code{proportional}.  Dividend size at the given \code{time} is
#'   then expected to be equal to \code{fixed + proportional * S / S0}
#' @param S Stock prices
#' @param h Default intensities
#' @param t_final Time beyond which to ignore dividends
#' @inheritParams blackscholes
#' @family Dividends
#' @return Sum of dividends, at each grid node
#' @export time_adj_dividends
time_adj_dividends = function(relevant_divs, t_final, r, h, S, S0)
{
  ## Returns a vector with one dividend sum per entry in S
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

#' Find the sum of time-adjusted dividend values and adjust grid prices according
#'  to their size in the given interval
#'
#' Analyze \code{dividends} to find ones paid in the interval \code{(t,t+dt]}. Form
#'  present value as of time \code{t} for them, and then use spline interpolation
#'  to adjust instrument values accordingly.
#'
#' @inheritParams take_implicit_timestep
#' @inheritParams blackscholes
#' @param S0 Time zero price of the base equity
#' @param dt Interval to end of timestep
#' @param h Default intensities
#' @param grid_values A \code{matrix} with one row for each level of \code{S} and one
#'  column per set of \code{S}-associated instrument values
#' @param dividends A \code{data.frame} with columns \code{time}, \code{fixed},
#'   and \code{proportional}.  Dividend size at the given \code{time} is
#'   then expected to be equal to \code{fixed + proportional * S / S0}
#' @return An object like \code{grid_values} with entries modified according to the dividends
#' @import futile.logger
#' @family Dividends
#' @export adjust_for_dividends
adjust_for_dividends = function(grid_values, t, dt, r, h, S, S0, dividends)
{
  ## Grid values are expected to come from f[k,m-1,]
  flog.info("adjust_for_dividends() grid_values %s by %s, S length %s",
            size_in_dimension(grid_values,1), size_in_dimension(grid_values,2), length(S),
            name='ragtop.cashflows.dividends')
  num_layers = size_in_dimension(grid_values,2)
  num_grid_vals_per_layer = size_in_dimension(grid_values,1)
  num_underlying_vals = size_in_dimension(S,1)
  if ( num_grid_vals_per_layer!=num_underlying_vals ) {
    stop(paste("grid values length",num_grid_vals_per_layer,
               "must match underlying prices length",num_underlying_vals))
  }
  if (is.blank(dividends)) {
    flog.debug("No discrete dividends to treat",
               name='ragtop.cashflows.dividends')
  } else {
    # Divs are present. Sum any relevant ones
    div_sum = 0 * S
    included_ix = (dividends$time > t)  & (dividends$time <= t + dt)
    relevant_divs = dividends[included_ix,c('time', 'fixed', 'proportional')]
    if (nrow(relevant_divs) > 0) {
      flog.info("Found %s dividends from t=%s to t+dt=%s",
                nrow(relevant_divs), t, t + dt,
                name='ragtop.cashflows.dividends')
      div_sum = time_adj_dividends(relevant_divs, t + dt, r, h, S, S0)
    } else {
      flog.info("No discrete dividends to treat in time interval (%s, %s]",
                 t, t+dt,
                name='ragtop.cashflows.dividends')
    }
    if (any(div_sum != 0)) {
      if (is.null(dim(grid_values)) || length(dim(grid_values))==1 || dim(grid_values)[2]==1) {
        flog.info("Grid values for div adjustment have just one dimension: length(dim(grid_values))=%s, dim(grid_values)=%s",
                  length(dim(grid_values)), dim(grid_values),
                  name='ragtop.cashflows.dividends')
        grid_values = as.matrix(shift_for_dividends(grid_values, S, div_sum), ncol=1)
      } else {
        flog.info("Grid values for div adjustment have %s layers",
                  num_layers,
                  name='ragtop.cashflows.dividends')
        grid_values = apply(grid_values, 2, shift_for_dividends, S, div_sum) # MARGIN=2 means to apply column-wise
      }
    }
  }
  grid_values
}

#' Present value of past coupons paid
#'
#' Present value as of time \code{t} for coupons paid since the \code{model_t}
#'
#' @param model_t The payment time beyond which coupons will be included in this computation
#' @param t The time toward which all coupons should be present valued
#' @param coupons_df A data.frame of details for each coupon.  It should have the
#'   columns \code{payment_time} and \code{payment_size}.
#' @param discount_factor_fcn A function specifying how the contract says future coupons should be discounted for this instrument in case the acceleration clause is triggered
#' @family Bond Coupons
value_from_prior_coupons = function(t, coupons_df, discount_factor_fcn, model_t=0)
{
  coups = coupons_df[(coupons_df$payment_time<=t) & (coupons_df$payment_time>model_t),]
  disc_factors = discount_factor_fcn(coups$payment_time, t)
  pvs = coups$payment_size / disc_factors
  sum(pvs)
}

#' Present value of coupons according to an acceleration schedule
#'
#' Compute "present" value as of time t for coupons that
#' would otherwise have been paid up to time acceleration_t, in the
#' case of accelerated coupon provisions for forced conversions (or
#' sometimes even unforced ones).
#' @inheritParams value_from_prior_coupons
#' @param acceleration_t Time limit, up to which coupons will be accelerated
#' @param discount_factor_fcn A function specifying how the contract says future coupons should be discounted for this instrument in case the acceleration clause is triggered
#' @param coupons_df A data.frame of details for each coupon.  It should have the
#'   columns \code{payment_time} and \code{payment_size}.
#' @family Bond Coupons
#' @family Bond Coupon Acceleration
#' @export accelerated_coupon_value
accelerated_coupon_value = function(t, coupons_df, discount_factor_fcn,
                                    acceleration_t=Inf)
{
  coups = coupons_df[(coupons_df$payment_time<=acceleration_t) & (coupons_df$payment_time>t),]
  disc_factors = discount_factor_fcn(t, coups$payment_time)
  pvs = coups$payment_size / disc_factors
  sum(pvs)
}

#' Present value of coupons according to an acceleration schedule
#'
#' Compute "present" value as of time t for coupons that
#' would otherwise have been paid up to time \code{acceleration_t}, in the
#' case of accelerated coupon provisions for forced conversions (or
#' sometimes even unforced ones).
#' @inheritParams value_from_prior_coupons
#' @inheritParams accelerated_coupon_value
#' @param discount_factor_fcn A function specifying how future cashflows should generally be discounted for this instrument
#' @param acceleration_t The maximum time up to which future coupons will be counted for acceleration, passed on to \code{\link{accelerated_coupon_value}}
#' @param acceleration_discount_factor_fcn A function specifying how future coupons should be discounted for this instrument under coupon acceleration conditions
#' @param accelerate_future_coupons If \code{TRUE}, future coupons will be accelerated on exercise to pad present value
#' @param model_t Model timestamp passed to \code{\link{value_from_prior_coupons}}
#' @return A scalar equal to the present value
#' @family Bond Coupons
#' @family Bond Coupon Acceleration
#' @export coupon_value_at_exercise
coupon_value_at_exercise = function(t, coupons_df, discount_factor_fcn, model_t=0,
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
