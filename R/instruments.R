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

## Define classes to represent various financial instruments
## to grid pricing schemes.
## Our main motive for using objects is that derivatives, by
## definition, suffer interdependencies between instrument values
## and so it is nearly impossible to make pricing models
## for complex instruments without tracking state.
##
## With apologies to all the world's wonderful functional
## programmers.

## We are currently using reference classes for instruments
## though R6 classes are a tempting alternative.  Calling
## syntax is close or identical so the change should be
## easy.

## Logging
require(futile.logger)

# TODO: Allow initialization of some instruments from QuantLib classes

#' Representation of financial instrument amenable to grid pricing schemes
#'
#' Our basic instrument defines a tenor/maturity, a method
#' to provide values in case of default, and a method
#' to correct instrument prices in light of exercise
#' decisions.
#'
#' @field maturity The tenor, expiration date or terminal date by which the value of this security will be certain.
#' @field last_computed_grid The most recently computed set of values from a grid pricing scheme.  Used internally for pricing chains of derivatives.
#' @export GridPricedInstrument
#' @exportClass GridPricedInstrument
GridPricedInstrument = setRefClass(
  "GridPricedInstrument",
  fields = list(maturity = "numeric",
                last_computed_grid = "numeric"),
  methods = list(
    recovery_fcn = function(v,S,t,...) {
      "Return recovery value, given non-default values {v} at time {t}.  Subclasses may be more elaborate, this method simply returns 0.0."
      0.0
    },
    optionality_fcn = function(v,...) {
      # To be overridden by subclasses
      "Return a terminal value, or a version of {v} at time {t} corrected for any optionality conditions."
      last_computed_grid <<- as.vector(v)
      v
    }
  )
)

#' An option contract with call or put terms
#'
#' @field strike A decision price for the contract
#' @field callput Either \code{1} for a call or \code{-1} for a put
#' @export EquityOption
#' @exportClass EquityOption
EquityOption = setRefClass(
  "EquityOption",
  contains = "GridPricedInstrument",
  fields = list(strike = "numeric",
                callput = "numeric")
)

#' A standard option contract
#'
#' At maturity, the call option holder will "exercise", i.e. choose stock, with value \code{S}, if the
#' stock price is above the strike \code{K}, paying \code{K} to the option issuer,
#' realizing value \code{S-K}.  The put option holder will exercise, receiving \code{K} while surrendering
#' stock worth \code{S}, if the stock price is below \code{K}.
#'
#' Therefore the value at maturity is equal to \code{max(0,callput*(S-K))}
#' @export EuropeanOption
#' @exportClass EuropeanOption
EuropeanOption = setRefClass(
  "EuropeanOption",
  contains = "EquityOption",
  methods = list(
    optionality_fcn = function(v,S,t,...) {
      "Return {v} up to maturity time.  Return exercise value after that time."
      if (t >= maturity) {
        v = callput * (S - strike)
      }
      v[v < 0] = 0
      last_computed_grid <<- as.vector(v)
      v
    }
  )
)

#' A standard option contract allowing for \emph{early} exercise at the choice of the option holder
#' @export AmericanOption
#' @exportClass AmericanOption
AmericanOption = setRefClass(
  "AmericanOption",
  contains = "EquityOption",
  methods = list(
    optionality_fcn = function(v,S,t,...) {
      "Return the greater of hold value {v} or early exercise value at each stock price level in {S} up to maturity time.  Return exercise value after that time."
      if (t < maturity) {
        exer_value = callput * (S - strike)
        exer_value[exer_value < 0] = 0
        exercise_ix = (v < exer_value)
        v[exercise_ix] = exer_value[exercise_ix]
      } else {
        v = callput * (S - strike)
      }
      v[v < 0] = 0
      last_computed_grid <<- as.vector(v)
      v
    }
  )
)


#' A simple contract paying the \code{notional} amount at the \code{maturity}
#'
#' @field notional The amount that will be paid at \code{maturity}, conditional on survival
#' @field recovery_rate The proportion of notional that would be expected to be paid to bond holders after bankruptcy court proceedings
#' @field discount_factor_fcn A function specifying how cashflows should generally be discounted for this instrument
#' @export ZeroCouponBond
#' @exportClass ZeroCouponBond
ZeroCouponBond = setRefClass(
  "ZeroCouponBond",
  contains="GridPricedInstrument",
  fields=list(notional="numeric",
              recovery_rate="numeric",
              discount_factor_fcn="function"),
  methods = list(
    optionality_fcn = function(v,S,t,...) {
      "Return the notional value in the shape of {S} at any time on or after maturity, otherwise just return {v}"
      if (t >= maturity) {
        v = 0.0*(v+S) + notional
      }
      v[v < 0] = 0
      last_computed_grid <<- as.vector(v)
      v
    }
  )
)

#' Standard corporate or government bond
#'
#' @field coupons A data.frame of details for each coupon.  It should have the
#'   columns \code{payment_time} and \code{payment_size}.
#' @export CouponBond
#' @exportClass CouponBond
CouponBond = setRefClass(
  "CouponBond",
  contains="ZeroCouponBond",
  fields=list(coupons="data.frame"),
  methods = list(
    accumulate_coupon_values_before = function(t,discount_factor_fcn=discount_factor_fcn) {
      "Compute the sum of coupon present values as of {t} according to {discount_factor_fcn}"
      paid_coupon_ix = (coupons$payment_time<=t)
      paid_coupons = coupons[paid_coupon_ix,]
      accumulation_proportions = discount_factor_fcn(t,paid_coupons$payment_time)
      sum(accumulation_proportions*paid_coupons$payment_size)
    },
    critical_times = function() {
      "Important times in the life of this instrument for simulation and grid solvers"
      ctimes = maturity
      if (nrow(coupons) > 0) {
        ctimes = c(ctimes, coupons$payment_time)
      }
      ctimes
    }
  )
)


#' Callable (and putable) corporate or government bond.
#'
#' When a bond is emph{callable}, the issuer may choose to pay the call
#' price to the bond holder and end the life of the contract.
#'
#' When a bond is emph{putable}, the bond holder may choose to force the
#' issuer pay the put price to the bond holder thus ending the life of the contract.
#'
#' @field calls A data.frame of details for each call.  It should have the columns \code{call_price} and \code{effective_time}.
#' @field puts A data.frame of details for each put.  It should have the columns \code{put_price} and \code{effective_time}.
#' @concept callable
#' @concept putable
#' @concept bond
#' @export CallableBond
#' @exportClass CallableBond
CallableBond = setRefClass(
  "CallableBond",
  contains="CouponBond",
  fields=list(calls="data.frame", puts="data.frame"),
  methods=list(
    critical_times = function() {
      "Important times in the life of this instrument for simulation and grid solvers"
      ctimes = maturity
      if (nrow(coupons) > 0) {  # Could just call superclass critical_times()
        ctimes = c(ctimes, coupons$payment_time)
      }
      if (nrow(calls) > 0) {
        ctimes = c(ctimes, calls$effective_time)
      }
      if (nrow(puts) > 0) {
        ctimes = c(ctimes, puts$effective_time)
      }
      ctimes
    }
  )
)


#' Convertible bond with exercise into stock
#'
#' @field conversion_ratio The number of shares, per bond, that result from exercise
#' @field dividend_ceiling The level of dividend protection (if any) specified in terms and conditions
#' @export ConvertibleBond
#' @exportClass ConvertibleBond
ConvertibleBond = setRefClass(
  "ConvertibleBond",
  contains="CallableBond",
  fields=list(conversion_ratio="numeric",
              dividend_ceiling="numeric"),
  methods=list(
    optionality_fcn = function(v,S,t,discount_factor_fcn=discount_factor_fcn,...) {
      "Return the greater of hold value {v} or conversion value at each stock price level in {S}"
      if (t > maturity) {
        v = 0.0*(v+S) + notional
      } else {
        flog.debug("Optionality at t=%s", t)
        exercise_values = S * conversion_ratio
        accumulated_past_coupons = accumulate_coupon_values_before(t,discount_factor_fcn=discount_factor_fcn)
        total_early_exercise_value = exercise_values + accumulated_past_coupons
        flog.debug("exercise_values %s total_early_exercise_value %s",
                  toString(exercise_values), toString(total_early_exercise_value))
        total_early_exercise_value[total_early_exercise_value < 0] = 0
        exercise_ix = (v < total_early_exercise_value)
        v[exercise_ix] = total_early_exercise_value[exercise_ix]
        last_computed_grid <<- as.vector(v)
      }
      v
    }
  )
)
