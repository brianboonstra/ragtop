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

## We are currently using reference classes for instruments,
## though R6 classes are a tempting alternative.  Calling
## syntax is close or identical so a change to R6 should be
## easy.

## Logging
library(futile.logger)

# TODO: Allow initialization of some instruments from QuantLib classes
# TODO: Handle daycount conventions and dates, perhaps using Quantlib

#' Representation of financial instrument amenable to grid pricing schemes
#'
#' Our basic instrument defines a tenor/maturity, a method
#' to provide values in case of default, and a method
#' to correct instrument prices in light of exercise
#' decisions.
#'
#' @field maturity The tenor, expiration date or terminal date by which the value of this security will be certain.
#' @field last_computed_grid The most recently computed set of values from a grid pricing scheme.  Used internally for pricing chains of derivatives.
#' @field name A mnemonic name for the instrument, not used by ragtop
#' @export GridPricedInstrument
#' @import methods
#' @exportClass GridPricedInstrument
GridPricedInstrument = setRefClass(
  "GridPricedInstrument",
  fields = list(maturity = "numeric",
                last_computed_grid = "numeric",
                name = "character"),
  methods = list(
    recovery_fcn = function(v,S,t,...) {
      "Return recovery value, given non-default values {v} at time {t}.  Subclasses may be more elaborate, this method simply returns 0.0."
      0.0
    },
    optionality_fcn = function(v,...) {
      # To be overridden by subclasses, which must update last_computed_grid
      "Return a version of {v} at time {t} corrected for any optionality conditions."
      last_computed_grid <<- as.vector(v)
      v
    },
    terminal_values = function(v,...) {
      # To be overridden by subclasses
      "Return a terminal value. defaults to simply calling {optionality_fcn}."
      optionality_fcn(v,...)
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
  fields = list(discount_factor_fcn="function"),
  methods = list(
    recovery_fcn = function(v,S,t,discount_factor_fctn=discount_factor_fcn,...) {
      "Return 0 for calls and discounted future payout for puts."
      recovery = 0
      if (-1==callput) {
        df = discount_factor_fctn(maturity, t)
        recovery = strike * df
        if (recovery>2*max(S)) {
          flog.warn("Unexpected large recovery value %s, strike %s, df %s, max(S) is %s",
                    recovery, strike, df, max(S),
                    name='ragtop.instruments.recovery.european')
        }
      }
      recovery
    },
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
    recovery_fcn = function(v,S,t,...) {
      "Return 0 for calls and discounted future payout for puts."
      recovery = 0
      if (-1==callput) {
        recovery = strike
      }
      recovery
    },
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
        flog.info("Timestep t=%s for %s is at or beyond maturity.  Using notional %s.",
                  t, name, notional, name="ragtop.instruments.optionality.zcb")
      }
      v[v < 0] = 0
      last_computed_grid <<- as.vector(v)
      v
    },
    recovery_fcn = function(v,S,t,...) {
      "Return recovery rate times notional, except where it exceeds bond value on the S grid."
      recovery = 0.0*(v+S)
      if (t < maturity && (!is.blank(recovery_rate))) {
        recovery = 0.0 * v + recovery_rate * notional
        recovery[v < recovery] = v[v < recovery]
        recovery = recovery + 0 * S  # ensure proper shape
      }
      recovery
    }
  )
)

#' Standard corporate or government bond
#'
#' A coupon bond is treated here as the entire collection of cashflows. In particular,
#'  coupons are included in the package even after they have been paid, accruing
#'  at the risk-free rate.
#'
#' @field coupons A data.frame of details for each coupon.  It should have the
#'   columns \code{payment_time} and \code{payment_size}.
#' @export CouponBond
#' @exportClass CouponBond
CouponBond = setRefClass(
  "CouponBond",
  contains="ZeroCouponBond",
  fields=list(coupons="data.frame", last_computed_cash="numeric"),
  methods = list(
    accumulate_coupon_values_before = function(t, discount_factor_fctn=discount_factor_fcn) {
      "Compute the sum of coupon present values as of {t} according to {discount_factor_fctn}"
      ac = 0
      paid_coupon_ix = (coupons$payment_time<=t+maturity*TIME_RESOLUTION_FACTOR)
      if (sum(paid_coupon_ix) >= 1) {
        paid_coupons = coupons[paid_coupon_ix,]
        accumulation_proportions = sapply(paid_coupons$payment_time,
                                          function(ct) discount_factor_fctn(t, ct))
        ac = sum(paid_coupons$payment_size/accumulation_proportions)
        flog.info("Accumulated coupon values at %s are: %s", t, ac,
                  name="ragtop.instruments.cashflows.bond")
      }
      ac
    },
    total_coupon_values_between = function(small_t, big_t, discount_factor_fctn=discount_factor_fcn) {
      "Compute the sum (as of {big_t}) of present values of coupons paid between small_t and big_t"
      ac = 0
      paid_coupon_ix = ( coupons$payment_time>small_t & coupons$payment_time<=big_t )
      if (sum(paid_coupon_ix) >= 1) {
        paid_coupons = coupons[paid_coupon_ix,]
        accumulation_proportions = sapply(paid_coupons$payment_time,
                                          function(ct) discount_factor_fctn(big_t, ct))
        ac = sum(accumulation_proportions*paid_coupons$payment_size)
        flog.info("Totaled coupon present values in interval (%s,%s] are: %s", small_t, big_t, ac,
                  name="ragtop.instruments.cashflows.bond")
      }
      ac
    },
    critical_times = function() {
      "Important times in the life of this instrument for simulation and grid solvers"
      ctimes = maturity
      if (nrow(coupons) > 0) {
        ctimes = c(ctimes, coupons$payment_time)
      }
      ctimes = sort(ctimes)
      ctimes
    },
    update_cashflows = function(small_t, big_t, discount_factor_fctn=discount_factor_fcn, include_notional=TRUE, ...) {
      "Update last_computed_cash and return cashflow information for the given time period, valued at big_t"
      if (is.blank(last_computed_cash)) {
        last_computed_cash <<- 0
      }
      df = discount_factor_fctn(big_t, small_t)
      cashflows = total_coupon_values_between(small_t, big_t, discount_factor_fctn=discount_factor_fctn)
      if (include_notional && maturity>small_t && maturity<=big_t) {
        ndf = notional * discount_factor_fctn(maturity, big_t)
        flog.info("Totaled %s cashflows in interval (%s,%s]  discounted to t=%s were %s.  This also included notional %s, discounted to %s.",
                  name, small_t, big_t, small_t, cashflows + ndf, notional, ndf,
                  name="ragtop.instruments.cashflows.bond")
        cashflows = cashflows + ndf
      } else {
        flog.info("Totaled %s cashflows in interval (%s,%s]  discounted to t=%s were %s. No notional amount included.",
                  name, small_t, big_t, big_t, cashflows,
                  name="ragtop.instruments.cashflows.bond")
      }
      new_cash = df*(last_computed_cash - cashflows)
      flog.info("Updating %s last_computed_cash (was: %s) discounted by %s to %s and subtracted cashflows %s discounted to %s to get %s.",
                name, last_computed_cash, df, df*new_cash,
                cashflows, df*cashflows, new_cash, name="ragtop.instruments.cashflows.bond")
      last_computed_cash <<- new_cash
      cashflows
    },
    optionality_fcn = function(v,S,t,discount_factor_fctn=discount_factor_fcn,...) {
      "Return the greater of hold value {v} or exercise value at each stock price level in {S}.  If the given date is beyond maturity, return value at maturity."
      if (t >= maturity) {
        accumulated_past_coupons = accumulate_coupon_values_before(maturity, discount_factor_fctn=discount_factor_fctn)
        last_computed_cash <<- notional + accumulated_past_coupons
        flog.info("Timestep t=%s for %s is at or beyond maturity.  Setting last_computed_cash to notional %s plus accumulated coupon value %s.",
                  t, name, notional, accumulated_past_coupons, name="ragtop.instruments.optionality.bond")
        v = 0.0*(v+S)
      } else {
        v[v < 0] = 0
        last_computed_grid <<- as.vector(v)
      }
      v
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
              dividend_ceiling="numeric",
              last_computed_exercise_value="numeric",
              last_computed_exercise_decision="logical",
              last_used_S="numeric",
              last_used_t="numeric"),
  methods=list(
    compute_exercise_decision = function(v,S,t,discount_factor_fctn=discount_factor_fcn,...) {
      flog.info("Evaluating exercise decisions of convertible bond %s at t=%s", name, t,
                name="ragtop.instruments.exercise.convertible")
      # Because grid values for bonds represent existing bond plus all past
      #  coupons, exercise values must have those accrued coupons added for
      #  a fair comparison
      exercise_values = S * conversion_ratio
      if (t >= maturity) {
        if (all(v <= 0)) {  # Must be initialization of grid
          total_early_exercise_value = exercise_values
          v = notional + 0*S
          flog.info("Since instrument %s t=%s was at or beyond the maturity %s, initialized total_early_exercise_value with S * conversion_ratio",
                    name, t, maturity,
                    name="ragtop.instruments.exercise.convertible")
        }
      } else {
        accumulated_past_coupons = accumulate_coupon_values_before(t, discount_factor_fctn=discount_factor_fctn)
        total_early_exercise_value = exercise_values + accumulated_past_coupons
      }
      flog.debug("exercise_values %s total_early_exercise_value %s",
                 toString(exercise_values), toString(total_early_exercise_value))
      total_early_exercise_value[total_early_exercise_value < 0] = 0  # Should just be a safety measure
      flog.info("Differences between %s grid value and exercise value range from %s to %s, averaging %s",
                name, min(v - total_early_exercise_value), max(v - total_early_exercise_value),
                mean(v - total_early_exercise_value),
                name="ragtop.instruments.exercise.convertible")
      exercise_ix = as.vector(v < total_early_exercise_value)
      if (sum(exercise_ix)>0) {
        flog.info("%s at %s has %s of %s nodes that appear exercisable, for value (unadjusted by coupons) ranging from %s to %s",
                  name, t, sum(exercise_ix), length(S),
                  min(exercise_values[exercise_ix]), max(exercise_values[exercise_ix]),
                  name="ragtop.instruments.exercise.convertible")
      }
      if (sum(exercise_ix) < length(S)) {
        flog.info("%s at %s has %s of %s nodes that do not appear exercisable, for value (including past coupons) ranging from %s to %s",
                  name, t, length(S) - sum(exercise_ix), length(S),
                  min(v[!exercise_ix]), max(v[!exercise_ix]),
                  name="ragtop.instruments.exercise.convertible")
      }
      last_computed_exercise_decision <<- exercise_ix
      last_computed_exercise_value <<- total_early_exercise_value
      last_used_S <<- S
      last_used_t <<- t
      last_computed_exercise_value
    },
    exercise_decision = function(v,S,t,discount_factor_fctn=discount_factor_fcn,...) {
      "Find indexes where hold value {v} will be inferior to conversion value at each stock price level in {S}, adjusted to include all past coupons"
      # Memoize for efficiency
      if (is.blank(last_used_S) || anyNA(last_used_S) || is.blank(last_used_t) || any(S!=last_used_S) || t!=last_used_t) {
        compute_exercise_decision(v, S, t, discount_factor_fctn=discount_factor_fctn, ...)
      } else {
        flog.info("Reusing previously computed %s exercise decisions at t=%s",
                  name, t, name="ragtop.instruments.exercise.convertible")
      }
      list(exercise_indexes=last_computed_exercise_decision,
           exercise_values=last_computed_exercise_value)
    },
    update_cashflows = function(small_t, big_t, discount_factor_fctn=discount_factor_fcn, ...) {
      exercise_free_cashflows = callSuper(small_t, big_t, discount_factor_fctn=discount_factor_fctn, include_notional=FALSE, ...)
      if (is.blank(last_computed_grid) || is.blank(last_computed_exercise_decision) ) {
        cashflows = exercise_free_cashflows
        flog.info("Since our last computed %s exercise decision or grid was blank, setting cashflows from %s to %s the same as a straight bond, to %s",
                  name, small_t, big_t, exercise_free_cashflows,
                  name="ragtop.instruments.cashflows.convertible")
      } else {
        cashflows = 0 * last_computed_grid + exercise_free_cashflows
        flog.info("Exercise-free cashflows for %s were %s in time interval (%s,%s].  Zeroing them out for %s of %s grid entry cash flows due to conversion.",
                  name, mean(exercise_free_cashflows), small_t, big_t,
                  sum(last_computed_exercise_decision), length(last_computed_grid),
                  name="ragtop.instruments.cashflows.convertible")
        cashflows[last_computed_exercise_decision] = 0
      }
      cashflows
    },
    optionality_fcn = function(v,S,t,discount_factor_fctn=discount_factor_fcn,...) {
      "Return the greater of hold value {v} or conversion value at each stock price level in {S}, adjusted to include all past coupons"
      exer = exercise_decision(v,S,t,discount_factor_fctn=discount_factor_fctn,...)
      ev = as.vector(v)
      flog.info("Last computed exer values for %s have %s cases of early exercise",
                name, sum(exer$exercise_indexes), name="ragtop.instruments.exercise.convertible")
      ev[exer$exercise_indexes] = exer$exercise_values[exer$exercise_indexes]
      last_computed_grid <<- ev
      last_computed_grid
    },
    terminal_values = function(v,S,t,discount_factor_fctn=discount_factor_fcn,...) {
      accumulated_past_coupons = accumulate_coupon_values_before(t, discount_factor_fctn=discount_factor_fctn)
      total_bond_values = notional + 0 * S
      exercise_values = S * conversion_ratio
      total_early_exercise_value = exercise_values
      terminal = total_bond_values
      exer_ix = (terminal<total_early_exercise_value)
      some_exercisable = any(exer_ix)
      if (some_exercisable) {
        num_exercised = sum(exer_ix)
        cases_with_exer = total_early_exercise_value[exer_ix]
        flog.info("Terminal values for %s have %s of %s indicating exercise, exceeding notional of %s, ranging in value (not past coupons of %s) from %s to %s",
                  name, num_exercised, length(S), notional, accumulated_past_coupons,
                  min(cases_with_exer), max(cases_with_exer),
                  name="ragtop.instruments.exercise.convertible")
        terminal[exer_ix] = total_early_exercise_value[exer_ix]
      } else {
        num_exercised = 0
        flog.info("Terminal values for %s have no nodes indicating exercise, all values set to notional, with past coupons were %s + %s = %s",
                  name, notional, accumulated_past_coupons, notional + accumulated_past_coupons,
                  name="ragtop.instruments.exercise.convertible")
      }
      last_used_S <<- S
      last_used_t <<- t
      terminal
    }
  )
)
