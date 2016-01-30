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
require(futile.logger)

# TODO: Allow initialization of some instruments from QuantLib classes

## Our basic instrument defines a tenor/maturity, a method
##  to provide values in case of default, and a method
##  to correct instrument prices in light of exercise
##  decisions.
GridPricedInstrument = setRefClass(
  "GridPricedInstrument",
  fields = list(maturity = "numeric",
                last_computed_grid = "numeric"),
  methods = list(
    recovery_fcn = function(v,S,t,...) {
      0.0
    },
    optionality_fcn = function(v,...) {
      last_computed_grid <<- v
      v
    }
  )
)

EquityOption = setRefClass(
  "EquityOption",
  contains = "GridPricedInstrument",
  fields = list(strike = "numeric",
                callput = "numeric")
)

AmericanOption = setRefClass(
  "AmericanOption",
  contains = "EquityOption",
  methods = list(
    optionality_fcn = function(v,S,t,...) {
      if (t < maturity) {
        exer_value = callput * (S - strike)
        exer_value[exer_value < 0] = 0
        exercise_ix = (v < exer_value)
        v[exercise_ix] = exer_value[exercise_ix]
      } else {
        v = callput * (S - strike)
      }
      v[v < 0] = 0
      last_computed_grid <<- v
      v
    }
  )
)

ZeroCouponBond = setRefClass(
  "ZeroCouponBond",
  contains="GridPricedInstrument",
  fields=list(notional="numeric",
              recovery_rate="numeric",
              discount_factor_fcn="function"),
  methods = list(
    optionality_fcn = function(v,S,t,...) {
      if (t >= maturity) {
        v = 0.0*(v+S) + notional
      }
      v[v < 0] = 0
      last_computed_grid <<- v
      v
    }
  )
)

CouponBond = setRefClass(
  "CouponBond",
  contains="ZeroCouponBond",
  fields=list(coupons="data.frame"),
  methods = list(
    accumulate_coupon_values_before = function(t,discount_factor_fcn=discount_factor_fcn) {
      paid_coupon_ix = (coupons$payment_time<=t)
      paid_coupons = coupons[paid_coupon_ix,]
      accumulation_proportions = discount_factor_fcn(t,paid_coupons$payment_time)
      sum(accumulation_proportions*paid_coupons$payment_size)
    }
  )
)

CallableBond = setRefClass(
  "CallableBond",
  contains="CouponBond",
  fields=list(calls="data.frame")
)

ConvertibleBond = setRefClass(
  "ConvertibleBond",
  contains="CallableBond",
  fields=list(conversion_ratio="numeric",
              dividend_ceiling="numeric"),
  methods=list(
    optionality_fcn = function(v,S,t,discount_factor_fcn=discount_factor_fcn,...) {
      flog.debug("Optionality at t=%s", t)
      exercise_values = S * conversion_ratio
      accumulated_past_coupons = accumulate_coupon_values_before(t,discount_factor_fcn=discount_factor_fcn)
      total_early_exercise_value = exercise_values + accumulated_past_coupons
      flog.debug("exercise_values %s total_early_exercise_value %s",
                toString(exercise_values), toString(total_early_exercise_value))
      total_early_exercise_value[total_early_exercise_value < 0] = 0
      exercise_ix = (v < total_early_exercise_value)
      v[exercise_ix] = total_early_exercise_value[exercise_ix]
      last_computed_grid <<- v
      v
    }
  )
)
