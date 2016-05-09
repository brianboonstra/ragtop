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
library(stats)
library(futile.logger)

#' Implied volatility
#'
#' Find default-free volatility based on known interest rates and hazard rates, using
#'   a given option price.
#'
#' @param option_price Present option value
#' @param callput 1 for calls, -1 for puts
#' @param S0 initial underlying price
#' @param K strike
#' @param r risk-free interest rate
#' @param time Time from \code{0} until expiration
#' @param default_intensity hazard rate of underlying default
#' @param divrate A continuous rate for dividends and other cashflows such as foreign interest rates
#' @param borrow_cost A continuous rate for stock borrow costs
#' @param dividends A \code{data.frame} with columns \code{time}, \code{fixed},
#'   and \code{proportional}.  Dividend size at the given \code{time} is
#'   then expected to be equal to \code{fixed + proportional * S / S0}.  Fixed
#'   dividends will be converted to proprtional for purposes of this algorithm.
#' @return A scalar volatility
#' @export implied_volatility
implied_volatility = function(option_price, callput, S0, K, r, time,
                              default_intensity=0, divrate=0, borrow_cost=0,
                              dividends=NULL,
                                            relative_tolerance=1.e-6,
                                            max.iter=100,
                              max_vola=4.00)
{
  min_vola = 1.e-12
  min_value = blackscholes(callput, S0, K, r, time, vola=min_vola,
                           default_intensity=default_intensity,
                           divrate=divrate, borrow_cost=borrow_cost,
                           dividends=dividends)$Price
  max_value = blackscholes(callput, S0, K, r, time, vola=max_vola,
                           default_intensity=default_intensity,
                           divrate=divrate, borrow_cost=borrow_cost,
                           dividends=dividends)$Price
  if (option_price<min_value) {
    flog.warn("The provided option price %s is so low that no positive volatility can explain it.  Minimum values would be %s",
              option_price, min_value,
              name='ragtop.calibration.implied_volatility')
    return(NA)
  } else {
    flog.debug("Option price price %s exceeds the minimum price %s so a solution exists",
               option_price, min_value,
               name='ragtop.calibration.implied_volatility')
  }
  if (option_price>max_value) {
    flog.warn("The provided option price %s is so high that it exceeds the maximum volatility specified %s, at which the price is %s",
              option_price, max_vola, max_value,
              name='ragtop.calibration.implied_volatility')
    return(NA)
  } else {
    flog.debug("Option price price %s less than the maximum price %s so a solution exists",
               option_price, max_value,
               name='ragtop.calibration.implied_volatility')
  }
  vola = 0.5
  bs_vals = blackscholes(callput, S0, K, r, time, vola=vola,
                         default_intensity=default_intensity,
                         divrate=divrate, borrow_cost=borrow_cost,
                         dividends=dividends)
  flog.debug("Defaultable volatility initial test used %s against %s at err %s vega %s",
             vola, option_price, bs_vals$Price - option_price, bs_vals$Vega,
             name='ragtop.calibration.implied_volatility')
  if (bs_vals$Price > option_price) {
    max_vola = vola
  } else if (bs_vals$Price < option_price) {
    min_vola = vola
  }
  iter = 0
  while ((abs(bs_vals$Price/option_price-1.0)>relative_tolerance) && (iter<max.iter)) {
    vola_adj = -(bs_vals$Price-option_price)/bs_vals$Vega  # Newton's method
    vola = vola + vola_adj
    if (vola>max_vola || vola<min_vola) {
      flog.debug("The implied volatility solver needed to take a bisection step",
                name='ragtop.calibration.implied_volatility')
      vola = 0.5 * (max_vola + min_vola)
    }
    bs_vals = blackscholes(callput, S0, K, r, time, vola=vola,
                           default_intensity=default_intensity,
                           divrate=divrate, borrow_cost=borrow_cost,
                           dividends=dividends)
    if (bs_vals$Price > option_price) {
      max_vola = vola
    } else if (bs_vals$Price < option_price) {
      min_vola = vola
    }
    flog.debug("Defaultable volatility testing %s at err %s vega %s",
               vola, bs_vals$Price - option_price, bs_vals$Vega,
               name='ragtop.calibration.implied_volatility')
    iter = iter + 1
  }
  vola
}

#' Implied volatilities
#'
#' Find default-free volatilities based on known interest rates and hazard rates, using
#'   a given option price.
#'
#' @param option_price Present option values
#' @param callput 1 for calls, -1 for puts
#' @param S0 initial underlying price
#' @param K strike
#' @param r risk-free interest rate
#' @param time Time from \code{0} until expiration
#' @param default_intensity hazard rate of underlying default
#' @param divrate A continuous rate for dividends and other cashflows such as foreign interest rates
#' @param borrow_cost A continuous rate for stock borrow costs
#' @param dividends A \code{data.frame} with columns \code{time}, \code{fixed},
#'   and \code{proportional}.  Dividend size at the given \code{time} is
#'   then expected to be equal to \code{fixed + proportional * S / S0}.  Fixed
#'   dividends will be converted to proprtional for purposes of this algorithm.
#' @return Scalar volatilities
#' @export implied_volatilities
implied_volatilities = Vectorize(implied_volatility)

#' Find default-free volatility
#'
#' Find default-free volatility based on known interest rates and hazard rates, using
#'   and at-the-money put option at the given tenor to set the standard price.
#'
#' @param default_free_vola BlackScholes volatility of an option with no default assumption
#' @param time Time to expiration of associated option contracts
#' @param const_short_rate A constant to use for the instantaneous interest rate in case \code{discount_factor_fcn}
#'  is not given
#' @param const_default_intensity A constant to use for the instantaneous default intensity in case \code{survival_probability_fcn}
#'  is not given
#' @return A scalar volatility
#' @export calibrate_defaultable_volatility
calibrate_defaultable_volatility = function(default_free_vola, time,
                                          const_short_rate=0, const_default_intensity=0,
                                          discount_factor_fcn = function(T, t, ...){exp(-const_short_rate*(T-t))},
                                          survival_probability_fcn = function(T, t, ...){exp(-const_default_intensity*(T-t))},
                                          dividends=NULL,
                                          borrow_cost=0.0,
                                          dividend_rate=0.0,
                                          relative_tolerance=1.e-6,
                                          max.iter=100)
{
  default_free_price = black_scholes_on_term_structures(-1, 1, 1, time,
                                                        const_volatility=default_free_vola,
                                                        discount_factor_fcn = discount_factor_fcn,
                                                        survival_probability_fcn = function(T,t){1},
                                                        dividends=dividends,
                                                        borrow_cost=borrow_cost,
                                                        dividend_rate=dividend_rate)$Price
  min_value = black_scholes_on_term_structures(-1, 1, 1, time,
                                                        const_volatility=1.e-10,
                                                        discount_factor_fcn = discount_factor_fcn,
                                                        survival_probability_fcn = survival_probability_fcn,
                                                        dividends=dividends,
                                                        borrow_cost=borrow_cost,
                                                        dividend_rate=dividend_rate)$Price
  if (default_free_price<=min_value) {
    flog.warn("The provided survival probability is so low that it alone implies more than the provided default-free volatility")
    return(NA)
  } else {
    flog.debug("Default free price %s exceeds the minimum price %s so a solution exists",
               default_free_price, min_value)
  }
  vola = 0.9*default_free_vola  # Reasonable starting guess without doing much math
  bs_vals = black_scholes_on_term_structures(-1, 1, 1, time,
                                                       const_volatility=vola,
                                                       discount_factor_fcn = discount_factor_fcn,
                                                       survival_probability_fcn = survival_probability_fcn,
                                                       dividends=dividends,
                                                       borrow_cost=borrow_cost,
                                                       dividend_rate=dividend_rate)
  flog.debug("Defaultable volatility testing %s against %s at err %s vega %s",
             vola, default_free_price, bs_vals$Price - default_free_price, bs_vals$Vega,
             name='ragtop.calibration.calibrate_defaultable_volatility')
  iter = 0
  while ((abs(bs_vals$Price/default_free_price-1.0)>relative_tolerance) && (iter<max.iter)) {
    vola_adj = -(bs_vals$Price-default_free_price)/bs_vals$Vega  # Newton's method
    vola = vola + vola_adj
    if (vola<0) {
      flog.warn("The provided survival probability is so low that we cannot find the provided default-free volatility")
      return(NA)
    }
    bs_vals = black_scholes_on_term_structures(-1, 1, 1, time,
                                               const_volatility=vola,
                                               discount_factor_fcn = discount_factor_fcn,
                                               survival_probability_fcn = survival_probability_fcn,
                                               dividends=dividends,
                                               borrow_cost=borrow_cost,
                                               dividend_rate=dividend_rate)
    flog.debug("Defaultable volatility testing %s at err %s vega %s",
               vola, bs_vals$Price - default_free_price, bs_vals$Vega,
               name='ragtop.calibration.calibrate_defaultable_volatility')
    iter = iter + 1
  }
  vola
}

#' Find defaultable volatility
#'
#' Find defaultable volatility based on known interest rates and hazard rates, using
#'   an at-the-money put option at the given tenor to set the standard price.
#'
#' @param defaultable_volatility Volatility of default-free process
#' @param time Time to expiration of associated option contracts
#' @param const_short_rate A constant to use for the instantaneous interest rate in case \code{discount_factor_fcn}
#'  is not given
#' @param const_default_intensity A constant to use for the instantaneous default intensity in case \code{survival_probability_fcn}
#'  is not given
#' @return A scalar defaultable volatility of an option
#' @export equivalent_bs_vola_to_defaultable
equivalent_bs_vola_to_defaultable = function(defaultable_volatility, time,
                                            const_short_rate=0, const_default_intensity=0,
                                            discount_factor_fcn = function(T, t, ...){exp(-const_short_rate*(T-t))},
                                            survival_probability_fcn = function(T, t, ...){exp(-const_default_intensity*(T-t))},
                                            dividends=NULL,
                                            borrow_cost=0.0,
                                            dividend_rate=0.0,
                                            relative_tolerance=1.e-6,
                                            max.iter=100)
{
  defaultable_price = black_scholes_on_term_structures(-1, 1, 1, time,
                                                        const_volatility=defaultable_volatility,
                                                        discount_factor_fcn = discount_factor_fcn,
                                                        survival_probability_fcn = survival_probability_fcn,
                                                        dividends=dividends,
                                                        borrow_cost=borrow_cost,
                                                        dividend_rate=dividend_rate)$Price
  vola = defaultable_volatility  # Reasonable starting guess without doing much math
  bs_vals = black_scholes_on_term_structures(-1, 1, 1, time,
                                             const_volatility=vola,
                                             discount_factor_fcn = discount_factor_fcn,
                                             survival_probability_fcn = function(T,t){1},
                                             dividends=dividends,
                                             borrow_cost=borrow_cost,
                                             dividend_rate=dividend_rate)
  flog.debug("Default free volatility testing %s against %s at err %s vega %s",
             vola, default_free_price, bs_vals$Price - defaultable_price, bs_vals$Vega,
             name='ragtop.calibration.calibrate_defaultable_volatility')
  iter = 0
  while ((abs(bs_vals$Price/defaultable_price-1.0)>relative_tolerance) && (iter<max.iter)) {
    vola_adj = -(bs_vals$Price-defaultable_price)/bs_vals$Vega  # Newton's method
    vola = vola + vola_adj
    bs_vals = black_scholes_on_term_structures(-1, 1, 1, time,
                                               const_volatility=vola,
                                               discount_factor_fcn = discount_factor_fcn,
                                               survival_probability_fcn = function(T,t){1},
                                               dividends=dividends,
                                               borrow_cost=borrow_cost,
                                               dividend_rate=dividend_rate)
    flog.debug("Defaultable volatility testing %s at err %s vega %s",
               vola, bs_vals$Price - default_free_price, bs_vals$Vega,
               name='ragtop.calibration.calibrate_defaultable_volatility')
    iter = iter + 1
  }
  vola
}
