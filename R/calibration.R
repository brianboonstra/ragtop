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

#' Implied volatility of european-exercise option under Black-Scholes or a jump-process extension
#'
#' Find default-free volatility (not necessarily just Black-Scholes) based on
#'  known interest rates and hazard rates, using a given option price.
#'
#' To get a straight Black-Scholes implied volatility, simply call this function with
#'  \code{const_default_intensity} set to zero (the default).
#'
#' @param option_price Present option value
#' @param callput 1 for calls, -1 for puts
#' @param S0 initial underlying price
#' @param K strike
#' @param r risk-free interest rate
#' @param time Time from \code{0} until expiration
#' @param const_default_intensity hazard rate of underlying default
#' @param divrate A continuous rate for dividends and other cashflows such as foreign interest rates
#' @param borrow_cost A continuous rate for stock borrow costs
#' @param dividends A \code{data.frame} with columns \code{time}, \code{fixed},
#'   and \code{proportional}.  Dividend size at the given \code{time} is
#'   then expected to be equal to \code{fixed + proportional * S / S0}.  Fixed
#'   dividends will be converted to proportional for purposes of this algorithm.  To handle
#'   truly fixed dividends, see \code{\link{implied_jump_process_volatility}}
#' @return A scalar volatility
#' @keywords Black-Scholes
#' @family Implied Volatilities
#' @family Equity Independent Default Intensity
#' @family European Options
#' @export implied_volatility
implied_volatility = function(option_price, callput, S0, K, r, time,
                              const_default_intensity=0, divrate=0, borrow_cost=0,
                              dividends=NULL,
                              relative_tolerance=1.e-6,
                              max.iter=100,
                              max_vola=4.00)
{
  min_vola = 1.e-12
  min_value = blackscholes(callput, S0, K, r, time, vola=min_vola,
                           default_intensity=const_default_intensity,
                           divrate=divrate, borrow_cost=borrow_cost,
                           dividends=dividends)$Price
  max_value = blackscholes(callput, S0, K, r, time, vola=max_vola,
                           default_intensity=const_default_intensity,
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
                         default_intensity=const_default_intensity,
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
                           default_intensity=const_default_intensity,
                           divrate=divrate, borrow_cost=borrow_cost,
                           dividends=dividends)
    if (bs_vals$Price > option_price) {
      max_vola = vola
    } else if (bs_vals$Price < option_price) {
      min_vola = vola
    }
    flog.debug("Defaultable volatility blackscholes() testing %s at err %s vega %s",
               vola, bs_vals$Price - option_price, bs_vals$Vega,
               name='ragtop.calibration.implied_volatility')
    iter = iter + 1
  }
  vola
}

#' Implied volatilities of european-exercise options under Black-Scholes or a jump-process extension
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
#' @family Implied Volatilities
#' @family European Options
#' @family Equity Independent Default Intensity
#' @export implied_volatilities
implied_volatilities = Vectorize(implied_volatility)

#' Find jump process volatility with a given default risk from a straight Black-Scholes volatility
#'
#' Find default-free volatility (i.e. volatility of a Wiener process with a
#'   companion jump process to default) based on known interest rates and hazard rates, using
#'   and at-the-money put option at the given tenor to set the standard price.
#'
#' @param bs_vola BlackScholes volatility of an option with no default assumption
#' @param time Time to expiration of associated option contracts
#' @param const_short_rate A constant to use for the instantaneous interest rate in case \code{discount_factor_fcn}
#'  is not given
#' @param const_default_intensity A constant to use for the instantaneous default intensity in case \code{survival_probability_fcn}
#'  is not given
#' @param relative_tolerance Relative tolerance in instrument price defining the root-finder halting condition
#' @param max.iter Maximum number of root-finder iterations allowed
#' @return A scalar volatility
#' @family Implied Volatilities
#' @family Equity Independent Default Intensity
#' @export equivalent_jump_vola_to_bs
equivalent_jump_vola_to_bs = function(bs_vola, time,
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
                                                        const_volatility=bs_vola,
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
               default_free_price, min_value, name='ragtop.calibration.equivalent_jump_vola_to_bs')
  }
  vola = 0.9*default_free_vola  # Reasonable starting guess without doing much math
  bs_vals = black_scholes_on_term_structures(-1, 1, 1, time,
                                                       const_volatility=vola,
                                                       discount_factor_fcn = discount_factor_fcn,
                                                       survival_probability_fcn = survival_probability_fcn,
                                                       dividends=dividends,
                                                       borrow_cost=borrow_cost,
                                                       dividend_rate=dividend_rate)
  flog.debug("Defaultable volatility black_scholes_on_term_structures() initially testing %s against %s at err %s vega %s",
             vola, default_free_price, bs_vals$Price - default_free_price, bs_vals$Vega,
             name='ragtop.calibration.equivalent_jump_vola_to_bs')
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
    flog.debug("Defaultable volatility testing black_scholes_on_term_structures(PUT,...) %s at err %s vega %s",
               vola, bs_vals$Price - default_free_price, bs_vals$Vega,
               name='ragtop.calibration.equivalent_jump_vola_to_bs')
    iter = iter + 1
  }
  vola
}

#' Find straight Black-Scholes volatility equivalent to jump process with a given default risk
#'
#' Find Black-Scholes volatility based on known
#'   interest rates and hazard rates, using
#'   an at-the-money put option at the given tenor to set the standard price.
#'
#' @param jump_process_vola Volatility of default-free process
#' @param time Time to expiration of associated option contracts
#' @param const_short_rate A constant to use for the instantaneous interest rate in case \code{discount_factor_fcn}
#'  is not given
#' @param const_default_intensity A constant to use for the instantaneous default intensity in case \code{survival_probability_fcn}
#'  is not given
#' @param relative_tolerance Relative tolerance in instrument price defining the root-finder halting condition
#' @param max.iter Maximum number of root-finder iterations allowed
#' @return A scalar defaultable volatility of an option
#' @family Implied Volatilities
#' @family Equity Independent Default Intensity
#' @export equivalent_bs_vola_to_jump
equivalent_bs_vola_to_jump = function(jump_process_vola, time,
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
                                                        const_volatility=jump_process_vola,
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
  flog.debug("Equiv BS  volatility black_scholes_on_term_structures() testing %s against %s at err %s vega %s",
             vola, default_free_price, bs_vals$Price - defaultable_price, bs_vals$Vega,
             name='ragtop.calibration.equivalent_bs_vola_to_jump')
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
    flog.debug("Equiv BS volatility testing black_scholes_on_term_structures(PUT,...) %s at err %s vega %s",
               vola, bs_vals$Price - default_free_price, bs_vals$Vega,
               name='ragtop.calibration.equivalent_bs_vola_to_jump')
    iter = iter + 1
  }
  vola
}



#' Find the implied volatility of a european-exercise option with term structures
#'
#' Use the Black-Scholes formula to generate European option values and run them
#'   through Newton's method until a constant volatility matching the provided
#'   option price has been found.
#'
#' Differs from \code{implied_volatility} by calling \code{black_scholes_on_term_structures} for
#'   pricing, thereby allowing term structures of rates nontrivial \code{survival_probability_fcn}
#'
#' @seealso \code{\link{implied_volatility}} for simpler cases with constant
#'   parameters, \code{\link{black_scholes_on_term_structures}} for the underlying
#'   pricing algorithm
#' @inheritParams form_present_value_grid
#' @param option_price Option price to match
#' @param callput 1 for calls, -1 for puts
#' @param S0 initial underlying price
#' @param K strike
#' @param time Time to expiration
#' @param ... Further arguments to be passed on to \code{black_scholes_on_term_structures}
#' @param starting_volatility_estimate The Newton method's original guess
#' @param survival_probability_fcn (Implied argument) A function for probability
#'   of survival, with arguments \code{T}, \code{t} and \code{T>t}.  E.g. with
#'   a constant volatility \eqn{s} this takes the form \eqn{(T-t)s^2}.
#' @param dividends A \code{data.frame} with columns \code{time}, \code{fixed},
#'   and \code{proportional}.  Dividend size at the given \code{time} is
#'   then expected to be equal to \code{fixed + proportional * S / S0}.  Fixed
#'   dividends will be converted to proportional for purposes of this algorithm.
#' @param relative_tolerance Relative tolerance in instrument price defining the root-finder halting condition
#' @param max.iter Maximum number of root-finder iterations allowed
#' @param max_vola Maximum volatility to try
#' @return Estimated volatility
#' @family Implied Volatilities
#' @family Equity Independent Default Intensity
#' @family European Options
#' @export implied_volatility_with_term_struct
implied_volatility_with_term_struct = function(option_price, callput, S0, K, time,
                                       ...,
                                       starting_volatility_estimate=0.5,
                                       relative_tolerance=1.e-6,
                                       max.iter=100,
                                       max_vola=4.00)
{
  compute_price = function(v) {
    black_scholes_on_term_structures(callput, S0, K, time,
                                     const_volatility=v, # Redundant, but left for clarity
                                     variance_cumulation_fcn = function(T, t){v^2*(T-t)},
                                     ...)
  }
  min_vola = 1.e-12
  min_value = compute_price(min_vola)$Price
  max_value = compute_price(max_vola)$Price
  if (option_price<min_value) {
    flog.warn("The provided option price %s is so low that no positive volatility can explain it.  Minimum values would be %s",
              option_price, min_value,
              name='ragtop.calibration.implied_volatility_with_term_struct')
    return(NA)
  } else {
    flog.debug("Option price price %s exceeds the minimum price %s so a solution exists",
               option_price, min_value,
               name='ragtop.calibration.implied_volatility_with_term_struct')
  }
  if (option_price>max_value) {
    flog.warn("The provided option price %s is so high that it exceeds the maximum volatility specified %s, at which the price is %s",
              option_price, max_vola, max_value,
              name='ragtop.calibration.implied_volatility_with_term_struct')
    return(NA)
  } else {
    flog.debug("Option price price %s less than the maximum price %s so a solution exists",
               option_price, max_value,
               name='ragtop.calibration.implied_volatility_with_term_struct')
  }
  vola = starting_volatility_estimate
  bs_vals = compute_price(vola)
  flog.debug("Impvol term struct initial test used black_scholes_on_term_structures(%s) against %s at err %s vega %s",
             vola, option_price, bs_vals$Price - option_price, bs_vals$Vega,
             name='ragtop.calibration.implied_volatility_with_term_struct')
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
                 name='ragtop.calibration.implied_volatility_with_term_struct')
      vola = 0.5 * (max_vola + min_vola)
    }
    bs_vals = compute_price(vola)
    if (bs_vals$Price > option_price) {
      max_vola = vola
    } else if (bs_vals$Price < option_price) {
      min_vola = vola
    }
    flog.debug("Impvol term struct testing black_scholes_on_term_structures(%s) at err %s vega %s",
               vola, bs_vals$Price - option_price, bs_vals$Vega,
               name='ragtop.calibration.implied_volatility_with_term_struct')
    iter = iter + 1
  }
  vola
}


#' Implied volatility of an american option with equity-independent term structures
#'
#' Use the grid solver to generate american option values under
#'   a default model with survival probabilities not linked to equity prices. and run them
#'   through a bisective root search method until a constant volatility matching the provided
#'   option price has been found.
#'
#' @seealso \code{\link{implied_volatility_with_term_struct}} for implied volatility
#'   of European options under the same conditions, \code{\link{american}} for the
#'   underlying pricing algorithm
#' @inheritParams form_present_value_grid
#' @inheritParams american
#' @param option_price Option price to match
#' @param ... Additional arguments to be passed on to \code{\link{implied_volatility_with_term_struct}}
#'   and \code{\link{american}}
#' @param survival_probability_fcn (Implied argument) A function for probability of survival, with
#'   arguments \code{T}, \code{t} and \code{T>t}.  E.g. with
#'   a constant volatility \eqn{s} this takes the form \eqn{(T-t)s^2}. Should be matched
#'   to \code{default_intensity_fcn}
#' @param default_intensity_fcn A function for computing default intensity
#'   occurring at a given time, dependent on time and stock price, with
#'   arguments \code{t}, \code{S}.  Should be matched to \code{survival_probability_fcn}
#' @param dividends A \code{data.frame} with columns \code{time}, \code{fixed},
#'   and \code{proportional}.  Dividend size at the given \code{time} is
#'   then expected to be equal to \code{fixed + proportional * S / S0}.  Fixed
#'   dividends will be converted to proportional for purposes of this algorithm.
#' @param relative_tolerance Relative tolerance in instrument price defining the root-finder halting condition
#' @param max.iter Maximum number of root-finder iterations allowed
#' @param max_vola Maximum volatility to try
#' @return Estimated volatility
#' @keywords calibration
#' @concept implied volatility
#' @family Implied Volatilities
#' @family Equity Independent Default Intensity
#' @family American Exercise Equity Options
#' @examples
#' american_implied_volatility(25,CALL,S0=100,K=100,time=2.2, const_short_rate=0.03)
#' df250 =  function(t) ( exp(-0.02*t)*exp(-0.03*max(0,t-1.0))) # Simple term structure
#' df25 = function(T,t){df250(T)/df250(t)} # Relative discount factors
#' american_implied_volatility(25,-1,100,100,2.2,discount_factor_fcn=df25)
#' @export american_implied_volatility
american_implied_volatility = function(option_price, callput, S0, K, time,
                                       const_default_intensity = 0,
                                       survival_probability_fcn = function(T, t, ...){exp(-const_default_intensity*(T-t))},
                                       default_intensity_fcn = function(t, S, ...){const_default_intensity+0.0*S},
                                       ...,
                                       num_time_steps=30,
                                       structure_constant=2.0,
                                       std_devs_width=5.0,
                                       relative_tolerance=1.e-4,
                                       max.iter=100,
                                       max_vola=4.00)
{
  compute_amer_cv_price = function(v) {
    vcum_f = function(T, t){v^2*(T-t)}
    px = american(callput=callput, S0=S0, K=K, time=time,
             survival_probability_fcn = survival_probability_fcn,
             default_intensity_fcn = default_intensity_fcn,
             const_volatility=v,
             variance_cumulation_fcn = vcum_f,
             ...,
             num_time_steps=num_time_steps,
             structure_constant=structure_constant,
             std_devs_width=std_devs_width)
    flog.debug("Implied vol got %s = american(callput=%s, S0=%s, K=%s, time=%s, const_volatility=%s, variance_cumulation_fcn=function(T, t){%s^2*(T-t)},...)",
               px, callput, S0, K, time, v, v,
               name='ragtop.calibration.american_implied_volatility')

    px
  }
  min_vola = 1.e-6
  min_value = compute_amer_cv_price(min_vola)
  max_value = compute_amer_cv_price(max_vola)
  if (option_price<min_value) {
    flog.warn("The provided option price %s is so low that no positive volatility can explain it.  Minimum values would be %s",
              option_price, min_value,
              name='ragtop.calibration.american_implied_volatility')
    return(NA)
  } else {
    flog.debug("Option price price %s exceeds the minimum price %s so a solution exists",
               option_price, min_value,
               name='ragtop.calibration.american_implied_volatility')
  }
  if (option_price>max_value) {
    flog.warn("The provided option price %s is so high that it exceeds the maximum volatility specified %s, at which the price is %s",
              option_price, max_vola, max_value,
              name='ragtop.calibration.american_implied_volatility')
    return(NA)
  } else {
    flog.debug("Option price price %s less than the maximum price %s so a solution exists",
               option_price, max_value,
               name='ragtop.calibration.american_implied_volatility')
  }
  # Start with a European volatility estimate.  We expect this to
  #   be higher than the American volatility due to the value of early exercise. This
  #   is particularly useful because it constrains us away from really testing
  #   high volatilities where computaiton expenses and errors get large (though
  #   we did test max_vol once up above).
  starting_volatility_estimate = implied_volatility_with_term_struct(
    option_price, callput, S0, K, time,
    survival_probability_fcn=survival_probability_fcn,
    ...,
    relative_tolerance=relative_tolerance,
    max.iter=max.iter,
    max_vola=max_vola)
  vola = starting_volatility_estimate
  grid_val = compute_amer_cv_price(vola)
  flog.debug("Amer impvol initial test used american(%s) against %s at err %s",
             vola, option_price, grid_val - option_price,
             name='ragtop.calibration.american_implied_volatility')
  if (grid_val > option_price) {
    max_vola = vola
  } else if (grid_val < option_price) {
    flog.warn("Unexpected case of European option appearing to have lower implied vol than American")
    min_vola = vola
  }
  # Run the bisection
  iter = 0
  while ((abs(grid_val/option_price-1.0)>relative_tolerance) && (iter<max.iter)) {
    vola = 0.5 * (max_vola + min_vola)
    grid_val = compute_amer_cv_price(vola)
    if (grid_val > option_price) {
      max_vola = vola
    } else if (grid_val < option_price) {
      min_vola = vola
    }
    flog.debug("Amer impvol test used american(%s) at err %s ",
               vola, grid_val - option_price,
               name='ragtop.calibration.american_implied_volatility')
    iter = iter + 1
  }
  vola
}


#' Implied volatility of any instrument
#'
#' Use the grid solver to generate instrument prices via \code{find_present_value} and run them
#'   through a bisective root search method until a constant volatility matching the provided
#'   instrument price has been found.
#'
#' Unlike \code{american_implied_volatility}, this routine allows for any legal
#'  term structures and equity-linked default intensities.  For that reason, it eschews
#'  the control variate tricks that make \code{american_implied_volatility} so much faster.
#'
#' Note that equity-linked default intensities can result in instrument prices that
#'   are not monotonic in volatility.  This bisective root finder will find a solution
#'   but not necessarily any particular one.
#'
#' @seealso \code{\link{find_present_value}} for the underlying
#'   pricing algorithm, \code{\link{implied_volatility_with_term_struct}} for European options
#'   without equity dependence of default intensity, \code{\link{american_implied_volatility}} for the same on American options
#' @inheritParams find_present_value
#' @inheritParams american
#' @param instrument_price Target price for root finder
#' @param starting_volatility_estimate Bisection method original guess
#' @param ... Additional arguments to be passed on to \code{\link{find_present_value}}
#' @param relative_tolerance Relative tolerance in instrument price defining the root-finder halting condition
#' @param max.iter Maximum number of root-finder iterations allowed
#' @param max_vola Maximum volatility to try
#' @return A list of present values, with the same names as \code{instruments}
#' @keywords calibration
#' @concept implied volatility
#' @family Implied Volatilities
#' @family Equity Dependent Default Intensity
#' @examples
#' implied_jump_process_volatility(25, AmericanOption(maturity=1.1, strike=100, callput=-1), S0=100, num_time_steps=50, relative_tolerance=1.e-3)
#'
#' @export implied_jump_process_volatility
implied_jump_process_volatility = function(instrument_price, instrument, variance_cumulation_fcn=NULL,
                                         ...,
                                         starting_volatility_estimate=0.85,
                                         relative_tolerance=5.e-3,
                                         max.iter=100,
                                         max_vola=4.00)
{
  instrument_list = list(impl_vol_tgt=instrument)
  compute_price = function(v) {
    find_present_value(instruments=instrument_list,
                       variance_cumulation_fcn = function(T, t){v^2*(T-t)},
                       ...)[[1]]
  }
  min_vola = 1.e-12
  min_value = compute_price(min_vola)
  max_value = compute_price(max_vola)
  if (instrument_price<min_value) {
    flog.warn("The provided instrument price %s is so low that no positive volatility can explain it.  Minimum values would be %s",
              instrument_price, min_value,
              name='ragtop.calibration.implied_jump_process_volatility')
    return(NA)
  } else {
    flog.debug("Instrument price price %s exceeds the minimum price %s so a solution exists",
               instrument_price, min_value,
               name='ragtop.calibration.implied_jump_process_volatility')
  }
  if (instrument_price>max_value) {
    flog.warn("The provided option price %s is so high that it exceeds the maximum volatility specified %s, at which the price is %s",
              instrument_price, max_vola, max_value,
              name='ragtop.calibration.implied_jump_process_volatility')
    return(NA)
  } else {
    flog.debug("Option price price %s less than the maximum price %s so a solution exists",
               instrument_price, max_value,
               name='ragtop.calibration.implied_jump_process_volatility')
  }
  vola = starting_volatility_estimate
  grid_val = compute_price(vola)
  flog.debug("Defaultable volatility initial test used %s against %s at err %s",
             vola, instrument_price, grid_val - instrument_price,
             name='ragtop.calibration.implied_jump_process_volatility')
  if (grid_val > instrument_price) {
    max_vola = vola
  } else if (grid_val < instrument_price) {
    flog.warn("Unexpected case of European option appearing to have lower implied vol than American")
    min_vola = vola
  }
  # Run the bisection
  iter = 0
  while ((abs(grid_val/instrument_price-1.0)>relative_tolerance) && (iter<max.iter)) {
    vola = 0.5 * (max_vola + min_vola)
    grid_val = compute_price(vola)
    if (grid_val > instrument_price) {
      max_vola = vola
    } else if (grid_val < instrument_price) {
      min_vola = vola
    }
    flog.debug("Defaultable volatility testing %s at err %s ",
               vola, grid_val - instrument_price,
               name='ragtop.calibration.implied_jump_process_volatility')
    iter = iter + 1
  }
  ## TODO: Return min_vola or max_vola if one of them was closer
  vola
}

#' Fit piecewise constant volatilities to a set of equity options
#'
#' Given a set of equity options with increasing tenors, along with
#' target prices for those options, and a set of equity-lined default
#' SDE parameters, fit a vector of piecewise constant
#' volatilities and an associated cumulative variance function
#' to them.
#'
#' By default, the fitting happens in implied Black-Scholes volatility
#' space for better normalization.  That is to say, the fitting does pricing
#' using the \emph{full} SDE and PDE solver via \code{\link{find_present_value}}, but
#' judges fit quality on the basis of running resulting prices through a
#' nonlinear transformation that just
#' happens to come from the straight Black-Scholes model.
#'
#' @inheritParams find_present_value
#' @param eq_options A list of options to find prices for.  Each must have fields \code{callput},
#'                   \code{maturity}, and \code{strike}.  This list must be in strictly increasing order of maturity.
#' @param mid_prices Prices to match
#' @param relative_spread_tolerance Tolerance multiplier on bid-ask spreads taken from vol normalization
#' @param initial_vols_guess Initial set of volatilities to try in the root finder
#' @param spreads Spreads within which any match is tolerable
#' @param survival_probability_fcn A function for probability of survival, with
#'   arguments \code{T}, \code{t} and \code{T>t}.  E.g. with
#'   a constant volatility \eqn{s} this takes the form \eqn{(T-t)s^2}.  This argument is
#'   only used in normalization of prices to vols for root finder tolerance, and
#'   is therefore entirely optional
#' @param force_same_grid Price all options on the same grid, rather than having smaller timestep sizes for earlier maturities
#' @param use_impvol Judge fit quality on implied vol distance rather than price distance
#' @param ... Futher arguments to \code{\link{find_present_value}}
#' @return A list with two elements, \code{volatilities} and \code{cumulation_function}.  The \code{cumulation_function} will
#'   be a 2-parameter function giving cumulated variances, as created by code{\link{variance_cumulation_from_vols}}
#' @keywords calibration
#' @keywords Black-Scholes
#' @concept implied volatility
#' @family Implied Volatilities
#' @family Equity Dependent Default Intensity
#'
#' @export fit_variance_cumulation
fit_variance_cumulation = function(S0, eq_options, mid_prices, spreads=NULL,
                                   initial_vols_guess=0.55 + 0*mid_prices,
                                   use_impvol=TRUE,
                                   relative_spread_tolerance=0.01,
                                   force_same_grid=FALSE,
                                   num_time_steps=100,
                                   const_short_rate=0, const_default_intensity=0,
                                   discount_factor_fcn = function(T, t, ...){exp(-const_short_rate*(T-t))},
                                   survival_probability_fcn = function(T, t, ...){exp(-const_default_intensity*(T-t))},
                                   default_intensity_fcn = function(t, S, ...){const_default_intensity+0.0*S},
                                   dividends=NULL,
                                   borrow_cost=0.0,
                                   dividend_rate=0.0,
                                   ...)
{
  N = length(eq_options)
  if (N != length(mid_prices) || (0==N)) {
    stop("Number of prices to match must agree with number of options (",N,") in the calibration of variance")
  } else {
    flog.info("fit_variance_cumulation on %s options...", N,
              name='ragtop.calibration.fit_variance_cumulation')
  }
  maturities = as.numeric(lapply(eq_options, function(x){x$maturity}))
  compute_var_cum_f = function(vs) {
    variance_cumulation_from_vols(data.frame(volatility=as.numeric(vs), time=maturities))
  }
  compute_bsimpvol = function(eq_opt, tgt) {
    iv = implied_volatility_with_term_struct(tgt, callput=eq_opt$callput,
                                        S0=S0, K=eq_opt$strike, time=eq_opt$maturity,
                                        discount_factor_fcn=discount_factor_fcn,
                                        survival_probability_fcn=survival_probability_fcn,
                                        dividends=dividends, borrow_cost=borrow_cost,
                                        dividend_rate=dividend_rate)
    flog.debug("BS impvol of %s maturity %s at price %s is %s",
               eq_opt$name, eq_opt$maturity, tgt, iv,
               name='ragtop.calibration.fit_variance_cumulation')
    iv
  }
  override_Tmax = NA
  solver_tolerance = max(0.005, 0.005 * abs(mid_prices))
  vols = as.numeric(initial_vols_guess)
  if (force_same_grid) {
    override_Tmax = eq_options[[length(eq_options)]]$maturity
  }
  impvols = NA * mid_prices
  if (use_impvol) {
    for (i in 1:N) {
      eq_opt = eq_options[[i]]
      impvols[i] = compute_bsimpvol(eq_opt, mid_prices[[i]])
      flog.info("Normalized t=%s price %s to BS impvol %s",
                eq_opt$maturity, mid_prices[[i]], impvols[i],
                name='ragtop.calibration.fit_variance_cumulation')
    }
    if (is.blank(spreads)) {
      solver_tolerance = pmax(0.0001, 0.03 * abs(impvols))
      flog.info("No spreads to set tolerances in impvol terms, just chose values ranging from %s to %s",
                min(solver_tolerance), max(solver_tolerance),
                name='ragtop.calibration.fit_variance_cumulation')
    } else {
      for (i in 1:N) {
        bidvol = NA
        askvol = NA
        bidvol = try(compute_bsimpvol(eq_opt, mid_prices[[i]] - 0.5*spreads[[i]]), silent=TRUE)
        askvol = try(compute_bsimpvol(eq_opt, mid_prices[[i]] + 0.5*spreads[[i]]), silent=TRUE)
        if (is.blank(bidvol) || is.blank(askvol)) {
          flog.info("Either bid vol %s or ask vol %s was not computable from %s and %s",
                    bidvol, askvol, mid_prices[[i]] - 0.5*spreads[[i]], mid_prices[[i]] + 0.5*spreads[[i]],
                    name='ragtop.calibration.fit_variance_cumulation')
        } else {
          solver_tolerance[i] = relative_spread_tolerance * (askvol-bidvol)
        }
      }
      flog.info("Using spreads to set tolerances in impvol terms, ranging from %s to %s",
                min(solver_tolerance), max(solver_tolerance),
                name='ragtop.calibration.fit_variance_cumulation')
    }
  } else if (!(is.blank(spreads))) {
    solver_tolerance = pmax(0.005 * abs(mid_prices), spreads)
    flog.info("Using spreads to set tolerances, ranging from %s to %s",
              min(solver_tolerance), max(solver_tolerance),
              name='ragtop.calibration.fit_variance_cumulation')
  }
  last_cumul_variance = 0
  for (i in 1:N) {
    eq_opt = eq_options[[i]]
    flog.info("Using instrument %s to fit vol at t=%s",
              eq_opt$name, eq_opt$maturity,
              name='ragtop.calibration.fit_variance_cumulation')
    distfunc = function(v) {
      vols[i] = v
      cumul_func = compute_var_cum_f(vols)
      flog.debug("Term struct solver for instrument %s will test vola %s giving cumulative variance %s",
                i, v, cumul_func(eq_opt$maturity),
                name='ragtop.calibration.fit_variance_cumulation')
      computed_price = find_present_value(S0=S0, num_time_steps=num_time_steps,
                                          override_Tmax=override_Tmax,
                                          instruments=list(tsopt=eq_opt),
                                          variance_cumulation_fcn=cumul_func,
                                          ...)$tsopt
      if (use_impvol) {
        civ = compute_bsimpvol(eq_opt, computed_price)
        distance = civ - impvols[[i]]
        flog.info("Term struct solver for instrument %s tested vola %s and found price %s for an impvol of %s which is distance %s from %s",
                  i, v, computed_price, civ, distance, impvols[[i]],
                  name='ragtop.calibration.fit_variance_cumulation')
      } else {
        distance = computed_price - mid_prices[[i]]
        flog.info("Term struct solver for instrument %s tested vola %s and found price %s which is distance %s from %s",
                  i, v, computed_price, civ, distance, mid_prices[[i]],
                  name='ragtop.calibration.fit_variance_cumulation')
      }
      distance
    }
    brent_tol = solver_tolerance[[i]]
    done = FALSE
    # uniroot interval extension fails.  Extend ourselves
    test_vol = max(sqrt(last_cumul_variance/eq_opt$maturity)*1.02, initial_vols_guess[[i]])
    first_distance = distfunc(test_vol)
    if (first_distance>0) {
      max_vol = test_vol
      f.upper = first_distance
      min_vol = max_vol
      iter = 0
      next_distance = first_distance
      f.lower = next_distance
      while (iter<20 && (next_distance*first_distance>0) && abs(next_distance) > brent_tol) {
        iter = iter + 1
        min_vol = max(min_vol/1.1, sqrt(last_cumul_variance/eq_opt$maturity)*1.02)
        next_distance = distfunc(min_vol)
        f.lower = next_distance
      }
      if (abs(next_distance)<=brent_tol) {
        found_v = min_vol
        done = TRUE
      }
      if (iter==20) {
        stop("Could not fit term structure for t=",eq_opt$maturity,", number ",i,". Could not find min")
      }
    } else {
      min_vol = test_vol
      f.lower = first_distance
      max_vol = min_vol
      iter = 0
      next_distance = first_distance
      f.upper = next_distance
      while ((iter<20) && next_distance*first_distance>0 && (abs(next_distance)> brent_tol)) {
        iter = iter + 1
        max_vol = max_vol * 1.1
        next_distance = distfunc(max_vol)
        f.upper = next_distance
      }
      if (abs(next_distance)<=brent_tol) {
        found_v = max_vol
        done = TRUE
      }
      if (iter==20) {
        stop("Could not fit term structure for t=",eq_opt$maturity,", number ",i,". Could not find max")
      }
    }
    if (!done) {
      # Run Brent's method
      flog.info("Min possible volatility allowed for %s is %s, max is %s.  Will search with tolerance %s",
                 eq_opt$maturity, min_vol, max_vol, brent_tol,
                 name='ragtop.calibration.fit_variance_cumulation')
      orig_search_interval = c(min_vol, max_vol)
      solv = uniroot(distfunc, interval=orig_search_interval, extendInt="no",
                     f.lower=f.lower, f.upper=f.upper,
                     tol=brent_tol, maxiter=20)
      found_v = unlist(solv['root'])[[1]]
      done = TRUE
    }
    vols[i] = found_v
    flog.info("Term struct solver for instrument %s concluded with vola %s",
              i, vols[i],
              name='ragtop.calibration.fit_variance_cumulation')
    last_cumul_variance = vols[i]^2 * eq_opt$maturity
  }
  cfunc = compute_var_cum_f(vols)
  list(volatilities=vols, cumulation_function=cfunc)
}
