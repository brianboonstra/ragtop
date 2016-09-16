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

#' Form instrument objects for vanilla options
#'
#' Form a list twice as long as the longest of the arguments
#' \code{callput}, \code{K}, \code{time} whose
#' first half consists of AmericanOption objects and
#' second half consists of EuropeanOption objects
#' having the same exercise specification
#' @param callput 1 for calls, -1 for puts
#' @param K strike
#' @param time Time from \code{0} until expiration
#' @family American Exercise Equity Options
control_variate_pairs = function(callput, K, time)
{
  # Vectorize
  callput = 0 * K + 0 * time + 1 * callput
  K = 1 * K + 0 * time + 0 * callput
  time = 0 * K + 1 * time + 0 * callput
  instr_details = data.frame(callput=callput, K=K, time=time)
  make_amer_instrument = function(row) {
    AmericanOption(maturity=row['time'], strike=row['K'], callput=row['callput'],
                        name=paste0('A',row['K'],'_',as.integer(365*row['time']),'_',row['callput']+1))
  }
  make_euro_instrument = function(row) {
    EuropeanOption(maturity=row['time'], strike=row['K'], callput=row['callput'],
                        name=paste0('E',row['K'],'_',as.integer(365*row['time']),'_',row['callput']+1))
  }
  amer_instruments = apply(instr_details, 1, make_amer_instrument)
  euro_instruments = apply(instr_details, 1, make_euro_instrument)
  as.list(c(amer_instruments, euro_instruments))
}


#' Price one or more american-exercise options
#'
#' Use a control-variate scheme to simultaneously estimate the present
#'  values of a collection of one or more American-exercise options under
#'  a default model with survival probabilities not linked to equity prices.
#'
#'  The scheme
#'  uses find_present_value() to price the options and their European-exercise
#'  equivalents.  It then compares the latter to black-scholes formula output
#'  and uses the results as an error correction on the prices of the
#'  American-exercise options.
#'
#' @inheritParams form_present_value_grid
#' @inheritParams construct_implicit_grid_structure
#' @param callput 1 for calls, -1 for puts (may be a vector of the same)
#' @param S0 initial underlying price
#' @param K strike (may be a vector)
#' @param time Time from \code{0} until expiration (may be a vector)
#' @param num_time_steps Number of steps to use in the grid solver.  Can usually be
#'   set quite low due to the control variate scheme.
#' @param discount_factor_fcn A function for computing present values to
#'   time \code{t} of various cashflows occurring during this timestep, with
#'   arguments \code{T}, \code{t}
#' @param survival_probability_fcn (Implied argument) A function for probability of survival, with
#'   arguments \code{T}, \code{t} and \code{T>t}.  E.g. with
#'   a constant volatility \eqn{s} this takes the form \eqn{(T-t)s^2}. Should be matched to \code{default_intensity_fcn}
#' @param default_intensity_fcn A function for computing default intensity
#'   occurring at a given time, dependent on time and stock price, with
#'   arguments \code{t}, \code{S}.  Should be matched to \code{survival_probability_fcn}
#' @param ... Further arguments passed on to \code{\link{find_present_value}}
#' @return A vector of estimated option present values
#' @family Equity Independent Default Intensity
#' @family American Exercise Equity Options
#' @examples
#' american(PUT, S0=100, K=110, time=0.77, const_short_rate = 0.06,
#'          const_volatility=0.20, num_time_steps=200)
#' american(callput=-1, S0=100, K=90, time=1, const_short_rate=0.025,
#'          variance_cumulation_fcn = function(T, t) {  # Term structure of vola
#'              0.45 ^ 2 * (T - t) + 0.15^2 * max(0, T-0.25)
#'          })
#' @import futile.logger
#' @export american
american = function(callput, S0, K, time,
                    const_short_rate=0, const_default_intensity=0,
                    discount_factor_fcn = function(T, t, ...){exp(-const_short_rate*(T-t))},
                    survival_probability_fcn = function(T, t, ...){exp(-const_default_intensity*(T-t))},
                    default_intensity_fcn = function(t, S, ...){const_default_intensity+0.0*S},
                    ...,
                    num_time_steps=100,
                    structure_constant=2.0,
                    std_devs_width=5.0)
{
  exact_euro_prices = black_scholes_on_term_structures(callput, S0, K, time,
                                                       discount_factor_fcn=discount_factor_fcn,
                                                       survival_probability_fcn=survival_probability_fcn, ...)$P
  flog.debug("american(callput=%s, S0=%s, K=%s, time=%s) on %s options. European BS formula vals from %s to %s",
             callput, S0, K, time, length(exact_euro_prices), min(exact_euro_prices), max(exact_euro_prices),
             name="ragtop.american_options")
  control_variate_instruments = control_variate_pairs(callput, K, time)
  N = length(control_variate_instruments) / 2
  # Get a list of grid-estimated prices
  all_prices = find_present_value(S0=S0, num_time_steps=num_time_steps,
                                  instruments=control_variate_instruments,
                                  const_short_rate=const_short_rate,
                                  discount_factor_fcn=discount_factor_fcn,
                                  default_intensity_fcn=default_intensity_fcn,
                                  ...,
                                  structure_constant=structure_constant,
                                  std_devs_width=std_devs_width)
  grid_euro_prices = unlist(all_prices[(N+1):(2*N)])
  grid_amer_prices = unlist(all_prices[1:N])
  # Find amounts we need to correct by
  european_errs = exact_euro_prices - grid_euro_prices
  flog.debug("american() correcting american option prices by errors on european exercise options ranging from %s to %s",
             min(european_errs), max(european_errs),
             name="ragtop.american_options")
  # Include corrections
  cv_amer_prices = grid_amer_prices + european_errs
  # Make sure no correction took us below early exercise value
  callput = 0 * K + 0 * time + 1 * callput  # Vectorize so indexing works below
  K = 1 * K + 0 * time + 0 * callput
  early_exer_values = callput*(S0-K)
  cv_amer_prices[cv_amer_prices<early_exer_values] = early_exer_values[cv_amer_prices<early_exer_values]
  cv_amer_prices
}
