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


# TODO: Should we stop using t and T as variables since they could
#       be confused with the transpose function?


#' Structure of implicit numerical integration grid
#'
#' Infer a reasonable structure for our implicit grid solver based
#' on the voltime, structure constant, and requested grid width
#' in standard deviations.
#'
#' Generally speaking pricing will be good to about 10bp of
#' relative accuracy when the ratio of timesteps to voltime
#' (in annualized units) is over 200.
#'
#' Cases with pathologically low volatility may go awry (in the sense of
#'  yielding ultimately inaccurate PDE solutions), as the
#'  {structure_constant} will force a step in {z} space much bigger
#'  than the width in standard deviations.
#'
#' @param tenors Tenors of instruments to be treated on this grid
#' @param M Minimum number of timesteps on this grid
#' @param S0 An initial stock price, for setting grid scale
#' @param K An instrument reference stock price, for setting grid scale
#' @param c A continuous stock drift rate
#' @param min_z_width Minimum grid width, in log space
#' @param sigma Volatility of diffusion process (without jumps to default)
#' @param structure_constant The maximum ratio between time intervals \code{dt}
#'  and the square of space intervals \code{dz^2}
#' @param std_devs_width The number of standard deviations, in \code{sigma * sqrt(T)}
#'  units, to incorporate into the grid
#' @family Implicit Grid Solver
#'
#' @return A list with elements \describe{
#'   \item{\code{T}}{The maximum time for this grid}
#'   \item{\code{dt}}{Largest permissible timestep size}
#'   \item{\code{dz}}{Distance between space grid points}
#'   \item{\code{z0}}{Center of space grid}
#'   \item{\code{z_width}}{Width in \eqn{z} space}
#'   \item{\code{half_N}}{A misnomer, actually \eqn{(N-1)/2}}
#'   \item{\code{N}}{The number of space points}
#'   \item{\code{z}}{Locations of space points}
#' }
#' @import futile.logger
construct_implicit_grid_structure = function(tenors, M, S0, K, c, sigma, structure_constant, std_devs_width,
                                             min_z_width=0)
{
  T = max(tenors)
  flog.info("construct_implicit_grid_structure(%s, %s, %s, %s, %s, %s, %s, %s)",
            T, M, S0, K, c, sigma, structure_constant, std_devs_width,
            name='ragtop.implicit.setup'
  )
  dt = T / M
  z0 = log(S0 / K) + (c - 0.5 * sigma^2) * T
  dz = sqrt(dt / structure_constant)
  z_width = std_devs_width * sigma * sqrt(T)
  if (dz*4>z_width) {
    flog.warn("Volatility is tiny in comparison to timestep size.  Very few space steps for the underlying will be used, and answers should be treated skeptically.",
              name='ragtop.implicit.setup.width')
  }
  half_N = as.integer(ceiling(z_width / dz))
  N = 2 * half_N + 1
  z = z0 + dz * (-half_N:half_N)
  flog.info("Grid structure with %s timesteps to time %s at vol %s has %s space steps to (%s sdevs widths, z width %s) of size %s",
            M, T, sigma, N, std_devs_width, z_width, dz,
            name='ragtop.implicit.setup')
  list(
    T = T, dt = dt, dz = dz, z0 = z0,
    z_width = z_width, half_N = half_N, N = N, z = z
  )
}




#' Matrix entries for implicit numerical differentiation using Neumann boundary conditions
#'
#' @param drift Vector of drift rate of underlying equity grid points, including
#'   induced drift from default intensity
#' @param sigma Volatility of diffusion process (without jumps to default)
#' @param structure_constant The ratio between time interval \code{dt}
#'  and the square of space interval \code{dz^2}
#' @return A list with elements \code{super}, \code{diag} and \code{sub}
#'   containing the superdiagonal, diagonal and subdiagonal of the implicit
#'   timestep differencing matrix
#' @import futile.logger
#' @export construct_tridiagonals
construct_tridiagonals = function(sigma, structure_constant, drift)
{
  N = length(drift)
  Nm1 = N - 1
  # Entries in the center of the matrix come from finite differencing
  #  free of boundary conditions
  diag = rep(1.0 + sigma ^ 2 * structure_constant, N)
  subdiag = -0.5 * (sigma ^ 2 * structure_constant - drift[2:N])
  superdiag = -0.5 * (sigma ^ 2 * structure_constant + drift[1:Nm1])
  # Entries on the low-z/low-S boundary (index 1)
  neumann_drift_low = drift[1]
  diag[1] = 1. + neumann_drift_low
  superdiag[1] = -neumann_drift_low
  # Entries on the high-z/high-S boundary (index N for diag and Nm1 for subdiag)
  neumann_drift_high = drift[N]
  diag[N] = 1. + neumann_drift_high
  subdiag[Nm1] = -neumann_drift_high
  bad_ix = (abs(diag[2:Nm1])-abs(subdiag[1:(Nm1 - 1)])-abs(superdiag[2:Nm1])<0)
  if (any(bad_ix)) {
    flog.info("Implicit routine encountered potentially noninvertible matrix.  Bad indexes: %s",
              toString(which(bad_ix)),
              name='ragtop.implicit.timestep.construct_tridiagonals')
  }
  list(super = superdiag, diag = diag, sub = subdiag)
}


#' Backwardate grid values one timestep
#'
#' Take one timestep of an implicit solver for a given instrument
#'
#' @param t Time after this timestep has been taken
#' @param S Underlying equity values for the grid
#' @param survival_probabilities Vector of probabilities of survival
#'  for each space grid node
#' @param tridiag_matrix_entries Diagonal, superdiagonal and subdiagonal
#'  of tridiagonal matrix from the numerical integrator
#' @param full_discount_factor A discount factor for the transform from
#'  grid values to actual derivative prices
#' @param instr_name Name of instrument to use in log messages
#' @param local_discount_factor A discount factor to apply to
#'  recovery values
#' @param discount_factor_fcn A function for computing present values to
#'   time \code{t} of various cashflows occurring during this timestep, with
#'   arguments \code{T}, \code{t}
#' @param prev_grid_values A vector of space grid values from the
#'  previously calculated timestep
#' @param instrument If not NULL/NA,  must have a \code{recovery_fcn} and
#'  an \code{optionality_fcn} though those properties are
#'  themselves allowed to be NA.
#' @param dividends A \code{data.frame} with columns \code{time}, \code{fixed},
#'   and \code{proportional}.  Dividend size at the given \code{time} is
#'   then expected to be equal to \code{fixed + proportional * S / S0}
#' @import limSolve
#' @return Grid values for the instrument after taking the implicit timestep
#' @family Implicit Grid Solver
#' @export take_implicit_timestep
take_implicit_timestep = function(t, S, full_discount_factor,
                                  local_discount_factor,
                                  discount_factor_fcn,
                                  prev_grid_values, survival_probabilities,
                                  tridiag_matrix_entries,
                                  instrument=NULL,
                                  dividends=NULL,
                                  instr_name="this instrument")
{
  # N.B.When accessing optionality_fcn() or recovery_fcn() methods, we adjust
  #  putative grid values by the discount factor to obtain non-transformed
  #  prices suitable for comparison to strikes, etc.

  # The value of holding this security at time t, assuming it will survive
  # to t+dt just comes from inverting our finite difference matrix
  hold_cond_on_surv = limSolve::Solve.tridiag(tridiag_matrix_entries$sub,
                                              tridiag_matrix_entries$diag,
                                              tridiag_matrix_entries$super,
                                              prev_grid_values)
  # limSolve::Solve.tridiag() uses LAPACK DGTSV Gaussian elimination
  #  with partial pivoting rather than the naive tridiagonal algorithm

  # We assume our derivative can have no negative values, so
  # we floor it at zero
  hold_cond_on_surv[hold_cond_on_surv < 0.0] = 0.0
  # If it will have value in case of default, work out what that value is
  if (is.blank(instrument) || is.blank(instrument$recovery_fcn)) {
    recovery_values = 0.0
  } else {
    recovery_at_t = instrument$recovery_fcn(v=hold_cond_on_surv/full_discount_factor,
                                            S=S, t=t,
                                            discount_factor_fctn=discount_factor_fcn)
    recovery_values = full_discount_factor * recovery_at_t
  }
  # The overall value of the derivative, assuming both parties want to
  # keep it in existence, comes from the hold value conditional on
  # survival times the survival likelihood, plus the recovery
  # value times default likelihood
  survival_value = survival_probabilities * hold_cond_on_surv
  default_value = (1. - survival_probabilities) * local_discount_factor * recovery_values
  hold_value = survival_value + default_value
  flog.info("Timestep of %s to t=%s has hold conditional on survival averaging %s, recovery values averaging %s and survival probabilities averaging %s for an average value all-in of %s",
            instr_name, t, mean(hold_cond_on_surv), mean(recovery_values),
            mean(survival_probabilities), mean(hold_value),
            name='ragtop.implicit.timestep')
  # If optionality is in play, the derivative value could be altered. This
  # can depend on _other_ derivative values that may be in play, in which
  # the instrument's optionality_fcn should handle the dependencies
  # generally by having all instruments involved be reference classes
  # (RC) which are stateful.
  # Another use of this optionality_fcn() is to have reasonable prices
  # at timesteps occurring beyond the tenor of this layer, so that
  # for example a bond is just set to the notional value corrected
  # for time value of money when beyond maturity
  if (is.blank(instrument) || is.blank(instrument$optionality_fcn)) {
    new_value = hold_value
  } else {
    undiscounted_value = instrument$optionality_fcn(hold_value/full_discount_factor, S, t,
                                                    discount_factor_fctn=discount_factor_fcn)
    new_value = full_discount_factor * undiscounted_value
  }
  new_value
}


#' Take an implicit timestep for all the given instruments
#'
#' Backwardate grid values for all the given instruments from a set of grid
#'   values matched to time \code{t+dt} to form a new set of grid value as
#'   of time \code{t}.
#'
#' @inheritParams take_implicit_timestep
#' @param S0 Time zero price of the base equity
#' @param instruments Instruments corresponding to layers of the value grid in \code{prev_grid_values}
#' @param dt Interval to the end of this timestep
#' @param z Space grid value morphable to stock prices using \code{stock_level_fcn}
#' @param stock_level_fcn A function for changing space grid value to stock
#'   prices, with arguments \code{z} and \code{t}
#' @param default_intensity_fcn A function for computing default intensity
#'   occurring during this timestep, dependent on time and stock price, with
#'   arguments \code{t}, \code{S}.
#' @param variance_cumulation_fcn A function for computing total stock variance
#'   occurring during this timestep, with arguments \code{T}, \code{t}.  E.g. with
#'   a constant volatility \eqn{s} this takes the form \eqn{(T-t)s^2}.
#' @param prev_grid_values A matrix with one column for each
#'   instrument and one row for each of the \eqn{N} values of \code{z}
#' @return  Grid values after applying an implicit timestep
#' @family Implicit Grid Solver
timestep_instruments = function(z, prev_grid_values,
                                t, dt, S0,
                                instruments,
                                stock_level_fcn,
                                discount_factor_fcn,
                                default_intensity_fcn,
                                variance_cumulation_fcn,
                                dividends = NULL)
{
  full_discount_factor = discount_factor_fcn(t,0)
  prev_timestep_discount_factor = discount_factor_fcn(t+dt,0)
  local_discount_factor = discount_factor_fcn(t+dt,t)
  r = -log(local_discount_factor)/dt
  S = stock_level_fcn(z, t)
  h = default_intensity_fcn(t, S)
  if (length(h) != length(S)) {
    stop("Default intensity function must return an array of the same size as S")
  }
  sigma = sqrt(variance_cumulation_fcn(t+dt, t)/dt)
  div_adj_grid_values = adjust_for_dividends(prev_grid_values,
                                             t, dt, r, h, S, S0,
                                             dividends = dividends)
  survival_probabilities = exp(-h * dt)
  flog.info("Default intensity at %s averages %s.  Short rate averages %s.  Local volatility is %s. Surv probs range from %s to %s with avg %s.",
            t, mean(h), mean(r), sigma, min(survival_probabilities), max(survival_probabilities), mean(survival_probabilities),
            name='ragtop.implicit.timestep')
  dz = diff(z)[1] # We are assuming a regular grid in z space
  structure_constant = dt/dz^2
  flog.info("Structure constant at %s for dt=%s is %s.", t, dt, structure_constant,
            name='ragtop.implicit.timestep')
  matrix_entries = construct_tridiagonals(sigma, structure_constant, drift=h*dt/dz)
  for (k in (1:length(instruments))) {
    instrument = instruments[[k]]
    instr_name = names(instruments)[[k]]
    prev_instr_grid_values = div_adj_grid_values[,k]
    if (instrument$maturity >= t) {
      if (all(is.na(prev_instr_grid_values))) {
        instr_grid_vals = full_discount_factor * instrument$optionality_fcn(0, S, t)
      } else {
        # Update hold value for cashflows
        instr_methods = instrument$getRefClass()$methods()
        flog.info("Now timestepping %s from %s to %s on N=%s grid, current mean value on grid is %s for a mean price of %s",
                  instr_name, t+dt, t, length(prev_instr_grid_values),
                  mean(prev_instr_grid_values), mean(prev_instr_grid_values)/(full_discount_factor*local_discount_factor),
                  name='ragtop.implicit.timestep')
        if ("update_cashflows" %in% instr_methods) {
          cash_increase = instrument$update_cashflows(t, t+dt,
                                                      discount_factor_fctn=discount_factor_fcn)
          grid_cash_increase =  cash_increase * prev_timestep_discount_factor
          flog.info("Instrument %s has cashflows up to %s in interval (%s,%s].  Increasing grid values prior to taking the timestep by the corresponding change-of-variables (x %s) amount up to %s.",
                    instr_name, max(cash_increase), t, t+dt, prev_timestep_discount_factor, max(grid_cash_increase),
                    name='ragtop.implicit.timestep')
          prev_instr_grid_values = prev_instr_grid_values + grid_cash_increase
        }
        instr_grid_vals = take_implicit_timestep(t, S, full_discount_factor,
                                                 local_discount_factor,
                                                 discount_factor_fcn,
                                                 prev_instr_grid_values,
                                                 survival_probabilities,
                                                 matrix_entries,
                                                 instrument = instrument,
                                                 dividends = dividends,
                                                 instr_name=instr_name)
        flog.info("Done timestepping %s from %s to %s on N=%s grid, new mean value on grid is %s for a mean price of %s",
                  instr_name, t+dt, t, length(prev_instr_grid_values),
                  mean(instr_grid_vals), mean(instr_grid_vals)/full_discount_factor,
                  name='ragtop.implicit.timestep')
      }
    } else {
      flog.info("Instrument %s has maturity %s prior to interval [%s,%s].  Skipping any grid computation.",
                instr_name, instrument$maturity, t, t+dt,
                name='ragtop.implicit.timestep')
      instr_grid_vals = NA*(prev_instr_grid_values)
    }
    div_adj_grid_values[,k] = instr_grid_vals
  }
  div_adj_grid_values
}


#' A time grid with extra times inserted for coupons, calls and puts
#'
#' At its base, this function chooses a time grid with \code{1+min_num_time_steps}
#'   elements from 0 to \code{Tmax}.  Any coupon, call, or put times occurring in
#'   one of the supplied instruments are also inserted.
#'
#' @param min_num_time_steps The minimum number of timesteps the output vector should have
#' @param Tmax The maximum time on the grid
#' @param instruments A set of instruments whose maturity and terms
#'   and conditions can introduce extra timesteps.  Each will be queried for the output of
#'   a \code{critical_times} function.
#' @return A vector of times at which the grid should have nodes
#' @family Implicit Grid Solver
infer_conforming_time_grid = function(min_num_time_steps, Tmax, instruments=NULL)
{
  time_grid = seq(from=0, to=Tmax, length.out=1+min_num_time_steps)
  if (!is.blank(instruments)) {
    # Cycle through instruments, include their maturities and critical times
    for (k in (1:length(instruments))) {
      instrument = instruments[[k]]
      time_grid = c(time_grid, instrument$maturity)
      instr_name = names(instruments)[[k]]
      flog.info("Checking instrument %s of maturity %s for terms and conditions requiring time grid entries",
                instr_name, instrument$maturity,
                name='ragtop.implicit.setup')
      instr_methods = instrument$getRefClass()$methods()
      if ("critical_times" %in% instr_methods) {
        inst_crit = instrument$critical_times()
        time_grid = c(time_grid, inst_crit)
        flog.info("Relevant terms and conditions for instrument %s in time grid were: %s",
                  instr_name, dput(inst_crit),
                  name='ragtop.implicit.setup')
      } else {
        flog.info("No relevant terms and conditions for instrument %s in time grid",
                  instr_name,
                  name='ragtop.implicit.setup')
      }
    }
  }
  time_grid = unique(time_grid)
  time_grid = time_grid[time_grid>=0 & time_grid<=Tmax]
  # Clear out any extraneous timesteps.  If we have any such,
  #  then the unique digitized grid would be shorter, so that is
  #  how we check
  digitized_grid = unique(signif(time_grid,
                                 digits=TIME_RESOLUTION_SIGNIF_DIGITS))
  if (length(digitized_grid)<length(time_grid)) {
    # Shorter.  Keep the endpoints but digitize everything in between
    betw = unique(signif(time_grid[time_grid>0 & time_grid<Tmax],
                         digits=TIME_RESOLUTION_SIGNIF_DIGITS))
    flog.info("Some of the %s requested timesteps are very close to each other.  Combined them to create a grid with only %s timesteps.",
              length(time_grid), 2+length(betw),
              name='ragtop.implicit.setup')
    time_grid = unique(c(0,betw,Tmax))
  }
  time_grid = sort(time_grid)
  flog.info("%s time steps requested. Instrument terms and conditions bring the total number to %s up to max t=%s",
            min_num_time_steps, length(time_grid)-1, time_grid[length(time_grid)],
            name='ragtop.implicit.setup')
  time_grid
}

#' Numerically integrate the pricing differential equation
#'
#' Use an implicit integration scheme to numerically integrate
#' the pricing differential equation for each of the given instruments,
#' backwardating from time \code{Tmax} to time 0.
#' @inheritParams timestep_instruments
#' @param Tmax The maximum time on the grid, from which
#'   all backwardation steps will take place.
#' @param min_num_time_steps The minimum number of timesteps used.  Calls,
#'   puts and coupons may result in extra timesteps taken.
#' @param instruments A list of instruments to be priced.  Each
#'   one must have a \code{strike} and a \code{optionality_fcn}, as
#'   with \code{\link{GridPricedInstrument}} and its subclasses.
#'
#' @return A grid of present values of derivative prices, adapted to \code{z} at
#'   each timestep.  Time zero value will appear in the first index.
#' @family Implicit Grid Solver
#' @export integrate_pde
integrate_pde <- function(z, min_num_time_steps, S0, Tmax, instruments,
                            stock_level_fcn,
                            discount_factor_fcn,
                            default_intensity_fcn,
                            variance_cumulation_fcn,
                            dividends=NULL)
{
  time_pts = infer_conforming_time_grid(min_num_time_steps, Tmax, instruments=instruments)
  num_time_pts = length(time_pts)
  num_time_steps = num_time_pts-1
  num_space_pts = length(z)
  num_instruments = length(instruments)
  grid = array(data=NA, dim=c(num_time_pts, num_space_pts, num_instruments))
  S_final = stock_level_fcn(z, Tmax)
  df_final = discount_factor_fcn(Tmax, 0)
  flog.info("Discount factor to Tmax=%s is %s", Tmax, df_final, name='ragtop.implicit.setup')
  # Set the initial condition for the PDE from instrument values at
  #  maximum time, discounted by df_final according to our change of variables
  for (k in (1:length(instruments))) {
    instrument = instruments[[k]]
    instr_name = names(instruments)[[k]]
    if (instrument$maturity>=Tmax) {
      undisc_terminal_vals = instrument$terminal_values(0*S_final, S=S_final,
                                                        t=instrument$maturity,
                                                        discount_factor_fctn=discount_factor_fcn)
      if (any(is.na(undisc_terminal_vals))) {  # Make sure nobody gave us a crazy instrument object
        flog.error("Bug in instrument class.  Terminal values at Tmax=%s for %s returned NA in %s of %s cases",
                   Tmax, instr_name, sum(is.na(undisc_terminal_vals)),
                   length(undisc_terminal_vals),
                   name='ragtop.implicit.setup')
        stop(paste("Invalid terminal values for", instr_name))
      }
      grid[num_time_pts,,k] = df_final * undisc_terminal_vals
      flog.info("Terminal values at Tmax=%s for %s avg is %s after present valuing using discount factor %s",
                Tmax, instr_name, mean(grid[num_time_pts,,k]), df_final,
                name='ragtop.implicit.setup')
    } else {
      flog.info("Terminal values of %s are being set to NA because its maturity %s is strictly less than Tmax of %s",
                instr_name, instrument$maturity, Tmax,
                name='ragtop.implicit.setup')
      grid[num_time_pts,,k] = NA
    }
    instrument$last_computed_grid = NA * S_final
  }
  grid = iterate_grid_from_timestep(num_time_steps, time_pts, z, S0, instruments,
                             stock_level_fcn=stock_level_fcn,
                             discount_factor_fcn=discount_factor_fcn,
                             default_intensity_fcn=default_intensity_fcn,
                             variance_cumulation_fcn=variance_cumulation_fcn,
                             dividends = dividends,
                             grid = grid)
  # Return grid values at all t.  The present values will be on the
  # grid at index 1
  grid
}

#' Iterate over a set of timesteps to integrate the pricing differential equation
#'
#' Timestep an implicit integration scheme to numerically integrate
#' the pricing differential equation for each of the given instruments,
#' backwardating from time \code{Tmax} to time 0.
#' @inheritParams timestep_instruments
#' @inheritParams integrate_pde
#' @param grid An optional grid into which results at each timestep will
#'   be written.  Its size should be at least
#'   \code{(1+starting_time_step, length(z), length(instruments))}
#' @param starting_time_step The index into time_pts of the first timestep
#'  to be emplyed.  This must be no larger than the length of time_pts, minus one
#' @param time_pts Time nodes to be treated on the grid
#' @param original_grid_values Grid values to timestep from
#' @param instruments A list of instruments to be priced.  Each
#'   one must have a \code{strike} and a \code{optionality_fcn}, as
#'   with \code{\link{GridPricedInstrument}} and its subclasses.
#'
#' @return Either a populated grid of present values of derivative prices, or a matrix
#'   of values at the first time point, adapted to \code{z} at
#'   each timestep.  Time zero value will appear in the first index of any grid.
#' @family Implicit Grid Solver
#' @export iterate_grid_from_timestep
iterate_grid_from_timestep = function(starting_time_step, time_pts, z, S0, instruments,
                                      stock_level_fcn,
                                      discount_factor_fcn,
                                      default_intensity_fcn,
                                      variance_cumulation_fcn,
                                      dividends = NULL,
                                      grid = NULL,
                                      original_grid_values=as.matrix(grid[1+starting_time_step,,]))
{
  prev_grid_values = original_grid_values
  new_grid_values = original_grid_values
  # Take starting_time_step timesteps to integrate
  for (m in (starting_time_step:1)) {
    t = time_pts[[m]]
    dt = time_pts[[m+1]] - time_pts[[m]]
    new_grid_values = timestep_instruments(
      z, prev_grid_values,
      t, dt, S0,
      instruments,
      stock_level_fcn=stock_level_fcn,
      discount_factor_fcn=discount_factor_fcn,
      default_intensity_fcn=default_intensity_fcn,
      variance_cumulation_fcn=variance_cumulation_fcn,
      dividends=dividends)
    if (!is.null(grid)) {
      # Save values at this timestep if necessary
      grid[m,,] = new_grid_values
    }
    prev_grid_values = new_grid_values
  }
  if (is.null(grid)) {
    iteration_result = new_grid_values
  } else {
    iteration_result = grid
  }
}

#' Use a model to estimate the present value of financial derivatives on a grid of initial underlying values
#'
#' Use a finite difference scheme to form estimates of present values for a variety
#'  of stock prices on a grid of initial underlying prices, determined by constructing
#'  a logarithmic equivalent conforming to the grid parameters \code{structure_constant}
#'  and \code{structure_constant}
#'
#' If any instrument in the \code{instruments} has a strike, then the grid will be
#'  normalized to the last such instrument's strike.
#'
#' @inheritParams construct_implicit_grid_structure
#' @inheritParams timestep_instruments
#' @param num_time_steps Minimum number of time steps in the grid
#' @param const_volatility A constant to use for volatility in case \code{variance_cumulation_fcn}
#'  is not given
#' @param const_short_rate A constant to use for the instantaneous interest rate in case \code{discount_factor_fcn}
#'  is not given
#' @param const_default_intensity A constant to use for the instantaneous default intensity in case \code{default_intensity_fcn}
#'  is not given
#' @param grid_center A reasonable central value for the grid, defaults to S0 or an instrument strike
#' @param borrow_cost Stock borrow cost, affecting the drift rate
#' @param dividend_rate Continuous dividend rate, affecting the drift rate
#' @param override_Tmax A different maximum time on the grid to enforce
#' @param instruments A list of instruments to be priced.  Each
#'   one must have a \code{strike} and a \code{optionality_fcn}, as
#'   with \code{\link{GridPricedInstrument}} and its subclasses.
#' @family Equity Dependent Default Intensity
#' @family Implicit Grid Solver
#'
#' @export form_present_value_grid
form_present_value_grid = function(S0, num_time_steps, instruments,
                              const_volatility=0.5, const_short_rate=0,
                              const_default_intensity=0, override_Tmax=NA,
                              discount_factor_fcn = function(T, t, ...){exp(-const_short_rate*(T-t))},
                              default_intensity_fcn = function(t, S, ...){const_default_intensity+0.0*S},
                              variance_cumulation_fcn = function(T, t){const_volatility^2*(T-t)},
                              dividends=NULL,
                              borrow_cost=0.0,
                              dividend_rate=0.0,
                              structure_constant=2.0,
                              std_devs_width=3.0,
                              grid_center=NA)
{
  if (is.blank(override_Tmax)) {
    Tmax = 0
  } else {
    Tmax = override_Tmax
  }
  K = S0
  # Get a strike and maximum maturity from the instruments
  for (k in (1:length(instruments))) {
    instrument = instruments[[k]]
    instr_name = names(instruments)[[k]]
    instr_fields = instrument$getRefClass()$fields()
    if ("strike" %in% instr_fields && instrument$strike > 0) {
      # Set the target strike for grid structure based on the
      # strike of the *last* instrument in our list
      K = instrument$strike
    }
    if (instrument$maturity > Tmax && is.blank(override_Tmax)) {
      Tmax = instrument$maturity
    }
    flog.info("Instrument %s: %s", k, instr_name,
              name='ragtop.implicit.setup')
  }
  if (!is.blank(grid_center)) {
    flog.info("Grid center forced to K=%s", grid_center, name='ragtop.implicit.setup')
    K = grid_center
  }
  if (Tmax <= 0) {
    stop("Cannot compute present value when no instrument maturity is positive")
  } else {
    flog.info("Max maturity: %s", Tmax,
              name='ragtop.implicit.setup')
  }
  sigma = sqrt(variance_cumulation_fcn(Tmax, 0) / Tmax)
  r = -log(discount_factor_fcn(Tmax, t=0) ) / Tmax
  c = r - dividend_rate - borrow_cost
  flog.info("Constant sigma equivalent sigma=%s and rate equivalent r=%s, implies c=%s",
            sigma, r, c,
            name='ragtop.implicit.setup')
  grid_structure = construct_implicit_grid_structure(Tmax, num_time_steps, S0, K,
                                           c, sigma,
                                           structure_constant=structure_constant,
                                           std_devs_width=std_devs_width)
  stock_level_fcn = function(z, t) {
    S_levs = K * exp(z - (c - 0.5 * sigma^2) * (Tmax-t))
    flog.info("z = %s stock levels at %s from min(S)=%s to max(S)=%s",
              length(S_levs), t, min(S_levs), max(S_levs),
              name='ragtop.implicit.setup')
    flog.debug("z %s", deparse(z, width.cutoff=150),
               name='ragtop.implicit.setup')
    flog.debug("S %s", deparse(S_levs, width.cutoff=150),
               name='ragtop.implicit.setup')
    S_levs
  }
  grid = integrate_pde(grid_structure$z,
                       num_time_steps,
                       S0, Tmax, instruments,
                        stock_level_fcn,
                        discount_factor_fcn,
                        default_intensity_fcn,
                        variance_cumulation_fcn,
                        dividends=dividends)
  flog.info("Completed PDE integration",
            name='ragtop.implicit.form_present_value_grid')
  present_value_grid = cbind(as.matrix(grid[1,,]),
                             matrix(stock_level_fcn(grid_structure$z,0), ncol=1))
  if (is.null(names(instruments))) {
    colnames(present_value_grid)[length(instruments)+1] = "Underlying"
  } else {
    colnames(present_value_grid) = c(names(instruments), "Underlying")
  }
  present_value_grid
}


#' Use a model to estimate the present value of financial derivatives
#'
#' Use a finite difference scheme to form estimates of present values for a variety
#'  of stock prices.  Once the grid has been created, interpolate to obtain the
#'  value of each instrument at the present stock price \code{S0}
#'
#' @inheritParams form_present_value_grid
#' @inheritParams construct_implicit_grid_structure
#' @return A list of present values, with the same names as \code{instruments}
#' @family Equity Dependent Default Intensity
#' @family Implicit Grid Solver
#'
#' @export find_present_value
find_present_value = function(S0, num_time_steps, instruments,
                                   const_volatility=0.5, const_short_rate=0, const_default_intensity=0, override_Tmax=NA,
                                   discount_factor_fcn = function(T, t, ...){exp(-const_short_rate*(T-t))},
                                   default_intensity_fcn = function(t, S, ...){const_default_intensity+0.0*S},
                                   variance_cumulation_fcn = function(T, t){const_volatility^2*(T-t)},
                                   dividends=NULL,
                                   borrow_cost=0.0,
                                   dividend_rate=0.0,
                                   structure_constant=2.0,
                                   std_devs_width=3.0)
{
  if (is.blank(dividends)) {
    divs_descr = "NULL"
  } else {
    divs_descr = nrow(dividends)
  }
  flog.info("find_present_value(S0=%s, num_time_steps=%s, <%s instruments>, <divs: %s>, structure_constant=%s, std_devs_width=%s)",
            S0, num_time_steps, length(instruments), divs_descr, structure_constant, std_devs_width,
            name='ragtop.implicit.find_present_value')
  if (is.null(names(instruments))) {
    names(instruments) = lapply(instruments,
                                function(inst) {inst$name})
    names(instruments)[is.na(names(instruments))] = seq(1,sum(is.na(names(instruments))))
  }
  present_value_grid = form_present_value_grid(S0=S0, num_time_steps=num_time_steps, instruments=instruments,
                                               const_volatility=const_volatility, const_short_rate=const_short_rate,
                                               const_default_intensity=const_default_intensity, override_Tmax=override_Tmax,
                                               discount_factor_fcn = discount_factor_fcn,
                                               default_intensity_fcn = default_intensity_fcn,
                                               variance_cumulation_fcn = variance_cumulation_fcn,
                                               dividends=dividends,
                                               borrow_cost=borrow_cost,
                                               dividend_rate=dividend_rate,
                                               structure_constant=structure_constant,
                                               std_devs_width=std_devs_width)
  present_values = list()
  for (k in (1:length(instruments))) {
    instrument = instruments[[k]]
    instr_name = names(instruments)[[k]]
    present_values[[instr_name]] = stats::spline(
      x=present_value_grid[,"Underlying"],
      y=present_value_grid[,instr_name],
      xout=S0)$y
  }
  present_values
}
