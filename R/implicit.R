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
require(futile.logger)

stock_level = function(z, tau, K, c, sigma)
{
  ## Stock prices as a function of log-prices z
  K * exp(z - (c - 0.5 * sigma^2) * tau)
}


hazard_alpha = function(hazard_scale, S0, T, M, hzd, min_hzd, hzd_power)
{
  ## Drift rates in stock prices must incorporate hazard rates
  ## in such a way that other default-sensitive instruments
  ## like bonds are properly priced
  ## Here we work out that dependency as a function of stock
  ## prices and their corresponding hazard rates, creating an
  ## array of shifts in the constant term of the hazard rate
  ## given a sequence of timesteps

  al = (hzd - min_hzd) * array(1,dim = M + 1) / S0 ^ hzd_power
  al[1] = 0
  last_index = 2
  for (ix in 1:length(hazard_scale)) {
    hazard_scale_point = hazard_scale[ix]
    next_index = hazard_scale_point$timestep_number + 1
    al[last_index:next_index] = al[last_index:next_index] * hazard_scale_point$scale
  }
  rev(al) # Reverse the array to go in timestep backwardation order
}


implicit_grid_structure = function(tenors, M, S0, K, c, sigma, structure_constant, std_devs)
{
  ## Infer a reasonable structure for our implicit grid solver based
  ## on the voltime, structure constant, and requested grid width
  ## in standard deviations

  T = max(tenors)
  dt = T / M
  z0 = log(S0 / K) + (c - 0.5 * sigma^2) * T
  dz = sqrt(dt / structure_constant)
  z_width = std_devs * sigma * sqrt(T)
  half_N = as.integer(ceiling(z_width / dz))
  N = 2 * half_N + 1
  z = z0 + dz * (-half_N:half_N)
  flog.info("Grid structure with %s timesteps to time %s at vol %s has %s space steps to %s sdevs widths",
            M, T, sigma, N, std_devs)
  list(
    T = T, dt = dt, dz = dz, z0 = z0, z_width = z_width, half_N = half_N, N =
      N, z = z
  )
}



tridiagonals = function(S_hzd, sigma, structure_constant, a_dt_dzinv)
{
  ## Matrix entries for implicit numerical differentiation
  ## using Neumann boundary conditions
  ## a_dt_dzinv = the drift alpha at this timestep multiplied by
  ##                the timestep size dt and divided by the grid spacing dz
  ## S_hzd = S to the hazard power (S**hzd_power)
  N = length(S_hzd)
  K = N - 1
  # Entries in the center of the matrix come from finite differencing
  #  free of boundary conditions
  diag = rep(1.0 + sigma ^ 2 * structure_constant, N)
  subdiag = -0.5 * (sigma ^ 2 * structure_constant - S_hzd[2:N] * a_dt_dzinv)
  superdiag = -0.5 * (sigma ^ 2 * structure_constant + S_hzd[1:K] * a_dt_dzinv)
  # Entries on the low-z/low-S boundary (index 1)
  neumann_drift_low = a_dt_dzinv * S_hzd[1]
  diag[1] = 1. + neumann_drift_low
  superdiag[1] = -neumann_drift_low
  # Entries on the high-z/high-S boundary (index N for diag and K for subdiag)
  neumann_drift_high = a_dt_dzinv * S_hzd[N]
  diag[N] = 1. + neumann_drift_high
  subdiag[K] = -neumann_drift_high
  bad_ix = (abs(diag[2:K])-abs(subdiag[1:(K - 1)])-abs(superdiag[2:K])<0)
  if (any(bad_ix)) {
    warning(paste0("Implicit routine encounterned potentially noninvertible matrix.  Bad indexes: ",
                   toString(which(bad_ix))))
  }
  list(super = superdiag, diag = diag, sub = subdiag)
}


take_implicit_timestep = function(t, dt, dz, r, h, S0, alph,
                                  discount_factor,
                                  S, S_hzd,
                                  prev_grid_values, survival_probabilities,
                                  tridiag_matrix_entries,
                                  instrument=NULL,
                                  dividends=NULL)
{
  ## Take one timestep of an implicit solver for a given instrument
  ## The instrument, if not NULL,  must have a 'recovery_fcn' and
  ## an 'optionality_fcn' though those properties are themselves
  ## allowed to be NULL.

  prev_grid_values = adjusted_for_dividends(
    grid_values = prev_grid_values,
    t = t, dt = dt,
    r = r, h = h, S = S, S0 = S0,
    dividends = dividends
  )
  a_dt_dzinv  = alph * dt / dz
  # The value of holding this security at time t, assuming it will survive
  # to t+dt just comes from inverting our finite difference matrix
  hold_cond_on_surv = limSolve::Solve.tridiag(tridiag_matrix_entries$sub,
                                              tridiag_matrix_entries$diag,
                                              tridiag_matrix_entries$super,
                                              prev_grid_values)
  # We assume our derivative can have no negative values, so
  # we floor it at zero
  hold_cond_on_surv[hold_cond_on_surv < 0.0] = 0.0
  # If it will have value in case of default, work out what that value is
  if (is.null(instrument) || is.null(instrument$recovery_fcn)) {
    recovery_values = 0.0
  } else {
    recovery_values = instrument$recovery_fcn(S, t, hold_cond_on_surv)
  }
  # The overall value of the derivative, assuming both parties want to
  # keep it in existence, comes from the hold value conditional on
  # survival times the appropriate likelihood, plus the recovery
  # value times default likelihood
  hold_value = (survival_probabilities * hold_cond_on_surv +
                  (1. - survival_probabilities) * discount_factor * recovery_values)
  # If optionality is in play, the derivative value could be altered. This
  # can depend on _other_ derivative values that may be in play, in which
  # the instrument's optionality_fcn should handle the dependencies
  # generally by having all instruments involved be reference classes
  # (RC) which are stateful.
  # Another use of this optionality_fcn() is to have reasonable prices
  # at timesteps occurring beyond the tenor of this layer, so that
  # for example a bond is just set to the notional value corrected
  # for time value of money
  if (is.null(instrument) || is.null(instrument$optionality_fcn)) {
    new_value = hold_value
  } else {
    new_value = instrument$optionality_fcn(hold_value, S, t)
  }
  new_value
}

timestep_instruments = function(prev_grid_values,
                                t, dt, dz, r, h, z, S0, alph,
                                K, c, sigma, tau, min_hzd,
                                discount_factor, structure_constant,
                                hzd_power,
                                instruments,
                                dividends = NULL)
{
  ## Take an implicit timestep for all the given instruments, under the
  ## assumption that prev_grid_values is a matrix with one row for
  ## each instrument and one column for each of the N values of z
  ## alph is the constant terms of the hazard drift rate at this timestep
  ## Each instrument, if not NULL,  must have a 'recovery_fcn' and
  ## an 'optionality_fcn' though those properties are themselves
  ## allowed to be NULL.

  grid_values = prev_grid_values
  a_dt_dzinv = alph * dt / dz
  S = stock_level(z, tau, K, c, sigma)
  S_hzd = S^hzd_power
  matrix_entries = tridiagonals(S_hzd, sigma, structure_constant, a_dt_dzinv)
  h = min_hzd + alph * S_hzd
  survival_probabilities = exp(-h * dt)
  for (k in (1:length(instruments))) {
    instrument = instruments[k]
    prev_instr_grid_values = prev_grid_values[k,]
    flog.info("Now timestepping %s on N=%s grid", instrument,
              length(prev_instr_grid_values))
    instr_grid_vals = take_implicit_timestep(t, dt, dz, r, h, S0, alph,
                                             discount_factor,
                                             S, S_hzd,
                                             prev_instr_grid_values, survival_probabilities,
                                             matrix_entries,
                                             instrument=instrument,
                                             dividends = dividends)
    grid_values[k,] = instr_grid_vals
  }
  grid_values
}
