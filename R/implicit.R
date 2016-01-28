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
require(stats)
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

