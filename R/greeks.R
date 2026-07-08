## ragtop -- convertibles pricing in R
##
## Copyright (C) 2026  Brian Boonstra <ragtop@boonstra.org>
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

#' Greeks available for computation, presently  "delta", "gamma", "theta", "vega",
#' "rates_dv01", and "credit_dv01".  Note in particular that we use the term
#' \code{rates_dv01} rather than \code{rho}
#' @family Greeks
#' @export
GREEK_NAMES = c("delta", "gamma", "theta", "vega", "rates_dv01", "credit_dv01")

#' Substitute discount factor function with altered rates
#'
#' In order to test the effects of changes in risk-free rates to instrument
#' prices, we need a way to take a given curve and alter it slightly, with
#' bumped rates.  This routine returns a new discount factor function with
#' the given bump included.
#'
#' @param discount_factor_fcn The unbumped discount factor function, a baseline
#' @param rate_bump The change in rates \eqn{r} to apply, effected through multiplying discount factors by \eqn{e^{-r(T-t)}}
#' @return A discount factor function, suitable for use in pricing routines
#' @family Greeks
#'
#' @noRd
bump_discount_factor_fcn = function(discount_factor_fcn, rate_bump)
{
  #  force() evals the original function right away to avoid circular reference
  force(discount_factor_fcn); force(rate_bump)
  function(T, t, ...) discount_factor_fcn(T, t, ...) * exp(-rate_bump * (T - t))
}

#' Substitute default intensity function with altered hazard rate
#'
#' In order to test the effects of changes in default intensity to instrument
#' prices, we need a way to take a given curve and alter it slightly, with
#' bumped hazard rates.  This routine returns a new default intensity function with
#' the given bump included.
#'
#' @param default_intensity_fcn The unbumped default intensity function, a baseline
#' @param hazard_bump The additive change in hazard rate \eqn{h} to apply
#' @return A default intensity function, suitable for use in pricing routines
#' @family Greeks
#'
#' @noRd
bump_default_intensity_fcn = function(default_intensity_fcn, hazard_bump)
{
  force(default_intensity_fcn); force(hazard_bump)
  function(t, S, ...) default_intensity_fcn(t, S, ...) + hazard_bump
}

#' Substitute variance cumulation function with altered volatility
#'
#' In order to test the effects of changes in volatility to instrument
#' prices, we need a way to take a given variance cumulation curve and alter it slightly, with
#' bumped volatility.  This routine returns a new variance cumulation function with
#' the given bump included.
#'
#' @param variance_cumulation_fcn The unbumped variance cumulation function, a baseline
#' @param vola_bump The additive change in volatility  to apply
#' @return A variance cumulation function, suitable for use in pricing routines
#' @family Greeks
#'
#' @noRd
bump_variance_cumulation_fcn = function(variance_cumulation_fcn, vola_bump)
{
  force(variance_cumulation_fcn); force(vola_bump)
  function(T, t) {
    dt = T - t
    base_var = variance_cumulation_fcn(T, t)
    sigma = ifelse(dt > 0, sqrt(base_var / dt), 0)
    (sigma + vola_bump)^2 * dt
  }
}


#' Reprice instruments with a bumped parameter
#'
#' Evaluate one repricing of the whole instrument set, catching all exceptions, with NA wherever a
#' price is missing or non-finite.  A thrown error fails every instrument at once.
#'
#' @param reprice_runner A function that reprices instruments, possibly unsafely (for example if a bump results in negative variance)
#' @param bump The change in model parameter to be applied
#' @param instrument_names The names of instruments being analyzed
#' @param greek The name of the greek being bumped
#' @return Instrument prices, with the given names
#' @family Greeks
#'
#' @keywords internal
safe_reprice = function(reprice_runner, bump, instrument_names, greek="Unspecified")
{
  vals = tryCatch(reprice_runner(bump),
                  error = function(e) {
                    flog.warn("Repricing for greek '%s' at bump %s failed: %s",
                              greek, bump, conditionMessage(e), name="ragtop.greeks")
                    stats::setNames(rep(NA_real_, length(instrument_names)), instrument_names)
                  })
  vals = as.numeric(vals[instrument_names])
  vals[!is.finite(vals)] = NA_real_
  stats::setNames(vals, instrument_names)
}

#' Find a greek by finite differencing
#'
#' Combine base values and their up/down bumped values into per-instrument
#'  derivative estimates.  Uses a central difference where possible, otherwise a
#'  one-sided difference, otherwise reports "failed".  All arguments are named
#'  vectors over the instruments; the result is vectorized the same way.
#' @param v0 Baseline prices
#' @param v_up Prices with an "up" bump
#' @param v_down Prices with a "down" bump
#' @param bump The size of the parameter bump
#' @return A list with \code{value} containing a named numeric vector of greek values found, and \code{method} recording central versus one-sided differences used
#' @family Greeks
#'
#' @keywords internal
greek_by_fd = function(v0, v_up, v_down, bump)
{
  up_ok = is.finite(v_up)
  down_ok = is.finite(v_down)
  v0_ok = is.finite(v0)
  central = up_ok & down_ok
  forward = up_ok & !down_ok & v0_ok
  backward = down_ok & !up_ok & v0_ok

  value = rep(NA_real_, length(v0))
  method = rep("failed", length(v0))
  value[central]  = (v_up[central] - v_down[central]) / (2 * bump)
  value[forward]  = (v_up[forward] - v0[forward]) / bump
  value[backward] = (v0[backward] - v_down[backward]) / bump
  method[central]  = "central"
  method[forward]  = "forward"
  method[backward] = "backward"

  list(value = stats::setNames(value, names(v0)), method = stats::setNames(method, names(v0)))
}


#' Find a greek for a set of instruments, using a sequence of possible bump sizes.
#'
#' This function uses finite differencing to compute greek values.  It attempts to trade off
#' errors and failures from (a) too-large bump size (giving rise to, e.g. negative variances) against (b)
#' numerical noise interfering with estimates from using a too-small bump size.
#'
#' Here, we walk a descending sequence of bump magnitudes for the whole instrument set.  At
#'  each, reprice up and down (either may fail per instrument) and form the best
#'  available per-instrument difference.  Smaller bumps are more accurate (ignoring
#'  underflow), so keep shrinking until every instrument has at least a one-sided
#'  estimate, or the smallest bump is reached.  Every instrument takes its value
#'  and method from this single final iteration.
#' @inheritParams safe_reprice
#' @param v0 Baseline prices
#' @param descending_bumps Potential bump changes in model parameter to be applied
#' @return A list with \code{value} containing a named numeric vector of greek values found, \code{bump} containing the bump size used, and \code{method} recording central versus one-sided differences used
#' @family Greeks
#'
#' @keywords internal
robust_greek = function(v0, reprice_runner, descending_bumps, instrument_names, greek="Unspecified")
{
  n = length(descending_bumps)
  est = list(value = stats::setNames(rep(NA_real_, length(instrument_names)), instrument_names),
             method = stats::setNames(rep("failed", length(instrument_names)), instrument_names),
             bump = NA_real_)
  for (i in seq_len(n)) {
    bump = descending_bumps[[i]]
    v_up = safe_reprice(reprice_runner, bump, instrument_names, greek)
    v_down = safe_reprice(reprice_runner, -bump, instrument_names, greek)
    est = greek_by_fd(v0, v_up, v_down, bump)
    est$bump = bump
    if (all(est$method != "failed")) break
    if (i < n) {
      failed = instrument_names[est$method == "failed"]
      flog.warn("Greek '%s' failed for %s at bump %s; trying a smaller bump.",
                greek, toString(failed), bump, name="ragtop.greeks")
    }
  }
  est
}


#' Representative constant-equivalent level of the model volatility, used to scale
#'  default bump sizes
#' @family Greeks
#'
#' @noRd
effective_volatility = function(variance_cumulation_fcn, Tmax) {
  sqrt(variance_cumulation_fcn(Tmax, 0) / Tmax)
}
#' Representative constant-equivalent level of the interest rates, used to scale
#'  default bump sizes
#' @family Greeks
#'
#' @noRd
effective_rate = function(discount_factor_fcn, Tmax) {
  -log(discount_factor_fcn(Tmax, 0)) / Tmax
}
#' Representative constant-equivalent level of the hazard rates, used to scale
#'  default bump sizes
#' @family Greeks
#'
#' @noRd
effective_hazard = function(default_intensity_fcn, S0) {
  mean(default_intensity_fcn(0, S0))
}

#' Reasonable bump sizes to try when finite-differencing for greeks
#'
#' A descending sequence of bump sizes derived from a representative \code{x}.
#'  The first bump is \code{x/divisor} (e.g. a \code{0.02} vola gives \code{0.02/40}), then
#'  \code{x/(2*divisor)}, \code{x/(3*divisor)}  When \code{x} is
#'  zero or non-finite there is no scale to work from, so we return an empty
#'  sequence.  A user-overridden bump size is required to actually receive a calculated
#'  greek in cases of zero or non-finite automatic bumps.
#'
#' @param x A baseline value for the parameter in question
#' @param divisor A reasonable fraction of the parameter to try bumping with
#' @param n The number of multiples of \code{divisor} to include in the bump sequence returned
#' @return An empty numeric vector when \code{x} is zero or non-finite, otherwise a descending sequence of fractions of \code{x}
#' @family Greeks
#'
#' @keywords internal
construct_descending_bumps = function(x, divisor=40, n=3) {
  if (!is.finite(x) || x == 0) return(numeric(0))
  base = abs(x) / divisor
  base / seq_len(n)
}

#' Choose reasonable bump sizes to try, or accept user overrides,  when finite-differencing for greeks
#'
#' Choose bumps for calculating for one greek.  If user_bump is NULL then auto-scale from the
#'  representative \code{level} and retry with progressively smaller bumps on failure.
#'  An explicit user bump results in only a single bump if nonzero, or means
#'  "skip this greek" if zero.
#'
#' @inheritDotParams construct_descending_bumps x divisor n
#' @param user_bump User-specified bump size to use (normally left NULL for automatic operation)
#' @param level Baseline level of the model parameter to be bumped
#' @return An empty numeric vector when the resolved \code{level} or \code{user_bump} is zero or non-finite, otherwise a descending sequence of bump sizes
#' @family Greeks
#'
#' @keywords internal
resolve_bumps = function(user_bump, level, ...) {
  if (is.null(user_bump)) return(construct_descending_bumps(level, ...))
  if (!is.finite(user_bump) || user_bump == 0) return(numeric(0))
  abs(user_bump)
}


#' Compute Delta and Gamma greeks by examining a completed pricing grid
#'
#' Since our grid pricing scheme scales its node equity values according to the initial
#' equity price \code{S0}, the entire scheme is almost surely locally linear in \code{S0} (the
#' exception being when a small change in \code{S0} causes some node to cross the early
#' exercise boundary).
#' This \emph{prevents} multiple pricing calls from making an accurate estimation of
#' Gamma entirely, and interferes with the accuracy of Delta as well.
#' Happily, it is easy to estimate Delta and Gamma, with high accuracy, from the
#' pricing grid used by the scheme, so we do that here.
#'
#' We fit a cubic interpolating spline through the grid's \code{(Underlying, price)}
#' nodes for each instrument and read Delta and Gamma directly off the spline as
#' its first and second derivatives at \code{S0}.  Because the spline derivative is
#' analytic, there is no finite-difference step and hence no bump size to choose.
#'
#' @param present_value_grid A grid of already-computed security prices, with final dimension corresponding to instrument names
#' @param S0 The equity price at start-of-computation time \eqn{t}
#' @param instrument_names The names of instruments being analyzed
#' @return A list with \code{delta} and \code{gamma}, each themselves a list with \code{instrument_names} as names and greek values as entries
#' @family Greeks
#'
#' @keywords internal
grid_delta_gamma <- function(present_value_grid, S0, instrument_names)
{
  underlying <- present_value_grid[, "Underlying"]
  # Rows: derivative order 1 (delta) then 2 (gamma); columns: instruments
  dg <- vapply(instrument_names, function(nm) {
    y <- present_value_grid[, nm]
    ok <- is.finite(underlying) & is.finite(y)
    interp <- stats::splinefun(x = underlying[ok], y = y[ok])
    c(interp(S0, deriv = 1), interp(S0, deriv = 2))
  }, FUN.VALUE = numeric(2))
  list(delta = dg[1, ], gamma = dg[2, ])
}

#' Find present value and selected greeks for given set of instruments
#'
#' Price the given instruments using \code{\link{find_present_value}}, and also compute the specified greeks
#' via finite-difference calculations.  The finite differences are higher-order
#' 2-sided estimates where feasible, otherwise they are 1-sided estimates.  Bump
#' sizes are automatically chosen by \code{\link{resolve_bumps}}, unless the corresponding
#' \code{dS}, \code{dr}, \code{dvola} or \code{dhazard} override has been given a
#' non-\code{NULL} value.
#'
#' If interest rates, default intensity, or volatility is zero then the corresponding
#' greek will \emph{not} be computed without a specifically-chosen user override
#' for \code{dr}, \code{dvola} or \code{dhazard}.
#'
#' Bump sizes ultimately used, and central/one-sided estimate status, are available as as \code{attributes()} of the result
#'
#' @inheritParams find_present_value
#' @inheritDotParams find_present_value
#' @param greeks Sequence of names of greeks to be estimated (see
#'  \code{\link{GREEK_NAMES}} for the available values).  Delta and Gamma are
#'  essentially zero-cost to include since they are read directly off the pricing
#'  grid as spline derivatives.  All other greeks require repeated calls to the
#'  pricing engine and multiply computation expense accordingly.
#' @param dvola User-specified change in volatility to use for Vega estimation by
#'  multiple pricing calls.  Usually left NULL, for automatic selection.
#' @param dr User-specified change in risk-free rates to use for Rates DV01 estimation
#'  by multiple pricing calls.  Usually left NULL, for automatic selection.
#' @param dhazard User-specified change in default intensity to use for Credit DV01
#'  estimation by multiple pricing calls.  Usually left NULL, for automatic selection.
#' @return A list, with at least one entry for \code{price} resulting from
#'  \code{find_present_value()}.  For each greek requested, the list will have another entry.
#' @examples
#' S0=241.80
#' varc = variance_cumulation_from_vols(
#'   data.frame(
#'     time=c(0.13, 0.38, 0.63, 0.72, 1.72),
#'     volatility=c(0.48, 0.46, 0.46, 0.45, 0.45)
#'   ))
#' disct_fcn = spot_to_df_fcn(
#'   data.frame(time=c(1, 5, 10, 15),rate=c(0.01, 0.02, 0.03, 0.05))
#' )
#' cb = ConvertibleBond(
#' maturity=2.87, conversion_ratio=2.7788, notional=1000,
#' coupons=data.frame(
#' payment_time=seq(2.8, 0, by=-0.25),
#'     payment_size=1000*0.0025/4),
#'   discount_factor_fcn = disct_fcn
#' )
#' price_and_greeks = find_greeks(
#'   greeks=c("delta", "vega", "credit_dv01"),
#'   S0=S0,
#' instruments=list(CB=cb),
#'   num_time_steps=55,
#'   default_intensity_fcn = function(t, S, ...){0.03 + 0.01 * (S0/S)^1.5},
#'   discount_factor_fcn = disct_fcn,
#'   variance_cumulation_fcn=varc)
#' round(unlist(price_and_greeks),4)
#' @family Greeks
#'
#' @export
find_greeks = function(S0, num_time_steps, instruments,
                       greeks = c("delta", "gamma"),
                       const_volatility=0.5, const_short_rate=0, const_default_intensity=0,
                       discount_factor_fcn = function(T, t, ...){exp(-const_short_rate*(T-t))},
                       default_intensity_fcn = function(t, S, ...){const_default_intensity+0.0*S},
                       variance_cumulation_fcn = function(T, t){const_volatility^2*(T-t)},
                       dvola=NULL, dr=NULL, dhazard=NULL,
                       ...)
{
  greeks = tolower(greeks)
  unknown = setdiff(greeks, GREEK_NAMES)
  if (length(unknown) == 1) {
    stop("Unknown greek requested: ", toString(unknown),
         ".  Valid choices are: ", toString(GREEK_NAMES))
  }
  if (length(unknown) > 1) {
    stop("Unknown greeks requested: ", toString(unknown),
         ".  Valid choices are: ", toString(GREEK_NAMES))
  }

  # Arguments we needn't alter
  arg_dots = list(...)

  # Name the instruments; greeks are reported for every one of them.
  named = instruments
  if (is.null(names(named))) names(named) = vapply(named, function(i) i$name, character(1))
  instrument_names = names(named)
  na_per_instrument = function() stats::setNames(rep(NA_real_, length(instrument_names)), instrument_names)

  # Time horizon used to scale auto bump sizes
  Tmax = if (!is.null(arg_dots$override_Tmax) && is.finite(arg_dots$override_Tmax)) arg_dots$override_Tmax
         else max(vapply(named, function(i) i$maturity, numeric(1)))

  # Common arguments for a single unbumped pricing pass on the shared grid.
  base_args = c(list(S0=S0, num_time_steps=num_time_steps, instruments=named,
                     discount_factor_fcn=discount_factor_fcn,
                     default_intensity_fcn=default_intensity_fcn,
                     variance_cumulation_fcn=variance_cumulation_fcn),
                arg_dots)

  # Price all instruments (optionally overriding a term-structure function for
  #  a bump) on their shared grid
  reprice_with = function(overrides=list()) {
    pv = do.call(find_present_value, utils::modifyList(base_args, overrides))
    vapply(instrument_names, function(nm) {
      x = pv[[nm]]
      if (length(x) != 1) NA_real_ else as.numeric(x)
    }, numeric(1))
  }

  # Baseline prices plus a grid we can get delta, gamma and theta from
  base_grid = do.call(form_present_value_grid, base_args)
  v0 = vapply(instrument_names, function(nm) {
    stats::spline(x=base_grid[, "Underlying"], y=base_grid[, nm], xout=S0)$y
  }, numeric(1))


  results = list(price = v0)
  if (("delta" %in% greeks)||("gamma" %in% greeks)) {
    gdg = grid_delta_gamma(base_grid, S0, instrument_names)
    results$delta = gdg$delta
    results$gamma = gdg$gamma
  }
  if ("theta" %in% greeks) results$theta = na_per_instrument()  # TODO from base grid + adjacent time row

  provenance = list()

  # Compute one greek over its descending bumps, recording the bump used
  #  and each instrument's difference method (central / one-sided / failed).
  run_bump_greek = function(greek, bumps, make_bumped, fcn_arg) {
    if (!(greek %in% greeks)) return(invisible())
    if (length(bumps) == 0) {
      flog.warn("Greek '%s' was requested but no usable bump size was determined (auto-sizing needs a non-zero, finite input level; otherwise pass an explicit positive bump); skipping it.",
                greek, name="ragtop.greeks")
      results[[greek]] <<- na_per_instrument()
      return(invisible())
    }
    reprice = function(b) reprice_with(stats::setNames(list(make_bumped(b)), fcn_arg))
    est = robust_greek(v0, reprice, bumps, instrument_names, greek)
    results[[greek]] <<- est$value
    provenance[[greek]] <<- list(bump=est$bump, method=est$method)
    one_sided = instrument_names[est$method %in% c("forward", "backward")]
    if (length(one_sided) > 0) {
      flog.warn("Greek '%s' computed with a one-sided difference at bump %s for: %s; the other bump direction failed.",
                greek, est$bump, toString(one_sided), name="ragtop.greeks")
    }
    failed = instrument_names[est$method == "failed"]
    if (length(failed) > 0) {
      flog.warn("Greek '%s' could not be computed for: %s; all bump sizes failed.",
                greek, toString(failed), name="ragtop.greeks")
    }
  }

  run_bump_greek("vega",
                 resolve_bumps(dvola, effective_volatility(variance_cumulation_fcn, Tmax)),
                 function(b) bump_variance_cumulation_fcn(variance_cumulation_fcn, b),
                 "variance_cumulation_fcn")
  run_bump_greek("rates_dv01",
                 resolve_bumps(dr, effective_rate(discount_factor_fcn, Tmax)),
                 function(b) bump_discount_factor_fcn(discount_factor_fcn, b),
                 "discount_factor_fcn")
  run_bump_greek("credit_dv01",
                 resolve_bumps(dhazard, effective_hazard(default_intensity_fcn, S0)),
                 function(b) bump_default_intensity_fcn(default_intensity_fcn, b),
                 "default_intensity_fcn")

  # Record bump size and difference method per greek, for reproducibility.
  attr(results, "bump_provenance") = provenance
  results
}

