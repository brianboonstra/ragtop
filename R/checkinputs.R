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

#' Check that a discount factor function is internally consistent
#'
#' Ensure discount factor function satisfies the
#' time-reversal identity \deqn{df(t,T) = 1/df(T,t)}.  Emits a warning,
#' whenever that identity is violated, returns  TRUE if no problems are found.
#'
#' Discount factors that have underflowed to zero or overflowed to
#' infinity carry no usable information about the reversed factor, so
#' such pairs are silently ignored rather than flagged.
#'
#' @param discount_factor_fcn A function for computing present values, with
#'   arguments \code{T}, \code{t}, returning the discount factor from
#'   time \code{T} back to time \code{t}
#' @param time_pts A vector of time points at which to test the
#'   consistency of \code{discount_factor_fcn}
#' @param tolerance Relative tolerance within which the reversal
#'   identity is considered to hold
#' @return \code{TRUE} if the discount factor function appears
#'   consistent at every tested pair, otherwise \code{FALSE}
#' @family Input Checking
#' @export
check_discount_factor_fcn = function(discount_factor_fcn, time_pts,
                                     tolerance=1e-6)
{
  time_pts = sort(unique(time_pts))
  num_pts = length(time_pts)
  if (num_pts < 2) {
    return(TRUE)
  }
  # Test every point against the central point only, keeping the work O(N).
  # Evaluate pointwise so we tolerate scalar-only discount factor functions.
  mid = time_pts[ceiling(num_pts / 2)]
  df_fwd = vapply(time_pts, function(tp) discount_factor_fcn(tp, mid), numeric(1))   # discount from each point back to mid
  df_rev = vapply(time_pts, function(tp) discount_factor_fcn(mid, tp), numeric(1))   # discount from mid back to each point
  products = df_fwd * df_rev
  # Underflowed or overflowed factors hold no usable information, so ignore them
  usable = is.finite(df_fwd) & is.finite(df_rev) & (df_fwd != 0) & (df_rev != 0)
  bad = usable & (abs(products - 1.0) > tolerance)
  if (any(bad)) {
    flog.warn("Discount factor function violates df(t,T)=1/df(T,t) at mid=%s for t/T=%s: df(T,t)=%s, df(t,T)=%s, product=%s",
              mid, toString(time_pts[bad]), toString(df_fwd[bad]),
              toString(df_rev[bad]), toString(products[bad]),
              name='ragtop.checkinputs.check_discount_factor_fcn')
  }
  !any(bad)
}

#' Check that a variance cumulation function is nondecreasing
#'
#' Total  variance accumulated from a fixed origin can never
#' decrease as the horizon moves forward in time.  This routine
#' evaluates the supplied \code{variance_cumulation_fcn} on the ordered
#' \code{time_pts}.  Emits a warning whenever the cumulated variance fails to
#' be nondecreasing.
#'
#' @param variance_cumulation_fcn A function for computing total stock
#'   variance, with arguments \code{T}, \code{t}, returning the variance
#'   cumulated from time \code{t} to time \code{T}
#' @param time_pts A vector of time points at which to test that
#'   \code{variance_cumulation_fcn} is nondecreasing
#' @return \code{TRUE} if the cumulated variance is nondecreasing across
#'   \code{time_pts}, otherwise \code{FALSE}
#' @family Input Checking
#' @export
check_variance_cumulation_fcn = function(variance_cumulation_fcn, time_pts)
{
  time_pts = sort(unique(time_pts))
  num_pts = length(time_pts)
  if (num_pts < 2) {
    return(TRUE)
  }
  # Variance cumulated from the earliest point out to each time point.
  # Evaluate pointwise so we tolerate scalar-only cumulation functions.
  origin = time_pts[1]
  cumulated = vapply(time_pts, function(tp) variance_cumulation_fcn(tp, origin), numeric(1))
  decreases = diff(cumulated) < 0
  if (any(decreases)) {
    flog.warn("Variance cumulation function is decreasing from origin=%s across time points %s with cumulated variances %s",
              origin, toString(time_pts), toString(cumulated),
              name='ragtop.checkinputs.check_variance_cumulation_fcn')
  }
  !any(decreases)
}

#' Wrap a default intensity function so it warns at most once about negatives
#'
#' The default intensity (hazard rate) must be nonnegative.  A solver
#' such as \code{\link{integrate_pde}} calls the intensity function many
#' times, on various stock-price grids and at various times, so a naive
#' check would emit a flood of redundant warnings.  This routine returns
#' a function that behaves exactly like \code{default_intensity_fcn} but
#' latches a flag in its enclosing environment so that, across all of its
#' invocations, at most one warning about negative intensities is
#' emitted in the usual \code{futile.logger} manner.
#'
#' Create a fresh wrapper for each top-level pricing call so that the
#' "already warned" state is scoped to that call and resets on the next.
#'
#' @param default_intensity_fcn A function for computing default
#'   intensity, with arguments \code{t}, \code{S}
#' @return A function with the same interface as
#'   \code{default_intensity_fcn} that warns at most once about negative
#'   default intensities
#' @family Input Checking
#' @export
warn_once_negative_default_intensity = function(default_intensity_fcn)
{
  # Force the argument now so the returned closure captures the original
  #  function rather than a promise that could resolve to itself when the
  #  caller rebinds default_intensity_fcn to this wrapper
  force(default_intensity_fcn)
  warned = FALSE
  function(t, S, ...) {
    intensities = default_intensity_fcn(t, S, ...)
    if (!warned && any(intensities < 0, na.rm=TRUE)) {
      flog.warn("Default intensity function returned negative values, e.g. %s at t=%s.  Hazard rates should be nonnegative.",
                min(intensities, na.rm=TRUE), t,
                name='ragtop.checkinputs.warn_once_negative_default_intensity')
      warned <<- TRUE
    }
    intensities
  }
}

#' Check that a survival probability function is well-behaved
#'
#' Survival probabilities accumulated from a fixed origin can never
#' increase as the horizon moves forward in time, and must always lie
#' in the closed interval \eqn{[0,1]}.  This routine evaluates the
#' supplied \code{survival_probability_fcn} on the ordered
#' \code{time_pts}.  Emits a warning whenever the probabilities are
#' increasing or fall outside \eqn{[0,1]}.
#'
#' @param survival_probability_fcn A function for computing survival
#'   probabilities, with arguments \code{T}, \code{t}, returning the
#'   probability of survival from time \code{t} to time \code{T}
#' @param time_pts A vector of time points at which to test the
#'   behavior of \code{survival_probability_fcn}
#' @return \code{TRUE} if the survival probabilities are nonincreasing
#'   and within \eqn{[0,1]} across \code{time_pts}, otherwise \code{FALSE}
#' @family Input Checking
#' @export
check_survival_probability_fcn = function(survival_probability_fcn, time_pts)
{
  time_pts = sort(unique(time_pts))
  num_pts = length(time_pts)
  if (num_pts < 2) {
    return(TRUE)
  }
  # Survival probability from the earliest point out to each time point.
  # Evaluate pointwise so we tolerate scalar-only survival probability functions.
  origin = time_pts[1]
  probs = vapply(time_pts, function(tp) survival_probability_fcn(tp, origin), numeric(1))
  out_of_range = (probs < 0) | (probs > 1)
  increases = diff(probs) > 0
  ok = !any(out_of_range) & !any(increases)
  if (!ok) {
    flog.warn("Survival probability function is ill-behaved from origin=%s across time points %s with survival probabilities %s",
              origin, toString(time_pts), toString(probs),
              name='ragtop.checkinputs.check_survival_probability_fcn')
  }
  ok
}
