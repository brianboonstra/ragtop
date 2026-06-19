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

#' Solve the implicit-timestep linear system
#'
#' Solve the linear system \eqn{A x = rhs} arising from one implicit
#' finite-difference timestep, where \eqn{A} is the matrix assembled by
#' \code{\link{construct_tridiagonals}} and \eqn{b} is the vector of
#' price space-grid values from the previously calculated timestep (which
#' corresponds to one step further in our future).
#'
#' @param tridiag_matrix_entries Diagonal, superdiagonal and subdiagonal
#'  of the matrix from the numerical integrator, as a list with elements
#'  \code{sub}, \code{diag} and \code{super} (see
#'  \code{\link{construct_tridiagonals}})
#' @param rhs A vector forming the right-hand side of the linear system,
#'  typically the space grid values from the previously calculated timestep
#' @return The solution vector \eqn{x} of the linear system
#' @param method Character string selecting the linear-system solver.  One of
#'  \code{'lapack'} (use \code{limSolve::Solve.tridiag()}, LAPACK DGTSV),
#'  \code{'banded'} (assemble a sparse banded matrix and solve with the
#'  \pkg{Matrix} package) or \code{'thomas'} (a pure-R Thomas tridiagonal
#'  solve, always available).  Any other value (including the default
#'  \code{'auto'}) selects the first available method in that same logical
#'  order.
#' @details
#' By default (\code{method='auto'}) the solver prefers
#' \code{limSolve::Solve.tridiag()}, which calls LAPACK DGTSV Gaussian
#' elimination with partial pivoting, then falls back to a sparse banded
#' solve from the \pkg{Matrix} package, and finally to a pure-R Thomas
#' algorithm that requires no additional packages.  Passing an explicit
#' \code{method} forces that solver; forcing \code{'lapack'} or
#' \code{'banded'} requires the corresponding package to be installed.
#' @family Implicit Grid Solver
#' @export
pde_matrix_solve = function(tridiag_matrix_entries, rhs, method='auto')
{
  have_lapack = requireNamespace("limSolve", quietly = TRUE)
  have_banded = requireNamespace("Matrix", quietly = TRUE)

  chosen = switch(method,
                  lapack = "lapack",
                  banded = "banded",
                  thomas = "thomas",
                  if (have_lapack) "lapack"
                  else if (have_banded) "banded"
                  else "thomas")

  # rhs may be either a single vector or a matrix whose columns are several
  #  right-hand sides (one per instrument).  We preserve that shape on output.
  rhs_is_matrix = is.matrix(rhs)

  if (chosen == "lapack") {
    if (!have_lapack)
      stop("ragtop::pde_matrix_solve(method='lapack') requires the 'limSolve' package ",
           "to be installed.")
    # limSolve::Solve.tridiag() uses LAPACK DGTSV Gaussian elimination
    #  with partial pivoting rather than the naive tridiagonal algorithm.
    #  It always returns a matrix; coerce to a plain numeric vector when the
    #  right-hand side was a single vector so all three methods agree in shape.
    sol = limSolve::Solve.tridiag(tridiag_matrix_entries$sub,
                                  tridiag_matrix_entries$diag,
                                  tridiag_matrix_entries$super,
                                  rhs)
    if (rhs_is_matrix) as.matrix(sol) else as.numeric(sol)
  } else if (chosen == "banded") {
    if (!have_banded)
      stop("ragtop::pde_matrix_solve(method='banded') requires the 'Matrix' package ",
           "to be installed.")
    n = length(tridiag_matrix_entries$diag)
    A = Matrix::bandSparse(n, n, k = c(-1, 0, 1),
                            diagonals = list(
                              tridiag_matrix_entries$sub,
                              tridiag_matrix_entries$diag,
                              tridiag_matrix_entries$super)
                           )
    sol = Matrix::solve(A, rhs)
    if (rhs_is_matrix) as.matrix(sol) else as.numeric(sol)
  } else {
    # Pure R, Thomas algorithm.  Works column-by-column so it handles both a
    #  single vector and a matrix of right-hand sides.
    n = length(tridiag_matrix_entries$diag)
    aa = tridiag_matrix_entries$sub
    bb = tridiag_matrix_entries$diag
    cc = tridiag_matrix_entries$super # 2-letter var names to avoid namespace trouble

    # The forward elimination of the superdiagonal (g) and the pivots (mm)
    #  depend only on the matrix, so compute them once for all columns
    g = numeric(n)
    mm = numeric(n)
    mm[1] = bb[1]
    g[1] = cc[1] / bb[1]
    for (i in 2:n) {
      mm[i] = bb[i] - aa[i-1] * g[i-1]
      if (i < n) g[i] = cc[i] / mm[i]
    }

    B = if (rhs_is_matrix) rhs else matrix(rhs, ncol = 1L)
    X = matrix(0.0, nrow = n, ncol = ncol(B))
    for (col in seq_len(ncol(B))) {
      r = numeric(n)
      r[1] = B[1, col] / bb[1]
      # Forward
      for (i in 2:n) {
        r[i] = (B[i, col] - aa[i-1] * r[i-1]) / mm[i]
      }
      # Back substitution
      X[n, col] = r[n]
      for (i in (n-1):1) {
        X[i, col] = r[i] - g[i] * X[i+1, col]
      }
    }
    if (rhs_is_matrix) X else as.numeric(X)
  }
}
