## Tests for pde_matrix_solve() and its three solver backends

# Tridiagonal to dense matrix, for independent reference solution via base::solve().
dense_from_tridiag = function(ent) {
  n = length(ent$diag)
  A = diag(ent$diag, n, n)
  for (i in seq_len(n - 1)) {
    A[i + 1, i] = ent$sub[i]
    A[i, i + 1] = ent$super[i]
  }
  A
}

# A diagonally-dominant (hence invertible) tridiagonal system
make_system = function(n, seed = 1) {
  set.seed(seed)
  ent = list(sub   = runif(n - 1, -0.3, -0.05),
             diag  = runif(n,      1.5,  2.0),
             super = runif(n - 1, -0.3, -0.05))
  rhs = rnorm(n)
  list(ent = ent, rhs = rhs)
}

test_that("Each tridiag matrix solver method matches an independent dense solve", {
  s = make_system(40)
  ref = as.numeric(solve(dense_from_tridiag(s$ent), s$rhs))

  skip_if_not_installed("limSolve")
  expect_equal(pde_matrix_solve(s$ent, s$rhs, method = "lapack"), ref)
})

test_that("Banded tridiag matrix solver method matches dense solve", {
  s = make_system(40)
  ref = as.numeric(solve(dense_from_tridiag(s$ent), s$rhs))

  skip_if_not_installed("Matrix")
  expect_equal(pde_matrix_solve(s$ent, s$rhs, method = "banded"), ref)
})

test_that("Thomas tridiag matrix solver method matches dense solve", {
  s = make_system(40)
  ref = as.numeric(solve(dense_from_tridiag(s$ent), s$rhs))

  # The pure-R Thomas solver needs no extra packages
  expect_equal(pde_matrix_solve(s$ent, s$rhs, method = "thomas"), ref)
})

test_that("All available tridiag matrix solver methods agree", {
  s = make_system(75, seed = 7)
  x_thomas = pde_matrix_solve(s$ent, s$rhs, method = "thomas")

  if (requireNamespace("limSolve", quietly = TRUE)) {
    expect_equal(pde_matrix_solve(s$ent, s$rhs, method = "lapack"), x_thomas)
  }
  if (requireNamespace("Matrix", quietly = TRUE)) {
    expect_equal(pde_matrix_solve(s$ent, s$rhs, method = "banded"), x_thomas)
  }
})


test_that("Tridiag matrix solver works on entries from construct_tridiagonals", {
  drift = 0.0011 * (2:21)
  ent = construct_tridiagonals(0.5, 0.05, drift)
  rhs = seq(0.6, 1.0, length.out = length(ent$diag))
  ref = as.numeric(solve(dense_from_tridiag(ent), rhs))

  expect_equal(pde_matrix_solve(ent, rhs, method = "thomas"), ref)
  if (requireNamespace("limSolve", quietly = TRUE)) {
    expect_equal(pde_matrix_solve(ent, rhs, method = "lapack"), ref)
  }
  if (requireNamespace("Matrix", quietly = TRUE)) {
    expect_equal(pde_matrix_solve(ent, rhs, method = "banded"), ref)
  }
})

test_that("Tridiag matrix solver methods handle a matrix of right-hand sides (multiple instruments)", {
  s = make_system(35, seed = 23)
  A = dense_from_tridiag(s$ent)
  # Three distinct right-hand sides as columns of a matrix
  RHS = cbind(s$rhs, rev(s$rhs), seq_along(s$rhs) / length(s$rhs))
  ref = solve(A, RHS)

  for (m in c("lapack", "banded", "thomas")) {
    if (m == "lapack" && !requireNamespace("limSolve", quietly = TRUE)) next
    if (m == "banded" && !requireNamespace("Matrix",   quietly = TRUE)) next
    sol = pde_matrix_solve(s$ent, RHS, method = m)
    expect_true(is.matrix(sol), info = paste("method =", m))
    expect_equal(dim(sol), dim(RHS), info = paste("method =", m))
    expect_equal(unname(sol), unname(ref), info = paste("method =", m))
  }
})

test_that("Tridiag matrix solver methods on a vector right-hand side returns a plain numeric vector", {
  s = make_system(20, seed = 31)
  for (m in c("lapack", "banded", "thomas")) {
    if (m == "lapack" && !requireNamespace("limSolve", quietly = TRUE)) next
    if (m == "banded" && !requireNamespace("Matrix",   quietly = TRUE)) next
    sol = pde_matrix_solve(s$ent, s$rhs, method = m)
    expect_false(is.matrix(sol), info = paste("method =", m))
    expect_type(sol, "double")
  }
})

test_that("Forcing an unavailable tridiag matrix solver method raises an informative error", {
  s = make_system(10)
  # These branches only execute in environments lacking the package; when the
  #  package is installed the forced method simply succeeds.
  if (!requireNamespace("limSolve", quietly = TRUE)) {
    expect_error(pde_matrix_solve(s$ent, s$rhs, method = "lapack"), "limSolve")
  } else {
    expect_silent(pde_matrix_solve(s$ent, s$rhs, method = "lapack"))
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    expect_error(pde_matrix_solve(s$ent, s$rhs, method = "banded"), "Matrix")
  } else {
    expect_silent(pde_matrix_solve(s$ent, s$rhs, method = "banded"))
  }
})
