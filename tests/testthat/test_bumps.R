flog.threshold(WARN, name="ragtop")
flog.threshold(WARN)

test_that("Discount factor bump applies a parallel rate shift", {
  df = function(T, t, ...) exp(-0.03 * (T - t))
  up = bump_discount_factor_fcn(df, 0.0001)
  down = bump_discount_factor_fcn(df, -0.0001)
  # Bumping the rate up by dr multiplies the discount factor by exp(-dr*(T-t))
  expect_equal(up(2, 0), df(2, 0) * exp(-0.0001 * 2))
  expect_equal(down(2, 0), df(2, 0) * exp(0.0001 * 2))
  # A rate increase lowers the discount factor
  expect_lt(up(2, 0), df(2, 0))
  expect_gt(down(2, 0), df(2, 0))
})

test_that("Discount factor bump is vectorized and respects t arg", {
  df = function(T, t, ...) exp(-0.03 * (T - t))
  up = bump_discount_factor_fcn(df, 0.0001)
  Ts = c(1, 2, 5)
  expect_equal(up(Ts, 0), df(Ts, 0) * exp(-0.0001 * Ts))
  expect_equal(up(2, 1), df(2, 1) * exp(-0.0001 * 1))
})

test_that("Default intensity bump shifts hazard additively", {
  h = function(t, S, ...) 0.02 + 0.0 * S
  up = bump_default_intensity_fcn(h, 0.0002)
  down = bump_default_intensity_fcn(h, -0.0002)
  expect_equal(up(1, 50), 0.0202)
  expect_equal(down(1, 50), 0.0198)
  # Shift is applied across a vector of stock prices
  expect_equal(up(1, c(10, 50, 90)), rep(0.0202, 3))
})

test_that("Variance cumulation bump shifts effective volatility", {
  vc = function(T, t) 0.25^2 * (T - t)
  up = bump_variance_cumulation_fcn(vc, 0.01)
  down = bump_variance_cumulation_fcn(vc, -0.01)
  # sigma = sqrt(var/dt) = 0.25, bumped to 0.26/0.24, recumulated over dt=2
  expect_equal(up(2, 0), (0.26)^2 * 2)
  expect_equal(down(2, 0), (0.24)^2 * 2)
  # Bumping vol up increases cumulated variance
  expect_gt(up(2, 0), vc(2, 0))
  expect_lt(down(2, 0), vc(2, 0))
})

test_that("Variance cumulation bump guards the zero-interval case", {
  vc = function(T, t) 0.25^2 * (T - t)
  up = bump_variance_cumulation_fcn(vc, 0.01)
  # dt == 0 must not divide by zero; cumulated variance over no time is zero
  expect_equal(up(1, 1), 0)
})

test_that("Bump wrappers force their argument (no lazy self-reference)", {
  df = function(T, t, ...) exp(-0.03 * (T - t))
  wrapped = bump_discount_factor_fcn(df, 0.0001)
  # Rebinding df after wrapping must not change the captured function
  df = function(T, t, ...) 999
  expect_equal(wrapped(2, 0), exp(-0.03 * 2) * exp(-0.0001 * 2))
})

test_that("Zero bump leaves each function unchanged", {
  df = function(T, t, ...) exp(-0.03 * (T - t))
  h = function(t, S, ...) 0.02 + 0.0 * S
  vc = function(T, t) 0.25^2 * (T - t)
  expect_equal(bump_discount_factor_fcn(df, 0)(2, 0), df(2, 0))
  expect_equal(bump_default_intensity_fcn(h, 0)(1, 50), h(1, 50))
  expect_equal(bump_variance_cumulation_fcn(vc, 0)(2, 0), vc(2, 0))
})
