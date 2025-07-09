library(ragtop)
library(futile.logger)
context("Black-Scholes Modeling")

flog.threshold(WARN, name="ragtop")
flog.threshold(WARN)

test_that("Black-Scholes reproduces equity", {
  expect_equal(blackscholes(callput=1, S0=100, K=0, r=0,
                            default_intensity=0.07,
                            time=4.44, vola=0.5)$Price,
               100,
               tolerance=1.e-6)
  expect_equal(blackscholes(callput=1, S0=100, K=0, r=0.05,
                            default_intensity=0.07,
                            time=4.44, vola=0.5)$Price,
               100,
               tolerance=1.e-6)
})

a_proportional_div = data.frame(time=0.5, fixed=0, proportional=20)
test_that("Black-Scholes pricing", {
  expect_equal(blackscholes(callput=-1, S0=100, K=90, r=0,
                            default_intensity=0, time=1, vola=0.5),
               list(Price = 14.160002955196, Delta = -0.322499382043183, Vega = 35.8771182787839),
               tolerance=1.e-6)
  expect_equal(blackscholes(callput=-1, S0=100, K=90, r=0,
                            default_intensity=0.07, time=1, vola=0.5),
               list(Price = 17.2748000186925, Delta = -0.274012904361218, Vega = 33.3080388264046),
               tolerance=1.e-6)
  expect_equal(blackscholes(callput=-1, S0=100, K=90, r=0,
                            default_intensity=0.07, time=1, vola=0.5, dividends=a_proportional_div),
               list(Price = 24.2940532694562, Delta = -0.438633801658832, Vega = 31.5370543258318),
               tolerance=1.e-6)
  expect_equal(blackscholes(callput=-1, S0=100, K=90, r=0.02,
                            default_intensity=0.07, time=1.8, vola=0.45, dividends=a_proportional_div),
               list(Price = 27.0951778269375, Delta = -0.352261169093692, Vega = 39.9472917093042),
               tolerance=1.e-6)
})

test_that("Black-Scholes with term structures", {
  expect_equal(blackscholes(callput=-1, S0=100, K=90, r=0.03,
                            default_intensity=0.07, time=1, vola=0.45),
               black_scholes_on_term_structures(callput=-1, S0=100, K=90, time=1,
                                                discount_factor_fcn = function(T, t, ...) {
                                                  exp(-0.03 * (T - t))
                                                },
                                                survival_probability_fcn = function(T, t, ...) {
                                                  exp(-0.07 * (T - t))
                                                },
                                                variance_cumulation_fcn = function(T, t) {
                                                  0.45 ^ 2 * (T - t)
                                                }),
               tolerance=1.e-6)
})
