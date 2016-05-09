library(ragtop)
library(futile.logger)
context("Calibration")

flog.threshold(WARN, name="ragtop")
flog.threshold(WARN)

test_that("Black-Scholes volatility calibration", {
  expect_equal(calibrate_defaultable_volatility(0.75, 20, const_default_intensity=0.0, max.iter=20),
               0.75,
               tolerance=1.e-5)
  expect_equal(calibrate_defaultable_volatility(0.75, 20, const_default_intensity=0.07, max.iter=20),
               0.5649847,
               tolerance=1.e-5)
})

