library(ragtop)
library(futile.logger)
context("Calibration")

flog.threshold(WARN, name="ragtop")
flog.threshold(WARN)

vc = variance_cumulation_from_vols(data.frame(time=c(0.1,2,3), volatility=c(0.2,0.5,1.2)))
test_that("Variance cumulation", {
  expect_equal(sqrt(vc(0.0999999999)/0.1), 0.2, tolerance=1.e-4)
  expect_equal(sqrt(vc(0.1000001)/0.1), 0.2, tolerance=1.e-4)
  expect_equal(sqrt(vc(100)/100), 1.93613, tolerance=1.e-3)
})
