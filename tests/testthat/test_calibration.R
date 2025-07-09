library(ragtop)
library(futile.logger)
context("Calibration")

flog.threshold(WARN, name="ragtop")
flog.threshold(ERROR, name='ragtop.implicit.setup.width')
flog.threshold(WARN)

test_that("Black-Scholes volatility calibration", {
  expect_equal(equivalent_jump_vola_to_bs(0.75, 20, const_default_intensity=0.0, max.iter=20),
               0.75,
               tolerance=1.e-5)
  expect_equal(equivalent_bs_vola_to_jump(0.75, 20, const_default_intensity=0.0, max.iter=20),
               0.75,
               tolerance=1.e-5)
  expect_equal(equivalent_jump_vola_to_bs(0.75, 20, const_default_intensity=0.07, max.iter=20),
               0.5649847,
               tolerance=1.e-5)
})
equivalent_jump_vola_to_bs
aiv = american_implied_volatility(25,CALL,S0=100,K=100,time=2.2, const_short_rate=0.03)
bsiv = implied_volatility(25, CALL, S0=100, K=100, r=0.03, time=2.2)
test_that("American implied vola same as european when no early exercise", {
  expect_equal(bsiv, aiv, tolerance=1.e-3)
})

amer_price = american(PUT, 33, 30, 0.77, const_short_rate=0.02, const_default_intensity=0.05, const_volatility=0.55, num_time_steps=25)
aiv = american_implied_volatility(amer_price,PUT,S0=33,K=30,time=0.77, const_short_rate=0.02, const_default_intensity=0.05, num_time_steps=25)[[1]]
test_that(paste("American implied vola round trip: known vola -> price -> implied vola =", aiv,"= 0.55"), {
  expect_equal(0.55, aiv, tolerance=4.e-4)
})
