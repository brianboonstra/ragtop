library(ragtop)
library(futile.logger)
context("Payoffs")

flog.threshold(WARN, name="ragtop")
flog.threshold(WARN)

pct3 = function(T,t=0) {
  exp(-0.03*(T-t))
}

amer_put = AmericanOption$new(maturity=5, strike=11, callput=-1)
test_that("American Put Payoff", {
  expect_equal(amer_put$optionality_fcn(c(1, 15, 16, 17, -0.1), c(2,3,13,14,19), 1)[1], 9)
  expect_equal(amer_put$optionality_fcn(c(1, 15, 16, 17, -0.1), c(2,3,13,14,19), 1)[2], 15)
  expect_equal(amer_put$optionality_fcn(c(1, 15, 16, 17, -0.1), c(2,3,13,14,19), 1)[4], 17)
  expect_equal(amer_put$optionality_fcn(c(1, 15, 16, 17, -0.1), c(2,3,13,14,19), 1)[5], 0)
})


coups = data.frame(payment_time=c(0.2, 0.7, 1.2, 1.7),
                   payment_size=c(30,30,30,30))
conv_bond = ConvertibleBond$new(coupons=coups,
                                conversion_ratio=0,
                                maturity=10,
                                notional=0,
                                name='CB')
test_that("Convertible Bond Coupon Discounting", {
  expect_equal(conv_bond$optionality_fcn(c(0,0), c(0,0.1), 0.9, discount_factor_fctn=pct3)[1],
               30 * (exp((0.9-0.7)*0.03) + exp((0.9-0.2)*0.03)), tolerance=1e-3)
})
