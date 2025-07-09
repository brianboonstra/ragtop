library(ragtop)
library(futile.logger)
context("Equity Options")

flog.threshold(WARN)
flog.threshold(WARN, name="ragtop")
flog.threshold(ERROR, name='ragtop.implicit.setup.width')

pct2 = function(T,t) {exp(-0.02*(T-t))}

euro_zero_strike_call = EuropeanOption(maturity=3.53, strike=0, callput=1, name="CallStrike0")
test_that("European options equivalent to equity, trivial cases", {
  expect_equal(100,
               find_present_value(S0=100,instruments=list(call=euro_zero_strike_call),
                                  num_time_steps=3, const_volatility = 0.00001, const_short_rate=0)$call,
               tolerance=1.e-4)
  expect_equal(100,
               find_present_value(S0=100,instruments=list(call=euro_zero_strike_call),
                                  num_time_steps=3, const_volatility = 0.00001, const_short_rate=0.03434)$call,
               tolerance=1.e-4)
})

test_that("European options equivalent to equity, nontrivial default intensity", {
  expect_equal(100,
               find_present_value(S0=100, instruments=list(call=euro_zero_strike_call),
                                  num_time_steps=200, const_default_intensity=0.07,
                                  const_volatility = 0.4, const_short_rate=0.05, std_devs_width=4)$call,
               tolerance=1.e-1)
})

test_that("European options equivalent to equity, nontrivial volatility", {
  expect_equal(34.4,
               find_present_value(S0=34.4,instruments=list(call=euro_zero_strike_call),
                                  num_time_steps=200, const_volatility = 0.44,
                                  const_short_rate=0, std_devs_width = 4)$call,
               tolerance=1.e-1)
  expect_equal(222.2,
               find_present_value(S0=222.2,instruments=list(call=euro_zero_strike_call),
                                  num_time_steps=200, const_volatility = 0.44,
                                  const_short_rate=0.033, std_devs_width=4)$call,
               tolerance=1.e-1)
})

euro_simple_call = EuropeanOption$new(maturity=1, strike=100, callput=1)
euro_put = EuropeanOption$new(maturity=1, strike=90, callput=-1, discount_factor_fcn=function(T,t){1})
grid_european_prices = find_present_value(S0=100, instruments=list(call=euro_simple_call, put=euro_put),
                                          const_default_intensity = 0.07, num_time_steps=50, std_devs_width=3)
exact_european_prices = list(call=blackscholes(callput=1, S0=100, K=100, r=0, default_intensity =0.07, time=1, vola=0.5)$Price,
                             put=blackscholes(callput=-1, S0=100, K=90, r=0, default_intensity =0.07, time=1, vola=0.5)$Price)
test_that("European options correctly priced", {
  expect_equal(exact_european_prices, grid_european_prices, tolerance=1.e-2)
})

a_proportional_div = data.frame(time=0.5, fixed=0, proportional=20)
grid_european_1div_prices = find_present_value(S0=100, instruments=list(call=euro_simple_call, put=euro_put),
                                               const_default_intensity = 0.07, num_time_steps=50, std_devs_width=3,
                                               dividends=a_proportional_div)
exact_european_1div_prices = list(call=blackscholes(callput=1, S0=100, K=100, r=0,
                                                    default_intensity =0.07, time=1, vola=0.5,
                                                    dividends=a_proportional_div)$Price,
                                  put=blackscholes(callput=-1, S0=100, K=90, r=0,
                                                   default_intensity =0.07, time=1, vola=0.5,
                                                   dividends=a_proportional_div)$Price)
test_that("European options with a stock dividend correctly priced", {
  expect_equal(exact_european_1div_prices, grid_european_1div_prices, tolerance=1.e-1)
})


amer_put_price_20k_steps = 11.6570723  # 20,000 steps in a Leisen-Reimer tree, good to about 1 part in 10^5, 0.001 in this case.  The early exercise premium is about 1.626.
amer_put = AmericanOption(maturity=1, strike=110, callput=-1)
grid_amer_price = find_present_value(S0=100, instruments=list(amer_put=amer_put),
                                     const_short_rate = 0.06,
                                     const_volatility = 0.20,
                                     num_time_steps=200, std_devs_width=5)
test_that("American options correctly priced", {
  expect_equal(amer_put_price_20k_steps, grid_amer_price$amer_put, tolerance=1.e-1)
})
cv_amer_price = as.numeric(american(PUT, S0=100, K=amer_put$strike,
                                    time=amer_put$maturity,
                                    const_short_rate = 0.06, const_volatility=0.20,
                                    num_time_steps=200))
test_that("American options in control variate scheme", {
  expect_equal(amer_put_price_20k_steps, cv_amer_price, tolerance=2.e-2)
})

long_term_ITM_put = EuropeanOption(maturity=3.53, strike=200, callput=PUT,
                                  discount_factor_fcn=pct2, name='Put200')
default_intensity_fcn= function(t, S, ...){h=0.05;p=1;0.95*h+0.05*h*(100/S)^p}
price = find_present_value(S0=100, num_time_steps=250, instruments=list(p200=long_term_ITM_put),
                         const_volatility=0.5, discount_factor_fcn=pct2,
                         default_intensity_fcn=default_intensity_fcn ,
                         structure_constant=2.0,
                         std_devs_width=3.0)$p200
test_that("A power function works", {
  expect_equal(109.2, price, tolerance=0.1)
})
