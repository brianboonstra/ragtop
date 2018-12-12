library(ragtop)
library(futile.logger)
context("Convertibles")

flog.threshold(WARN, name="ragtop")
flog.threshold(WARN)
flog.threshold(ERROR, name='ragtop.implicit.setup.width')

pct4 = function(T,t=0) {
  exp(-0.04*(T-t))
}
pct0 = function(T,t=0) {
  exp(-0.0*(T-t))
}

zccb0 = ConvertibleBond(conversion_ratio=0, maturity=1.3, notional=100,
                        discount_factor_fcn=pct0, name='ConvertibleBondZeroCoup')
zcb = ZeroCouponBond(maturity=zccb0$maturity, notional=zccb0$notional,
                     discount_factor_fcn=pct0, name='EquivalentStraightZCB')
gvs_zc_const_h = form_present_value_grid(100, 2,
                                      instruments=list(zcb=zcb, convertible_bond=zccb0),
                                      const_default_intensity=0.06,
                                      const_short_rate = 0)
gvs_zc_const_r = form_present_value_grid(100, 2,
                                         instruments=list(zcb=zcb, convertible_bond=zccb0),
                                         const_short_rate=0.11)
gvs_zc_const_rh = form_present_value_grid(100, 2,
                                          instruments=list(zcb=zcb, convertible_bond=zccb0),
                                          const_default_intensity=0.05,
                                          const_short_rate=0.11)
differing_zc_const_h = gvs_zc_const_h[,'zcb'] - gvs_zc_const_h[,'convertible_bond']
differing_zc_const_r = gvs_zc_const_r[,'zcb']  - gvs_zc_const_r[,'convertible_bond']
differing_zc_const_rh = gvs_zc_const_rh[,'zcb']  - gvs_zc_const_rh[,'convertible_bond']
test_that("Straight Zero Coupon Bond Equivalents", {
  expect_equal(max(abs(differing_zc_const_h)), 0.0, tolerance=1.e-6)
  expect_equal(max(abs(differing_zc_const_r)), 0.0, tolerance=1.e-6)
  expect_equal(max(abs(differing_zc_const_rh)), 0.0, tolerance=1.e-6)
})


cb0 = ConvertibleBond(conversion_ratio=0, maturity=1, notional=100,
                      coupons=data.frame(payment_time=0.875, payment_size=1),
                      discount_factor_fcn=pct0, name='ConvertibleBond')
split_coup = ZeroCouponBond(maturity=0.875, notional=1, discount_factor_fcn=pct0)
no_coup = ZeroCouponBond(maturity=1, notional=100)
gvs_0_rh = form_present_value_grid(100, 2,
                                      instruments=list(zcb=no_coup, coupon=split_coup, convertible_bond=cb0),
                                      const_default_intensity=0,
                                      const_short_rate = 0)
gvs_const_h = form_present_value_grid(100, 2,
                                      instruments=list(zcb=no_coup, coupon=split_coup, convertible_bond=cb0),
                                      const_default_intensity=0.05862091,
                                      const_short_rate = 0)
gvs_const_r = form_present_value_grid(100, 2,
                                      instruments=list(zcb=no_coup, coupon=split_coup, convertible_bond=cb0),
                                      const_short_rate=0.11)
gvs_const_rh = form_present_value_grid(100, 2,
                                       instruments=list(zcb=no_coup, coupon=split_coup, convertible_bond=cb0),
                                       const_default_intensity=0.05,
                                       const_short_rate=0.11)
differing_const_h = gvs_const_h[,'zcb'] + gvs_const_h[,'coupon'] - gvs_const_h[,'convertible_bond']
differing_const_r = gvs_const_r[,'zcb'] + gvs_const_r[,'coupon'] - gvs_const_r[,'convertible_bond']
differing_0_rh = gvs_0_rh[,'zcb'] + gvs_0_rh[,'coupon'] - gvs_0_rh[,'convertible_bond']
differing_const_rh = gvs_const_rh[,'zcb'] + gvs_const_rh[,'coupon'] - gvs_const_rh[,'convertible_bond']
test_that("Straight Coupon Bond Equivalents", {
  expect_equal(max(abs(differing_const_h)), 0.0, tolerance=1.e-6)
  expect_equal(max(abs(differing_const_r)), 0.0, tolerance=1.e-6)
  expect_equal(max(abs(differing_const_rh)), 0.0, tolerance=1.e-6)
})


cbeq = ConvertibleBond(conversion_ratio=1, maturity=1.5, notional=0,
                       discount_factor_fcn=pct0, name='ConvertibleBond')
euro_zero_strike_call = EuropeanOption(maturity=cbeq$maturity, strike=0, callput=1, name="CallStrike0")
S0 = 100
equity_equiv_prices = find_present_value(S0=S0,
                                         instruments=list(call=euro_zero_strike_call, convertible_bond=cbeq),
                                         num_time_steps=200, const_default_intensity=0.07,
                                         const_volatility = 0.4, const_short_rate=0.03, std_devs_width=4)
test_that("Trivial equity-equivalent", {
  expect_equal(S0, equity_equiv_prices$call, tolerance=1.e-1)
  expect_equal(S0, equity_equiv_prices$convertible_bond, tolerance=1.e-1)
})




cbopt = ConvertibleBond(conversion_ratio=5./3., maturity=1.5, notional=100,
                       discount_factor_fcn=pct0, name='ConvertibleBond')
S0 = 0.95 * cbopt$notional/cbopt$conversion_ratio
euro_simple_call = EuropeanOption$new(maturity=cbopt$maturity,
                                      strike=cbopt$notional/cbopt$conversion_ratio,
                                      callput=1, name="EquivCall")
opt_equiv_prices = find_present_value(S0=S0, instruments=list(call=euro_simple_call, convertible_bond=cbopt),
                                          num_time_steps=200, const_default_intensity=0.0,
                                          const_volatility = 0.4, const_short_rate=0.0, std_devs_width=4)
test_that("Option-equivalent convertible", {
  expect_equal(cbopt$conversion_ratio * opt_equiv_prices$call+cbopt$notional, opt_equiv_prices$convertible_bond, tolerance=1.e-2)
})




cb_zero_recov = ConvertibleBond(maturity=2.87,
                     conversion_ratio=2.7788,
                     notional=1000,
                     coupons=data.frame(payment_time=seq(2.8,0, by=-0.25),
                                        payment_size=1000*0.0025/4),
                     discount_factor_fcn = pct4,
                     name='NoRecov');
cb_40_recov = ConvertibleBond(maturity=2.87,
                      conversion_ratio=2.7788,
                      notional=1000,
                      recovery_rate=0.4,
                      coupons=data.frame(payment_time=seq(2.8,0, by=-0.25),
                                         payment_size=1000*0.0025/4),
                      discount_factor_fcn = pct4,
                      name='StdRecov');
vs_with_and_without = find_present_value(S0=252,
                                         instruments=list(NoRecov=cb_zero_recov, StdRecov=cb_40_recov),
                                         num_time_steps=25,
                                         const_default_intensity=0.1,
                                         const_short_rate = 0.04,
                                         const_volatility=0.5)
test_that("Recovery rates matter", {
  expect_gt(vs_with_and_without$StdRecov, 1.05 * vs_with_and_without$NoRecov)
})
