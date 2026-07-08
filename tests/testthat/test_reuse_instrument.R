flog.threshold(WARN, name="ragtop")
flog.threshold(WARN)
flog.threshold(ERROR, name='ragtop.implicit.setup.width')

pct4 = function(T,t=0) {
  exp(-0.04*(T-t))
}
pct0 = function(T,t=0) {
  exp(-0.0*(T-t))
}


# Repricing a single convertible bond (e.g. for bump-and-reprice greeks) must
#  match pricing fresh objects, because integrate_pde resets per-run caches.
make_reuse_cb = function() ConvertibleBond(maturity=2.87,
                                           conversion_ratio=2.7788,
                                           notional=1000,
                                           recovery_rate=0.4,
                                           coupons=data.frame(payment_time=seq(2.8, 0, by=-0.25),
                                                              payment_size=1000*0.0025/4),
                                           discount_factor_fcn=pct4,
                                           name='ReuseCB')
reuse_base_fresh = find_present_value(S0=252, instruments=list(ReuseCB=make_reuse_cb()),
                                      num_time_steps=25, const_default_intensity=0.1,
                                      const_short_rate=0.04, const_volatility=0.5)$ReuseCB
reuse_bump_fresh = find_present_value(S0=252, instruments=list(ReuseCB=make_reuse_cb()),
                                      num_time_steps=25, const_default_intensity=0.1,
                                      const_short_rate=0.04, const_volatility=0.6)$ReuseCB
reuse_cb = make_reuse_cb()
reuse_base = find_present_value(S0=252, instruments=list(ReuseCB=reuse_cb),
                                num_time_steps=25, const_default_intensity=0.1,
                                const_short_rate=0.04, const_volatility=0.5)$ReuseCB
reuse_bump = find_present_value(S0=252, instruments=list(ReuseCB=reuse_cb),
                                num_time_steps=25, const_default_intensity=0.1,
                                const_short_rate=0.04, const_volatility=0.6)$ReuseCB
test_that("Reusing a convertible bond object matches fresh objects", {
  # Cache reset in integrate_pde makes each valuation self-contained
  expect_equal(reuse_base, reuse_base_fresh)
  expect_equal(reuse_bump, reuse_bump_fresh)
  # Sanity: the vol bump actually changed the price, so this is a real test
  expect_false(isTRUE(all.equal(reuse_base, reuse_bump)))
})

test_that("reset_caches clears per-run state across the inheritance chain", {
  cb = make_reuse_cb()
  cb$last_computed_grid = c(1, 2, 3)
  cb$last_computed_cash = 50
  cb$last_computed_exercise_value = c(9, 9)
  cb$last_computed_exercise_decision = c(TRUE, FALSE)
  cb$last_used_S = c(80, 120)
  cb$last_used_t = 0.5
  cb$reset_caches()
  expect_length(cb$last_computed_grid, 0)
  expect_length(cb$last_computed_cash, 0)
  expect_length(cb$last_computed_exercise_value, 0)
  expect_length(cb$last_computed_exercise_decision, 0)
  expect_length(cb$last_used_S, 0)
  expect_length(cb$last_used_t, 0)
})

make_reuse_coupon = function() CouponBond(maturity=2.87,
                                          notional=1000,
                                          recovery_rate=0.4,
                                          coupons=data.frame(payment_time=seq(2.8, 0, by=-0.25),
                                                             payment_size=1000*0.0025/4),
                                          discount_factor_fcn=pct4,
                                          name='ReuseCoupon')
coupon_base_fresh = find_present_value(S0=252, instruments=list(ReuseCoupon=make_reuse_coupon()),
                                       num_time_steps=25, const_default_intensity=0.1,
                                       const_short_rate=0.04, const_volatility=0.5)$ReuseCoupon
coupon_bump_fresh = find_present_value(S0=252, instruments=list(ReuseCoupon=make_reuse_coupon()),
                                       num_time_steps=25, const_default_intensity=0.1,
                                       const_short_rate=0.08, const_volatility=0.5)$ReuseCoupon
reuse_coupon = make_reuse_coupon()
coupon_base = find_present_value(S0=252, instruments=list(ReuseCoupon=reuse_coupon),
                                 num_time_steps=25, const_default_intensity=0.1,
                                 const_short_rate=0.04, const_volatility=0.5)$ReuseCoupon
coupon_bump = find_present_value(S0=252, instruments=list(ReuseCoupon=reuse_coupon),
                                 num_time_steps=25, const_default_intensity=0.1,
                                 const_short_rate=0.08, const_volatility=0.5)$ReuseCoupon
test_that("Reusing a coupon bond object matches fresh objects", {
  expect_equal(coupon_base, coupon_base_fresh)
  expect_equal(coupon_bump, coupon_bump_fresh)
  # Sanity: the rate bump actually moved the price, so this is a real test
  expect_false(isTRUE(all.equal(coupon_base, coupon_bump)))
})

test_that("reset_caches clears coupon bond per-run state", {
  bond = make_reuse_coupon()
  bond$last_computed_grid = c(1, 2, 3)
  bond$last_computed_cash = 50
  bond$reset_caches()
  expect_length(bond$last_computed_grid, 0)
  expect_length(bond$last_computed_cash, 0)
})
