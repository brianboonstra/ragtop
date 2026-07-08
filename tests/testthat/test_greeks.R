flog.threshold(WARN, name = "ragtop")
flog.threshold(WARN)
flog.threshold(ERROR, name = 'ragtop.implicit.setup.width')


N_Steps = 100
dS = 0.1
dVola = 0.005
dr = 0.0001
dh = 0.0001


pct4 = function(T, t = 0) {
  exp(-0.04 * (T - t))
}
pct0 = function(T, t = 0) {
  exp(-0.0 * (T - t))
}


cbeq = ConvertibleBond(
  conversion_ratio = 1,
  maturity = 1.5,
  notional = 0,
  discount_factor_fcn = pct0,
  name = 'ConvertibleBond'
)


cbopt = ConvertibleBond(
  conversion_ratio = 5. / 3.,
  maturity = 1.5,
  notional = 100,
  discount_factor_fcn = pct0,
  name = 'OptionEquiv'
)

coupon_bond = CouponBond(
  maturity = 2.87,
  notional = 1000,
  recovery_rate = 0.4,
  coupons = data.frame(
    payment_time = seq(2.8, 0, by = -0.25),
    payment_size = 1000 * 0.0025 /
      4
  ),
  discount_factor_fcn = pct4,
  name = 'StraightBond'
)

cb_40_recov = ConvertibleBond(
  maturity = 2.87,
  conversion_ratio = 2.7788,
  notional = 1000,
  recovery_rate = 0.4,
  coupons = data.frame(
    payment_time = seq(2.8, 0, by = -0.25),
    payment_size = 1000 * 0.0025 /
      4
  ),
  discount_factor_fcn = pct4,
  name = 'NontrivCB'
)



package_greeks = find_greeks(
  greeks=GREEK_NAMES,
  S0 = 252,
  instruments = list(
    EquityEquiv = cbeq,
    OptionEquiv = cbopt,
    NontrivCB = cb_40_recov,
    StraightBond = coupon_bond
  ),
  num_time_steps = N_Steps,
  const_default_intensity = 0.1,
  const_short_rate = 0.04,
  const_volatility = 0.5
)




test_that("Greek computations have the right baseline price", {
  skip_on_cran()
  baseline = unlist(
    find_present_value(
      S0 = 252,
      instruments = list(
        EquityEquiv = cbeq,
        OptionEquiv = cbopt,
        NontrivCB = cb_40_recov,
        StraightBond = coupon_bond
      ),
      num_time_steps = N_Steps,
      const_default_intensity = 0.1,
      const_short_rate = 0.04,
      const_volatility = 0.5
    )
  )


  for (nm in names(baseline)) {
    expect_equal(as.double(package_greeks$price[nm]),
                 as.double(baseline[nm]),
                 info = paste("Instrument:", nm)
    )
  }
})



test_that("Vega computations are working", {
  up  = unlist(find_present_value(
    S0 = 252,
    instruments = list(
      EquityEquiv = cbeq,
      OptionEquiv = cbopt,
      NontrivCB = cb_40_recov,
      StraightBond = coupon_bond
    ),
    num_time_steps = N_Steps,
    const_default_intensity = 0.1,
    const_short_rate = 0.04,
    const_volatility = 0.5 + dVola
  ))
  down  = unlist(find_present_value(
    S0 = 252,
    instruments = list(
      EquityEquiv = cbeq,
      OptionEquiv = cbopt,
      NontrivCB = cb_40_recov,
      StraightBond = coupon_bond
    ),
    num_time_steps = N_Steps,
    const_default_intensity = 0.1,
    const_short_rate = 0.04,
    const_volatility = 0.5 - dVola
  ))

  by_repricing = (up - down) / (2 * dVola)
  for (nm in names(up)) {
    expect_equal(
      as.double(package_greeks$vega[nm]),
      as.double(by_repricing[nm]),
      tolerance = 0.01,
      info = paste("Instrument:", nm)
    )
  }
})


test_that("Rho/Rates DV01 computations are working", {
  skip_on_cran()
  up = unlist(
    find_present_value(
      S0 = 252,
      instruments = list(
        EquityEquiv = cbeq,
        OptionEquiv = cbopt,
        NontrivCB = cb_40_recov,
        StraightBond = coupon_bond
      ),
      num_time_steps = N_Steps,
      const_default_intensity = 0.1,
      const_short_rate = 0.04 + dr,
      const_volatility = 0.5
    )
  )
  down = unlist(
    find_present_value(
      S0 = 252,
      instruments = list(
        EquityEquiv = cbeq,
        OptionEquiv = cbopt,
        NontrivCB = cb_40_recov,
        StraightBond = coupon_bond
      ),
      num_time_steps = N_Steps,
      const_default_intensity = 0.1,
      const_short_rate = 0.04 - dr,
      const_volatility = 0.5
    )
  )


  by_repricing = (up - down) / (2 * dr)
  for (nm in names(up)) {
    expect_equal(
      as.double(package_greeks$rates_dv01[nm]),
      as.double(by_repricing[nm]),
      tolerance = 0.05,
      info = paste("Instrument:", nm)
    )
  }
})

test_that("Credit DV01 computations are working", {
  skip_on_cran()
  up = unlist(find_present_value(
    S0 = 252,
    instruments = list(
      EquityEquiv = cbeq,
      OptionEquiv = cbopt,
      NontrivCB = cb_40_recov,
      StraightBond = coupon_bond
    ),
    num_time_steps = N_Steps,
    const_default_intensity = 0.1 + dh,
    const_short_rate = 0.04,
    const_volatility = 0.5
  ))
  down = unlist(find_present_value(
    S0 = 252,
    instruments = list(
      EquityEquiv = cbeq,
      OptionEquiv = cbopt,
      NontrivCB = cb_40_recov,
      StraightBond = coupon_bond
    ),
    num_time_steps = N_Steps,
    const_default_intensity = 0.1 - dh,
    const_short_rate = 0.04,
    const_volatility = 0.5
  ))


  by_repricing = (up - down) / (2 * dh)
  for (nm in names(up)) {
    expect_equal(
      as.double(package_greeks$credit_dv01[nm]),
      as.double(by_repricing[nm]),
      tolerance = 0.05,
      info = paste("Instrument:", nm)
    )
  }
})



