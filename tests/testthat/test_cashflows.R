library(ragtop)
library(futile.logger)
context("Cash flow treatments")

flog.threshold(WARN, name="ragtop")
flog.threshold(WARN)

test_that("Spline interpolation on nontrivial vector", {
  expect_equal(shift_for_dividends(sin(0.3*(2:8)), 20+(2:8), 0.4+0.05*(2:8))[1], 0.448, tolerance=1.e-3)
  expect_equal(shift_for_dividends(sin(0.3*(2:8)), 20+(2:8), 0.4+0.05*(2:8))[7], 0.830, tolerance=1.e-3)
})

test_that("Spline interpolation on trivial vector", {
  expect_equal(shift_for_dividends(sin(0.3*(2:3)), 20+(2:3), 0.4+0.05*(2:3))[1], 0.455, tolerance=1.e-3)
  expect_equal(shift_for_dividends(sin(0.3*(2:3)), 20+(2:3), 0.4+0.05*(2:3))[2], 0.663, tolerance=1.e-3)
})

div_near_1_yr = data.frame(
  time=c(0.2, 1.01),
  fixed=c(0.001,0.003),
  proportional=c(100,0.5))
mult_divs = data.frame(
  time=c(0.2, 1.01, 1.02, 1.03, 1.04, 5),
  fixed=c(0.001,0.003, 0.005, 0.007, 0.009, 0.5),
  proportional=c(100,0.5, 0.2, 0.1, 0.15, 100))
nontriv_layer_1 = 0.1*(20:28)
nontriv_layer_2 = 2 + 0.1*(20:28)

nontriv_grid = matrix(c(0.1*(20:28), 2+0.1*(20:28)),ncol=2)

test_that("Dividend discounting and summing", {
  expect_equal(time_adj_dividends(mult_divs[2,], 1.2, 0, 0, 11:16, 12),
               c(0.461333333333333, 0.503, 0.544666666666667, 0.586333333333333,
                 0.628, 0.669666666666667),
               tolerance=1.e-6)
  expect_equal(time_adj_dividends(mult_divs[2:3,], 1.2, 0, 0, 11:16, 12),
               time_adj_dividends(mult_divs[2,], 1.2, 0, 0, 11:16, 12)+
                 time_adj_dividends(mult_divs[3,], 1.2, 0, 0, 11:16, 12),
               tolerance=1.e-6)
  expect_equal(time_adj_dividends(mult_divs[2:3,], 1.2, 0.2, 0, 11:16, 12),
               time_adj_dividends(mult_divs[2,], 1.2, 0.2, 0, 11:16, 12)+
                 time_adj_dividends(mult_divs[3,], 1.2, 0.2, 0, 11:16, 12),
               tolerance=1.e-6)
  expect_equal(time_adj_dividends(mult_divs[2:3,], 1.2, 0.2, 0.3, 11:16, 12),
               time_adj_dividends(mult_divs[2,], 1.2, 0.2, 0.3, 11:16, 12)+
                 time_adj_dividends(mult_divs[3,], 1.2, 0.2, 0.3, 11:16, 12),
               tolerance=1.e-6)
  expect_equal(time_adj_dividends(mult_divs[2:3,], 1.2, 0.2, 0.3, 12:13, 12),
               time_adj_dividends(mult_divs[2,], 1.2, 0.2, 0.3, 12:13, 12)+
                 time_adj_dividends(mult_divs[3,], 1.2, 0.2, 0.3, 12:13, 12),
               tolerance=1.e-6)
})

test_that("Basic dividend adjustments no discounting", {
  expect_equal(adjust_for_dividends(
    grid_values=3:7, t=0.9, dt=0.3, r=0, h=0, S=250*(23:27), S0=250*25,
    data.frame(time=1.01,fixed=200,proportional=0))[,1],
    c(2.2, 3.2, 4.2, 5.2, 6.2), tolerance=1.e-3)
  expect_equal(adjust_for_dividends(
    grid_values=3:7, t=0.9, dt=0.3, r=0, h=0, S=250*(23:27), S0=250*25,
    data.frame(time=1.01,fixed=0, proportional=200))[,1],
    c(2.264, 3.232, 4.2, 5.168, 6.136), tolerance=1.e-3)
})

test_that("Basic dividend adjustments with interest rates", {
  expect_equal(adjust_for_dividends(
    grid_values=3:7, t=0.9, dt=0.3, r=0.2, h=0, S=250*(23:27), S0=250*25,
    data.frame(time=1.01,fixed=200, proportional=0))[,1],
    c(2.1690150136972, 3.1690150136972, 4.1690150136972, 5.1690150136972,
      6.1690150136972), tolerance=1.e-3)
  expect_equal(adjust_for_dividends(
    grid_values=3:7, t=0.9, dt=0.3, r=0.2, h=0, S=250*(23:27), S0=250*25,
    data.frame(time=1.01,fixed=0, proportional=200))[,1],
    c(2.23549381260142, 3.20225441314931, 4.1690150136972, 5.13577561424509,
      6.10253621479298), tolerance=1.e-3)
})

test_that("Basic dividend adjustments with hazard rates", {
  expect_equal(adjust_for_dividends(
    grid_values=3:7, t=0.9, dt=0.3, r=0, h=0.2+0*(3:7), S=250*(23:27), S0=250*25,
    data.frame(time=1.01,fixed=200,proportional=0))[,1],
    c(2.1690150136972, 3.1690150136972, 4.1690150136972, 5.1690150136972,
      6.1690150136972), tolerance=1.e-5)
  expect_equal(adjust_for_dividends(
    grid_values=3:7, t=0.9, dt=0.3, r=0, h=0.2+0*(3:7), S=250*(23:27), S0=250*25,
    data.frame(time=1.01,fixed=0, proportional=200))[,1],
    c(2.23549381260142, 3.20225441314931, 4.1690150136972, 5.13577561424509,
      6.10253621479298), tolerance=1.e-5)
  expect_equal(adjust_for_dividends(
    grid_values=3:7, t=0.9, dt=0.3, r=0, h=0.05+0.04*(3:7), S=250*(23:27), S0=250*25,
    data.frame(time=1.01,fixed=0, proportional=200))[,1],
    c(2.23983910202986, 3.20073725569112, 4.16108303910297, 5.12087030001549,
      6.08009282445678), tolerance=1.e-5)
  expect_equal(adjust_for_dividends(
    grid_values=3:7, t=0.9, dt=0.3, r=0.3, h=0.05+0.04*(3:7), S=250*(23:27), S0=250*25,
    data.frame(time=1.01,fixed=0, proportional=200))[,1],
    c(2.19525124852845, 3.15385585170476, 4.11187568490703, 5.06930423502015,
      6.02613492358666), tolerance=1.e-5)
})


test_that("Multilayer grid, one div", {
  expect_equal(adjust_for_dividends(
    grid_values=nontriv_grid, t=0.9, dt=0.3, r=0, h=0, S=250*(20:28), S0=250*25,
    data.frame(time=1.01,fixed=200,proportional=0))[,1],
    adjust_for_dividends(
      grid_values=nontriv_layer_1, t=0.9, dt=0.3, r=0, h=0, S=250*(20:28), S0=250*25,
      data.frame(time=1.01,fixed=200,proportional=0))[,1],
    tolerance=1.e-5)
  expect_equal(adjust_for_dividends(
    grid_values=nontriv_grid, t=0.9, dt=0.3, r=0, h=0, S=250*(20:28), S0=250*25,
    data.frame(time=1.01,fixed=200,proportional=0))[,2],
    adjust_for_dividends(
      grid_values=nontriv_layer_2, t=0.9, dt=0.3, r=0, h=0, S=250*(20:28), S0=250*25,
      data.frame(time=1.01,fixed=200,proportional=0))[,1],
    tolerance=1.e-5)
  expect_equal(adjust_for_dividends(
    grid_values=nontriv_grid, t=0.9, dt=0.3, r=0, h=0.006*(20:28), S=250*(20:28), S0=250*25,
    data.frame(time=1.01,fixed=200,proportional=0))[,2],
    adjust_for_dividends(
      grid_values=nontriv_layer_2, t=0.9, dt=0.3, r=0, h=0.006*(20:28), S=250*(20:28), S0=250*25,
      data.frame(time=1.01,fixed=200,proportional=0))[,1],
    tolerance=1.e-5)
})

test_that("Multilayer grid, multiple divs", {
  expect_equal(adjust_for_dividends(
    grid_values=nontriv_grid, t=0.9, dt=0.3, r=0, h=0, S=250*(20:28), S0=250*25, mult_divs)[,1],
    adjust_for_dividends(
      grid_values=nontriv_layer_1, t=0.9, dt=0.3, r=0, h=0, S=250*(20:28), S0=250*25,mult_divs)[,1],
    tolerance=1.e-5)
  expect_equal(adjust_for_dividends(
    grid_values=nontriv_grid, t=0.9, dt=0.3, r=0, h=0, S=250*(20:28), S0=250*25,mult_divs)[,2],
    adjust_for_dividends(
      grid_values=nontriv_layer_2, t=0.9, dt=0.3, r=0, h=0, S=250*(20:28), S0=250*25,mult_divs)[,1],
    tolerance=1.e-5)
  expect_equal(adjust_for_dividends(
    grid_values=nontriv_grid, t=0.9, dt=0.3, r=0, h=0.006*(20:28), S=250*(20:28), S0=250*25,mult_divs)[,2],
    adjust_for_dividends(
      grid_values=nontriv_layer_2, t=0.9, dt=0.3, r=0, h=0.006*(20:28), S=250*(20:28), S0=250*25,mult_divs)[,1],
    tolerance=1.e-5)
})

test_that("Dividend adjustments working", {
  expect_equal(adjust_for_dividends(50*(20:28), 0.9, 0.3, 0.05, 0.09, (20:28), 25.5, div_near_1_yr)[,1],
               c(979.359595601146, 1028.29833225701, 1077.23349674515, 1126.16507994684,
                 1175.09307272263, 1224.01746591225, 1272.93825033462, 1321.85541678775,
                 1370.76895604874),
               tolerance=1.e-3)
  expect_equal(adjust_for_dividends(20*(20:28), 1.9, 0.3, 0.05, 0.09, (20:28), 25.5, div_near_1_yr)[2],
               420)
})
