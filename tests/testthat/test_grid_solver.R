library(ragtop)
library(futile.logger)
context("Grid solver")

flog.threshold(WARN, name="ragtop")
flog.threshold(WARN)

test_that("Inferred grid parameters are correct", {
  expect_equal(construct_implicit_grid_structure(5, 100, 33, 33, 0.01, 0.5, 0.25, 2)$T, 5)
  expect_equal(construct_implicit_grid_structure(5, 100, 33, 33, 0.01, 0.5, 0.25, 2)$dt, .0500)
  expect_equal(construct_implicit_grid_structure(5, 100, 33, 33, 0.01, 0.5, 0.25, 2)$dz, 0.447213595499958)
  expect_equal(construct_implicit_grid_structure(5, 100, 33, 33, 0.01, 0.5, 0.25, 2)$z0, -0.575)
  expect_equal(construct_implicit_grid_structure(5, 100, 33, 33, 0.01, 0.5, 0.25, 2)$z_width, 2.23606797749979)
  expect_equal(construct_implicit_grid_structure(5, 100, 33, 33, 0.01, 0.5, 0.25, 2)$N, 11)
})

test_that("Tridiagonal matrix entries", {
  expect_equal(construct_tridiagonals(0., 0.05, 0+0.0*(2:6))$sub, c(0, 0, 0, 0), tolerance=1.e-4)
  expect_equal(construct_tridiagonals(0., 0.05, 0+0.0*(2:6))$diag, c(1,1,1,1,1), tolerance=1.e-4)
  expect_equal(construct_tridiagonals(0., 0.05, 0+0.0*(2:6))$super, c(0, 0, 0, 0), tolerance=1.e-4)
  expect_equal(construct_tridiagonals(0., 0.05, 0+0.0011*(2:6))$diag, c(1.0022, 1, 1, 1, 1.0066), tolerance=1.e-4)
  expect_equal(construct_tridiagonals(0., 0.05, 0+0.0011*(2:6))$super, c(-0.0022, -0.00165, -0.0022, -0.00275), tolerance=1.e-4)
  expect_equal(construct_tridiagonals(0., 0.05, 0+0.0011*(2:6))$sub, c(0.00165, 0.0022, 0.00275, -0.0066), tolerance=1.e-4)
  expect_equal(construct_tridiagonals(0.5, 0.05, 0+0.0*(2:6))$sub, c(-0.00625, -0.00625, -0.00625, 0), tolerance=1.e-4)
  expect_equal(construct_tridiagonals(0.5, 0.05, 0+0.0*(2:6))$diag, c(1, 1.0125, 1.0125, 1.0125, 1), tolerance=1.e-4)
  expect_equal(construct_tridiagonals(0.5, 0.05, 0+0.0011*(2:6))$diag, c(1.0022, 1.0125, 1.0125, 1.0125, 1.0066), tolerance=1.e-4)
  expect_equal(construct_tridiagonals(0.5, 0.05, 0+0.0011*(2:6))$sub, c(-0.0046, -0.00405, -0.0035, -0.0066), tolerance=1.e-4)
  expect_equal(construct_tridiagonals(0.5, 0.05, 0+0.0011*(2:6))$super, c(-0.0022, -0.0079, -0.00845, -0.009), tolerance=1.e-4)
})

v_1 = c(0.95,0.9, 0.8, 0.6)
v_2 = 2*v_1
v_both = matrix(c(v_1,v_2), ncol=2)
test_that("Single-instrument timestep", {
  expect_equal(as.vector( # as.vector because get a matrix with 1 column back
    take_implicit_timestep(1, 1, c(0.95,0.9, 0.8), 1,
                           list(super=c(0,0), diag=c(1,1,1), sub=c(0,0)))),
    c(0.95, 0.9, 0.8),
    tolerance=1.e-6)
  expect_equal(as.vector( # as.vector because get a matrix with 1 column back
    take_implicit_timestep(1, 1, c(0.95,0.9, 0.8), 0.8,
                           list(super=c(0,0), diag=c(1,1,1), sub=c(0,0)))),
    c(0.76, 0.72, 0.64),
    tolerance=1.e-6)
  expect_equal(as.vector( # as.vector because get a matrix with 1 column back
    take_implicit_timestep(1, 0.8, c(0.95,0.9, 0.8), 1,
                           list(super=c(0,0), diag=c(1,1,1), sub=c(0,0)))),
    c(0.95, 0.9, 0.8),
    tolerance=1.e-6)
  expect_equal(as.vector( # as.vector because get a matrix with 1 column back
    take_implicit_timestep(1, 1, v_1, 1,
                           list(super=c(-0.15,0,0), diag=c(1.15,1,1,1.1), sub=c(0,0,-0.1)))
  ),
  c(0.943478260869565, 0.9, 0.8, 0.618181818181818),
  tolerance=1.e-6)
  expect_equal(as.vector( # as.vector because get a matrix with 1 column back
    take_implicit_timestep(1, 1, v_1, 1,
                           list(super=c(-0.15,-0.15,-0.15), diag=c(1.15,1.3,1.3,1.1), sub=c(-0.15,-0.15,-0.1)))
  ),
  c(0.942454838513472, 0.892153761936622, 0.789544431603917, 0.617231311963992),
  tolerance=1.e-6)
  expect_equal(as.vector( # as.vector because get a matrix with 1 column back
    take_implicit_timestep(1, 1, v_1, 0.8,
                           list(super=c(-0.15,-0.15,-0.15), diag=c(1.15,1.3,1.3,1.1), sub=c(-0.15,-0.15,-0.1)))
  ),
  c(0.753963870810778, 0.713723009549298, 0.631635545283134, 0.493785049571194),
  tolerance=1.e-6)
  expect_equal(as.vector( # as.vector because get a matrix with 1 column back
    take_implicit_timestep(1, 1, v_2, 0.8,
                           list(super=c(-0.15,-0.15,-0.15), diag=c(1.15,1.3,1.3,1.1), sub=c(-0.15,-0.15,-0.1)))
  ),
  2*c(0.753963870810778, 0.713723009549298, 0.631635545283134, 0.493785049571194),
  tolerance=1.e-6)
  expect_equal(
    take_implicit_timestep(1, 1, v_both, 0.8,
                           list(super=c(-0.15,-0.15,-0.15), diag=c(1.15,1.3,1.3,1.1), sub=c(-0.15,-0.15,-0.1))),
  structure(c(0.753963870810778, 0.713723009549298, 0.631635545283134,
              0.493785049571194, 1.50792774162156, 1.4274460190986, 1.26327109056627,
              0.987570099142388), .Dim = c(4L, 2L)),
  tolerance=1.e-6)
})

very_short_bond = ZeroCouponBond(notional=1.0, maturity=0.125, discount_factor_fcn=function(...){NULL})
zcb_prices_h = find_present_value(100, 1, instruments=list(zcb=very_short_bond), const_default_intensity=0.4103464)
zcb_prices_r = find_present_value(100, 1, instruments=list(zcb=very_short_bond), const_short_rate=0.4103464)
test_that("Trivial 1 timestep solutions", {
  expect_equal(max(zcb_prices_h), 0.95, tolerance=1.e-6)
  expect_equal(min(zcb_prices_h), max(zcb_prices_h), tolerance=1.e-6)
  expect_equal(max(zcb_prices_r), 0.95, tolerance=1.e-6)
  expect_equal(min(zcb_prices_r), max(zcb_prices_r), tolerance=1.e-6)
})
