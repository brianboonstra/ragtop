.onLoad = function(libname, pkgname) {
  library(futile.logger)
  flog.threshold(WARN, name="ragtop")
  # The following fire warnings during calibration routines and I do not know how
  #  to temporarily disable them and then restore the user state
  flog.threshold(ERROR, name='ragtop.implicit.timestep.construct_tridiagonals')
  flog.threshold(ERROR, name='ragtop.calibration.implied_volatility_with_term_struct')
  flog.threshold(ERROR, name='ragtop.calibration.implied_volatility.lowprice')
  invisible()
}


.onAttach = function(libname, pkgname) {
  packageStartupMessage("Welcome to ragtop.  Logging can be enabled with commands such as\n  flog.threshold(INFO, name='ragtop.calibration')")
}
