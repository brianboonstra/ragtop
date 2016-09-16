#' Setup standard logging
#'
#' @import futile.logger
.onLoad = function(libname, pkgname) {
  futile.logger::flog.threshold(futile.logger::WARN, name="ragtop")
  # The following fire warnings during calibration routines and I do not know how
  #  to temporarily disable them and then restore the user state
  futile.logger::flog.threshold(futile.logger::ERROR, name='ragtop.implicit.timestep.construct_tridiagonals')
  futile.logger::flog.threshold(futile.logger::ERROR, name='ragtop.calibration.implied_volatility_with_term_struct')
  futile.logger::flog.threshold(futile.logger::ERROR, name='ragtop.calibration.implied_volatility.lowprice')
  futile.logger::flog.threshold(futile.logger::ERROR, name='ragtop.implicit.setup.width')
  invisible()
}


.onAttach = function(libname, pkgname) {
  packageStartupMessage("Welcome to ragtop.  Logging can be enabled with commands such as\n  futile.logger::flog.threshold(futile.logger::INFO, name='ragtop.calibration')")
}
