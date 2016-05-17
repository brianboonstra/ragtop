.onLoad = function(libname, pkgname) {
  library(futile.logger)
  flog.threshold(WARN, name="ragtop")
  invisible()
}


.onAttach = function(libname, pkgname) {
  packageStartupMessage("Welcome to ragtop.  Logging can be enabled with commands such as\n  flog.threshold(INFO, name='ragtop.calibration')")
}
