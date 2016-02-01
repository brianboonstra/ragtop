library(futile.logger)

#' Constant to define when times are considered so close to each other
#'  that they should be treated as simultaneous
TIME_RESOLUTION_SIGNIF_DIGITS = 7
TIME_RESOLUTION_FACTOR = 10^(-TIME_RESOLUTION_SIGNIF_DIGITS)

log_layout_fcn = function(level, msg, ...)
{
  if (length(list(...)) > 0) {
    parsed <- lapply(list(...), function(x) ifelse(is.null(x),
                                                   "NULL", x))
    msg <- do.call(sprintf, c(msg, parsed))
  }
  sprintf("%s %s\n", names(level), msg)
}

flog.layout(log_layout_fcn, name='ragtop')
flog.layout(log_layout_fcn)

size_in_dimension = function(x,d=1)
{
  n = dim(x)[d]
  if (is.null(n)) {
    n = nrow(x)
  }
  if (is.null(n)) {
    n = length(x)
  }
  n
}
