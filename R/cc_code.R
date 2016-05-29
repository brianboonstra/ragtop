## ragtop -- convertibles pricing in R
##
## Copyrights by the authors, Creative Commons Licenses
##
## This module contains snippets of open-source code from other
## sources that are subject to the Creative Commons license
##  http://creativecommons.org/licenses/by-sa/3.0/

## Modified from:
## http://stackoverflow.com/questions/19655579/a-function-that-returns-true-on-na-null-nan-in-r
## Copyright user Backlin


#' Return TRUE if the argument is empty, NULL or NA
#'
#' @param x Argument to test
#' @param false.triggers Whether to allow nonempty vectors of all FALSE to trigger this condition
#' @export is.blank
is.blank <- function(x, false.triggers=FALSE){
  if(is.function(x) || typeof(x)=="S4") return(FALSE) # Some of the tests below trigger
  # warnings when used on functions
  return(
    is.null(x) ||                # Actually this line is unnecessary since
      length(x) == 0 ||            # length(NULL) = 0, but I like to be clear
      all(is.na(x)) ||
      all(x=="") ||
      (false.triggers && all(!x))
  )
}
