## Submission comments
Version raised to 1.1, still no known bugs in 3 years.  

One minor enhancement: added a helper function to use daycount-derived coupon dates as output in the list object from BondValuation package.  

This does not introduce BondValuation as a dependency, since the only actual use of BondValuation is a call in the vignette.

## Test environments
* local MacOS 10.14.4 install, R 3.5.1
* Ubuntu 14.04.5 LTS (on travis-ci), R 3.6.0
* win-builder, R unstable 2019-05-16 r76519

## R CMD check results
There were no ERRORs, WARNINGs or NOTEs. 


## Downstream dependencies
There are currently no downstream dependencies for this package
