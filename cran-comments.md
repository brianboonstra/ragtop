## Submission comments
Version raised to 1.2.1:

ragtop version 1.2.1 fixes a bug in bond past-coupon evaluation.  It also adds vectorization and input checking.

Details:
  - Coupon accumulation in previous versions was discounting rather than accumulating.  Many thanks fto Max Muecke (author of the treasury package) for finding the bug!
  - The new checkinputs.R has functions for ensuring the sanity of user-supplied rates, vol, and hazard functions.  These are now used to log at the warning level when integrate_pde is called with, for example, negative interest rates. 
  - The vignette now has guards around the BondValuation package, bypassing errors on some platforms.
  - documentation formatting fixes

## Test environments
* local MacOS 26.5.1 install, R 4.4.2
* Github CI MacOS 14.7.6, R 4.5.1
* Github CI Ubuntu-latest 24.04.2 LTS (devel, release, oldrel-1), R 4.5.1
* Github CI Windows Server 2022 x64 (build 20348), R 4.5.1

## R CMD check results
There were no ERRORs, WARNINGs or NOTEs. 


## Downstream dependencies
There are currently no downstream dependencies for this package
