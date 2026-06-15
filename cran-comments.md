## Submission comments
Version raised to 1.2.1:

ragtop version 1.2.1 fixes a bug in bond past-coupon evaluation.  It also adds vectorization and input checking.

Details:
  - Coupon accumulation in previous versions was discounting rather than accumulating.  Many thanks to Max Muecke (author of the treasury package) for finding the bug!
  - Argument vectorization added (thanks again to Max) in a couple key places
  - The new checkinputs.R has functions for ensuring the sanity of user-supplied rates, vol, and hazard functions.  These are now used to log at the warning level when integrate_pde is called with, for example, negative interest rates. 
  - The vignette now has guards around the BondValuation package, bypassing errors on some platforms.
  - Documentation formatting fixes

## Test environments
* local MacOS 26.5.1 install, R 4.4.2
* local MacOS 26.5.1 install, R 4.6.0
* Github CI MacOS 15.7.7, R 4.6.0
* Github CI Ubuntu-latest 24.04.04 LTS (release), R 4.6.0
* Github CI Ubuntu-latest 24.04.04 LTS (oldrel-1), R 4.5.3
* Github CI Ubuntu-latest 24.04.04 LTS (devel), unstable 2026-06-13 r90149
* Github CI Windows Server 2025 x64 (20260608.135.2), R 4.6.0
* CRAN Win-builder Windows Server 2022 x64 (build 20348, x86_64-w64-mingw32) unstable 2026-06-14 r90150 ucrt

## R CMD check results
There were no ERRORs, WARNINGs or NOTEs. 


## Downstream dependencies
There are currently no downstream dependencies for this package
