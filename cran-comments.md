## Submission comments
Version raised to 2.0.0:

ragtop version 2.0.0 allows for calculating greeks, and better cache control.  

Details:
  - New reset_caches() allows for resetting price grid caches.  This makes it possible to price an instrument multiple times without recreating the instrument object.
  - Calculation of greeks, from the pricing grid where possible, and from "bumping" input values where not.
  - Fixed the futile.logger layout so logger names no longer leak into formatted log messages.

## Test environments
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
