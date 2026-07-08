## Submission comments
Version raised to 2.0.0:

ragtop version 2.0.0 allows for calculating greeks, and better cache control.  

Details:
  - New reset_caches() allows for resetting price grid caches.  This makes it possible to price an instrument multiple times without recreating the instrument object.
  - Calculation of greeks, from the pricing grid where possible, and from "bumping" input values where not.
  - Fixed the futile.logger layout so logger names no longer leak into formatted log messages.

## Test environments
* local MacOS 26.5.1 install, R 4.6.0
* Github CI MacOS 26.4 (20260630.0213.1), R 4.6.1 (2026-06-24)
* Github CI Ubuntu-latest 24.04.04 LTS (release), R 4.6.1
* Github CI Ubuntu-latest 24.04.04 LTS (oldrel-1), R 4.5.3
* Github CI Ubuntu-latest 24.04.04 LTS (devel), unstable 2026-06-21 r90185
* Github CI Windows Server 2025 x64 (build 26100), R 4.6.1 (2026-06-24 ucrt)
* CRAN Win-builder Windows Server 2022 x64 (build 20348) unstable 2026-07-07 r90210 ucrt

## R CMD check results
There were no ERRORs, WARNINGs or NOTEs. 


## Downstream dependencies
There are currently no downstream dependencies for this package
