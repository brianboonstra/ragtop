## Submission comments
Version raised to 1.3.0:

ragtop version 1.2.1 deprecates limSolve to a suggested package, providing substitutes when unavailable.

Details:
  - The `limSolve` package has gone through some periods of unreliability.
  - But, `limSolve` provides a convenient gateway to LAPACK DGTSV Gaussian elimination with partial pivoting.
  - We use `limSolve::Solve.tridiag()` guarded by a namespace check, otherwise fall back to banded methods from `Matrix` or Thomas' algorithm in pure R.
  - New unit tests check on matrix solvers.

## Test environments
* local MacOS 26.5.1 install, R 4.6.0
* Github CI MacOS 15.7.7, R 4.6.0
* Github CI Ubuntu-latest 24.04.04 LTS (release), R 4.6.0
* Github CI Ubuntu-latest 24.04.04 LTS (oldrel-1), R 4.5.3
* Github CI Ubuntu-latest 24.04.04 LTS (devel), unstable 2026-06-17 r90169
* Github CI Windows Server 2025 x64 (20260614.141.1), R 4.6.0
* CRAN Win-builder Windows Server 2022 x64 (build 20348, x86_64-w64-mingw32) unstable 2026-06-18 r90173 ucrt

## R CMD check results
One NOTE from win-builder:
* checking CRAN incoming feasibility ... [11s] NOTE
Maintainer: 'Brian K. Boonstra <ragtop@boonstra.org>'

Days since last update: 4

I did not expect to be updating so soon, but I received email early this morning from CRAN that my dependency limSolve is failing again.  This update eliminates the problem by deprecating that package to a suggestion.

There were no ERRORs, WARNINGs or other NOTEs. 


## Downstream dependencies
There are currently no downstream dependencies for this package
