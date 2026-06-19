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
* Github CI Ubuntu-latest 24.04.04 LTS (devel), unstable 2026-06-13 r90149
* Github CI Windows Server 2025 x64 (20260608.135.2), R 4.6.0
* CRAN Win-builder Windows Server 2022 x64 (build 20348, x86_64-w64-mingw32) unstable 2026-06-14 r90150 ucrt

## R CMD check results
There were no ERRORs, WARNINGs or NOTEs. 


## Downstream dependencies
There are currently no downstream dependencies for this package
