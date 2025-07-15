## Submission comments
Version raised to 1.2.0:

ragtop version 1.2.0 recovers from deprecation / archiving which happened due to CI finding test failures in ragtop's limSolve dependency.  limSolve has now resolved those problems.

Details:
  - Quandl was purchased by NASDAQ, term structure queries updated to use treasury package instead
  - bugfix in dependency limSolve allows us to recover from deprecation
  - documentation formatting fixes

## Test environments
* local MacOS 15.5 24F74 install, R 4.4.2
* Github CI MacOS 14.7.6, R 4.5.1
* Github CI Ubuntu-latest 24.04.2 LTS (devel, release, oldrel-1), R 4.5.1
* Github CI Windows Server 2022 x64 (build 20348), R 4.5.1

## R CMD check results
There were no ERRORs, WARNINGs or NOTEs. 


## Downstream dependencies
There are currently no downstream dependencies for this package
