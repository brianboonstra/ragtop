# ragtop 1.3.1 (CRAN)

Version 1.3.1 uses `limSolve` for matrix operations when available, otherwise falls back to `Matrix::bandSparse` solvers.  As a result, `limSolve` is demoted from a dependency to a suggested package.

The file `matrix.R` is added to handle this, along with new regression tests.  Vignette dependency on `BondValuation` is now guarded.

# ragtop 1.2.1 (CRAN)

Version 1.2.1 incorporates an important bugfix thanks to M Muecke, handling accumulation of previously-paid coupons.

Version 1.2.1 also add some vectorization and sanity checks with warnings, as well as other minor fixes suggested by Muecke.

# ragtop 1.2.0 (CRAN)

Version 1.2.0 recovers from deprecation / archiving which happened due to CI finding test failures in its limSolve dependency.  limSolve has now resolved those problems.

Version 1.2.0 also switches from Quandl to the treasury package for optional US interest rate downloads, and improves documentation formatting.

# ragtop 1.1.1 (CRAN)

Version 1.1.1 fixes extraneous parameter inheritance in documentation.

# ragtop 1.1.0 (CRAN)

Version 1.1.0 adds the `detail_from_AnnivDates()` function to make it simple for bond definitions to take advantage of the excellent daycount convention treatments coded into the `BondValuation` package.  A few typos in documentation were also fixed.

# ragtop 1.0.0 (CRAN)

After 2 stable years, `ragtop` is advancing to version 1.0.  The following minor changes were made since the 0.5 release:

* Added a `recovery_fcn` to the bond objects that reads the previously unused `recovery_rate` variable.  The implicit finite difference solver is able to use `recovery_fcn` to set default-conditional values on the grid.
* Documentation now includes more hyperlinks.
* Added a `NEWS.md` file to track changes to the package.

# ragtop 0.5 (CRAN)

The initial version of the package, released in 2016.
