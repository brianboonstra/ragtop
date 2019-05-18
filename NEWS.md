# ragtop 1.1.0 (CRAN)

Version 1.1.0 adds the `detail_from_AnnivDates()` function to make it simple for bond definitions to take advantage of the excellent daycount convention treatments coded into the `BondValuation` package.  A few typos in documentation were also fixed.

# ragtop 1.0.0 (CRAN)

After 2 stable years, `ragtop` is advancing to version 1.0.  The following minor changes were made since the 0.5 release:

* Added a `recovery_fcn` to the bond objects that reads the previously unused `recovery_rate` variable.  The implicit finite difference solver is able to use `recovery_fcn` to set default-conditional values on the grid.
* Documentation now includes more hyperlinks.
* Added a `NEWS.md` file to track changes to the package.

# ragtop 0.5 (CRAN)

The initial version of the package, released in 2016.
