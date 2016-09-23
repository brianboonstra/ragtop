<!-- README.md is generated from README.Rmd. Please edit that file -->
Description And Installation
============================

**ragtop** prices equity derivatives are priced on variants of the famous Black-Scholes model, with special attention paid to the case of American and European exercise options and to convertible bonds. To install the development version, use the command

``` r
devtools::install_github('brianboonstra/ragtop')
```

Usage
=====

Basic Usage
-----------

You can price american and european exercise options, either individually, or in groups. In the simplest case that looks like this for European exercise

``` r
blackscholes(c(CALL, PUT), S0=100, K=c(100,110), time=0.77, r = 0.06, vola=0.20)
#> $Price
#> [1] 9.326839 9.963285
#> 
#> $Delta
#> [1]  0.6372053 -0.5761608
#> 
#> $Vega
#> [1] 32.91568 34.36717
```

and like this for American exercise

``` r
american(PUT, S0=100, K=c(100,110), time=0.77, const_short_rate = 0.06, const_volatility=0.20)
#> A100_281_0 A110_281_0 
#>    5.24386   11.27715
```

### Including Term Structures

There are zillions of implementations of the Black-Scholes formula out there, and quite a few simple trees as well. One thing that makes **ragtop** a bit more useful than most other packages is that it treats dividends and term structures without too much pain. Assume we have some nontrivial term structures and dividends

``` r
## Dividends
divs = data.frame(time=seq(from=0.11, to=2, by=0.25),
                  fixed=seq(1.5, 1, length.out=8),
                  proportional = seq(1, 1.5, length.out=8))

## Interest rates
disct_fcn = ragtop::spot_to_df_fcn(data.frame(time=c(1, 5, 10), 
                                              rate=c(0.01, 0.02, 0.035)))

## Default intensity
disc_factor_fcn = function(T, t, ...) {
  exp(-0.03 * (T - t)) }
surv_prob_fcn = function(T, t, ...) {
  exp(-0.07 * (T - t)) }

## Variance cumulation / volatility term structure
vc = variance_cumulation_from_vols(
   data.frame(time=c(0.1,2,3),
              volatility=c(0.2,0.5,1.2)))
paste0("Cumulated variance to 18 months is ", vc(1.5, 0))
[1] "Cumulated variance to 18 months is 0.369473684210526"
```

then we can price vanilla options

``` r
black_scholes_on_term_structures(
   callput=TSLAMarket$options[500,'callput'], 
   S0=TSLAMarket$S0, 
   K=TSLAMarket$options[500,'K'], 
   discount_factor_fcn=disct_fcn, 
   time=TSLAMarket$options[500,'time'], 
   variance_cumulation_fcn=vc,
   dividends=divs)
$Price
[1] 62.55998

$Delta
[1] 0.7977684

$Vega
[1] 52.21925
```

American exercise options

``` r
american(
    callput = TSLAMarket$options[400,'callput'], 
    S0 = TSLAMarket$S0, 
    K=TSLAMarket$options[400,'K'], 
    discount_factor_fcn=disct_fcn, 
    time = TSLAMarket$options[400,'time'],
    survival_probability_fcn=surv_prob_fcn,
    variance_cumulation_fcn=vc,
    dividends=divs)
A360_137_2 
  2.894296 
```

We can also find volatilities of European exercise options

``` r
implied_volatility_with_term_struct(
    option_price=19, callput = PUT, 
    S0 = 185.17,K=182.50, 
    discount_factor_fcn=disct_fcn, 
    time = 1.12,
    survival_probability_fcn=surv_prob_fcn,
    dividends=divs)
[1] 0.1133976
```

as well as American exercise options

``` r
american_implied_volatility(
    option_price=19, callput = PUT, 
    S0 = 185.17,K=182.50, 
    discount_factor_fcn=disct_fcn, 
    time = 1.12,
    survival_probability_fcn=surv_prob_fcn,
    dividends=divs)
[1] 0.113407
```

More Sophisticated Calibration
------------------------------

You can also find more complete calibration routines in **ragtop**. See the vignette or the documentation for *fit\_variance\_cumulation* and *fit\_to\_option\_market*.

Technical Documentation
=======================

The source for the technical paper is in this repository. You can also find the pdf [here](http://thureoscapital.com/ragtop.pdf)

[![Travis-CI Build Status](https://travis-ci.org/brianboonstra/ragtop.svg?branch=master)](https://travis-ci.org/brianboonstra/ragtop)
