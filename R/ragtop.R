#' Pricing schemes for derivatives using equity-linked default intensity
#'
#' Using numerical integration, we price convertible bonds, straight
#'  bonds, equity options and various other derivatives consistently using a
#'  jump-diffusion model in which default intensity can vary with equity price
#'  in a user-specified deterministic manner.
#'
#' We apply the stochastic model
#' \deqn{dS/S=(r+h-q) dt + \sigma dZ - dJ}
#' where \eqn{r} and \eqn{q} play their usual roles, \eqn{h} is a
#' deterministic function of stock price and time, and \eqn{J}
#' is a Poisson jump process adapted to the \emph{default intensity}
#' or \emph{hazard rate} \eqn{h}.  This model is a \emph{jump-diffusion}
#' extension of Black-Scholes, with the jump process \emph{J} representing
#' default, compensated by extra drift in the equity at rate \eqn{h}.
#'
#' Volatilities, default intensities and risk-free rates may all be represented
#'   with arbitrary term structures.  Default intensity term structures may also
#'   take the underlying equity price into account.
#'
#' Pricing in the standard Black-Scholes model is a special case with default
#'  intensity set to zero.  Therefore this package \emph{also} serves to price securities
#'  in the standard Black-Scholes model, while still allowing risk-free rates and volatilites
#'  have nontrivial term structures.
#'
#' @section Important Features:
#'
#' \describe{
#'  \item{\emph{Black-Scholes} }{The standard model is automatically supported as a special case, but also has optimized routines}
#'  \item{\emph{Term Structures} }{The package allows for any kind of instrument to be priced with time-varying rates, volatility and default intensity}
#'  \item{\emph{Dividends} }{Allows for discrete dividends in an arbitrary combination of fixed and proportional amounts.  The difference between fixed and proprtional can be up to 10 percent in implied volatility terms.}
#'  \item{\emph{Calibration} }{Model calibration routines are included}
#'  \item{\emph{Bankruptcy Realism} }{A parsimonious deterministic model of default intensity gives rich behavior and conforms reasonably well to observed market data}
#'  \item{\emph{Algorithm Parameters} }{Default parameters for the algorithm work well for a very wide variety of pricing and implied volatility scenarios}
#'  }
#'
#'
#' @examples
#' ## Vanilla European exercise
#' blackscholes(callput=-1, S0=100, K=90, r=0.03, time=1, vola=0.5)
#' blackscholes(PUT, S0=100, K=90, r=0.03, time=1, vola=0.5,
#'              default_intensity=0.07, borrow_cost=0.005)
#' ## With a term structure of volatility
#' \dontrun{
#' black_scholes_on_term_structures(callput=-1, S0=100, K=90, time=1,
#'                                  const_short_rate=0.025,
#'                                  variance_cumulation_fcn = function(T, t) {
#'                                    0.45 ^ 2 * (T - t) + 0.15^2 * max(0, T-0.25)
#'                                  })
#' }
#'
#' ## Vanilla American exercise
#' \dontrun{
#' american(PUT, S0=100, K=110, time=0.77, const_short_rate = 0.06,
#'          const_volatility=0.20, num_time_steps=200)
#' }
#' ## With a term structure of volatility
#' \dontrun{
#' american(callput=-1, S0=100, K=90, time=1, const_short_rate=0.025,
#'          variance_cumulation_fcn = function(T, t) {
#'              0.45 ^ 2 * (T - t) + 0.15^2 * max(0, T-0.25)
#'          })
#' }
#' ## With discrete dividends, combined fixed and proportional
#' divs = data.frame(time=seq(from=0.11, to=2, by=0.25),
#'                   fixed=seq(1.5, 1, length.out=8),
#'                   proportional = seq(1, 1.5, length.out=8))
#' \dontrun{
#' american(callput=-1, S0=100, K=90, time=1, const_short_rate=0.025,
#'          const_volatility=0.20, dividends=divs)
#' }
#'
#' ## American Exercise Implied Volatility
#' american_implied_volatility(25,CALL,S0=100,K=100,time=2.2, const_short_rate=0.03)
#' df250 =  function(t) ( exp(-0.02*t)*exp(-0.03*max(0,t-1.0))) # Simple term structure
#' df25 = function(T,t){df250(T)/df250(t)} # Relative discount factors
#' \dontrun{
#' american_implied_volatility(25,-1,100,100,2.2,discount_factor_fcn=df25)
#' }
#'
#' ## Convertible Bond
#' ## Not Run
#' pct4 = function(T,t=0) { exp(-0.04*(T-t)) }
#' cb = ConvertibleBond(conversion_ratio=3.5, maturity=1.5, notional=100,
#'                      discount_factor_fcn=pct4, name='Convertible')
#' S0 = 10; p = 6.0; h = 0.10
#' h_fcn = function(t, S, ...){0.9 * h + 0.1 * h * (S0/S)^p }  # Intensity linked to equity price
#' \dontrun{
#' find_present_value(S0=S0, instruments=list(Convertible=cb), num_time_steps=250,
#'                    default_intensity_fcn=h_fcn,
#'                    const_volatility = 0.4, discount_factor_fcn=pct4,
#'                    std_devs_width=5)
#' }
#'
#' ## Fitting Term Structure of Volatility
#' ## Not Run
#' opts = list(m1=AmericanOption(callput=-1, strike=9.9, maturity=1/12, name="m1"),
#'             m2=AmericanOption(callput=-1, strike=9.8, maturity=1/6, name="m2"))
#' \dontrun{
#' vfit = fit_variance_cumulation(S0, opts, c(0.6, 0.8), default_intensity_fcn=h_fcn)
#' print(vfit$volatilities)
#' }
#'
#' @keywords internal
#' @import futile.logger
"_PACKAGE"
