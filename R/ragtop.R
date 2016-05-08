#' Pricing schemes for convertible bonds
#'
#' Using numerical integration, we price convertible bonds, straight
#'  bonds, equity options and various other derivatives consistently using a jump-diffusion model.
#'
#' Price derivatives using the stochastic model
#' \deqn{dS/S=(r+h-q) dt + \sigma dZ - dJ}
#' where \eqn{r} and \eqn{q} play their usual roles, \eqn{h} is a
#' deterministic function of stock price and time, and \eqn{J}
#' is a Poisson jump process adapted to the \emph{default intensity}
#' or \emph{hazard rate} \eqn{h}.  This model is a \emph{jump-diffusion}
#' extension of Black-Scholes, with the jump process \emph{J} representing
#' default, compensated by extra drift in the equity at rate \eqn{h}.
#' @docType package
#' @name ragtop
NULL
#> NULL
