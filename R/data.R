#' Market information snapshot for TSLA options
#'
#' A dataset containing option contract details and a
#'  snapshot of market prices
#'  for Tesla Motors (TSLA) equity options, interest rates and an equity price.
#'
#' The TSLAMarket list contains three elements:
#' \describe{
#' \item{\code{S0}}{The stock price as of snapshot time}
#' \item{\code{risk_free_rates}}{The spot risk-free rate curve as of snapshot time}
#' \item{\code{options}}{A data frame with details of the options market}
#' }
#' @docType data
#' @usage data(TSLAMarket)
#' @format A list with these prices and rates
#' @export
"TSLAMarket"
