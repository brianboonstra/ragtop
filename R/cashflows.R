shift_for_dividends = function(grid_values_before_shift, stock_prices, div_sum)
{
  ## We need to shift grid values to account for the dividend paid
  ## In theory a linear interpolation is fine but in practice
  ## a smoother interpolator gives better results with little
  ## computational penalty, so we use splines
  shifted_prices = stock_prices - div_sum
  # Interpolate and extrapolate using a natural cubic spline
  #TODO: We can actually let grid_values_before_shift be a matrix
  # with columns for multiple layers, but we would then want to apply()
  # the spline interpolation over each column
  spl = stats::spline(
    x = stock_prices, y = grid_values_before_shift,
    xout = shifted_prices,
    method = "natural"
  ) # natural method allows linear extrapolation
  new_grid_values = spl$y
  new_grid_values
}



adjusted_for_dividends = function(grid_values, t, dt, r, h, S, S0, dividends)
{
  ## Adjust grid values according to any dividends paid during this timestep
  ## The dividends are expected to be represented by a data frame with
  ## columns 'time', 'fixed' and 'proportional' where the proportional amount
  ## is a multiplier of the stock price to form conditional dividend size
  ## Grid values are expected to come from f[k,m-1,]
  if (length(grid_values) != length(S)) {
    stop(paste("grid values length",length(grid_values),"must match underlying prices length",length(S)))
  }
  if (!is.null(dividends)) {
    # Divs are present. Sum any relevant ones
    div_sum = 0 * S
    included_ix = (dividends$time > t)  & (dividends$time <= t + dt)
    relevant_divs = dividends[included_ix,c('time', 'fixed', 'proportional')]
    if (nrow(relevant_divs) > 0) {
      flog.info("Found %s dividends from t=%s to t+dt=%s",
                nrow(relevant_divs), t, t + dt)
      inv_discount_factor = exp(+(t + dt - relevant_divs$time) * (r + h)) # positive/inverse because we are carrying forward
      fixed_contrib = inv_discount_factor * relevant_divs$fixed
      proportional_contrib = inv_discount_factor * relevant_divs$proportional
      div_amt_by_price = fixed_contrib + (proportional_contrib %o% S) /
        S0 # outer product %o%
      div_sum = colSums(div_amt_by_price)
    }
    if (any(div_sum != 0)) {
      grid_values = shift_for_dividends(grid_values, S, div_sum)
    }
  }
  grid_values
}
