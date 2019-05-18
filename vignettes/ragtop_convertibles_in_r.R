## ----global_options, include=FALSE---------------------------------------
library(ragtop)
library(futile.logger)
library(ggplot2)

flog.threshold(ERROR)
flog.threshold(ERROR, name='ragtop.implicit.timestep.construct_tridiagonals')
flog.threshold(ERROR, name='ragtop.calibration.implied_volatility.lowprice')
flog.threshold(ERROR, name='ragtop.calibration.implied_volatility_with_term_struct')
flog.threshold(ERROR, name='ragtop.implicit.setup.width')

knitr::opts_chunk$set(fig.width=6.5, fig.height=4, fig.path='Figs/',
                      echo=FALSE, warning=FALSE, message=FALSE, comment=FALSE)


## ----show_TSLA_S0, comment="", echo=TRUE---------------------------------
TSLAMarket$S0

## ----show_TSLA_rf, comment="", echo=TRUE---------------------------------
knitr::kable(TSLAMarket$risk_free_rates, digits=3, row.names = F)

## ---- results='asis', comment=""-----------------------------------------
knitr::kable(TSLAMarket$options[c(200,300,400,500,600, 800),], digits=3, row.names = F)

## ----blackscholes, comment=""--------------------------------------------
blackscholes(TSLAMarket$options[500,'callput'], 
             TSLAMarket$S0, 
             TSLAMarket$options[500,'K'], 
             0.005, 
             TSLAMarket$options[500,'time'], 
             0.50)

## ----implied_volatility, echo=TRUE, comment=""---------------------------
implied_volatility(option_price = TSLAMarket$options[400,'ask'], 
                   S0 = TSLAMarket$S0, 
                   callput = TSLAMarket$options[400,'callput'], 
                   K=TSLAMarket$options[400,'K'], 
                   r = 0.005, 
                   time = TSLAMarket$options[400,'time'])

## ----amer, echo=TRUE, comment=""-----------------------------------------
american(
       callput = TSLAMarket$options[400,'callput'], 
       S0 = TSLAMarket$S0, 
       K=TSLAMarket$options[400,'K'], 
       const_short_rate = 0.005, 
       time = TSLAMarket$options[400,'time'])

## ----american_implied_volatility, echo=TRUE, comment=""------------------
american_implied_volatility(option_price = TSLAMarket$options[400,'ask'], 
     S0 = TSLAMarket$S0, 
     callput = TSLAMarket$options[400,'callput'], 
     K=TSLAMarket$options[400,'K'], 
     const_short_rate = 0.005, 
     time = TSLAMarket$options[400,'time'])

## ----implied_volatility_def, comment="", echo=TRUE-----------------------
implied_volatility(option_price = 17, 
                   S0 = 250, callput = CALL,  K=245,
                   r = 0.005, time = 2,
                   const_default_intensity = 0.03)

## ----american_implied_volatility_def, comment="", echo=T-----------------
american_implied_volatility(option_price = 19.1, 
     S0 = 223.17, callput = PUT, K=220, 
     const_short_rate = 0.005, time = 1.45,
     const_default_intensity = 0.0200)

## ----ts_fcns, echo=TRUE, comment=""--------------------------------------
## Dividends
divs = data.frame(time=seq(from=0.11, to=2, by=0.25),
                  fixed=seq(1.5, 1, length.out=8),
                  proportional = seq(1, 1.5, length.out=8))

## Interest rates
disct_fcn = ragtop::spot_to_df_fcn(
  data.frame(time=c(1, 5, 10, 15),
             rate=c(0.01, 0.02, 0.03, 0.05))
)

## Default intensity
surv_prob_fcn = function(T, t, ...) {
  exp(-0.07 * (T - t)) }

## Variance cumulation / volatility term structure
vc = variance_cumulation_from_vols(
   data.frame(time=c(0.1,2,3),
              volatility=c(0.2,0.5,1.2)))
paste0("Cumulated variance to 18 months is ", vc(1.5, 0))

## ----blackscholes_ts, comment="", echo=TRUE------------------------------
black_scholes_on_term_structures(
   callput=TSLAMarket$options[500,'callput'], 
   S0=TSLAMarket$S0, 
   K=TSLAMarket$options[500,'K'], 
   discount_factor_fcn=disct_fcn, 
   time=TSLAMarket$options[500,'time'], 
   survival_probability_fcn=surv_prob_fcn,
   variance_cumulation_fcn=vc,
   dividends=divs)

## ----amer_ts, echo=TRUE, comment=""--------------------------------------
american(
    callput = TSLAMarket$options[400,'callput'], 
    S0 = TSLAMarket$S0, 
    K=TSLAMarket$options[400,'K'], 
    discount_factor_fcn=disct_fcn, 
    time = TSLAMarket$options[400,'time'],
    survival_probability_fcn=surv_prob_fcn,
    variance_cumulation_fcn=vc,
    dividends=divs)

## ----implied_bs_volatility_def_ts, comment="", echo =T-------------------
implied_volatility_with_term_struct(
   option_price = TSLAMarket$options[400,'ask'], 
   S0 = TSLAMarket$S0, 
   callput = TSLAMarket$options[400,'callput'], 
   K=TSLAMarket$options[400,'K'], 
   discount_factor_fcn=disct_fcn, 
   time = TSLAMarket$options[400,'time'],
   survival_probability_fcn=surv_prob_fcn,
   dividends=divs)

## ----amer_ts_iv, echo=TRUE, comment=""-----------------------------------
american_implied_volatility(
    option_price=TSLAMarket$options[400,'ask'], 
    callput = TSLAMarket$options[400,'callput'], 
    S0 = TSLAMarket$S0, 
    K=TSLAMarket$options[400,'K'], 
    discount_factor_fcn=disct_fcn, 
    time = TSLAMarket$options[400,'time'],
    survival_probability_fcn=surv_prob_fcn,
    dividends=divs)

## ----def_int_fav, echo=TRUE, comment=""----------------------------------
  def_ints_fcn = function(t, S, ...){
    0.09+0.01*(S0/S)^1.5
    }

## ----impl_curve_elems, echo=TRUE, comment=""-----------------------------
  options_df = TSLAMarket$options
  S0 = TSLAMarket$S0
  make_option = function(x) {
    if (x['callput']>0) cp='C' else cp='P'
    ragtop::AmericanOption(callput=x['callput'], strike=x['K'], maturity=x['time'],
                           name=paste(cp,x['K'],as.integer(100*x['time']), sep='_'))
  }
  atm_put_price = max(options_df$K[options_df$K<=S0])
  atm_put_ix = ((options_df$K==atm_put_price) & (options_df$callput==PUT)
                & (options_df$time>1/12))
  atm_puts = apply(options_df[atm_put_ix,], 1, make_option)

  atm_put_prices = options_df$mid[atm_put_ix]
  knitr::kable(options_df[atm_put_ix,], digits=3, row.names = F)

## ----varfit, echo=TRUE, comment=""---------------------------------------
vcm = fit_variance_cumulation(S0, eq_options=atm_puts,
       mid_prices=atm_put_prices,
       spreads=0.01*atm_put_prices,
       use_impvol=TRUE,
       discount_factor_fcn=disct_fcn, 
       default_intensity_fcn = def_ints_fcn,
       num_time_steps=100)
vcm$volatilities

## ----fit_elems, echo=FALSE, comment=""-----------------------------------
  h0 = 0.05
  fit_target_ix = ((options_df$K==210 | options_df$K==220 ) & (options_df$callput==PUT)
                & (options_df$time>6/12))
  fit_targets = apply(options_df[atm_put_ix,], 1, make_option)

  fit_target_prices = options_df$mid[atm_put_ix]
  knitr::kable(options_df[fit_target_ix,], digits=3, row.names = F)

## ----fit_pen, echo=TRUE, comment=""--------------------------------------
  h0 = 0.05
  fit_penalty = function(p, s) {
  def_intens_f = function(t,S,...) {h0 * (s + (1-s) * (S0/S)^p)}
  varnce = fit_variance_cumulation(
              S0, eq_options=atm_puts,
              mid_prices=atm_put_prices,
              spreads=0.01*atm_put_prices,
              use_impvol=TRUE,
              default_intensity_fcn = def_intens_f,
              discount_factor_fcn=disct_fcn, 
              num_time_steps=100)
  pvs_list = find_present_value(
              S0=S0, instruments=fit_targets,
              default_intensity_fcn=def_intens_f,
              variance_cumulation_fcn=varnce$cumulation_function,
              discount_factor_fcn=disct_fcn,
              num_time_steps=45)
  pvs = as.numeric(pvs_list)
  pensum = sum((fit_target_prices - pvs)^2)
  pensum
}
fit_penalty(1, 0.5)

## ----define_cb, echo=TRUE, comment=""------------------------------------
cb = ragtop::ConvertibleBond(
  maturity=2.87, conversion_ratio=2.7788, notional=1000,
  coupons=data.frame(payment_time=seq(2.8,0, by=-0.25),
                     payment_size=1000*0.0025/4),
  discount_factor_fcn = disct_fcn,
  name='CBond'
)

s = 0.75
h0 = 0.04
p = 2.5

## ----associated_f, echo=TRUE, comment=""---------------------------------
calibrated_intensity_f = function(t, S, ...){
  0.03+0.01*(S0/S)^1.5
}

calib_varnce = fit_variance_cumulation(
              S0, eq_options=atm_puts,
              mid_prices=atm_put_prices,
              spreads=0.01*atm_put_prices,
              use_impvol=TRUE,
              default_intensity_fcn = calibrated_intensity_f,
              discount_factor_fcn=disct_fcn, 
              num_time_steps=100)
calib_varnce$volatilities

## ----price_cb, echo=TRUE, comment=""-------------------------------------
cb_value = form_present_value_grid(
  S0=S0, grid_center=S0,
  instruments=list(Convertible=cb),
  num_time_steps=250,
  default_intensity_fcn=calibrated_intensity_f,
  discount_factor_fcn = disct_fcn,
  variance_cumulation_fcn=calib_varnce$cumulation_function)

## ----grid_cb_delta, echo=TRUE, comment=""--------------------------------
cbprices = ragtop::form_present_value_grid(
  S0=S0, grid_center=S0,
  instruments=list(Convertible=cb),
  num_time_steps=250,
  default_intensity_fcn=calibrated_intensity_f,
  discount_factor_fcn = disct_fcn,
  variance_cumulation_fcn=calib_varnce$cumulation_function,
  std_devs_width=5)
cbgrid = na.omit(as.data.frame(cbprices))

present_value_interp = splinefun(
  x=cbgrid[,"Underlying"],
  y=cbgrid[,"Convertible"])

delta = present_value_interp(S0, deriv=1)
delta

## ----grid_cb_plot, echo=TRUE, comment=""---------------------------------
present_value = present_value_interp(S0)
cbplot = ( ggplot(cbgrid,
                  aes(x=Underlying,y=Convertible)) +
             geom_line(size=1.2) +
             scale_x_continuous(limits=c(0,2.5*S0)) +
             scale_y_continuous(limits=c(0,2.5*cb$notional)) +
             geom_point(aes(x=S0,y=present_value), color="red") +
             labs(title="Convertible Bond Value")
)
cbplot

## ----twitter_cb_via_BondValuation, echo=TRUE, comment=""-----------------
twitter_bv = BondValuation::AnnivDates(
    Em=as.Date('2018-06-11'),    # Issue date
    Mat=as.Date('2024-06-15'), 
    CpY=2, 
    FIPD=as.Date('2018-12-15'),   # First coupon
    FIAD=as.Date('2018-06-15'),   # Beginning of first coupon accrual
    RV=1000,   # Notional
    Coup=0.25, 
    DCC=which(BondValuation::List.DCC$DCC.Name=='30/360'),  # 30/360 daycount convention
    EOM=0
  )

twitter_specs = ragtop::detail_from_AnnivDates(
    twitter_bv, 
    as_of=as.Date('2018-02-15')
  )

twtr_cb = ragtop::ConvertibleBond(
    maturity=twitter_specs$maturity, 
    conversion_ratio=17.5001, 
    notional=twitter_specs$notional,
    coupons=twitter_specs$coupons,
    discount_factor_fcn = disct_fcn,
    name='TwitterConvertWithGreenshoe'
  )

pvs = ragtop::find_present_value(
    S0=33.06,
    num_time_steps=200,
    instruments=list(TWTR=twtr_cb),
    const_volatility=0.47,
    const_default_intensity=0.01,
    discount_factor_fcn=disct_fcn,
  )

paste("Twitter bond value is", pvs$TWTR)

