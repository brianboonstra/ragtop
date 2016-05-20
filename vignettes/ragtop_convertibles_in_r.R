## ----global_options, include=FALSE---------------------------------------
library(ragtop)
library(futile.logger)
library(ggplot2)
library(reshape2)
library(stringr)
library(MASS)

flog.threshold(ERROR)
flog.threshold(ERROR, name='ragtop.implicit.timestep.construct_tridiagonals')
flog.threshold(ERROR, name='ragtop.calibration.implied_volatility.lowprice')
flog.threshold(ERROR, name='ragtop.calibration.implied_volatility_with_term_struct')
flog.threshold(ERROR, name='ragtop.implicit.setup.width')

knitr::opts_chunk$set(fig.width=6.5, fig.height=4, fig.path='Figs/',
                      echo=FALSE, warning=FALSE, message=FALSE, comment=FALSE)
# PALETTES: see http://www.r-bloggers.com/the-paul-tol-21-color-salute/
# FOCUS PALETTES
# Red as highlight
redfocus = c("#CB181D", "#252525", "#525252", "#737373", "#969696", "#BDBDBD", "#D9D9D9", "#F0F0F0")
 
# Green as highlight
greenfocus = c("#41AB5D", "#252525", "#525252", "#737373", "#969696", "#BDBDBD", "#D9D9D9", "#F0F0F0")
 
# Blue as highlight
bluefocus = c("#0033FF", "#252525", "#525252", "#737373", "#969696", "#BDBDBD", "#D9D9D9", "#F0F0F0")
 
# EQUAL WEIGHT
# Generated with rainbow(12, s = 0.6, v = 0.75)
rainbow12equal = c("#BF4D4D", "#BF864D", "#BFBF4D", "#86BF4D", "#4DBF4D", "#4DBF86", "#4DBFBF", "#4D86BF", "#4D4DBF", "#864DBF", "#BF4DBF", "#BF4D86")
rainbow10equal = c("#BF4D4D", "#BF914D", "#A8BF4D", "#63BF4D", "#4DBF7A", "#4DBFBF", "#4D7ABF", "#634DBF", "#A84DBF", "#BF4D91")
rainbow8equal = c("#BF4D4D", "#BFA34D", "#86BF4D", "#4DBF69", "#4DBFBF", "#4D69BF", "#864DBF", "#BF4DA3")
rainbow6equal = c("#BF4D4D", "#BFBF4D", "#4DBF4D", "#4DBFBF", "#4D4DBF", "#BF4DBF")
 
# Generated with package "gplots" function rich.colors(12)
rich12equal = c("#000040", "#000093", "#0020E9", "#0076FF", "#00B8C2", "#04E466", "#49FB25", "#E7FD09", "#FEEA02", "#FFC200", "#FF8500", "#FF3300")
rich10equal = c("#000041", "#0000A9", "#0049FF", "#00A4DE", "#03E070", "#5DFC21", "#F6F905", "#FFD701", "#FF9500", "#FF3300")
rich8equal = c("#000041", "#0000CB", "#0081FF", "#02DA81", "#80FE1A", "#FDEE02", "#FFAB00", "#FF3300")
rich6equal = c("#000043", "#0033FF", "#01CCA4", "#BAFF12", "#FFCC00", "#FF3300")
 
# Generated with package "fields" function tim.colors(12), which is said to emulate the default matlab colorset
tim12equal = c("#00008F", "#0000EA", "#0047FF", "#00A2FF", "#00FEFF", "#5AFFA5", "#B5FF4A", "#FFED00", "#FF9200", "#FF3700", "#DB0000", "#800000")
tim10equal = c("#00008F", "#0000FF", "#0070FF", "#00DFFF", "#50FFAF", "#BFFF40", "#FFCF00", "#FF6000", "#EF0000", "#800000")
tim8equal = c("#00008F", "#0020FF", "#00AFFF", "#40FFBF", "#CFFF30", "#FF9F00", "#FF1000", "#800000")
tim6equal = c("#00008F", "#005AFF", "#23FFDC", "#ECFF13", "#FF4A00", "#800000")
 
# Generated with sort(brewer.pal(8,"Dark2")) #Dark2, Set2
dark8equal = c("#1B9E77", "#666666", "#66A61E", "#7570B3", "#A6761D", "#D95F02", "#E6AB02", "#E7298A")
dark6equal = c("#1B9E77", "#66A61E", "#7570B3", "#D95F02", "#E6AB02", "#E7298A")
set8equal = c("#66C2A5", "#8DA0CB", "#A6D854", "#B3B3B3", "#E5C494", "#E78AC3", "#FC8D62", "#FFD92F")
set6equal = c("#66C2A5", "#8DA0CB", "#A6D854", "#E78AC3", "#FC8D62", "#FFD92F")
 


## ----show_TSLA_S0, echo=TRUE---------------------------------------------
TSLAMarket$S0

## ----show_TSLA_rf, echo=TRUE---------------------------------------------
TSLAMarket$risk_free_rates

## ---- results='asis'-----------------------------------------------------
knitr::kable(TSLAMarket$options[c(200,300,400,500,600, 800),], digits=3, row.names = F)

## ----blackscholes, echo=TRUE---------------------------------------------
blackscholes(TSLAMarket$options[500,'callput'], 
             TSLAMarket$S0, 
             TSLAMarket$options[500,'K'], 
             0.005, 
             TSLAMarket$options[500,'time'], 
             0.50)

## ----implied_volatility, echo=TRUE---------------------------------------
implied_volatility(option_price = TSLAMarket$options[400,'ask'], 
                   S0 = TSLAMarket$S0, 
                   callput = TSLAMarket$options[400,'callput'], 
                   K=TSLAMarket$options[400,'K'], 
                   r = 0.005, 
                   time = TSLAMarket$options[400,'time'])

## ----amer, echo=TRUE-----------------------------------------------------
american(
       callput = TSLAMarket$options[400,'callput'], 
       S0 = TSLAMarket$S0, 
       K=TSLAMarket$options[400,'K'], 
       const_short_rate = 0.005, 
       time = TSLAMarket$options[400,'time'])

## ----american_implied_volatility, echo=TRUE------------------------------
american_implied_volatility(option_price = TSLAMarket$options[400,'ask'], 
     S0 = TSLAMarket$S0, 
     callput = TSLAMarket$options[400,'callput'], 
     K=TSLAMarket$options[400,'K'], 
     const_short_rate = 0.005, 
     time = TSLAMarket$options[400,'time'])

## ----implied_volatility_def, echo=TRUE-----------------------------------
implied_volatility(option_price = TSLAMarket$options[400,'ask'], 
                   S0 = TSLAMarket$S0, 
                   callput = TSLAMarket$options[400,'callput'], 
                   K=TSLAMarket$options[400,'K'], 
                   r = 0.005, 
                   time = TSLAMarket$options[400,'time'],
                   const_default_intensity = 0.03)

## ----american_implied_volatility_def, echo=TRUE--------------------------
american_implied_volatility(option_price = TSLAMarket$options[400,'ask'], 
     S0 = TSLAMarket$S0, 
     callput = TSLAMarket$options[400,'callput'], 
     K=TSLAMarket$options[400,'K'], 
     const_short_rate = 0.005, 
     time = TSLAMarket$options[400,'time'],
     const_default_intensity = 0.0200)

