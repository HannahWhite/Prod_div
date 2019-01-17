
###############################################################################
####### Productivity-diversity model for entire dataset - area included #######
###############################################################################

## HAnnah White 14.11.2018

## This code carries out the productivity diversity relationships but includes an area term to account for the species-area relationship. The z value 0.08 was determined from an nls

#### Non-smoothed area with z parameter included



### 

## Prepare data

library(mgcv)
library(ggplot2)
library(cowplot)
library(visreg)
library(plyr)

### Read in data


vasc.env <- read.csv('G:/Postdoc Grassland Resilience/Species richness/vasc.env.csv', header = TRUE)

# Add land cover data

corine.agg <- read.csv('G:\\Postdoc Grassland Resilience\\LandCoverData\\corine.agg.csv', header = TRUE)
corine.agg$east <- corine.agg$east + 5000
corine.agg$north <- corine.agg$north + 5000

# merge dataframes
vasc.env <- merge(vasc.env, corine.agg, by = c('east', 'north', 'dom'), all.x = TRUE, all.y = TRUE)



## dominant habitat function
dominant <- function (x, threshold) {
  ifelse(max(x) >=threshold, names(which.max(x)), 'Heterogenous') 
}

vasc.env$dom50 <- apply(vasc.env[,52:66], 1, function(x) dominant(x, threshold = 0.5))

vasc.env$hetero.dom50 <- ifelse(vasc.env$dom50 == 'Heterogenous', 'H', 'N')

vasc.env$pasture.dom50 <- ifelse(vasc.env$dom50 == 'Pasture', 'P', 'N')


vasc.env$hetpast <- ifelse(vasc.env$dom50 == 'Heterogenous', 'H',
                           ifelse(vasc.env$dom50 == 'Pasture', 'P', 'N'))

vasc.env$non.past <- 1 - vasc.env$Pasture

vasc.env$non.past.area <- 100*vasc.env$non.past

vasc.env$area.0.08 <- vasc.env$non.past.area^0.08


rm(corine.agg)


###


m.full.area <- 
  gamm(nat.fres ~ s(evi.mean,  bs = 'cs', k = -1) + 
         s(evi.sd, bs = 'cs', k = -1) +
         s(inter.evi.sd, bs = 'cs', k = -1) +
         s(annual.evi.sd, bs = 'cs', k = -1) +
         te(east, north, bs = c('tp', 'tp')) +
         area.0.08, 
       correlation = corExp(form = ~ east + north), data = vasc.env, 
       method = 'ML', control = lmeControl(opt = 'optim'))

## AIC = 6799.005



m.full.exp <- 
  gamm(nat.fres ~ s(evi.mean,  bs = 'cs', k = -1) + 
         s(evi.sd, bs = 'cs', k = -1) +
         s(inter.evi.sd, bs = 'cs', k = -1) +
         s(annual.evi.sd, bs = 'cs', k = -1) +
         te(east, north, bs = c('tp', 'tp')), 
       correlation = corExp(form = ~ east + north), data = vasc.env, 
       method = 'ML', control = lmeControl(opt = 'optim'))


## AIC = 6843.558




m.full.dom3level.area <- gamm(nat.fres ~ s(evi.mean,  bs = 'cs', by = as.factor(hetpast), k = -1) + 
                               s(evi.sd, bs = 'cs', by = as.factor(hetpast), k = -1) +
                               s(inter.evi.sd, bs = 'cs', by = as.factor(hetpast), k = -1) +
                               s(annual.evi.sd, bs = 'cs', by = as.factor(hetpast), k = -1) +
                               te(east, north, bs = c('tp', 'tp')) + hetpast +
                               area.0.08, 
                             correlation = corExp(form = ~ east + north), data = vasc.env, 
                             method = 'ML', control = lmeControl(opt = 'optim'))


## AIC = 6788.297


#### Allow GAM to smooth area term



m.full.area2 <- 
  gamm(nat.fres ~ s(evi.mean,  bs = 'cs', k = -1) + 
         s(evi.sd, bs = 'cs', k = -1) +
         s(inter.evi.sd, bs = 'cs', k = -1) +
         s(annual.evi.sd, bs = 'cs', k = -1) +
         te(east, north, bs = c('tp', 'tp')) +
         s(non.past.area, bs = 'cs', k = -1), 
       correlation = corExp(form = ~ east + north), data = vasc.env, 
       method = 'ML', control = lmeControl(opt = 'optim'))

## AIC = 6733.796


m.full.area3 <-  # same as previous model but reduces number of knots
  gamm(nat.fres ~ s(evi.mean,  bs = 'cs', k = -1) + 
         s(evi.sd, bs = 'cs', k = -1) +
         s(inter.evi.sd, bs = 'cs', k = -1) +
         s(annual.evi.sd, bs = 'cs', k = -1) +
         te(east, north, bs = c('tp', 'tp')) +
         s(non.past.area, bs = 'cs', k = 3), 
       correlation = corExp(form = ~ east + north), data = vasc.env, 
       method = 'ML', control = lmeControl(opt = 'optim'))

## AIC = 6754.801



m.full.dom3level.area2 <- gamm(nat.fres ~ s(evi.mean,  bs = 'cs', by = as.factor(hetpast), k = -1) + 
                                s(evi.sd, bs = 'cs', by = as.factor(hetpast), k = -1) +
                                s(inter.evi.sd, bs = 'cs', by = as.factor(hetpast), k = -1) +
                                s(annual.evi.sd, bs = 'cs', by = as.factor(hetpast), k = -1) +
                                te(east, north, bs = c('tp', 'tp')) + hetpast +
                                s(non.past.area, bs ='cs', k = -1), # doesn't converge when k set to 3
                              correlation = corExp(form = ~ east + north), data = vasc.env, 
                              method = 'ML', control = lmeControl(opt = 'optim')) 
## AIC = 6745.789







## check order

m.full.dom3level.areatop <- gamm(nat.fres ~ s(non.past.area, bs = 'cs', k = -1) +
                                 s(evi.mean,  bs = 'cs', by = as.factor(hetpast), k = -1) + 
                                 s(evi.sd, bs = 'cs', by = as.factor(hetpast), k = -1) +
                                 s(inter.evi.sd, bs = 'cs', by = as.factor(hetpast), k = -1) +
                                 s(annual.evi.sd, bs = 'cs', by = as.factor(hetpast), k = -1) +
                                 te(east, north, bs = c('tp', 'tp')) + hetpast, 
                               correlation = corExp(form = ~ east + north), data = vasc.env, 
                               method = 'ML', control = lmeControl(opt = 'optim')) ## does not converge

## AIC = 6745.789


## Order covariates are in model does not make a difference

## null

m0 <- gam(nat.fres ~  1 + te(east, north, bs = c('tp', 'tp'), sp = m.full.dom3level.area2$sp[13:14]),
                             correlation = corExp(form = ~ east + north), data = vasc.env, method = 'ML')




## Contribution of each whilst maintaining original smooth - I think this might be impossible

m.nomean <- gamm(nat.fres ~ # nope this isn't maintaining smooth
                   s(evi.sd, bs = 'cs', by = as.factor(hetpast), k = -1) +
                   s(inter.evi.sd, bs = 'cs', by = as.factor(hetpast), k = -1) +
                   s(annual.evi.sd, bs = 'cs', by = as.factor(hetpast), k = -1) +
                   te(east, north, bs = c('tp', 'tp')) + hetpast +
                   s(non.past.area, bs ='cs', k = -1), 
                 correlation = corExp(form = ~ east + north), data = vasc.env, 
                 method = 'ML', control = lmeControl(opt = 'optim'), sp = m.full.dom3level.area2$gam$sp[4:16])
# AIc = 6739.612




#m.nomean <- gam(nat.fres ~  ## this is wrong as gam ignores the correlation structure
#                  s(evi.sd, bs = 'cs', by = as.factor(hetpast), k = -1, sp = m.full.dom3level.area2$sp[4:6]) +
#                  s(inter.evi.sd, bs = 'cs', by = as.factor(hetpast), k = -1, sp = m.full.dom3level.area2$sp[7:9]) +
#                  s(annual.evi.sd, bs = 'cs', by = as.factor(hetpast), k = -1, sp = m.full.dom3level.area2$sp[10:12]) +
#                  te(east, north, bs = c('tp', 'tp'), sp = m.full.dom3level.area2$sp[13:14]) + hetpast +
#                  s(non.past.area, bs ='cs', k = -1, sp = m.full.dom3level.area2$sp[15]), 
#                correlation = corExp(form = ~ east + north), data = vasc.env, 
#                method = 'ML')
# AIc = 6961.678


m.nomean.nonspec <- gamm(nat.fres ~ 
                           s(evi.sd, bs = 'cs', by = as.factor(hetpast), k = -1) +
                           s(inter.evi.sd, bs = 'cs', by = as.factor(hetpast), k = -1) +
                           s(annual.evi.sd, bs = 'cs', by = as.factor(hetpast), k = -1) +
                           te(east, north, bs = c('tp', 'tp')) + hetpast +
                           s(non.past.area, bs ='cs', k = -1), 
                         correlation = corExp(form = ~ east + north), data = vasc.env, 
                         method = 'ML', control = lmeControl(opt = 'optim'))
# AIc = 6739.612





(deviance(m.nomean$gam)-deviance(m.full.dom3level.area2$gam))/deviance(m0)




















