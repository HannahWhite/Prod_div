
###############################################################################
####### Productivity-diversity model for entire dataset - area included #######
###############################################################################

## Hannah White 14.11.2018

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


vasc.env <- read.csv('vascular_evi.csv', header = TRUE)


## dominant habitat function
dominant <- function (x, threshold) {
  ifelse(max(x) >=threshold, names(which.max(x)), 'Heterogenous') 
}

vasc.env$dom50 <- apply(vasc.env[,8:22], 1, function(x) dominant(x, threshold = 0.5))

vasc.env$hetero.dom50 <- ifelse(vasc.env$dom50 == 'Heterogenous', 'H', 'N')

vasc.env$pasture.dom50 <- ifelse(vasc.env$dom50 == 'Pasture', 'P', 'N')


vasc.env$hetpast <- ifelse(vasc.env$dom50 == 'Heterogenous', 'H',
                           ifelse(vasc.env$dom50 == 'Pasture', 'P', 'N'))

vasc.env$non.past <- 1 - vasc.env$Pasture

vasc.env$non.past.area <- 100*vasc.env$non.past



rm(corine.agg)




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



m.full.dom3level.area2 <- gamm(nat.fres ~ s(evi.mean,  bs = 'cs', by = as.factor(hetpast), k = -1) + 
                                s(evi.sd, bs = 'cs', by = as.factor(hetpast), k = -1) +
                                s(inter.evi.sd, bs = 'cs', by = as.factor(hetpast), k = -1) +
                                s(annual.evi.sd, bs = 'cs', by = as.factor(hetpast), k = -1) +
                                te(east, north, bs = c('tp', 'tp')) + hetpast +
                                s(non.past.area, bs ='cs', k = -1), # doesn't converge when k set to 3
                              correlation = corExp(form = ~ east + north), data = vasc.env, 
                              method = 'ML', control = lmeControl(opt = 'optim')) 
## AIC = 6745.789




























