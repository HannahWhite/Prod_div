#######################################################################################################
### Models for Productivity-Diversity with shrinkage and exponential correlation for spatial errors ###
#######################################################################################################

### Hannah White 24.09.2018

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


rm(corine.agg)




############################################################################
####### Productivity-diversity model for entire dataset - no habitat #######
############################################################################


m.full.exp <- 
  gamm(nat.fres ~ s(evi.mean,  bs = 'cs', k = -1) + 
         s(evi.sd, bs = 'cs', k = -1) +
         s(inter.evi.sd, bs = 'cs', k = -1) +
         s(annual.evi.sd, bs = 'cs', k = -1) +
         te(east, north, bs = c('tp', 'tp')), 
       correlation = corExp(form = ~ east + north), data = vasc.env, 
       method = 'ML', control = lmeControl(opt = 'optim'))
## AIC = 6843.558




# create a sequence of evi that spans evi.mean
maxevi<-max(vasc.env$evi.mean)
minevi<-min(vasc.env$evi.mean)
evi.seq<-seq(minevi, maxevi, length=300)
evi.seq<-data.frame(evi.mean=evi.seq, evi.sd = rep(0.11, 300), inter.evi.sd = rep(0.016, 300), annual.evi.sd = rep(0.097, 300),
                    east = rep(205096.4, 300), north = rep(239771.1, 300))

# predict only the evi.mean term (the sum of the
# term predictions and the intercept gives you the overall 
# prediction)

preds.mean<-predict(m.full.exp$gam, type="terms", newdata=evi.seq,
                    se.fit=TRUE)

# set up the temperature, the fit and the upper and lower
# confidence interval

evi<-evi.seq$evi.mean
fit<-preds.mean$fit[,1]
fit.up95<-fit-1.96*preds.mean$se.fit[,1]
fit.low95<-fit+1.96*preds.mean$se.fit[,1]


## plot in ggplot2 - this code can be adjusted for presentations

evi.df <- data.frame(evi = evi, fit = fit, fit.up95 = fit.up95, fit.low95 = fit.low95)

p1 <- ggplot(evi.df, aes(x = evi, y = fit)) + geom_line(aes(x=evi, y=fit), size = 1.5) + 
  geom_ribbon(aes(ymin = fit.low95, ymax = fit.up95, x = evi, fill = 'band', alpha = 0.3)) +
  xlab('Mean EVI') + ylab('Standardised Species Richness') + 
  scale_colour_manual('', values = 'black', guide = 'none') + scale_fill_manual('', values = 'black', guide = 'none') + 
  scale_alpha(guide = 'none')
p1 <- p1 + theme(axis.title = element_text(size = 18), axis.text = element_text(size = 18))

# get partial residuals
pred.orig <- predict(m.full.exp$gam, type = 'terms')
pr.full <- residuals(m.full.exp$gam, type = 'working') + pred.orig[,1]

orig.df <- data.frame(evi.mean = vasc.env$evi.mean, pr.full = as.numeric(pr.full))

p1 <- p1 + geom_point(data = orig.df, mapping = aes(x = evi.mean, y = pr.full))





#### Can now redo this for other predictors.....

### evi.sd

# create a sequence of evi that spans evi.mean
maxsd<-max(vasc.env$evi.sd)
minsd<-min(vasc.env$evi.sd)
sd.seq<-seq(minsd, maxsd, length=300)
sd.seq<-data.frame(evi.mean=rep(0.5, 300), evi.sd = sd.seq, inter.evi.sd = rep(0.016, 300), annual.evi.sd = rep(0.097, 300), east = rep(205096.4, 300), north = rep(239771.1, 300))

# predict only the evi.mean term (the sum of the
# term predictions and the intercept gives you the overall 
# prediction)

preds.sd<-predict(m.full.exp$gam, type="terms", newdata=sd.seq,
                  se.fit=TRUE)

# set up the sd, the fit and the upper and lower
# confidence interval

evi.sd<-sd.seq$evi.sd
fit<-preds.sd$fit[,2]
fit.up95<-fit-1.96*preds.sd$se.fit[,2]
fit.low95<-fit+1.96*preds.sd$se.fit[,2]





sd.df <- data.frame(evi.sd = evi.sd, fit = fit, fit.up95 = fit.up95, fit.low95 = fit.low95)

p1.sd <- ggplot(sd.df, aes(x = evi.sd, y = fit)) + geom_line(aes(x=evi.sd, y=fit, colour = 'darkorchid4'), size = 1.5) + 
  geom_ribbon(aes(ymin = fit.low95, ymax = fit.up95, x = evi.sd, fill = 'band', alpha = 0.3)) +
  xlab('Spatial SD of EVI') + ylab('Standardised Species Richness') +
  scale_colour_manual('', values = 'black', guide = 'none') + scale_fill_manual('', values = 'black', guide = 'none') + 
  scale_alpha(guide = 'none')
p1.sd <- p1.sd + theme(axis.title = element_text(size = 18), axis.text = element_text(size = 18))

# get partial residuals
pred.orig <- predict(m.full.exp$gam, type = 'terms')
pr.full <- residuals(m.full.exp$gam, type = 'working') + pred.orig[,2]

orig.df <- data.frame(evi.sd = vasc.env$evi.sd, pr.full = as.numeric(pr.full))

p1.sd <- p1.sd + geom_point(data = orig.df, mapping = aes(x = evi.sd, y = pr.full))





### inter.evi.sd

# create a sequence of evi that spans inter.evi.sd
maxsd<-max(vasc.env$inter.evi.sd)
minsd<-min(vasc.env$inter.evi.sd)
inter.sd.seq<-seq(minsd, maxsd, length=300)
inter.sd.seq<-data.frame(evi.mean=rep(0.5, 300), evi.sd = rep(0.11, 300), inter.evi.sd = inter.sd.seq, annual.evi.sd = rep(0.097, 300), east = rep(205096.4, 300), north = rep(239771.1, 300))

# predict only the inter.evi.sd term (the sum of the
# term predictions and the intercept gives you the overall 
# prediction)

preds.sd<-predict(m.full.exp$gam, type="terms", newdata=inter.sd.seq,
                  se.fit=TRUE)

# set up the sd, the fit and the upper and lower
# confidence interval

inter.evi.sd<-inter.sd.seq$inter.evi.sd
fit<-preds.sd$fit[,3]
fit.up95<-fit-1.96*preds.sd$se.fit[,3]
fit.low95<-fit+1.96*preds.sd$se.fit[,3]





inter.sd.df <- data.frame(inter.evi.sd = inter.evi.sd, fit = fit, fit.up95 = fit.up95, fit.low95 = fit.low95)

p1.intersd <- ggplot(inter.sd.df, aes(x = inter.evi.sd, y = fit)) + geom_line(aes(x=inter.evi.sd, y=fit, colour = 'darkorchid4'), size = 1.5) + 
  geom_ribbon(aes(ymin = fit.low95, ymax = fit.up95, x = inter.evi.sd, fill = 'band', alpha = 0.3)) +
  xlab('Interannual SD of EVI') + ylab('Standardised Species Richness') +
  scale_colour_manual('', values = 'black', guide = 'none') + scale_fill_manual('', values = 'black', guide = 'none') + 
  scale_alpha(guide = 'none')
p1.intersd <- p1.intersd + theme(axis.title = element_text(size = 18), axis.text = element_text(size = 18))

# get partial residuals
pred.orig <- predict(m.full.exp$gam, type = 'terms')
pr.full <- residuals(m.full.exp$gam, type = 'working') + pred.orig[,3]

orig.df <- data.frame(inter.evi.sd = vasc.env$inter.evi.sd, pr.full = as.numeric(pr.full))

p1.intersd <- p1.intersd + geom_point(data = orig.df, mapping = aes(x = inter.evi.sd, y = pr.full))



### annual sd


# create a sequence of evi that spans annual.evi.sd
maxsd<-max(vasc.env$annual.evi.sd)
minsd<-min(vasc.env$annual.evi.sd)
annual.sd.seq<-seq(minsd, maxsd, length=300)
annual.sd.seq<-data.frame(evi.mean=rep(0.5, 300), evi.sd = rep(0.11, 300), inter.evi.sd = rep(0.016, 300), annual.evi.sd = annual.sd.seq, east = rep(205096.4, 300), north = rep(239771.1, 300))

# predict only the annual.evi.sd term (the sum of the
# term predictions and the intercept gives you the overall 
# prediction)

preds.sd<-predict(m.full.exp$gam, type="terms", newdata=annual.sd.seq,
                  se.fit=TRUE)

# set up the sd, the fit and the upper and lower
# confidence interval

annual.evi.sd<-annual.sd.seq$annual.evi.sd
fit<-preds.sd$fit[,4]
fit.up95<-fit-1.96*preds.sd$se.fit[,4]
fit.low95<-fit+1.96*preds.sd$se.fit[,4]


annual.sd.df <- data.frame(annual.evi.sd = annual.evi.sd, fit = fit, fit.up95 = fit.up95, fit.low95 = fit.low95)

p1.annualsd <- ggplot(annual.sd.df, aes(x = annual.evi.sd, y = fit)) + geom_line(aes(x=annual.evi.sd, y=fit, colour = 'darkorchid4'), size = 1.5) + 
  geom_ribbon(aes(ymin = fit.low95, ymax = fit.up95, x = annual.evi.sd, fill = 'band', alpha = 0.3)) +
  xlab('Annual SD of EVI') + ylab('Standardised Species Richness') +
  scale_colour_manual('', values = 'black', guide = 'none') + scale_fill_manual('', values = 'black', guide = 'none') + 
  scale_alpha(guide = 'none')
p1.annualsd <- p1.annualsd + theme(axis.title = element_text(size = 18), axis.text = element_text(size = 18))

# get partial residuals
pred.orig <- predict(m.full.exp$gam, type = 'terms')
pr.full <- residuals(m.full.exp$gam, type = 'working') + pred.orig[,4]

orig.df <- data.frame(annual.evi.sd = vasc.env$annual.evi.sd, pr.full = as.numeric(pr.full))

p1.annualsd <- p1.annualsd + geom_point(data = orig.df, mapping = aes(x = annual.evi.sd, y = pr.full))



plot_grid(p1, p1.sd, p1.intersd, p1.annualsd, nrow = 2, ncol = 2)




################################################################################
####### Productivity-diversity model smoothed by three levels of habitat #######
################################################################################

m.full.dom3level.exp <- gamm(nat.fres ~ s(evi.mean,  bs = 'cs', by = as.factor(hetpast), k = -1) + 
                           s(evi.sd, bs = 'cs', by = as.factor(hetpast), k = -1) +
                           s(inter.evi.sd, bs = 'cs', by = as.factor(hetpast), k = -1) +
                           s(annual.evi.sd, bs = 'cs', by = as.factor(hetpast), k = -1) +
                           te(east, north, bs = c('tp', 'tp')) + hetpast, 
                         correlation = corExp(form = ~ east + north), data = vasc.env, 
                         method = 'ML', control = lmeControl(opt = 'optim'))

## AIC = 6817.637 Rs1 = 0.681

###########################################
###### Graphics for land cover model ######
###########################################

##### mean EVI
### simulate data with each land cover in turn to then combine on to same graphic

# Simulate data with heterogeneous habitat
# add dominant habitat column to evi.seq

# create a sequence of evi that spans evi.mean
maxevi<-max(vasc.env$evi.mean)
minevi<-min(vasc.env$evi.mean)
evi.seq<-seq(minevi, maxevi, length=300)
evi.seq<-data.frame(evi.mean=evi.seq, evi.sd = rep(0.11, 300), inter.evi.sd = rep(0.016, 300), annual.evi.sd = rep(0.097, 300), east = rep(205096.4, 300), north = rep(239771.1, 300))

evi.seq$hetpast <- rep('H', 300) # can be changed to N or P for other graphs


# predict only the evi.mean term (the sum of the
# term predictions and the intercept gives you the overall 
# prediction)

preds.mean.dom<-predict(m.full.dom3level.exp$gam, type="terms", newdata=evi.seq,
                        se.fit=TRUE)

# set up evi.mean, the fit and the upper and lower
# confidence interval

evi<-evi.seq$evi.mean
fit<-preds.mean.dom$fit[,2] ### eek what number should this be?
fit.up95<-fit-1.96*preds.mean.dom$se.fit[,2] 
fit.low95<-fit+1.96*preds.mean.dom$se.fit[,2]

evi.df <- data.frame(evi = evi, fit = fit, fit.up95 = fit.up95, fit.low95 = fit.low95)


# get partial residuals
pred.orig <- predict(m.full.dom3level.exp$gam, type = 'terms')
pr.3level.dom <- residuals(m.full.dom3level.exp$gam, type = 'working') + pred.orig[,2]

orig.df <- data.frame(evi.mean = vasc.env$evi.mean, pr.3level.dom = as.numeric(pr.3level.dom))
orig.df$dom <- vasc.env$hetpast


## Simulate data with 'other' as dominant habitat

evi.seq$hetpast <- rep('N', 300) # 

preds.mean.domN<-predict(m.full.dom3level.exp$gam, type="terms", newdata=evi.seq,
                         se.fit=TRUE)

# set up evi.mean, the fit and the upper and lower
# confidence interval

evi<-evi.seq$evi.mean
fit<-preds.mean.domN$fit[,3]
fit.up95N<-fit-1.96*preds.mean.domN$se.fit[,3] 
fit.low95N<-fit+1.96*preds.mean.domN$se.fit[,3]

evi.dfN <- data.frame(evi = evi, fit = fit, fit.up95 = fit.up95N, fit.low95 = fit.low95N)


# get partial residuals
pred.origN <- predict(m.full.dom3level.exp$gam, type = 'terms')
pr.3level.domN <- residuals(m.full.dom3level.exp$gam, type = 'working') + pred.origN[,3]

orig.dfN <- data.frame(evi.mean = vasc.env$evi.mean, pr.3level.domN = as.numeric(pr.3level.domN))
orig.dfN$dom <- vasc.env$hetpast


## Simulate data with pasture as dominant habitat

evi.seq$hetpast <- rep('P', 300)

preds.mean.domP <- predict(m.full.dom3level.exp$gam, type = 'terms', newdata = evi.seq, se.fit = TRUE)

# set up to calcualte confidence intervals

evi <- evi.seq$evi.mean
fit <- preds.mean.domP$fit[,4]
fit.up95P <- fit - 1.96*preds.mean.domP$se.fit[,4]
fit.low95P <- fit + 1.96*preds.mean.domP$se.fit[,4]

evi.dfP <- data.frame(evi = evi, fit = fit, fit.up95 = fit.up95P, fit.low95 = fit.low95P)

# get partial residuals
pred.origP <- predict(m.full.dom3level.exp$gam, type = 'terms')
pr.3level.domP <- residuals(m.full.dom3level.exp$gam, type = 'working') + pred.origP[,4]

orig.dfP <- data.frame(evi.mean = vasc.env$evi.mean, pr.3level.domP = as.numeric(pr.3level.domP))
orig.dfP$dom <- vasc.env$hetpast




### Build graphic that combines both heterogenous and homogenous landscapes

# add to residual dfs whether habitat is heterogenous or not

orig.df$dom <- vasc.env$hetpast
orig.dfN$dom <- vasc.env$hetpast
orig.dfP$dom <- vasc.env$hetpast


# extract range of EVIs where each dominant habitat is found

ve.H <- vasc.env[which(vasc.env$hetpast=='H'),]
ve.H <- droplevels(ve.H)
range(ve.H$evi.mean) # 0.30 0.58

ve.N <- vasc.env[which(vasc.env$hetpast=='N'),]
ve.N <- droplevels(ve.N)
range(ve.N$evi.mean) # 0.25 0.50

ve.P <- vasc.env[which(vasc.env$hetpast=='P'),]
ve.P <- droplevels(ve.P)
range(ve.P$evi.mean) # 0.44 0.62


evi.df.Hless <- evi.df[evi.df$evi >= 0.3 & evi.df$evi <= 0.58,]
evi.df.Hless <- droplevels(evi.df.Hless)

evi.df.Nless <- evi.dfN[evi.dfN$evi >= 0.24 & evi.dfN$evi <= 0.5,]
evi.df.Nless <- droplevels(evi.df.Nless)

evi.df.Pless <- evi.dfP[evi.dfP$evi >= 0.44 & evi.dfP$evi <= 0.62,]
evi.df.Pless <- droplevels(evi.df.Pless)



p.nice <- ggplot(evi.df.Hless, aes(x = evi, y = fit)) + geom_line(aes(x = evi, y = fit, colour = 'darkorchid4'), size = 1.5) +
  geom_ribbon(aes(ymin = fit.low95, ymax = fit.up95, x = evi, fill = 'band', alpha = 0.1), fill = 'darkorchid4') +
  xlab('Mean EVI') + ylab('Standardised Species Richness') +
  xlim(0.24, 0.625)


p.nice <- p.nice + geom_line(data = evi.df.Pless, aes(x = evi, y = fit, colour = 'forestgreen'), size = 1.5) +
  geom_ribbon(data = evi.df.Pless, aes(ymin = fit.low95, ymax = fit.up95, x = evi, fill = 'band', alpha = 0.1), fill = 'forestgreen')

p.nice <- p.nice + geom_line(data = evi.df.Nless, aes(x = evi, y = fit, colour = 'darkorange3'), size = 1.5) +
  geom_ribbon(data = evi.df.Nless, aes(ymin = fit.low95, ymax = fit.up95, x = evi, fill = 'band', alpha = 0.1), fill = 'darkorange3')

p.nice <- p.nice + geom_point(data = orig.df[orig.df$dom=='H',], mapping = aes(x = evi.mean, y = pr.3level.dom), colour = 'darkorchid4', size = 3) # plot ast so top layer


p.nice <- p.nice + scale_colour_manual(values = c('darkorange3'='darkorange3', 'darkorchid4'='darkorchid4', 'forestgreen' = 'forestgreen'),guide = guide_legend(title = NULL),
                                       labels = c('Other', 'Heterogeneous', 'Pasture')) + 
  scale_fill_manual(guide=  'none') + 
  scale_alpha(guide = 'none')

p.nice <- p.nice + geom_point(data = orig.dfN[orig.dfN$dom=='N',], mapping = aes(x = evi.mean, y = pr.3level.domN), colour = 'darkorange3', size = 3)
p.nice <- p.nice + geom_point(data = orig.dfP[orig.dfP$dom=='P',], mapping = aes(x = evi.mean, y = pr.3level.domP), colour = 'forestgreen', size = 3)

p.nice <- p.nice + theme(axis.title = element_text(size = 25), axis.text = element_text(size = 20), legend.position = 'none') 




##### EVI sd
### simulate data with each land cover in turn to then combine on to same graphic

# Simulate data with heterogeneous habitat
# add dominant habitat column to evi.seq

# create a sequence of evi that spans evi.mean
maxevi<-max(vasc.env$evi.sd)
minevi<-min(vasc.env$evi.sd)
evi.seq<-seq(minevi, maxevi, length=300)
evi.seq<-data.frame(evi.mean = rep(0.5, 300), evi.sd = evi.seq, inter.evi.sd = rep(0.016, 300), annual.evi.sd = rep(0.097, 300), east = rep(205096.4, 300), north = rep(239771.1, 300))

evi.seq$hetpast <- rep('H', 300) # can be changed to N or P for other graphs


# predict only the evi.sd term (the sum of the
# term predictions and the intercept gives you the overall 
# prediction)

preds.sd.dom<-predict(m.full.dom3level.exp$gam, type="terms", newdata=evi.seq,
                      se.fit=TRUE)

# set up evi.sd, the fit and the upper and lower
# confidence interval

evi<-evi.seq$evi.sd
fit<-preds.sd.dom$fit[,5] 
fit.up95<-fit-1.96*preds.sd.dom$se.fit[,5] 
fit.low95<-fit+1.96*preds.sd.dom$se.fit[,5]

evi.df <- data.frame(evi = evi, fit = fit, fit.up95 = fit.up95, fit.low95 = fit.low95)


# get partial residuals
pred.orig <- predict(m.full.dom3level.exp$gam, type = 'terms')
pr.3level.dom <- residuals(m.full.dom3level.exp$gam, type = 'working') + pred.orig[,5]

orig.df <- data.frame(evi.sd = vasc.env$evi.sd, pr.3level.dom = as.numeric(pr.3level.dom))
orig.df$dom <- vasc.env$hetpast


## Simulate data with 'other' as dominant habitat

evi.seq$hetpast <- rep('N', 300) # 

preds.sd.domN<-predict(m.full.dom3level.exp$gam, type="terms", newdata=evi.seq,
                       se.fit=TRUE)

# set up evi.sd, the fit and the upper and lower
# confidence interval

evi<-evi.seq$evi.sd
fit<-preds.sd.domN$fit[,6]
fit.up95N<-fit-1.96*preds.sd.domN$se.fit[,6] 
fit.low95N<-fit+1.96*preds.sd.domN$se.fit[,6]

evi.dfN <- data.frame(evi = evi, fit = fit, fit.up95 = fit.up95N, fit.low95 = fit.low95N)


# get partial residuals
pred.origN <- predict(m.full.dom3level.exp$gam, type = 'terms')
pr.3level.domN <- residuals(m.full.dom3level.exp$gam, type = 'working') + pred.origN[,6]


orig.dfN <- data.frame(evi.sd = vasc.env$evi.sd, pr.3level.domN = as.numeric(pr.3level.domN))
orig.dfN$dom <- vasc.env$hetpast


## Simulate data with pasture as dominant habitat

evi.seq$hetpast <- rep('P', 300)

preds.sd.domP <- predict(m.full.dom3level.exp$gam, type = 'terms', newdata = evi.seq, se.fit = TRUE)

# set up to calcualte confidence intervals

evi <- evi.seq$evi.sd
fit <- preds.sd.domP$fit[,7]
fit.up95P <- fit - 1.96*preds.sd.domP$se.fit[,7]
fit.low95P <- fit + 1.96*preds.sd.domP$se.fit[,7]

evi.dfP <- data.frame(evi = evi, fit = fit, fit.up95 = fit.up95P, fit.low95 = fit.low95P)

# get partial residuals
pred.origP <- predict(m.full.dom3level.exp$gam, type = 'terms')
pr.3level.domP <- residuals(m.full.dom3level.exp$gam, type = 'working') + pred.origP[,7]

orig.dfP <- data.frame(evi.sd = vasc.env$evi.sd, pr.3level.domP = as.numeric(pr.3level.domP))
orig.dfP$dom <- vasc.env$hetpast




### Build graphic that combines both heterogenous and homogenous landscapes

# add to residual dfs whether habitat is heterogenous or not

orig.df$dom <- vasc.env$hetpast
orig.dfN$dom <- vasc.env$hetpast
orig.dfP$dom <- vasc.env$hetpast


# extract range of EVIs where each dominant habitat is found

ve.H <- vasc.env[which(vasc.env$hetpast=='H'),]
ve.H <- droplevels(ve.H)
range(ve.H$evi.sd) # 0.078 0.25

ve.N <- vasc.env[which(vasc.env$hetpast=='N'),]
ve.N <- droplevels(ve.N)
range(ve.N$evi.sd) # 0.068 0.16

ve.P <- vasc.env[which(vasc.env$hetpast=='P'),]
ve.P <- droplevels(ve.P)
range(ve.P$evi.sd) # 0.076 0.25


evi.df.Hless <- evi.df[evi.df$evi >= 0.078 & evi.df$evi <= 0.25,]
evi.df.Hless <- droplevels(evi.df.Hless)

evi.df.Nless <- evi.dfN[evi.dfN$evi >= 0.068 & evi.dfN$evi <= 0.16,]
evi.df.Nless <- droplevels(evi.df.Nless)

evi.df.Pless <- evi.dfP[evi.dfP$evi >= 0.076 & evi.dfP$evi <= 0.25,]
evi.df.Pless <- droplevels(evi.df.Pless)



p.nice.sd <- ggplot(evi.df.Hless, aes(x = evi, y = fit)) + geom_line(aes(x = evi, y = fit, colour = 'darkorchid4'), size = 1.5) +
  geom_ribbon(aes(ymin = fit.low95, ymax = fit.up95, x = evi, fill = 'band', alpha = 0.1), fill = 'darkorchid4') +
  xlab('Spatial SD of EVI') + ylab('Standardised Species Richness') +
  xlim(0.068, 0.25)


p.nice.sd <- p.nice.sd + geom_line(data = evi.df.Pless, aes(x = evi, y = fit, colour = 'forestgreen'), size = 1.5) +
  geom_ribbon(data = evi.df.Pless, aes(ymin = fit.low95, ymax = fit.up95, x = evi, fill = 'band', alpha = 0.1), fill = 'forestgreen')

p.nice.sd <- p.nice.sd + geom_line(data = evi.df.Nless, aes(x = evi, y = fit, colour = 'darkorange3'), size = 1.5) +
  geom_ribbon(data = evi.df.Nless, aes(ymin = fit.low95, ymax = fit.up95, x = evi, fill = 'band', alpha = 0.1), fill = 'darkorange3')

p.nice.sd <- p.nice.sd + geom_point(data = orig.df[orig.df$dom=='H',], mapping = aes(x = evi.sd, y = pr.3level.dom), colour = 'darkorchid4', size = 3) # plot ast so top layer


p.nice.sd <- p.nice.sd + scale_colour_manual(values = c('darkorange3'='darkorange3', 'darkorchid4'='darkorchid4', 'forestgreen' = 'forestgreen'),guide = guide_legend(title = NULL),
                                             labels = c('Other', 'Heterogeneous', 'Pasture')) + 
  scale_fill_manual(guide=  'none') + 
  scale_alpha(guide = 'none')

p.nice.sd <- p.nice.sd + geom_point(data = orig.dfN[orig.dfN$dom=='N',], mapping = aes(x = evi.sd, y = pr.3level.domN), colour = 'darkorange3', size = 3)
p.nice.sd <- p.nice.sd + geom_point(data = orig.dfP[orig.dfP$dom=='P',], mapping = aes(x = evi.sd, y = pr.3level.domP), colour = 'forestgreen', size = 3)

p.nice.sd <- p.nice.sd + theme(axis.title = element_text(size = 25), axis.text = element_text(size = 20), legend.position = 'none') 




##### Interannual sd

### simulate data with each land cover in turn to then combine on to same graphic

# Simulate data with heterogeneous habitat
# add dominant habitat column to evi.seq

# create a sequence of evi that spans evi.mean
maxevi<-max(vasc.env$inter.evi.sd)
minevi<-min(vasc.env$inter.evi.sd)
evi.seq<-seq(minevi, maxevi, length=300)
evi.seq<-data.frame(evi.mean = rep(0.5, 300), evi.sd = rep(0.11, 300), inter.evi.sd = evi.seq, annual.evi.sd = rep(0.097, 300), east = rep(205096.4, 300), north = rep(239771.1, 300))

evi.seq$hetpast <- rep('H', 300) # can be changed to N or P for other graphs


# predict only the inter.evi.sd term (the sum of the
# term predictions and the intercept gives you the overall 
# prediction)

preds.inter.dom<-predict(m.full.dom3level.exp$gam, type="terms", newdata=evi.seq,
                         se.fit=TRUE)

# set up evi.sd, the fit and the upper and lower
# confidence interval

evi<-evi.seq$inter.evi.sd
fit<-preds.inter.dom$fit[,8] 
fit.up95<-fit-1.96*preds.inter.dom$se.fit[,8] 
fit.low95<-fit+1.96*preds.inter.dom$se.fit[,8]

evi.df <- data.frame(evi = evi, fit = fit, fit.up95 = fit.up95, fit.low95 = fit.low95)


# get partial residuals
pred.orig <- predict(m.full.dom3level.exp$gam, type = 'terms')
pr.3level.dom <- residuals(m.full.dom3level.exp$gam, type = 'working') + pred.orig[,8]

orig.df <- data.frame(inter.evi.sd = vasc.env$inter.evi.sd, pr.3level.dom = as.numeric(pr.3level.dom))
orig.df$dom <- vasc.env$hetpast


## Simulate data with 'other' as dominant habitat

evi.seq$hetpast <- rep('N', 300) # 

preds.inter.domN<-predict(m.full.dom3level.exp$gam, type="terms", newdata=evi.seq,
                          se.fit=TRUE)

# set up evi.sd, the fit and the upper and lower
# confidence interval

evi<-evi.seq$inter.evi.sd
fit<-preds.inter.domN$fit[,9]
fit.up95N<-fit-1.96*preds.inter.domN$se.fit[,9] 
fit.low95N<-fit+1.96*preds.inter.domN$se.fit[,9]

evi.dfN <- data.frame(evi = evi, fit = fit, fit.up95 = fit.up95N, fit.low95 = fit.low95N)


# get partial residuals
pred.origN <- predict(m.full.dom3level.exp$gam, type = 'terms')
pr.3level.domN <- residuals(m.full.dom3level.exp$gam, type = 'working') + pred.origN[,9]

orig.dfN <- data.frame(inter.evi.sd = vasc.env$inter.evi.sd, pr.3level.domN = as.numeric(pr.3level.domN))
orig.dfN$dom <- vasc.env$hetpast


## Simulate data with pasture as dominant habitat

evi.seq$hetpast <- rep('P', 300)

preds.inter.domP <- predict(m.full.dom3level.exp$gam, type = 'terms', newdata = evi.seq, se.fit = TRUE)

# set up to calcualte confidence intervals

evi <- evi.seq$inter.evi.sd
fit <- preds.inter.domP$fit[,10]
fit.up95P <- fit - 1.96*preds.inter.domP$se.fit[,10]
fit.low95P <- fit + 1.96*preds.inter.domP$se.fit[,10]

evi.dfP <- data.frame(evi = evi, fit = fit, fit.up95 = fit.up95P, fit.low95 = fit.low95P)

# get partial residuals
pred.origP <- predict(m.full.dom3level.exp$gam, type = 'terms')
pr.3level.domP <- residuals(m.full.dom3level.exp$gam, type = 'working') + pred.origP[,10]

orig.dfP <- data.frame(inter.evi.sd = vasc.env$inter.evi.sd, pr.3level.domP = as.numeric(pr.3level.domP))
orig.dfP$dom <- vasc.env$hetpast




### Build graphic that combines both heterogenous and homogenous landscapes

# add to residual dfs whether habitat is heterogenous or not

orig.df$dom <- vasc.env$hetpast
orig.dfN$dom <- vasc.env$hetpast
orig.dfP$dom <- vasc.env$hetpast


# extract range of EVIs where each dominant habitat is found

ve.H <- vasc.env[which(vasc.env$hetpast=='H'),]
ve.H <- droplevels(ve.H)
range(ve.H$inter.evi.sd) # 0.0063 0.026

ve.N <- vasc.env[which(vasc.env$hetpast=='N'),]
ve.N <- droplevels(ve.N)
range(ve.N$inter.evi.sd) # 0.0070 0.029

ve.P <- vasc.env[which(vasc.env$hetpast=='P'),]
ve.P <- droplevels(ve.P)
range(ve.P$inter.evi.sd) # 0.0098 0.029


evi.df.Hless <- evi.df[evi.df$evi >= 0.0063 & evi.df$evi <= 0.026,]
evi.df.Hless <- droplevels(evi.df.Hless)

evi.df.Nless <- evi.dfN[evi.dfN$evi >= 0.0070 & evi.dfN$evi <= 0.029,]
evi.df.Nless <- droplevels(evi.df.Nless)

evi.df.Pless <- evi.dfP[evi.dfP$evi >= 0.0098 & evi.dfP$evi <= 0.029,]
evi.df.Pless <- droplevels(evi.df.Pless)



p.nice.inter <- ggplot(evi.df.Hless, aes(x = evi, y = fit)) + geom_line(aes(x = evi, y = fit, colour = 'darkorchid4'), size = 1.5) +
  geom_ribbon(aes(ymin = fit.low95, ymax = fit.up95, x = evi, fill = 'band', alpha = 0.1), fill = 'darkorchid4') +
  xlab('Interannual SD of EVI') + ylab('Standardised Species Richness') +
  xlim(0.0063, 0.029)


p.nice.inter <- p.nice.inter + geom_line(data = evi.df.Pless, aes(x = evi, y = fit, colour = 'forestgreen'), size = 1.5) +
  geom_ribbon(data = evi.df.Pless, aes(ymin = fit.low95, ymax = fit.up95, x = evi, fill = 'band', alpha = 0.1), fill = 'forestgreen')

p.nice.inter <- p.nice.inter + geom_line(data = evi.df.Nless, aes(x = evi, y = fit, colour = 'darkorange3'), size = 1.5) +
  geom_ribbon(data = evi.df.Nless, aes(ymin = fit.low95, ymax = fit.up95, x = evi, fill = 'band', alpha = 0.1), fill = 'darkorange3')

p.nice.inter <- p.nice.inter + geom_point(data = orig.df[orig.df$dom=='H',], mapping = aes(x = inter.evi.sd, y = pr.3level.dom), colour = 'darkorchid4', size=  3) # plot ast so top layer


p.nice.inter <- p.nice.inter + scale_colour_manual(values = c('darkorange3'='darkorange3', 'darkorchid4'='darkorchid4', 'forestgreen' = 'forestgreen'),guide = guide_legend(title = NULL),
                                                   labels = c('Other', 'Heterogeneous', 'Pasture')) + 
  scale_fill_manual(guide=  'none') + 
  scale_alpha(guide = 'none')

p.nice.inter <- p.nice.inter + geom_point(data = orig.dfN[orig.dfN$dom=='N',], mapping = aes(x = inter.evi.sd, y = pr.3level.domN), colour = 'darkorange3', size = 3)
p.nice.inter <- p.nice.inter + geom_point(data = orig.dfP[orig.dfP$dom=='P',], mapping = aes(x = inter.evi.sd, y = pr.3level.domP), colour = 'forestgreen', size = 3)

p.nice.inter <- p.nice.inter + theme(axis.title = element_text(size = 25), axis.text = element_text(size = 20), legend.position ='none') 






