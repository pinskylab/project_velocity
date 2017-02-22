library(mgcv)
library(data.table) 
library(ggplot2)
library(reshape2)
library(lattice)
setwd('/Users/jamesmorley/Documents/project_velocity')
    
load('dat_selectedspp_Feb_1_2017.Rdata')# load species catch data
dat <- data.table(dat)

# Some qaqc regarding occurences of NA or zeros for wtcpue
nrow(dat[wtcpue==0]) # 197,220 rows have zero wtcpue
abc <- dat[wtcpue==0] # Most of these in SoGulf, where they include all species with each tow, but also a fair number in both NEFSC and ScotianShelf in all seasons
# It seems NEFSC got a better scale or protocol in 2001, as few 0 weights thereafter_same might have happened on Scotian in 95
nrow(dat[is.na(wtcpue) & !sppocean == "NO CATCH"]) #246 rows have no catch weight
abc <- dat[is.na(wtcpue) & !sppocean == "NO CATCH"] # spread out among the surveys mostly_mostly all in unique hauls (ie no entire hauls are unweighed)
# These NAs will be dropped from biomass models (of course), but all but one species have <9 occurences, so not really an issue
# oncorhynchus tshawytscha_Pac (Chinook salmon) was not weighed on 100 occasions_and it was caught 1195 times_seems unlikely that it's weight was negligible_these all occurred after 2005 in NWFSC_WCAnn
# Maybe a tagging study or something..
abc <- dat[sppocean == 'oncorhynchus tshawytscha_Pac' & is.na(wtcpue)] # caught mostly off west coast and in GOA, but others in Aleuts and EBS
# All chinook caught in WC_ANN 2005-2014 were not weighed_estimate weight from #caught? Need to redo that from script 1.....
rm(abc)
dat <- dat[!(wtcpue == 0 & region == 'DFO_SoGulf')] # the zeros in SoGulf are actual zeros (ie not just so few to not register) and thus an absence
  
dat$survey <- dat$region
 
# ==================================================================================================================================================

# CHOOSE A TEST SPECIES_first one with zeros for wtcpue, then one with NAs_just to make sure everything else adjusts for those situations
mydat <- dat[sppocean == 'squalus acanthias_Atl'] 
mysurveys<-unique(mydat$survey)
#Add the empty hauls for the species of interest
haulsMod <- haulsTest[haulsTest$survey %in% mysurveys,]# restrict the hauls file to surveys of interest for this species
# Add a presence column to mydat and drop redundant columns for a merge with haulsMod
mydat <- data.frame(haulid = mydat$haulid, sppocean = mydat$sppocean, Freq = mydat$Freq, wtcpue = mydat$wtcpue, pres = TRUE, stringsAsFactors = F)
mydat$logwtcpue <- log(mydat$wtcpue + 1) # for starters we'll keep the zero cpue's in
# mydat$logwtcpue <- log(mydat$wtcpue)
# mydat$logwtcpue[is.infinite(mydat$logwtcpue)] <- NA # So these are essentially dropped from the biomass models_not necessarily a good thing, esp. for smaller spp.....
mydatMod <- merge(haulsMod, mydat, by='haulid', all.x = T, sort=F)
rm(mydat, haulsMod)

mydatMod$pres[is.na(mydatMod$pres)] <- FALSE
mydatMod$presNum <- as.numeric(mydatMod$pres)
# NEED TO GET THE 'ZEROS' INTO THE wtcpue and logwtcpue SO THAT THIS CAN BE USED TO GENERATE AVERAGE ANNUAL Density VALUES
mydatMod$remove.biom.mod <- ifelse(is.na(mydatMod$wtcpue) & mydatMod$pres == T, T, F) # an NA for wtcpue, but the species was present, this will flag it to remove from the biomass model and annual biomass estimates
mydatMod$wtcpue[is.na(mydatMod$wtcpue)] <- 0 # Now some of the zeros are presences, but presumably very light_while other zeros are absences
# IF WE END UP CALCULATING ANNUAL CPUE WITH LOGGED HAUL CPUE VALUES, I'LL NEED TO MODIFY THE ABOVE CODE_currently all the true absences have 'logwtcpue' of NA
 
# If trimming out bad areas drops a species from a survey (ie dropping inshore), that survey should be dropped too....
# I should keep years that have zero in the consistent sampling area, but that caught them in one of the oddball areas_Not sure...
haulsTrim <- haulsTrim[!haulsTrim$remove.biom.mod == TRUE,] # remove any NA for wtcpue
ave.catch.wt <- as.data.frame(tapply(haulsTrim$wtcpue, list(haulsTrim$year, haulsTrim$survey), mean, na.rm=T))
ave.catch.wt$year <- row.names(ave.catch.wt); row.names(ave.catch.wt) <- NULL
avecatchyrreg <- melt(ave.catch.wt, id.vars='year', variable.name='survey', value.name='biomassmean')
 
#Also save average biomassmean for each region (across all years) to use in later predictions.
avecatchyrreg <- avecatchyrreg[!is.na(avecatchyrreg$biomassmean),]
avemeanbiomass <- aggregate(avecatchyrreg$biomassmean, by=list(avecatchyrreg$survey), FUN=mean)
colnames(avemeanbiomass) <- c('survey', 'avecpue')

avecatchyrreg$survey <- as.character(avecatchyrreg$survey)
avecatchyrreg$year <- as.numeric(avecatchyrreg$year)
spdata <- merge(mydatMod, avecatchyrreg, all=T, by=c('year', 'survey'), sort=F)
rm(ave.catch.wt, avecatchyrreg, haulsMod, mydat)

save(dat, haulsTest, spdata, file='uncertainty_Feb15_2017.RData') 

# ===========================================================
# Model fitting =============================================
# ===========================================================

# The first four are just seeing how the soda data relate to observed data
plot(bottemp~SBT.actual, data=spdata); summary(lm(bottemp~SBT.actual, data=spdata)) # r2 .77
plot(bottemp~SBT.actual, data=spdata[spdata$pres == TRUE,]); summary(lm(bottemp~SBT.actual, data=spdata[spdata$pres == TRUE,])) # r2 .51
plot(surftemp~SST.actual, data=spdata); summary(lm(surftemp~SST.actual, data=spdata)) # r2 .91
plot(surftemp~SST.actual, data=spdata[spdata$pres == TRUE,]); summary(lm(surftemp~SST.actual, data=spdata[spdata$pres == TRUE,])) # r2 .82
#The remaining are relating potentially correlated predictors
plot(SBT.max~SST.max, data=spdata); summary(lm(SBT.max~SST.max, data=spdata)) # r2 .65_but pretty wonky looking_may want to include both
plot(SBT.max~SST.max, data=spdata[spdata$pres == TRUE,]); summary(lm(SBT.max~SST.max, data=spdata[spdata$pres == TRUE,])) # r2 .39_same residual pattern as above
plot(SBT.min~SST.min, data=spdata); summary(lm(SBT.min~SST.min, data=spdata)) # r2 .90
plot(SBT.min~SST.min, data=spdata[spdata$pres == TRUE,]); summary(lm(SBT.min~SST.min, data=spdata[spdata$pres == TRUE,])) # r2 .90
plot(bottemp~surftemp, data=spdata); summary(lm(bottemp~surftemp, data=spdata)) # r2 .84
plot(bottemp~surftemp, data=spdata[spdata$pres == TRUE,]); summary(lm(bottemp~surftemp, data=spdata[spdata$pres == TRUE,])) # r2 .90
plot(SST.seasonal.mean~surftemp, data=spdata); summary(lm(SST.seasonal.mean~surftemp, data=spdata)) # .92
plot(SST.seasonal.mean~surftemp, data=spdata[spdata$pres == TRUE,]); summary(lm(SST.seasonal.mean~surftemp, data=spdata[spdata$pres == TRUE,])) # .84
plot(SBT.seasonal~bottemp, data=spdata); summary(lm(SBT.seasonal~bottemp, data=spdata)) # .79
plot(SBT.seasonal~bottemp, data=spdata[spdata$pres == TRUE,]); summary(lm(SBT.seasonal~bottemp, data=spdata[spdata$pres == TRUE,])) # .53

# observed temps match up pretty well with SODA estimates (very well at the surface)
# May want to include both SBT.max and SST.max as a lot of variation at higher temps
# Both SBT.min and SST.min highly correlated, should stick with just SBT.min
# Bottom temp and surface temp highly correlated with this data set
# Seasonal means highly correlated with observed for SST, should drop seasonal mean of SST
# Seasonal mean for BTemp also correlated, especially for P/A model, probably can drop

# consider some interactions with temperature and things like habitat.....also try GRAVEL instead of GRAINSIZE
# Can you even interact a continuous with a factor? 
ptm <- Sys.time()
mod1 <- gam(presNum~s(SBT.seasonal) + s(SST.seasonal.mean) + s(rugosity) + s(GRAINSIZE) + s(SBT.min) + s(SBT.max) + s(SST.max) + s(biomassmean, by=regionfact) + regionfact + habitatFact - 1, family=binomial, data=spdata.fluke, select=TRUE, gamma=5.5)
Sys.time() - ptm
ptm <- Sys.time()
mod1b <- gam(presNum~s(SBT.seasonal, k=7) + s(SST.seasonal.mean, k=7) + s(rugosity, k=7) + s(GRAINSIZE, k=7) + s(SBT.min, k=7) + s(SBT.max, k=7) + s(SST.max, k=7) + s(biomassmean, by=regionfact) + regionfact + habitatFact - 1, family=binomial, data=spdata.fluke, select=TRUE, gamma=5.5)
Sys.time() - ptm

# above gamma is set to approximate log(nrow(spdata))/2   which is a method I saw on a SN Wood PPT.
# 'select=T' allows model terms to be penalized to have no effect (ie df is negligible)_might be a 'model selection' strategy that can work in our framework_also seems to reduce curviness, but not as much as the 'gamma' option can
summary(mod1) # a number of terms maxed out there df (knots)  

pdf(width=9, height=7, file='/Users/jamesmorley/Desktop/gam_PA.pdf')
plot(mod1b, scale=0, shade=TRUE, all.terms=T, pages=1)
dev.off()

#par(mfrow=c(1,2))
plot(mod1, all.terms=T, se=FALSE, select=13)
# plot(mod1, all.terms=T, se=FALSE, select=10)
par(mfrow=c(2,2))
gam.check(mod1)
par(mfrow=c(1,1))
concurvity(mod1, full=T) # concurvity of each term with the entire remaining model examined
# 0 indicating no problem, and 1 indicating total lack of identifiability
# WORST is a very 'pessimistic' measure; OBSERVED may be 'over-optimistic'; ESTIMATE sounds like it may be the best to use...
# below models one of the predictors with the rest of the variables, looks at how well a predictor is already modeled by other variables (in this case %deviance is 89%!)
summary(gam(surftemp~s(bottemp) + s(rugosity) + s(GRAINSIZE) + s(SBT.min) + s(SBT.max) + s(SST.max) + s(biomassmean) + regionfact + habitatFact, data=spdata))
concurvity(mod1, full=F)[3] # pairwise comparisons of terms_just prints 'ESTIMATE' value

mod2 <- gam(logwtcpue~s(SBT.seasonal) + s(SST.seasonal.mean) + s(rugosity) + s(GRAINSIZE) + s(SBT.min) + s(SBT.max) + s(SST.max) + s(biomassmean, by=regionfact) + regionfact + habitatFact - 1, data=spdata.dog[spdata.dog$pres == TRUE & spdata.dog$remove.biom.mod == FALSE,], select=TRUE, gamma=4.9)
mod2b <- gam(logwtcpue~s(SBT.seasonal, k=7) + s(SST.seasonal.mean, k=7) + s(rugosity, k=7) + s(GRAINSIZE, k=7) + s(SBT.min, k=7) + s(SBT.max, k=7) + s(SST.max, k=7) + s(biomassmean, by=regionfact) + regionfact + habitatFact - 1, data=spdata.dog[spdata.dog$pres == TRUE & spdata.dog$remove.biom.mod == FALSE,], select=TRUE, gamma=4.9)
mod2c <- gam(logwtcpue~s(SBT.seasonal, k=7) + s(SST.seasonal.mean, k=7) + s(rugosity, k=7) + s(GRAINSIZE, k=7) + s(SBT.min, k=7) + s(SBT.max, k=7) + s(SST.max, k=7) + s(biomassmean, by=regionfact) + regionfact + habitatFact - 1, data=spdata.dog[spdata.dog$pres == TRUE & spdata.dog$remove.biom.mod == FALSE,], gamma=4.9)
mod2d <- gam(logwtcpue~s(SBT.seasonal) + s(SST.seasonal.mean) + s(rugosity) + s(GRAINSIZE) + s(SBT.min) + s(SBT.max) + s(SST.max) + s(biomassmean, by=regionfact) + regionfact + habitatFact - 1, data=spdata.dog[spdata.dog$pres == TRUE & spdata.dog$remove.biom.mod == FALSE,], gamma=4.9)
mod2e <- gam(logwtcpue~s(SBT.seasonal, k=7) + s(SST.seasonal.mean, k=7) + s(rugosity, k=7) + s(GRAINSIZE, k=7) + s(SBT.min, k=7) + s(SBT.max, k=7) + s(SST.max, k=7) + s(biomassmean, by=regionfact) + regionfact + habitatFact - 1, data=spdata.dog[spdata.dog$pres == TRUE & spdata.dog$remove.biom.mod == FALSE,])
mod2f <- gam(logwtcpue~s(SBT.seasonal, k=6) + s(SST.seasonal.mean, k=6) + s(rugosity, k=6) + s(GRAINSIZE, k=6) + s(SBT.min, k=6) + s(SBT.max, k=6) + s(SST.max, k=6) + s(biomassmean, by=regionfact) + regionfact + habitatFact - 1, data=spdata.dog[spdata.dog$pres == TRUE & spdata.dog$remove.biom.mod == FALSE,])

mod2bb <- gam(logwtcpue~s(SBT.seasonal, k=7) + s(SST.seasonal.mean, k=7) + s(rugosity, k=7) + s(GRAINSIZE, k=7) + s(SBT.min, k=7) + s(SBT.max, k=7) + s(SST.max, k=7) + (biomassmean, by=regionfact) + regionfact + habitatFact - 1, data=spdata.dog[spdata.dog$pres == TRUE & spdata.dog$remove.biom.mod == FALSE,], select=TRUE, gamma=4.9)

# all models make slightly different predictions
predict(mod2, type='response')[10000:10010] # though I don't think 'response' is necessary as we aren't using a link function here
predict(mod2b, type='response')[10000:10010] 
predict(mod2c, type='response')[10000:10010] 
predict(mod2d, type='response')[10000:10010] 
 
specs <- mod1c$smooth[[1]]# Gives info on each predictor (i.e. [[1]] is SBT.seasonal in this case)_is itself a list of 23 things
str(specs) # str provides a better display of the list items
mod1c$coefficients # give effect values for everything, including each segment of each term = k-1_doesn't matter the effective df, each term has k-1 values
  # however, this doesn't pertain to the value on the plots (i.e. y-axis value)_it's a spline so I think it relates more to the slope in that section or something
  # though it seems when a segment has little effect the value is very close to zero (may be positive or negative)
model.matrix(mod1c)

# ESTIMATING MODEL UNCERTAINTY ====================================================================================================
   
# Specify a simple one variable model to learn how this all works
mod1 <- gam(logwtcpue~s(SBT.seasonal, k=7) + s(biomassmean, by=regionfact) + regionfact , data=spdata.dog[spdata.dog$pres == TRUE & spdata.dog$remove.biom.mod == FALSE,], select=TRUE, gamma=4.9)
mod2 <- gam(presNum~s(SBT.seasonal, k=7) + s(biomassmean, by=regionfact) + regionfact , family=binomial, data=spdata.dog, select=TRUE, gamma=5.6)
 
pd <- data.frame(SBT.seasonal=c(-1:25), biomassmean=rnorm(27, mean=60, sd=25), regionfact=rep(c('DFO_ScotianShelf', 'DFO_SoGulf', 'NEFSC_NEUS', 'SCDNR_SEUS', 'VIMS_NEAMAP'), length.out=27))# simple prediction matrix
pd <- expand.grid(pd)# expand this to get some better diagnostic graphs below
Xp.1 <- predict(mod1, pd, type='lpmatrix') # a prediction matrix, needed for ultimately getting variance estimates for quantities derived from the model
Xp.2 <- predict(mod2, pd, type='lpmatrix') 
# matrix that gets multiplied by model coefficients to yields predictions at new values in 'pd'_assigns value for each section of each curve, for each row of newdata to predict
# The lpmatrix breaks down each curve sections influence on the prediction.
beta.1 <- coef(mod1) # model coeffients for each curve section
beta.2 <- coef(mod2) # model coeffients for each curve section
Vb.1 <- vcov(mod1) # the variance of each coefficient (each curve section) in diagonal, and covariation b/w all factor levels and the sections of curves 
Vb.2 <- vcov(mod2) # the variance of each coefficient (each curve section) in diagonal, and covariation b/w all factor levels and the sections of curves 
n <- 100
library(MASS)
br.1 <- mvrnorm(n, beta.1, Vb.1) # simulate n rep coef vectors (for each section of curve) from posterior
# each column of br is a replicate parameter vector drawn from posterior of distribution of the parameter
br.2 <- mvrnorm(n, beta.2, Vb.2)
ilink <- family(mod2)$linkinv
a.range <- matrix(nrow=nrow(pd), ncol=n) # for filling with a loop
for(i in 1:n){
  # pred.a gives a response level prediction (I believe b/c we use no link function....)
  pred.a1 <- Xp.1%*%br.1[i,] # combine mean predictions from Xp with resamples of parameter estimates_ the %*% is a matrix algebra operator I believe
  pred.a2 <- Xp.2%*%br.2[i,]
  a.range[,i] <- pred.a1 * ilink(pred.a2) # mod2 uses the logit link function_here the 'ilink' function transforms it back to the scale of the response variable (i.e. a 0-1 value)
  # the ilink function is not needed for biomass model (mod1) as it does not use the link function
} 
 
# predict(mod2, pd, type='response') # similar to the values reported in a.range
   
uncer <- data.frame(cbind(a.range, pd))
uncer.melt <- melt(uncer, id.vars=c('SBT.seasonal', 'biomassmean', 'regionfact'), variable.name='iteration', value.name='prediction')   
plot(prediction~SBT.seasonal, data=uncer.melt)
plot(prediction~biomassmean, data=uncer.melt)
levelplot(prediction ~ SBT.seasonal * biomassmean, data=uncer.melt)
 