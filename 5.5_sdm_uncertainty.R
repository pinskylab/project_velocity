library(mgcv)
library(MASS)
library(data.table) 
library(ggplot2)
library(reshape2)
library(lattice)  
library(dismo)
require(Hmisc)
# library(gbm)
# library(grDevices)
setwd('/Users/jamesmorley/Documents/project_velocity')

load('data/master_hauls_March7_2017.RData') # import master hauls file, and trimmed hauls file for calculating annual mean cpue
load('data/dat_selectedspp_Feb_1_2017.Rdata')# load species catch data
load('data/ProjectionBathGrid_Feb27_2017.RData') 

dat <- dat[!(dat$wtcpue == 0 & dat$region == 'DFO_SoGulf'),] # the zeros in SoGulf are actual zeros (ie not just a scale issue) and thus are true absences
# NOTE: All chinook salmon ('oncorhynchus tshawytscha_Pac') caught in WC_ANN 2005-2014 (n = 100, out of 1195 total catches in Pacific) were not weighed_probably a tagging study or something....
# Weight could potentially be estimated from #caught? Would need to redo that from script 1.....maybe if a major revision is required (as it is an important species)

dat$wtcpue[dat$wtcpue == 0] <- 0.00001 # 'zeros' in dat are species too light to register on scales_here a value just below the lowest non-zero value is assigned_for transforming data
dat$logwtcpue <- log(dat$wtcpue)
# trim columns that are already in master hauls file, which will be merged in below with the hauls data
dat <- data.frame(haulid = dat$haulid, sppocean = dat$sppocean, Freq = dat$Freq, wtcpue = dat$wtcpue, logwtcpue = dat$logwtcpue, presfit = TRUE, stringsAsFactors = F)

# ==================================================================================================================================================
# BELOW IS A LOOP FOR CONDUCTION OF AN UNCERTAINTY ANALYSIS FOR ANY SPECIES
# ==================================================================================================================================================
allspp <- c('triglops forficatus_Pac', 'sebastes saxicola_Pac', 'paralithodes camtschaticus_Pac', 'microstomus pacificus_Pac', 'pandalus borealis_Atl',
            'menippe mercenaria_Atl', 'lutjanus campechanus_Atl', 'loligo pealeii_Atl', 'peprilus triacanthus_Atl', 'rhizoprionodon terraenovae_Atl')
oceans <- c('Pac', 'Pac', 'Pac', 'Pac', 'Atl', 'Atl', 'Atl', 'Atl', 'Atl', 'Atl')
options(warn=1) # print warnings as they occur

n = rep(NA, length(allspp))
modeldiag = data.frame(sppocean=n, npres=n, ntot=n, smear=n, pred_obsMedian=n, pred_obsMean=n, pred_obsMedianNosmear=n, pred_obsMeanNosmear=n, thresh=n, auc=n, tss=n, tssmax=n, acc=n, accmax=n, sens=n, spec=n, kappa=n, kappamax=n, rpb=n, r2.biomass=n, r2.all=n, r2.pres.1deg=n, r2.abun.1deg=n, dev.pres=n, dev.biomass=n,dev.pres.null=n, dev.biomass.null=n, stringsAsFactors=FALSE)

for(i in 1:length(allspp)){ 
  sp <- allspp[i]
  ocean <- oceans[i]  
  mydat <- dat[dat$sppocean == sp,] # trim dat file to the species of interest
  # Make master file for GAMs
  haulsMod <- hauls[hauls$ocean == ocean,] # trim master hauls file to the ocean of interest 
  haulsMod <- merge(haulsMod, mydat, by='haulid', all.x = T, sort=F)   # Add empty hauls from the relavant ocean
  haulsMod$presfit[is.na(haulsMod$presfit)] <- FALSE
  haulsMod <- droplevels(haulsMod) # mainly to drop the west coast 'regionfact' levels
  
  # Create a conditional for surveys that never catch 'sp', and conform that region to one that catches them_this prevents the factor of that region alone from explaining low abundance
  # This basically assumes that catchability is the same in both regions_based on gear_which may be reasonable if in the one region they aren't there anyway_it assumes that a species' absence isn't due to any sampling artifact, which is also reasonable in most cases
  haulsMod$logwtcpue.pad <- haulsMod$logwtcpue # will have some zeros changed to -18 to allow abundance model fitting
  haulsMod$presfit.pad <- haulsMod$presfit # will have some FALSE changed to TRUE to allow inclusion in abundance model fitting
  
  npres <- table(haulsMod$regionfact[haulsMod$presfit])
  if(any(npres < 1)){
    mywarn <- paste('Zero presences for', sp, 'in', paste(names(npres)[npres<1], collapse=', '))
    warning(mywarn)
    regstofill <- names(npres)[as.numeric(npres) == 0]
    haulsMod$regionfact[haulsMod$region %in% regstofill] <- names(npres)[which.max(npres)] # in regions with no observations, replace the region ID with that from a region with observations. this prevents a low region coefficient from explaining the zero observations.

  # pick some absences and force them to low abundance presences for abundance model fitting
  for(j in 1:length(regstofill)){
    theseinds <- haulsMod$region == regstofill[j]
    fake0s <- sample(which(theseinds), size = round(0.1 * sum(theseinds)))
    haulsMod$logwtcpue.pad[fake0s] <- -23
    haulsMod$presfit.pad[fake0s] <- TRUE
    print(paste(regstofill[j], ': Added', length(fake0s), 'fake zeros'))
    }
  }
     
  haulsMod <- droplevels(haulsMod) # may not be necessary as gam may not plot empty factor levels
  # set 'gamma' penalty levels for gam to prevent overfitting_got this from a presentation by Simon Wood (https://people.maths.bris.ac.uk/~sw15190/mgcv/tampere/mgcv-advanced.pdf)
  gammaPA <- log(nrow(haulsMod)) / 2
  gammaBiom <- log(nrow(haulsMod[haulsMod$presfit.pad == TRUE,])) / 2
   
  print(paste('Presence model fitting for', sp))
  mod1 <- gam(presfit~s(SBT.seasonal) + s(SST.seasonal.mean) + s(SBT.min) + s(SBT.max) + s(SST.max) + s(rugosity) + s(GRAINSIZE) + regionfact + habitatFact - 1, family=binomial, data=haulsMod, select=TRUE, gamma=gammaPA)
  print(paste('Abundance model fitting for', sp))
  mod2 <- gam(logwtcpue.pad~s(SBT.seasonal) + s(SST.seasonal.mean) + s(SBT.min) + s(SBT.max) + s(SST.max) + s(rugosity) + s(GRAINSIZE) + regionfact + habitatFact - 1, data=haulsMod[haulsMod$presfit.pad == TRUE,], select=TRUE, gamma=gammaBiom)
  mod1null <- gam(presfit~s(rugosity) + s(GRAINSIZE) + regionfact + habitatFact - 1, family=binomial, data=haulsMod, select=TRUE, gamma=gammaPA)
  mod2null <- gam(logwtcpue.pad~s(rugosity) + s(GRAINSIZE) + regionfact + habitatFact - 1, data=haulsMod[haulsMod$presfit.pad == TRUE,], select=TRUE, gamma=gammaBiom)
  # 'select=T' allows model terms to be penalized to have no effect (ie df is negligible)_might be a 'model selection' strategy that can work in our framework_also seems to reduce curviness, but not as much as the 'gamma' option can
     
  pdf(width=9, height=9, file=paste('/figures/uncertainty/gamFit_', sp, '.pdf', sep=''))
  plot(mod1, scale=0, shade=TRUE, all.terms=T, pages=1); mtext(paste(sp," presence_absence"),outer=T,line=-2)
  plot(mod2, scale=0, shade=TRUE, all.terms=T, pages=1); mtext(paste(sp," log(cpue)"),outer=T,line=-2)
  
  # Model diagnostics from script 6
  preds1 <- mod1$fitted.values # vector of predicted/modeled values based on observed data 
  preds2 <- exp(predict(mod2, newdata = haulsMod, type='response', na.action='na.pass')) # abundance predictions_needs different approach to fit to whole data set
  smear <- mean(exp(mod2$residuals)) # smearing estimator for re-transformation bias (see Duan 1983, http://www.herc.research.va.gov/include/page.asp?ID=cost-regression)
  preds <- preds1*preds2*smear # adds the bias correction as well_just a scaler that makes values more similar to observed (sometimes....)
  preds.nosmear <- preds1*preds2 # without the smear value
  haulsMod$wtcpue[haulsMod$presfit == FALSE] <- 0 
  haulsMod$predictions <- preds
  haulsMod$predictions.nosmear <- preds.nosmear
  par(mfrow=c(2,2))
  plot(wtcpue~predictions, ylab='observed cpue', main='Obs. vs pred. (w/ smear)', data=haulsMod)
  mtext(paste('Median diff=', summary(haulsMod$predictions - haulsMod$wtcpue)[3], sep=''), side=3, line=-1)
  mtext(paste('Mean diff=', summary(haulsMod$predictions - haulsMod$wtcpue)[4], sep=''), side=3, line=-2.5)
  plot(wtcpue~predictions.nosmear, ylab='observed cpue', main='Obs. vs pred. (no smear)', data=haulsMod)
  mtext(paste('Median diff=', summary(haulsMod$predictions.nosmear - haulsMod$wtcpue)[3], sep=''), side=3, line=-1)
  mtext(paste('Mean diff=', summary(haulsMod$predictions.nosmear - haulsMod$wtcpue)[4], sep=''), side=3, line=-2.5)
  
  # fill in diagnostics
  modeldiag$sppocean[i] = sp
  modeldiag$npres[i] = sum(haulsMod$presfit)
  modeldiag$ntot[i] = dim(haulsMod)[1]
  modeldiag$smear[i] = smear
  
  modeldiag$pred_obsMedian[i] = summary(haulsMod$predictions - haulsMod$wtcpue)[3]
  modeldiag$pred_obsMean[i] = summary(haulsMod$predictions - haulsMod$wtcpue)[4]
  modeldiag$pred_obsMedianNosmear[i] = summary(haulsMod$predictions.nosmear - haulsMod$wtcpue)[3]
  modeldiag$pred_obsMeanNosmear[i] = summary(haulsMod$predictions.nosmear - haulsMod$wtcpue)[4]
    
  # evaluate model (dismo package)
  # pick a threshold for pres/abs model evaluation (where needed)
  e <- evaluate(p=as.vector(preds1[haulsMod$presfit]), a=as.vector(preds1[!haulsMod$presfit])) # computes some model diagnostics based on predicted values at obersved presences and absences
  # the 'evaluate' function creates a confusion matrix, which is basically like a contingency table
  modeldiag$thresh[i] <- threshold(e, stat='prevalence') # Chooses a cut off value to consider something present or absent, based on 'e'
  e.ind <- which(e@t == modeldiag$thresh[i]) # index for the chosen threshold_I don't fully understand what e@t is....some hidden feature within e that isn't shown when you enter 'e'
  conf <- as.data.frame(e@confusion) # confusion matrices (all thresholds)
    
  # pres/abs model diagnostics (no threshold needed)
  modeldiag$dev.pres[i] = summary(mod1)$dev.expl
  modeldiag$auc[i] <- e@auc
  modeldiag$tssmax[i] <- max(with(conf, (tp*tn - fn*fp)/((tp+fp)*(fn+tn))), na.rm=TRUE) # maximum TSS (any threshold)
  modeldiag$accmax[i] <- max(with(conf, (tp+tn)/(tp+fp+fn+tn)), na.rm=TRUE) # maximum overall accuracy
  modeldiag$kappamax[i] <- max(e@kappa, na.rm=TRUE) # maximum kappa
  modeldiag$rpb[i] <- cor(preds1, haulsMod$presfit) # point biserial correlation
  
  # abundance model diagnostics
  modeldiag$dev.biomass[i] = summary(mod2)$dev.expl
  # Make a temp dataframe to do this correlation, because the 'cor' function doesn't work with NAs present
  abc <- data.frame(cbind(prediction = as.vector(log(preds2[haulsMod$presfit])), observed = as.vector(haulsMod$logwtcpue[haulsMod$presfit])))
  abc <- abc[!is.na(abc$observed),]
  modeldiag$r2.biomass[i] = cor(abc$prediction, abc$observed)^2 # correlation of log(biomass) when present
  
  # full model diagnostics
  abc <- data.frame(cbind(prediction = as.vector(preds), observed = as.vector(haulsMod$wtcpue)))
  abc <- abc[!is.na(abc$observed),]
  modeldiag$r2.all[i] = cor(abc$prediction, abc$observed)^2  # overall biomass correlation
  rm(abc)
  
  #Compare to models without temperature to ultimately calculation %explained by temp terms
  modeldiag$dev.pres.null[i] = summary(mod1null)$dev.expl
  modeldiag$dev.biomass.null[i] = summary(mod2null)$dev.expl
  
  # true skill statistic, accuracy, kappa, and other stats that require a threshold
  modeldiag$tss[i] = 	with(conf[e.ind,], (tp*tn - fn*fp)/((tp+fp)*(fn+tn))) # TSS for chosen threshold
  modeldiag$acc[i] = 	with(conf[e.ind,], (tp+tn)/(tp+fp+fn+tn)) # overall accuracy
  modeldiag$sens[i] = with(conf[e.ind,], (tp)/(tp+fn)) # sensitivity: fraction of correctly predicted presences
  modeldiag$spec[i] = with(conf[e.ind,], (tn)/(tn+fp)) # specificity: fraction of correctly predicted absences
  modeldiag$kappa[i] = e@kappa[e.ind] # Cohen's kappa
    
  test<-cbind(haulsMod,preds1,preds)
  # originally Malin had 'cs1' for 'surveyfact'...not sure what that is about
  t1<-tapply(test$preds1,list(test$year,test$surveyfact),mean) #average predicted p(occur)
  t2<-tapply(test$presfit,list(test$year,test$surveyfact),mean) #proportion of hauls with presence
  t3<-tapply(test$preds,list(test$year,test$surveyfact),mean) #average predicted abundance
  t4<-tapply(test$wtcpue,list(test$year,test$surveyfact),mean) #average observed abundance
  
  presr2<-round(cor(stack(as.data.frame(t2))[,1],stack(as.data.frame(t1))[,1],use="p")^2,2)
  abunr2<-round(cor(stack(as.data.frame(t4))[,1],(stack(as.data.frame(t3))[,1]),use="p")^2,2)
  par(mfrow=c(1,2))
  plot(stack(as.data.frame(t2))[,1],stack(as.data.frame(t1))[,1],xlab="Proportion of hauls with species present (by year and survey)",ylab="Mean predicted probability of occurrence", cex=0.5,main=sp)
  mtext(paste("r^2 =",presr2))
  plot(stack(as.data.frame(t4))[,1],(stack(as.data.frame(t3))[,1]),xlab="Average log(WTCPUE) (by year and survey)",ylab="Average predicted log(WTCPUE)", cex=0.5,main=sp)
  mtext(paste("r^2 =",abunr2))
  par(mfrow=c(1,1))
  # Need to change the name, b/c this isn't looking at it spatially anymore, or I need to add code to do it spatially
  modeldiag$r2.pres.1deg[i]<-presr2
  modeldiag$r2.abun.1deg[i]<-abunr2

  # ESTIMATING MODEL UNCERTAINTY ====================================================================================================
  pd <- proj.grid
  pd$rugosity <- log(pd$rugosity + 1)
  # Below some arbitrary space fillers while I'm waiting for climate files to be completed
  pd$SBT.seasonal <- rnorm(nrow(pd), mean=mean(haulsMod$SBT.seasonal), sd=2)
  pd$SST.seasonal.mean <- rnorm(nrow(pd), mean=mean(haulsMod$SST.seasonal.mean), sd=3)
  pd$SBT.min <- rnorm(nrow(pd), mean=mean(haulsMod$SBT.min), sd=1)
  pd$SBT.max <- rnorm(nrow(pd), mean=mean(haulsMod$SBT.max), sd=2)
  pd$SST.max <- rnorm(nrow(pd), mean=mean(haulsMod$SST.max), sd=4)

  pd$regionfact <- names(npres)[which.max(npres)] 
  # the region of max presences isn't necessarily the region with highest cpue, due to different survey durations
  # biomass <- avemeanbiomass[grep(names(npres)[which.max(npres)], avemeanbiomass$survey),] # extract the surveys within the region of max presences
  # logmeanbiom <- log(1 + mean(biomass$avecpue))
  # pd$logmeanbiom <- logmeanbiom # so use the mean biomass from surveys within the region used for regionfact
  
  # Trim prediction file to east or west coast
  if(mean(haulsMod$longrid) > -100) {
    pd <- pd[pd$lonBathgrid > -100]
  }
  if(mean(haulsMod$longrid) < -100) {
    pd <- pd[pd$lonBathgrid < -100]
  }

  Xp.1 <- predict(mod1, pd, type='lpmatrix') # a prediction matrix, needed for ultimately getting variance estimates for quantities derived from the model
  Xp.2 <- predict(mod2, pd, type='lpmatrix') 
  # matrix that gets multiplied by model coefficients to yields predictions at new values in 'pd'_assigns value for each section of each curve, for each row of newdata to predict
  # The lpmatrix breaks down each curve sections influence on the prediction.
  beta.1 <- coef(mod1) # model coeffients for each curve section
  beta.2 <- coef(mod2) # model coeffients for each curve section

  Vb.1 <- vcov(mod1) # the variance of each coefficient (each curve section) in diagonal, and covariation b/w all factor levels and the sections of curves 
  Vb.2 <- vcov(mod2) # the variance of each coefficient (each curve section) in diagonal, and covariation b/w all factor levels and the sections of curves 
  n <- 100
   
  br.1 <- mvrnorm(n, beta.1, Vb.1) # simulate n rep coef vectors (for each section of curve) from posterior
  # each column of br is a replicate parameter vector drawn from posterior of distribution of the parameter
  br.2 <- mvrnorm(n, beta.2, Vb.2)
  ilink <- family(mod1)$linkinv
  a.range <- matrix(nrow=nrow(pd), ncol=n) # for filling with a loop
  for(j in 1:n){
    pred.a1 <- Xp.1%*%br.1[j,] # combine mean predictions from Xp with resamples of parameter estimates_ the %*% is a matrix algebra operator I believe
    pred.a2 <- Xp.2%*%br.2[j,] # may need to add the smear estimate here
    a.range[,j] <- ilink(pred.a1) * exp(pred.a2) * smear# mod1 uses the logit link function_here the 'ilink' function transforms it back to the scale of the response variable (i.e. a 0-1 value)
    # the ilink function is not needed for biomass model (mod2) as it does not use the link function (gaussian distribution)
  } 
  
  uncer <- apply(a.range, 1, var)
  uncer.biom <- apply(a.range, 1, mean)
  uncer.grid <- data.frame(cbind(variance = uncer, mean = uncer.biom, latBathgrid = pd$latBathgrid, lonBathgrid = pd$lonBathgrid))
  
  # Make a threshold abundance for decalring suitable habitat?
  hab.thresh <- max(uncer.grid$mean) * 0.1 # 10% of maximum predicted catch
  uncer.grid$suitable <- uncer.grid$mean > hab.thresh
  # nrow(uncer.grid[uncer.grid$mean >= hab.thresh,]) / nrow(uncer.grid) # percentage of grid cells that are suitable habitat
  
  map1 <- levelplot(variance ~ lonBathgrid * latBathgrid, data=uncer.grid)
  map2 <- levelplot(mean ~ lonBathgrid * latBathgrid, data=uncer.grid)
  # map3 <- levelplot(suitable ~ lonBathgrid * latBathgrid, data=uncer.grid)
  print(map1); print(map2)#; print(map3)
  
  # a.rangeRaster <- data.frame(cbind(lonBathgrid = pd$lonBathgrid, latBathgrid = pd$latBathgrid, a.range))
  # habitat <- rasterFromXYZ(a.rangeRaster, res=c(0.05,0.05), digits=9)
  # plot(habitat[[50]]) # plots one layer of the rasterbrick 
  
  # compute centroids for each model iteration
  a.range <- data.frame(cbind(a.range, lonBathgrid = pd$lonBathgrid, latBathgrid = pd$latBathgrid))
  centroid <- data.frame(model.iter=numeric(), latcentroid=numeric(), loncentroid=numeric()) # this will become the final output file
  for(j in 1:n){
  lat.centroid <- wtd.mean(x = a.range$latBathgrid, weights = a.range[,j], na.rm = TRUE)
  lon.centroid <- wtd.mean(x = a.range$lonBathgrid, weights = a.range[,j], na.rm = TRUE)
  centroid[j,] <- data.frame(model.iter=j, latcentroid=lat.centroid, loncentroid=lon.centroid)
  }
  plot(latBathgrid~lonBathgrid, cex=.1, col='gray', data=uncer.grid[uncer.grid$suitable == F,])
  points(latBathgrid~lonBathgrid, cex=.1, col='red', data=uncer.grid[uncer.grid$suitable == T,])
  points(mean(centroid$latcentroid)~mean(centroid$loncentroid), col='black', cex=2, pch=8)
   
  # GO THROUGH LOOP ONE MORE TIME CAREFULLY AND MAKE SURE ITS WORKING RIGHT!!!!_LABEL (and rescale levelplots) GRAPHS etc.
  #r2.biomass produced a couple NAs last run through
  # Rescale levelplots so a light gray is the lowest value
  # Add some way to see how realistic the predicted values are with observed values_this is done graphically but.......
  # How to get better predicted values......some sort of rescaling based on the distribution of observed vs. predicted? 
  # Make internal loop for going through each prediction grid (n=39)
  # For each run I really only need to save the centroid value and sum wtcpue predictions; plus this must be done at the regional level too for U.S. locations
  # Create a third plot showing historical average cpue (averaged over seasons) within the bathgrid cells? Just for a rough comparison....
  # Do a loop with more species and then do some plots looking at sample size effects on dev.explained, smear, etc.
  # Read up on the dismo diagnostics above
  dev.off()
}

# SCRIPT END


# =========================================================================================================
# below is some code on estimating mean annual cpue (which we no longer use as a predictor),
# also BRT and gam output that I saved in case it is useful down the road
# =========================================================================================================

# Calculate annual average cpue within consistent survey sampling areas
haulsAbun <- haulsTrim[haulsTrim$ocean == ocean,]
haulsAbun <- merge(haulsAbun, mydat, by='haulid', all.x = T, sort=F) 
haulsAbun$presfit[is.na(haulsAbun$presfit)] <- FALSE
haulsAbun <- haulsAbun[!(is.na(haulsAbun$wtcpue) & haulsAbun$presfit == TRUE),] # remove rows with NA for wtcpue_for most species this is < 5 rows_this may not be necessary for tapply to work below
haulsAbun$wtcpue[is.na(haulsAbun$wtcpue)] <- 0 # convert NAs to zeros to calculate average annual cpue by survey
ave.catch.wt <- as.data.frame(tapply(haulsAbun$wtcpue, list(haulsAbun$year, haulsAbun$survey), mean, na.rm=T))
ave.catch.wt$year <- as.numeric(row.names(ave.catch.wt)); row.names(ave.catch.wt) <- NULL
avecatchyrreg <- melt(ave.catch.wt, id.vars='year', variable.name='survey', value.name='biomassmean')
avecatchyrreg <- avecatchyrreg[!is.na(avecatchyrreg$biomassmean),]
avecatchyrreg$survey <- as.character(avecatchyrreg$survey) # the melt (or tapply) messes up column classes
# merge with master file, adds mean biomass predictor
haulsMod <- merge(haulsMod, avecatchyrreg, all=T, by=c('year', 'survey'), sort=F)
haulsMod$logmeanbiom <- log(haulsMod$biomassmean + 1)

#Save average biomassmean for each region (across all years) to use in later predictions.
avemeanbiomass <- aggregate(avecatchyrreg$biomassmean, by=list(avecatchyrreg$survey), FUN=mean)
colnames(avemeanbiomass) <- c('survey', 'avecpue')
rm(ave.catch.wt, avecatchyrreg, haulsAbun, mydat)

#===================================================================================
#Boosted Regression Tree models using cross validation 
# I used Elith et al. 2008_also the supplemental with that has step by steps
#===================================================================================
  BRTB1 <- gbm.step(data=haulsMod,
                              gbm.x = c('SBT.seasonal', 'SST.seasonal.mean', 'SBT.min', 'SBT.max', 'SST.max', 'rugosity', 'GRAINSIZE', 'logmeanbiom', 'regionfact', 'habitatFact'), # 27:28, 
                              gbm.y = 'presNum',
                              family = "bernoulli",
                              tree.complexity = 5,
                              learning.rate = 0.01,
                              bag.fraction = 0.5)
  # loss function should be binomial deviance; is more robust to false negatives, which are very common_not sure if this is default for 'bernoulli'
  # the graph produced shows changes in predictive deviance with the addition of new trees_so how does predictability change, if it hits a minimum then rises, that rise is showing the model is being overfit at too many trees
   
  BRT2 <- gbm.step(data=haulsMod[haulsMod$presNum.pad == 1,],
                   gbm.x = c('SBT.seasonal', 'SST.seasonal.mean', 'SBT.min', 'SBT.max', 'SST.max', 'rugosity', 'GRAINSIZE', 'logmeanbiom', 'regionfact', 'habitatFact'), # 27:28, 
                   gbm.y = 'logwtcpue.pad',
                   family = "gaussian",
                   tree.complexity = 5,
                   learning.rate = 0.01,
                   bag.fraction = 0.5)
    
  BRT2.int <- gbm.interactions(BRT2) # returns a list
  BRT2.int$rank.list # gives the five most important interactions
  BRT2.int$interactions # gives all the interactions scores
  gbm.perspec(BRT2,x=3, y=6, z.range=c(0,3))# Plot interactions_the numbers indicate variable number from 'rank.list'
  
  names(BRT1) # objects stored in the model_e.g. on next line
  BRT2$contributions # sums to 100
  head(BRT1$fitted) #predicted values for observed data
  head(BRT1$fitted.vars) # variance of fitted values on response scale
  gbm.plot(BRTB1) # accounts for the average effects of all other variables, for each variable plot_I think effects average to zero for each subplot
  gbm.plot.fits(BRTB1) # plots fitted values in relation to each of the predictors
   
# BRT notes for future work
# More useful code is in the tutorial from Elith (2008)_especially on preditions and calculating ideal tree size (which may be less useful for lots of species)...
    # in active papers folder_an appendix is also there that I should go through if we use these
# For Malin's paper he used learning rate of 0.001, tree complexity of 10, and a bag fraction of 75%. 
# It seems I may have to increase learning rate_as it takes many trees to reach a minimum of predictive error_learning rate relates to how much influence each tree has; or.....
   # increase treecomplexity would be better given our large sample sizes (try this first)
# Maybe increase the 'bag.fraction' in model a bit, includes more of the data for each tree (~0.6 to 0.7?)
# Need to use bootstrap method to calculate CI bands around BRT plots_unless new methods have arisen
    # Malin's paper used 500 bootstrap replicates
# Is the smear effect needed with BRT_we could try not log transforming the response variable....
    # need to compare predicted verse observed and exp(predicted) for logged response variable
# Can we just do a BRT with all the zeros included for cpue?.....probably not, didn't try
# USE A normal regression tree OF PREDICTOR VARIABLE(S) (in replace of the usual respone) to see how the habitat variables group together_as in De'ath (2000)
# Do we want to use the model simplification step for BRT_probably not as no gam analogue_supposed to be slow_function = gbm.simplify

  
# ====================================================================
# Eraseble code below involving gam output===================
# ====================================================================

gam.check(mod1)
concurvity(mod1, full=T) # concurvity of each term with the entire remaining model examined
# 0 indicating no problem, and 1 indicating total lack of identifiability
# WORST is a very 'pessimistic' measure; OBSERVED may be 'over-optimistic'; ESTIMATE sounds like it may be the best to use...
# below models one of the predictors with the rest of the variables, looks at how well a predictor is already modeled by other variables (in this case %deviance is 89%!)
concurvity(mod1, full=F)[3] # pairwise comparisons of terms_just prints 'ESTIMATE' value

specs <- mod1c$smooth[[1]]# Gives info on each predictor (i.e. [[1]] is SBT.seasonal in this case)_is itself a list of 23 things
str(specs) # str provides a better display of the list items
mod1$coefficients # give effect values for everything, including each segment of each term = k-1_doesn't matter the effective df, each term has k-1 values
  # however, this doesn't pertain to the value on the plots (i.e. y-axis value)_it's a spline so I think it relates more to the slope in that section or something
  # though it seems when a segment has little effect the value is very close to zero (may be positive or negative)
model.matrix(mod1c)

# zero out the region and biomass effects_I think it's important to do both_this is tedious for biomass due to nine curve sections_but this is the same for any species, regardless of estimated degrees of freedom
regionfact <- paste('regionfact', names(npres)[which.max(npres)], sep='')
beta.2[names(beta.2)==regionfact] <- 0
regionfactbiom <- paste('s(logmeanbiom):', regionfact, '.1', sep='')
beta.2[names(beta.2)==regionfactbiom] <- 0

rank.preds <- data.frame(index = 1:nrow(uncer.grid), pred=sort(uncer.grid$mean))
plot(pred~index, cex=.1, rank.preds)
perc.change <- vector(length=nrow(rank.preds))
for(j in 2:nrow(rank.preds)){
  val1 <- rank.preds$pred[j-1]
  val2 <- rank.preds$pred[j]
  perc.change[j] <- val2 - val1
}
rank.preds$change <- perc.change
rank.preds$percent.change <- (rank.preds$change/rank.preds$pred) * 100
plot(percent.change~index, cex=.1, rank.preds)  
