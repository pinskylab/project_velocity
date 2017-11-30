# First part of this script builds the GAM models, the second part extracts the median temperature preferences for each species.
# Third part of script gives predicted biomass at each temperature for all Atlantic species

setwd('~/Documents/range_projections/')
#setwd('/Users/jamesmorley/Documents/project_velocity')
modfolder <- 'CEmodels_tempPrefs/'
runname <- "fitallreg_2017" # use all regions in each fit that are from the same ocean

library(mgcv)
library(Hmisc)

load('data/master_hauls_March7_2017.RData') # import master hauls file
load('data/dat_selectedspp_Feb_1_2017.Rdata')# load species catch data
dat <- dat[!(dat$wtcpue == 0 & dat$region == 'DFO_SoGulf'),] # the zeros in SoGulf are actual zeros (ie not just a scale issue) and thus are true absences
dat$wtcpue[dat$wtcpue == 0] <- 0.00002 # 'zeros' in dat are now species too light to register on scales_here a value below the lowest non-zero value is assigned_for transforming data
dat$logwtcpue <- log(dat$wtcpue)

# trim out recent years for NEUS_Patrick using these as testing data
hauls <- hauls[!(hauls$region == 'NEFSC_NEUS' & hauls$year > 1989),]

# trim columns that are already in master hauls file, which will be merged in below with the hauls data
dat <- data.frame(haulid = dat$haulid, sppocean = dat$sppocean, Freq = dat$Freq, wtcpue = dat$wtcpue, logwtcpue = dat$logwtcpue, presfit = TRUE, stringsAsFactors = F)
dat <- dat[!is.na(dat$wtcpue),] # drop NAs as it creates errors for some species. May want to go back and manually do 'oncorhynchus tshawytscha_Pac' as almost 10% of presence records have NA for wtcpue (tagging study?)

# Create table to store model diagnostics
allspp = sort(unique(dat$sppocean))
allspp <- allspp[-174] # This species didn't fit as all records in NEUS in more recent years
n = rep(NA, length(allspp))
modeldiag = data.frame(sppocean=n, npres=n, fakepres=n, dev.pres=n, dev.biomass=n, stringsAsFactors=FALSE) # pred_obsMedianNosmear=n, pred_obsMeanNosmear=n, # tt is for training/testing model
 
pdf(file=paste("figures/CEmodelGAMsmooths/GAMs_TempPrefs",runname,".pdf",sep=""),width=10,height=10)

options(warn=1) # print warnings as they occur
allwarnings = NULL
print(paste(length(allspp), 'models to fit'))
 
for(i in 1:length(allspp)){
  fittrain = TRUE
  mygam1 <- mygam2 <- NULL 
   
  sp<-allspp[i]
  print(paste(i,sp, Sys.time()))
  mydat<-dat[dat$sppocean==sp,] 
  
  spdata <- merge(hauls, mydat, by='haulid', all.x = T, sort=F) # Add empty hauls
  myocean <- head(spdata$ocean[spdata$presfit == TRUE])[1] # identify if this is a Pacific or Atlantic species
  spdata <- spdata[spdata$ocean == myocean,] # trim master hauls file to the ocean of interest 
  spdata$presfit[is.na(spdata$presfit)] <- FALSE
  spdata <- droplevels(spdata) # drop the west or east coast 'regionfact' levels
  spdata <- spdata[!is.na(spdata$bottemp),]
  myregions<-unique(spdata$region) # Not sure if this needs to be where species present or whole ocean_or it may be able to be dropped entirely
  
  ##############################################
  # Set up data in regions with no presences
  ##############################################
  spdata$logwtcpue.pad <- spdata$logwtcpue # has some zeros changed to -18 to allow abundance model fitting
  spdata$presfit.pad <- spdata$presfit # has some FALSE changed to TRUE to allow abundance model fitting
  
  npres <- table(spdata$regionfact[spdata$presfit])
  if(any(npres < 1)){
    mywarn <- paste('Zero presences for', i, sp, 'in', paste(names(npres)[npres<1], collapse=', '))
    allwarnings <- c(allwarnings, mywarn)
    warning(mywarn)
    regstofill <- names(npres)[as.numeric(npres) == 0]
    spdata$regionfact[spdata$region %in% regstofill] <- names(npres)[which.max(npres)] # in regions with no observations, replace the region ID with that from a region with observations. this prevents a low region coefficient from explaining the zero observations.
    
    # if a region has no presences
    # pick some absences and force them to low abundance presences for abundance model fitting
    for(j in 1:length(regstofill)){
      theseinds <- spdata$region == regstofill[j]
      fake0s <- sample(which(theseinds), size = round(0.1 * min(sum(theseinds), sum(spdata$presfit))))# Uses either 10% of species total samples, or 10% of hauls in a region, whichever is smaller
      spdata$logwtcpue.pad[fake0s] <- -23 
      spdata$presfit.pad[fake0s] <- TRUE
      print(paste(regstofill[j], ': Added', length(fake0s), 'fake zeros'))
    }
  }
  spdata <- droplevels(spdata)
  
  ####################################################
  # Figure out which model formula given data
  ####################################################
  
  #Default models. Leave out region factor if necessary
  # since fitallreg, using all regions in an ocean
  if(length(levels(spdata$regionfact))==1){
    mypresmod<-formula(presfit ~ s(bottemp) + s(rugosity) + s(GRAINSIZE))
    myabunmod<-formula(logwtcpue.pad ~ s(bottemp) + s(rugosity) + s(GRAINSIZE))
  } else {
    mypresmod<-formula(presfit ~ s(bottemp) + s(rugosity) + s(GRAINSIZE) + regionfact)
    myabunmod<-formula(logwtcpue.pad ~ s(bottemp) + s(rugosity) + s(GRAINSIZE) + regionfact)
  }
  
  ####################################################
  #Fit models to All data (no test/training split)
  ####################################################
  
  gammaPA <- log(nrow(spdata)) / 2
  gammaAbun <- log(nrow(spdata[spdata$presfit.pad,])) / 2
  
  try2 <- tryCatch({
    mygam1<-gam(mypresmod,family="binomial",data=spdata, select=TRUE, gamma=gammaPA)
    mygam2<-gam(myabunmod,data=spdata[spdata$presfit.pad,], na.action='na.exclude', select=TRUE, gamma=gammaAbun) # only fit where spp is present

  }, error = function(e) {
    mywarn <- paste('Error in gam fitting for', i, sp, ':', e)
    assign('allwarnings', c(get('allwarnings', envir=.GlobalEnv), mywarn), envir=.GlobalEnv) # these assigns outside the local scope are poor form in R. But not sure how else to do it here...
    warning(mywarn)
  })
  
  ####################################################
  # Plot gam smooths to check for unrealistic out-of-range responses 
  ####################################################
  
  freq <- sum(spdata$presfit) # Number of true presences
  plot(mygam1,pages=1,scale=0,all.terms=TRUE, shade=T);mtext(paste(sp,"presence", "  ", freq, "true presences"),outer=T,line=-2)
  plot(mygam2,pages=1,scale=0,all.terms=TRUE, shade=T);mtext(paste(sp,"abundance", "  ", freq, "true presences"),outer=T,line=-2)
  
  # fill in diagnostics
  modeldiag$sppocean[i] = sp
  modeldiag$npres[i] = freq
  modeldiag$fakepres[i] = nrow(spdata[spdata$presfit.pad,]) - nrow(spdata[spdata$presfit,])

  # pres/abs model diagnostics (no threshold needed)
  modeldiag$dev.pres[i] = summary(mygam1)$dev.expl

  # abundance model diagnostics
  modeldiag$dev.biomass[i] = summary(mygam2)$dev.expl

  ####################################################
  #### Save models for later projections
  ####################################################
  
  mods = list(mygam1 = mygam1, mygam2 = mygam2)
  
  sp <- gsub('/', '', sp) # would mess up saving the file
  
  save(mods, myregions, file=paste(modfolder, 'TempPrefs_',runname, '_', sp, '.RData', sep='')) 
  
  # write these files each time through the loop so that we can watch progress
  save(modeldiag, file=paste("output/modeldiag_TempPrefs_",runname,".Rdata",sep=""))
  write.csv(modeldiag, file=paste("output/modeldiag_TempPrefs_",runname,".csv",sep=""))
  write.csv(allwarnings, file=paste('output/warnings_TempPrefs_', runname, '.csv', sep=''))
  
}
dev.off()



# ==================================================================================
# Develop temperature prefs for each species_for Patrick ==========================
# I basically run this as a separate script (i.e., start from scratch, except for the 'dat' file...I think) 
# ==================================================================================

load('output/modeldiag_TempPrefs_fitallreg_2017.Rdata')
allspp <- modeldiag$sppocean
load('data/master_hauls_March7_2017.RData') # reload master hauls file to get full time series to develop prediction matrix
 
sppPrefs <- data.frame(SBTmed=numeric(), SBTmode=numeric(), SBTwtdMean=numeric(), SBTmedPA=numeric(), SBTmodePA=numeric(), SBTwtdMeanPA=numeric(), PAdev=numeric(), CPUEdev=numeric())
pdf(width=8, height=11, file='figures/Species_temp_prefs_Part3.pdf')
par(mfrow=c(4,2), tcl=-.1, mgp=c(1.5,0.1,0), mar=c(2.7,3,1.5,.5), oma=c(.1,.1,.1,.1))  
for(i in 1:length(allspp)){
  spp <- allspp[i]
  print(paste(spp, ': number ', i, sep=''))
  # Load habitat models
  load(paste(modfolder, 'TempPrefs_fitallreg_2017_', spp, '.RData', sep=''))
  mygam1 <- mods[[1]]
  mygam2 <- mods[[2]]
  rm(mods)
  # Load observed catch data for species to develop the prediction matrix
  catch <- dat[dat$sppocean == spp,] # catch weight data
  tows <- merge(catch, hauls, by='haulid', all.x=T) # join the haul info for these positive hauls
  # Get the total hauls completed on the relevant coast to figure out the range of seasonal temperatures we want to predict over
  if(grepl('_Atl', spp)){
    all.tows <- hauls[hauls$lon > -100,]
  }
  if(grepl('_Pac', spp)){
    all.tows <- hauls[hauls$lon < -100,]
  }
  # What is the most common region this species is caught in
  regionFreqs <- data.frame(table(tows$region))
  regionFact <- as.character(regionFreqs$Var1[regionFreqs$Freq == max(regionFreqs$Freq)])[1] # the '[1]' at the end is in case of a tie
   
  # The following species wasn't present in initial fit for NEUS, but was still assigned that based on above criteria (lots of observations in later part of NEUS time series that was omitted from fitting)
  if(i==149){
    regionFact <- 'DFO_ScotianShelf'
  }
  
  # make a temp range for graphing
  all.tows <- all.tows[!is.na(all.tows$bottemp),]
  SBTrange <- seq(range(all.tows$bottemp)[1], range(all.tows$bottemp)[2], by=.05)
  
  predSBT <- data.frame(bottemp=SBTrange, rugosity=mean(tows$rugosity, na.rm=T), GRAINSIZE=mean(tows$GRAINSIZE, na.rm=T), regionfact=regionFact)
     
  # Now do predictions for both prediction data.frames, each using both gams and integrate each separately and do calculations as on 'Thermal.env.broad.R'_cumsum function
  predictSBT <- cbind(predSBT, PREDPRESENCE = predict(mygam1, newdata=predSBT, type="response"), PREDCPUE = exp(predict(mygam2, newdata=predSBT, type="response")))
  predictSBT$PREDICT <- predictSBT$PREDPRESENCE * predictSBT$PREDCPUE
  
  SBTmax <- predictSBT$bottemp[which.max(predictSBT$PREDICT)]
  SBTmaxPA <- predictSBT$bottemp[which.max(predictSBT$PREDPRESENCE)]
  
  PREDsumSBT <- sum(predictSBT$PREDICT)/2
  predictSBT$CUM <- cumsum(predictSBT$PREDICT) 
  predictSBT$diff <- abs(predictSBT$CUM - PREDsumSBT)
  SBTmed <- predictSBT$bottemp[predictSBT$diff == min(predictSBT$diff)][1]
  
  PREDsumSBTPA <- sum(predictSBT$PREDPRESENCE)/2
  predictSBT$CUMPA <- cumsum(predictSBT$PREDPRESENCE) 
  predictSBT$diffPA <- abs(predictSBT$CUMPA - PREDsumSBTPA)
  SBTmedPA <- predictSBT$bottemp[predictSBT$diffPA == min(predictSBT$diffPA)][1]
  
  SBTwtdMean <- wtd.mean(x=predictSBT$bottemp, weights=predictSBT$PREDICT)
  SBTwtdMeanPA <- wtd.mean(x=predictSBT$bottemp, weights=predictSBT$PREDPRESENCE)
  
  PAdev <- modeldiag$dev.pres[modeldiag$sppocean==spp]
  CPUEdev <- modeldiag$dev.biomass[modeldiag$sppocean==spp]
  
  plot(PREDICT~bottemp, ylab='Predicted CPUE', xlab='SBT', main=spp, cex.main=.9, data=predictSBT)
  abline(v=SBTmed, col="red") # identify where the 50% is
  abline(v=SBTwtdMean, col="blue")  
  abline(v=SBTmedPA, col="red", lty=3) # identify where the 50% is
  abline(v=SBTwtdMeanPA, col="blue", lty=3)  
  
  sppPrefs[i,] <- data.frame(SBTmed=SBTmed, SBTmode=SBTmax, SBTwtdMean=SBTwtdMean, SBTmedPA=SBTmedPA, SBTmodePA=SBTmaxPA, SBTwtdMeanPA=SBTwtdMeanPA, PAdev=PAdev, CPUEdev=CPUEdev)
}
dev.off()

# Add species names
sppPrefs$sppocean <- allspp
save(sppPrefs, file='Temperature_preferences_Part3.RData')



# ===========================================================================================
# Other request for predicted biomass at temperature_I left this to Atlantic species only ===
# ===========================================================================================

AtlRef <- grepl('_Atl', allspp)
allspp <- allspp[AtlRef]
temps <- seq(0, 28, by=.01)

tempPreds <- matrix(data=NA, nrow=length(temps), ncol=length(allspp))
tempPredsPA <- matrix(data=NA, nrow=length(temps), ncol=length(allspp))

for(i in 1:length(allspp)){
  spp <- allspp[i]
  print(paste(spp, ': number ', i, sep=''))
  # Load habitat models
  load(paste(modfolder, 'TempPrefs_fitallreg_2017_', spp, '.RData', sep=''))
  mygam1 <- mods[[1]]
  mygam2 <- mods[[2]]
  rm(mods)
  # Load observed catch data for species to develop the prediction matrix
  catch <- dat[dat$sppocean == spp,] # catch weight data
  tows <- merge(catch, hauls, by='haulid', all.x=T) # join the haul info for these positive hauls

  # What is the most common region this species is caught in
  regionFreqs <- data.frame(table(tows$region))
  regionFact <- as.character(regionFreqs$Var1[regionFreqs$Freq == max(regionFreqs$Freq)])[1] # the '[1]' at the end is in case of a tie
  # The following species wasn't present in initial fit for NEUS, but was still assigned that based on above criteria (lots of observations in later part of NEUS time series that was omitted from fitting)
  if(i==90){
    regionFact <- 'DFO_ScotianShelf'
  }
  
  predSBT <- data.frame(bottemp=temps, rugosity=mean(tows$rugosity, na.rm=T), GRAINSIZE=mean(tows$GRAINSIZE, na.rm=T), regionfact=regionFact)
  # Now do predictions for both prediction data.frames, each using both gams and integrate each separately and do calculations as on 'Thermal.env.broad.R'_cumsum function
  predictSBT <- cbind(predSBT, PREDPRESENCE = predict(mygam1, newdata=predSBT, type="response"), PREDCPUE = exp(predict(mygam2, newdata=predSBT, type="response")))
  predictSBT$PREDICT <- predictSBT$PREDPRESENCE * predictSBT$PREDCPUE
  
  tempPreds[,i] <- predictSBT$PREDICT
  tempPredsPA[,i] <- predictSBT$PREDPRESENCE
}

TempPreds <- data.frame(cbind(botTemp=temps, tempPreds))
TempPredsPA <- data.frame(cbind(botTemp=temps, tempPredsPA))
cols <- c('botTemp', allspp)
colnames(TempPreds) <- cols
colnames(TempPredsPA) <- cols

save(TempPreds, TempPredsPA, file='Temperature_preferences_fullMatrix.RData')
  
