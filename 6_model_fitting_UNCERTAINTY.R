## Set working directories
setwd('/Users/jamesmorley/Documents/project_velocity')
modfolder <- 'output/CEmodels_Uncertainty_2018/'
       
library(mgcv); library(dismo)
library(gbm); library(MASS) 
library(caret)
   
load('data/hauls_catch_Dec2017.RData') # hauls and catch data
#allspp = sort(unique(dat$sppocean)) # for reference to spelling
sp <- 'centropristis striata_Atl'
  
# CAN SKIP THE FOLLOWING BLOCK IF RAW DATA ALREADY CREATED ==================
 
mydat <- dat[dat$sppocean==sp,] 
spdata <- merge(hauls, mydat, by='haulid', all.x = T, sort=F) # Add haul data and empty hauls
myocean <- head(spdata$ocean[spdata$presfit == TRUE])[1] # identify if this is a Pacific or Atlantic species
spdata <- spdata[spdata$ocean == myocean,] # trim master hauls file to the ocean of interest 
spdata$presfit[is.na(spdata$presfit)] <- FALSE
spdata <- droplevels(spdata) # drop the west or east coast 'regionfact' levels
myregions <- table(spdata$region[spdata$presfit==TRUE])  

  ##############################################
  # Set up data in regions with no presences
  ##############################################
  spdata$logwtcpue.pad <- spdata$logwtcpue # will have some zeros changed to -18 to allow abundance model fitting
  spdata$presfit.pad <- spdata$presfit # will have some FALSE changed to TRUE to allow abundance model fitting
   
  npres <- table(spdata$regionfact[spdata$presfit])
  if(any(npres < 1)){
    mywarn <- paste('Zero presences for', sp, 'in', paste(names(npres)[npres<1], collapse=', '))
    warning(mywarn)
    regstofill <- names(npres)[as.numeric(npres) == 0]
    spdata$regionfact[spdata$region %in% regstofill] <- names(npres)[which.max(npres)][1] # in regions with no observations, replace the region ID with that from a region with observations. this prevents a low region coefficient from explaining the zero observations.
     
    # if a region has no presences pick some absences and force them to low abundance presences for abundance model fitting
    for(j in 1:length(regstofill)){
      theseinds <- spdata$region == regstofill[j]
      fake0s <- sample(which(theseinds), size = round(0.1 * min(sum(theseinds), sum(spdata$presfit))))# Uses either 10% of species total samples, or 10% of hauls in a region, whichever is smaller
      spdata$logwtcpue.pad[fake0s] <- -23 
      spdata$presfit.pad[fake0s] <- TRUE
      print(paste(regstofill[j], ': Added', length(fake0s), 'fake zeros'))
    }
  }
  spdata <- droplevels(spdata)
   
# FOR SPECIES WITH VERY FEW OBSERVATIONS IN A REGION, ALSO CONFORM THEM TO MOST ABUNDANT REGION  
# Note that this is a bit of a custom-by-species job_so otherwise leave muted
  # Did this for two regions for Market Squid 
  
      #regstofill <- names(npres)[as.numeric(npres) < 15]
      #spdata$regionfact[spdata$region %in% regstofill] <- names(npres)[which.max(npres)][1] # in regions with no observations, replace the region ID with that from a region with observations. this prevents a low region coefficient from explaining the zero observations.
  
  # Also MODIFIED THIS SO THAT SOME FAKE ZEROS CAN STILL BE ADDED TO THOSE REGIONS
  # I'm doing this loop manually, b/c ea. time want to make sure I'm not replacing any of the actual values
      #theseinds <- spdata$region == regstofill[2]
      #fake0s <- sample(which(theseinds), size = round(0.1 * min(sum(theseinds), sum(spdata$presfit))))# Uses either 10% of species total samples, or 10% of hauls in a region, whichever is smaller
      #summary(spdata$presfit[fake0s]) # if any are TRUE, repeat previous line
      #spdata$logwtcpue.pad[fake0s] <- -23 
      #spdata$presfit.pad[fake0s] <- TRUE
      #print(paste(regstofill[2], ': Added', length(fake0s), 'fake zeros'))
      #spdata <- droplevels(spdata)
  
  ####################################################
  #Set up data for training and testing to evaluate performance
  ####################################################
   
  #Subset training and testing data by year by indexing row numbers (use first 80% to predict last 20%)
  spdata <- spdata[order(spdata$year,spdata$month),]
  # indices for both pres and abs
  ninds <- table(spdata$regionfact) # number of entries per region (regions as set up for fitting)
  traininds <- NULL; testinds <- NULL
  for(j in 1:length(ninds)){ # loop through each region to get first 80% and last 20%
    traininds <- c(traininds, which(as.character(spdata$regionfact) == names(ninds)[j])[1:round(ninds[j]*0.8)])
    testinds <- c(testinds, which(as.character(spdata$regionfact) == names(ninds)[j])[(round(ninds[j]*0.8)+1):ninds[j]])
  }

  # indices for only where present (for the abundance model), including fake zeros
  trainindsp <- intersect(traininds, which(spdata$presfit.pad))
  testindsp <- intersect(testinds, which(spdata$presfit.pad))
    
  # test if we have enough presences in testing and training sets (at least one per region as set up for model fitting)
  # Drop a region from testing data set if it was not represented in training data_this would lead to major bias in testing the P/A models and leads to an error for abundance models
  nprestrain <- table(spdata$regionfact[trainindsp])
  nprestest <- table(spdata$regionfact[testindsp])
  if(any(nprestrain == 0)){
    regstodrop <- names(nprestrain)[as.numeric(nprestrain) == 0]
    for(j in 1:length(regstodrop)){
      testinds <- setdiff(testinds, which(as.character(spdata$regionfact) == regstodrop[j]))
    } 
    testindsp <- intersect(testinds, which(spdata$presfit.pad))
    mywarn <- paste('Zero training presences for', sp, 'in', paste(names(nprestrain)[nprestrain==0], collapse=', '))
    warning(mywarn)
    # regstofill <- names(nprestrain)[as.numeric(nprestrain) == 0]
  }
  if(any(nprestest == 0)){
    mywarn <- paste('Zero testing presences for', sp, 'in', paste(names(nprestest)[nprestest==0], collapse=', '))
    warning(mywarn)
  }

  # warn if too few presences overall
  if(length(trainindsp)<2){
    mywarn <- paste('Only', length(trainindsp), 'presence values in training dataset for', sp)
    warning(mywarn)
  }
  if(length(testindsp)<2){
    mywarn <- paste('Only', length(testindsp), 'presence values in testing dataset for', sp)
    warning(mywarn)
  }
  
  # Need to drop NA values first for the biomass models_For BRTs
  spdataB <- spdata[!is.na(spdata$logwtcpue.pad),]
  spdataC <- spdata[trainindsp,]
  spdataC <- spdataC[!is.na(spdataC$logwtcpue.pad),]
  
  # Save the raw data, for use in graphing on script #9
  filename <- paste(modfolder, 'RawData_', sp, '_uncertainty2018.RData', sep='')  
  save(spdata, spdataB, spdataC, myregions, traininds, testinds, trainindsp, testindsp, file=filename)
  #load(filename)

#####################################################################################
# Format and then run each niche model
# This is assuming every species I look at is going to be present in multiple regions
#####################################################################################
   
# Standard GAM from PloS ONE paper ==========================================================
  mypresmod <- formula(presfit ~ s(SBT.seasonal) + s(SST.seasonal.mean) + s(SBT.min) + s(SBT.max) + s(SST.max) + s(rugosity) + s(GRAINSIZE) + regionfact)
  myabunmod <- formula(logwtcpue.pad ~ s(SBT.seasonal) + s(SST.seasonal.mean) + s(SBT.min) + s(SBT.max) + s(SST.max) + s(rugosity) + s(GRAINSIZE) + regionfact)
  
  # set 'gamma' penalty levels for gam to prevent overfitting_got this from a presentation by Simon Wood (https://people.maths.bris.ac.uk/~sw15190/mgcv/tampere/mgcv-advanced.pdf)
  gammaPAtt <- log(nrow(spdata[traininds,])) / 2
  gammaAbuntt <- log(nrow(spdata[trainindsp,])) / 2
  gammaPA <- log(nrow(spdata)) / 2
  gammaAbun <- log(nrow(spdata[spdata$presfit.pad,])) / 2
  
  mygam1tt <- gam(mypresmod, family="binomial", data=spdata[traininds,], select=TRUE, gamma=gammaPAtt) 
  mygam2tt <- gam(myabunmod, data=spdata[trainindsp,], select=TRUE, gamma=gammaAbuntt) # only fit where species is present
  mygam1 <- gam(mypresmod,family="binomial", data=spdata, select=TRUE, gamma=gammaPA)
  mygam2 <- gam(myabunmod, data=spdata[spdata$presfit.pad,], na.action='na.exclude', select=TRUE, gamma=gammaAbun) # only fit where spp is present

  filename <- paste(modfolder, 'nicheMods_GAMs_', sp, '_uncertainty2018.RData', sep='')  
  save(myregions, mygam1tt, mygam2tt, mygam1, mygam2, file=filename)
   
  mod1 <- mygam1
  mod2 <- mygam2
    
  beta.1 <- coef(mod1) # a vector of model coeffients for each curve section
  beta.2 <- coef(mod2) 
  
  Vb.1 <- vcov(mod1) # the variance of each coefficient (each curve section) in diagonal, and covariation b/w all factor levels and the sections of curves 
  Vb.2 <- vcov(mod2) 
  n <- 100 # number of iterations
  
  br.1 <- mvrnorm(n, beta.1, Vb.1) # simulate n rep coef vectors (for each section of curve) from posterior
  # each column of br is a replicate parameter vector drawn from posterior of distribution of the parameter
  br.2 <- mvrnorm(n, beta.2, Vb.2)
  
  filename <- paste(modfolder, 'Cluster_nicheMods_GAM_', sp, '_uncer2018.RData', sep='')  
  save(myregions, br.1, br.2, mod1, mod2, file=filename)
  
  
# Now fit a basic GLM =========================================================================================================================
library(brglm) # some biased reduction for binomial data
   
  mypresmod <- formula(presfit ~ SBT.seasonal + I(SBT.seasonal^2) + SST.seasonal.mean + I(SST.seasonal.mean^2) + SBT.min + I(SBT.min^2) + SBT.max + I(SBT.max^2) + SST.max + I(SST.max^2) + rugosity  + I(rugosity^2) + GRAINSIZE + I(GRAINSIZE^2) + regionfact)
  myabunmod <- formula(logwtcpue.pad ~ SBT.seasonal + I(SBT.seasonal^2) + SST.seasonal.mean + I(SST.seasonal.mean^2) + SBT.min + I(SBT.min^2) + SBT.max + I(SBT.max^2) + SST.max + I(SST.max^2) + rugosity  + I(rugosity^2) + GRAINSIZE + I(GRAINSIZE^2) + regionfact)
 
  myglm1tt <- brglm(mypresmod, data=spdata[traininds,], family=binomial)
  myglm2tt <- glm(myabunmod, data=spdata[trainindsp,])
  myglm1 <- brglm(mypresmod, data=spdata, family=binomial)
  myglm2 <- glm(myabunmod, data=spdata[spdata$presfit.pad,])
   
  filename <- paste(modfolder, 'nicheMods_GLMs_', sp, '_uncertainty2018.RData', sep='')  
  save(myregions, myglm1tt, myglm2tt, myglm1, myglm2, file=filename)

  mod1 <- myglm1
  mod2 <- myglm2
  
  beta.1 <- coef(mod1) # a vector of model coeffients
  beta.2 <- coef(mod2) 
  
  Vb.1 <- vcov(mod1) # the variance of each coefficient (each curve section) in diagonal, and covariation b/w all factor levels and the sections of curves 
  Vb.2 <- vcov(mod2) 
  n <- 100 # number of iterations
  
  br.1 <- mvrnorm(n, beta.1, Vb.1) # simulate n rep coef vectors 
  # each column of br is a replicate parameter vector drawn from posterior of distribution of the parameter
  br.2 <- mvrnorm(n, beta.2, Vb.2)
  
  filename <- paste(modfolder, 'Cluster_nicheMods_GLM_', sp, '_uncer2018.RData', sep='')  
  save(myregions, br.1, br.2, mod1, mod2, file=filename)
   
# Now fit a Boosted regression model =========================================================================================================================
  # uses cross validation from Elith et al. 2008_also the supplemental with that has step by steps
     
  tree.cmx <- myBRT1_caret$bestTune$interaction.depth
  le.rt <- myBRT1_caret$bestTune$shrinkage
 
  myBRT1 <- gbm.step(data=spdata,
                    gbm.x = c('SBT.seasonal', 'SST.seasonal.mean', 'SBT.min', 'SBT.max', 'SST.max', 'rugosity', 'GRAINSIZE', 'regionfact'), 
                    gbm.y = 'presfit',
                    family = "bernoulli",
                    tree.complexity = tree.cmx, # increasing this can reduce the number of trees needed to reach minima
                    learning.rate = le.rt, # increase this, may reduce number of trees_indicates how much ea. tree contributes to model
                    bag.fraction = 0.6,
                    n.trees = 50, step.size=500)

  # the graph produced shows changes in predictive deviance with the addition of new trees_so how does predictability change, if it hits a minimum then rises, that rise is showing the model is being overfit at too many trees
  BRT1.int <- gbm.interactions(myBRT1) # returns a list
  #BRT1.int$rank.list # gives the five most important interactions
  #BRT1.int$interactions # gives all the interactions scores
  #names(myBRT1) # objects stored in the model_e.g. on next line
  #myBRT1$contributions # sums to 100
  #head(myBRT1$fitted) #predicted values for observed data
  #head(myBRT1$fitted.vars) # variance of fitted values on response scale
  #gbm.plot.fits(myBRT1) # plots fitted values in relation to each of the predictors

  myBRT1tt <- gbm.step(data=spdata[traininds,],
                     gbm.x = c('SBT.seasonal', 'SST.seasonal.mean', 'SBT.min', 'SBT.max', 'SST.max', 'rugosity', 'GRAINSIZE', 'regionfact'), 
                     gbm.y = 'presfit',
                     family = "bernoulli",
                     tree.complexity = tree.cmx, # increasing this can also reduce the number of trees needed to reach minima
                     learning.rate = le.rt, # increase this, may reduce number of trees_indicates how much ea. tree contributes to model
                     bag.fraction = 0.6,
                     n.trees = 50, step.size=500)
  
  tree.cmx <- myBRT2_caret$bestTune$interaction.depth
  le.rt <- myBRT2_caret$bestTune$shrinkage

  myBRT2 <- gbm.step(data=spdataB,
                   gbm.x = c('SBT.seasonal', 'SST.seasonal.mean', 'SBT.min', 'SBT.max', 'SST.max', 'rugosity', 'GRAINSIZE', 'regionfact'), 
                   gbm.y = 'logwtcpue.pad',
                   family = "gaussian",
                   tree.complexity = tree.cmx,
                   learning.rate = le.rt,
                   bag.fraction = 0.6,
                   n.trees = 50, step.size=100) 
  BRT2.int <- gbm.interactions(myBRT2) # returns a list

  myBRT2tt <- gbm.step(data=spdataC,
                     gbm.x = c('SBT.seasonal', 'SST.seasonal.mean', 'SBT.min', 'SBT.max', 'SST.max', 'rugosity', 'GRAINSIZE', 'regionfact'), 
                     gbm.y = 'logwtcpue.pad',
                     family = "gaussian",
                     tree.complexity = tree.cmx,
                     learning.rate = le.rt,
                     bag.fraction = 0.6,
                     n.trees=50, step.size=100) 
  
  filename <- paste(modfolder, 'nicheMods_BRTs_', sp, '_uncertainty2018.RData', sep='')  
  save(myBRT1, myBRT2, myBRT1tt, myBRT2tt, file=filename)

  
# ======================================================================================
# Summarizing the results from the Caret package  
# ======================================================================================
  
  library(caret)
  
  filename <- paste(modfolder, sp, '_optimizeBRT_Jul2018.RData', sep='')
  load(filename)
  myBRT1_caret
  myBRT2_caret
  names(myBRT2_caret) # objects stored in the model_e.g. on next line
  
  plot(myBRT1_caret)
  plot(myBRT2_caret)
  varImp(myBRT1_caret) # scales most important predictor to 100 (caret package)
  varImp(myBRT2_caret) # scales most important predictor to 100 (caret package)
  summary(myBRT1_caret$finalModel) # determines % model explanatory power attributable to each variable (gbm package). Also makes a plot.
  summary(myBRT2_caret$finalModel) # determines % model explanatory power attributable to each variable (gbm package). Also makes a plot.
  myBRT1_caret$bestTune #final model features
  myBRT2_caret$bestTune #final model features
  

  
  
# ======================================================================================
# Saved code for mixed FX and tweedie distributions
# ======================================================================================
 
#myglmRAND1 <- glmmPQL(mypresmod, random=list(regionfact=~1), family=binomial, data=spdata)
#myglmRAND2 <- glmmPQL(myabunmod, random=list(regionfact=~1), family=gaussian, data=spdata[spdata$presfit.pad,])
#ranef(myglmRAND1) # gives random effects values
  
# Saved code for a Tweedie GAM
  # NEED TO CREATE A NEW BIOMASS COLUMN WITH ALL THE 'ZERO' VALUES IN PLACE
  #spdata$wtcpue_tweed <- spdata$wtcpue
  #spdata$wtcpue_tweed[is.na(spdata$wtcpue_tweed)] <- 0
  #spdata$logwtcpue_tweed <- log(spdata$wtcpue_tweed + 1)
  
  #mytweedmod <- formula(wtcpue_tweed ~ s(SBT.seasonal) + s(SST.seasonal.mean) + s(SBT.min) + s(SBT.max) + s(SST.max) + s(rugosity) + s(GRAINSIZE) + regionfact)
  #gammaTweed <- log(nrow(spdata)) / 2
  #mygam <- gam(mytweedmod, family=tw(link="log"), data=spdata, select=TRUE, gamma=gammaTweed) 
  