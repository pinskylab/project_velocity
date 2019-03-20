## Set working directories
if(Sys.info()["nodename"] == "amphiprion.deenr.rutgers.edu"){ # This works for Jim's amphiprion directory too.
  setwd('~/Documents/range_projections/')
  #.libPaths(new='~/R/x86_64-redhat-linux-gnu-library/3.1/') # For Jim's amphiprion this needs to be muted
  modfolder <- 'CEmodels/'
} 
if(Sys.info()["user"] == "jamesmorley"){
  setwd('/Users/jamesmorley/Documents/project_velocity')
  modfolder <- 'output/CEmodels_Nov2017/'
}
      
# options(repos='http://cran.rstudio.com/') # to get packages installed on amphiprion
# install.packages("dismo")

# Loop through species and fit models.
library(mgcv);library(dismo)
runname <- "fitallreg_2017" # use all regions in each fit that are from the same ocean
 
load('data/master_hauls_March7_2017.RData') # import master hauls file
load('data/dat_selectedspp_Feb_1_2017.Rdata')# load species catch data
dat <- dat[!(dat$wtcpue == 0 & dat$region == 'DFO_SoGulf'),] # the zeros in SoGulf are actual zeros (ie not just a scale issue) and thus are true absences
dat$wtcpue[dat$wtcpue == 0] <- 0.0002 # 'zeros' in dat are now species too light to register on scales_here a value below the lowest non-zero value is assigned_for transforming data
dat$wtcpue[dat$sppocean=="calappa sulcata_Atl" & is.na(dat$wtcpue)] <- 0.13 # the median weight of this species when observed
dat$logwtcpue <- log(dat$wtcpue)

# drop dubious GMex observations of temperate species
drops <- c('alosa pseudoharengus_Gmex', 'brevoortia tyrannus_Gmex', 'clupea harengus_Gmex', 'dipturus laevis_Gmex', 'paralichthys dentatus_Gmex',
           'hippoglossoides platessoides_Gmex', 'pseudopleuronectes americanus_Gmex', 'scomber scombrus_Gmex', 'cryptacanthodes maculatus_Gmex',
           'echinarachnius parma_Gmex', 'illex illecebrosus_Gmex', 'melanostigma atlanticum_Gmex', 'menidia menidia_Gmex', 'ovalipes ocellatus_Gmex','placopecten magellanicus_Gmex')
drops <- gsub('_Gmex', '_Atl', drops)
dat <- dat[!(dat$region=='SEFSC_GOMexFall' & dat$sppocean %in% drops),]
dat <- dat[!(dat$region=='SEFSC_GOMexSummer' & dat$sppocean %in% drops),]

# trim columns that are already in master hauls file, which will be merged in below with the hauls data
dat <- data.frame(haulid = dat$haulid, sppocean = dat$sppocean, Freq = dat$Freq, wtcpue = dat$wtcpue, logwtcpue = dat$logwtcpue, presfit = TRUE, stringsAsFactors = F)
dat <- dat[!dat$sppocean=='NO CATCH',]
#dat <- dat[!is.na(dat$wtcpue),] # drop NAs as it creates errors for some species. May want to go back and manually do 'oncorhynchus tshawytscha_Pac' as almost 10% of presence records have NA for wtcpue (tagging study?)
# Found a species naming error_correcting here
dat$sppocean[dat$sppocean=='heppasteria phygiana_Atl'] <- 'hippasteria phrygiana_Atl' # no need to aggregate as they never overlapped w/n hauls (only found in Newfoundland and they called it one way or the other)
save(dat, hauls, file='data/hauls_catch_Dec2017.RData') # For distribution to other people
  
######################
# Start the big loop #
######################
      
# Create table to store model diagnostics
allspp = sort(unique(dat$sppocean))
#allspp <- gsub('_Gmex', '_Atl', drops) # taken from MS_figures script
#allspp <- sort(sample(allspp, 30))
n = rep(NA, length(allspp))
modeldiag = data.frame(sppocean=n, npres=n, fakepres=n, npres.tr=n, npres.te=n, ntot=n, thresh.kappa=n, thresh=n, thresh.trKap=n, thresh.tr=n, auc=n, auc.tt=n, tss=n, tss.tt=n, tssmax=n, tssmax.tt=n, acc=n, acc.tt=n, accmax=n, accmax.tt=n, sens=n, sens.tt=n, spec=n, spec.tt=n, kappa=n, kappa.tt=n, kappamax=n, kappamax.tt=n, rpb=n, prev.obs=n, prev.pred.prev=n, prev.pred.kap=n, smear=n, pred_obsMedian=n, pred_obsMean=n, r2.biomass=n, r2.biomass.tt=n, r2.all=n, r2.all.tt=n, r2.pres.surv_year=n, r2.pres.surv_year.kap=n, r2.abun.surv_year=n, r2.predPATT.surv_year=n, r2.predPATTkap.surv_year=n, r2.predTT.surv_year=n, dev.pres=n, dev.biomass=n, dev.pres.null=n, dev.biomass.null=n, stringsAsFactors=FALSE) # pred_obsMedianNosmear=n, pred_obsMeanNosmear=n, # tt is for training/testing model
# redos <- c(40, 159)# I deleted some of these already that are non-salvagable
# redos2 <- c(40, 115, 159, 198, 201, 335, 361, 558, 651) # this group adds 6 other species that wasn't able to calculate an r2 value, but not one of our criteria so not going to look into it now

#Open pdf to print figures 
pdf(file=paste("figures/CEmodelGAMsmooths/GAMs_Revised_Nov2017_GmexDrops_",runname,".pdf",sep=""),width=10,height=10)
 
options(warn=1) # print warnings as they occur
allwarnings = NULL
print(paste(length(allspp), 'models to fit'))
    
for(i in 1:length(allspp)){  
  fittrain = TRUE
  mygam1tt <- mygam2tt <- mygam1 <- mygam2 <- preds <- preds1 <- preds2 <- predstt <- preds1tt <- preds2tt <- NULL 
         
  sp<-allspp[i]
  print(paste(i,sp, Sys.time()))
  mydat<-dat[dat$sppocean==sp,] 
 
  ####################################################
  # Add records for when there was a haul but no fish
  ####################################################
   
  spdata <- merge(hauls, mydat, by='haulid', all.x = T, sort=F) # Add empty hauls
  myocean <- head(spdata$ocean[spdata$presfit == TRUE])[1] # identify if this is a Pacific or Atlantic species
  spdata <- spdata[spdata$ocean == myocean,] # trim master hauls file to the ocean of interest 
  spdata$presfit[is.na(spdata$presfit)] <- FALSE
  #spdata$logwtcpue[spdata$region == "SEFSC_GOMex" & spdata$presfit==TRUE] <- NA
  #spdata$presfit[spdata$region == "SEFSC_GOMex" & spdata$presfit==TRUE] <- FALSE
  # The three anchovy species were not identified in SEUS_so should drop that region for those species_otherwise is erroneous zeros
  if(grepl("anchoa", sp)){
    spdata <- spdata[!spdata$regionfact=='SCDNR_SEUS',]
  }
  spdata <- droplevels(spdata) # drop the west or east coast 'regionfact' levels
  
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
    spdata$regionfact[spdata$region %in% regstofill] <- names(npres)[which.max(npres)][1] # in regions with no observations, replace the region ID with that from a region with observations. this prevents a low region coefficient from explaining the zero observations.
     
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
  #Set up data for training and testing to evaluate performance
  ####################################################
   
  #Subset training and testing data by year by indexing row numbers (use first 80% to predict last 20%)
  spdata<-spdata[order(spdata$year,spdata$month),]
  # indices for both pres and abs
  ninds<-table(spdata$regionfact) # number of entries per region (regions as set up for fitting)
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
    mywarn <- paste('Zero training presences for', i, sp, 'in', paste(names(nprestrain)[nprestrain==0], collapse=', '))
    allwarnings <- c(allwarnings, mywarn)
    warning(mywarn)
    # regstofill <- names(nprestrain)[as.numeric(nprestrain) == 0]
  }
  if(any(nprestest == 0)){
    mywarn <- paste('Zero testing presences for', i, sp, 'in', paste(names(nprestest)[nprestest==0], collapse=', '))
    allwarnings <- c(allwarnings, mywarn)
    warning(mywarn)
  }

  
  # warn if too few presences overall
  if(length(trainindsp)<2){
    mywarn <- paste('Only', length(trainindsp), 'presence values in training dataset for', i, sp)
    allwarnings <- c(allwarnings, mywarn)
    warning(mywarn)
  }
  if(length(testindsp)<2){
    mywarn <- paste('Only', length(testindsp), 'presence values in testing dataset for', i, sp)
    allwarnings <- c(allwarnings, mywarn)
    warning(mywarn)
  }
  
  
  # make sure we have at least 6 unique levels for each variable (necessary to fit gam with 4 knots)
  # look at training presence indices, since the most constraining (for mygam2tt)
  levs <- apply(spdata[trainindsp,c('SBT.seasonal', 'SST.seasonal.mean', 'SBT.min', 'SBT.max', 'SST.max', 'rugosity', 'GRAINSIZE')], 2, FUN=function(x) length(unique(x)))
  if(any(levs < 6)){
    mywarn <- paste("Not enough (>=6) unique levels in training presence set for", i, sp, ". Won't fit training models")
    allwarnings <- c(allwarnings, mywarn)
    warning(mywarn)
    fittrain = FALSE
  }	

  ####################################################
  # Figure out which model formula given data
  ####################################################
     
  #Default models. Leave out region factor if necessary
  # since fitallreg, using all regions in an ocean
  if(length(levels(spdata$regionfact))==1){
    mypresmod<-formula(presfit ~ s(SBT.seasonal) + s(SST.seasonal.mean) + s(SBT.min) + s(SBT.max) + s(SST.max) + s(rugosity) + s(GRAINSIZE))
    myabunmod<-formula(logwtcpue.pad ~ s(SBT.seasonal) + s(SST.seasonal.mean) + s(SBT.min) + s(SBT.max) + s(SST.max) + s(rugosity) + s(GRAINSIZE))
    mynullpresmod<-formula(presfit ~ s(rugosity) + s(GRAINSIZE)) #Null model w/o temp
    mynullabunmod<-formula(logwtcpue.pad ~ s(rugosity) + s(GRAINSIZE)) #Null model w/o temp
  } else {
    mypresmod<-formula(presfit ~ s(SBT.seasonal) + s(SST.seasonal.mean) + s(SBT.min) + s(SBT.max) + s(SST.max) + s(rugosity) + s(GRAINSIZE) + regionfact)
    myabunmod<-formula(logwtcpue.pad ~ s(SBT.seasonal) + s(SST.seasonal.mean) + s(SBT.min) + s(SBT.max) + s(SST.max) + s(rugosity) + s(GRAINSIZE) + regionfact)
    mynullpresmod<-formula(presfit ~ s(rugosity) + s(GRAINSIZE) + regionfact) #Null model w/o temp
    mynullabunmod<-formula(logwtcpue.pad ~ s(rugosity) + s(GRAINSIZE) + regionfact) #Null model w/o temp
  }
  # For training GAM, as it differs from full model for a few species_otherwise these species get dropped from analysis
  abc <- data.frame(reg.pres=nprestrain>0)
  if(length(abc[abc$reg.pres==TRUE,])==1){
    mypresmodtt<-formula(presfit ~ s(SBT.seasonal) + s(SST.seasonal.mean) + s(SBT.min) + s(SBT.max) + s(SST.max) + s(rugosity) + s(GRAINSIZE))
    myabunmodtt<-formula(logwtcpue.pad ~ s(SBT.seasonal) + s(SST.seasonal.mean) + s(SBT.min) + s(SBT.max) + s(SST.max) + s(rugosity) + s(GRAINSIZE))
  } else {
    mypresmodtt<-formula(presfit ~ s(SBT.seasonal) + s(SST.seasonal.mean) + s(SBT.min) + s(SBT.max) + s(SST.max) + s(rugosity) + s(GRAINSIZE) + regionfact)
    myabunmodtt<-formula(logwtcpue.pad ~ s(SBT.seasonal) + s(SST.seasonal.mean) + s(SBT.min) + s(SBT.max) + s(SST.max) + s(rugosity) + s(GRAINSIZE) + regionfact)
  }
  
  ####################################	
  # Fit the training/testing models
  ####################################	
  
  # set 'gamma' penalty levels for gam to prevent overfitting_got this from a presentation by Simon Wood (https://people.maths.bris.ac.uk/~sw15190/mgcv/tampere/mgcv-advanced.pdf)
  gammaPA <- log(nrow(spdata[traininds,])) / 2
  gammaAbun <- log(nrow(spdata[trainindsp,])) / 2
 
    if(fittrain){
    try1 <- tryCatch({
      mygam1tt<-gam(mypresmodtt, family="binomial",data=spdata[traininds,], select=TRUE, gamma=gammaPA) 
      mygam2tt<-gam(myabunmodtt, data=spdata[trainindsp,], select=TRUE, gamma=gammaAbun) # only fit where species is present
    }, error = function(e) { # ignore warnings, since no function to catch them
      mywarn <- paste('Error in training gam fitting for', i, sp, ':', e)
      assign('allwarnings', c(get('allwarnings', envir=.GlobalEnv), mywarn), envir=.GlobalEnv) # these assign outside the local scope are poor form in R. But not sure how else to do it here...
      assign('fittrain', FALSE, envir=.GlobalEnv) # if we hit an error in predictions, we can't calculate performance stats
      warning(mywarn)
    })
  }
   
  ####################################################
  #Fit models to All data (no test/training split)
  ####################################################
  
  gammaPA <- log(nrow(spdata)) / 2
  gammaAbun <- log(nrow(spdata[spdata$presfit.pad,])) / 2
  
  try2 <- tryCatch({
    mygam1<-gam(mypresmod,family="binomial",data=spdata, select=TRUE, gamma=gammaPA)
    mygam2<-gam(myabunmod,data=spdata[spdata$presfit.pad,], na.action='na.exclude', select=TRUE, gamma=gammaAbun) # only fit where spp is present
    mygam1null<-gam(mynullpresmod,family="binomial",data=spdata, select=TRUE, gamma=gammaPA)
    mygam2null<-gam(mynullabunmod,data=spdata[spdata$presfit.pad,], select=TRUE, gamma=gammaAbun) # only fit where spp is present
    
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
    
  ####################################################
  # Compare predictions to observations to assess model performance
  ####################################################
       
  # For FULL model
  preds1 <- mygam1$fitted.values
  preds2 <- exp(predict(mygam2, newdata = spdata, type='response', na.action='na.pass')) # abundance predictions
  smear = mean(exp(mygam2$residuals)) # smearing estimator for re-transformation bias (see Duan 1983, http://www.herc.research.va.gov/include/page.asp?ID=cost-regression)
  preds <- preds1*preds2#*smear # adds the bias correction as well
  # preds.nosmear <- preds1*preds2 # without the smear value_just for comparison
  preds[preds<0] = 0
    
  # And for training/testing data set
  if(fittrain){ 
    try3 <- tryCatch({
      preds1tt <- predict(mygam1tt, newdata = spdata[testinds,], type="response") 
      preds2tt <- exp(predict(mygam2tt, newdata = spdata[testinds,], type='response'))
      #smeartt = mean(exp(mygam2tt$residuals)) # smearing estimator for re-transformation bias (see Duan 1983, http://www.herc.research.va.gov/include/page.asp?ID=cost-regression)
      predstt <- preds1tt*preds2tt#*smeartt
      predstt[predstt<0] = 0
    }, error = function(e) {
      assign('fittrain', FALSE, envir=.GlobalEnv) # if we hit an error in predictions, we can't calculate performance stats
      mywarn <- paste('Error in predicting to test data for', i, sp, ':', e)
      assign('allwarnings', c(get('allwarnings', envir=.GlobalEnv), mywarn), envir=.GlobalEnv) # these assigns outside the local scope are poor form in R. But not sure how else to do it here...
      warning(mywarn)
    })
  }  
  # Need the predictions for the training data to develope an appropriate threshold value for the testing data
  if(fittrain){ 
    try3 <- tryCatch({
      preds1tr <- mygam1tt$fitted.values 
    }, error = function(e) {
      assign('fittrain', FALSE, envir=.GlobalEnv) # if we hit an error in predictions, we can't calculate performance stats
      mywarn <- paste('Error in developing prediction values from training models', i, sp, ':', e)
      assign('allwarnings', c(get('allwarnings', envir=.GlobalEnv), mywarn), envir=.GlobalEnv) # these assigns outside the local scope are poor form in R. But not sure how else to do it here...
      warning(mywarn)
    })
  }
    
  # fill in diagnostics
  modeldiag$sppocean[i] = sp
  modeldiag$npres[i] = freq
  modeldiag$fakepres[i] = nrow(spdata[spdata$presfit.pad,]) - nrow(spdata[spdata$presfit,])
  if(fittrain){ # Number of true presences in the training and testing data sets
    modeldiag$npres.tr[i] = sum(spdata$presfit[traininds])
    modeldiag$npres.te[i] = sum(spdata$presfit[testinds])
  }
  modeldiag$ntot[i] = dim(spdata)[1] # total rows of data for full prediction model
  
  # evaluate model (use dismo package)
  # pick a threshold for pres/abs model evaluation (where needed)
  # evaluate function acts on a column of predicted values at true presences (p); and predicted values at true absences
  # evaluate provides the AUC_some summary stats_and one type of threshold value
  e <- evaluate(p=as.vector(preds1[spdata$presfit]), a=as.vector(preds1[!spdata$presfit]))
  #plot(e, 'ROC')
  #str(e)
  #e@t # how you look at the different items w/n 'e'
  modeldiag$thresh.kappa[i] <- threshold(e, stat='kappa')
  modeldiag$thresh[i] <- threshold(e, stat='prevalence')
  e.ind <- which(e@t == modeldiag$thresh[i]) # index for the chosen threshold's confusion matrix
  conf <- as.data.frame(e@confusion) # confusion matrices (all thresholds consecutively from e@t)
   
  # THIS DEVELOPES A THRESHOLD BASED ON THE TRAINING DATA_the last two lines in this block may not be needed anymore
  if(length(trainindsp)>0 & fittrain){ # need presences in the test dataset
    e.tr <- evaluate(p=as.vector(preds1tr[spdata$presfit[traininds]]), a=as.vector(preds1tr[!spdata$presfit[traininds]]))
    modeldiag$thresh.tr[i] <- threshold(e.tr, stat='prevalence') # for testing/training
    modeldiag$thresh.trKap[i] <- threshold(e.tr, stat='kappa') # for testing/training
    e.ind.tr <- which(e.tr@t == modeldiag$thresh.tr[i]) # index for the chosen threshold
    conf.tr <- as.data.frame(e.tr@confusion) # confusion matrices (all thresholds)
  }
   
  # Full model pres/abs diagnostics (no threshold needed)
  modeldiag$dev.pres[i] = summary(mygam1)$dev.expl
  modeldiag$auc[i] <- e@auc
  #modeldiag$tssmax[i] <- max(with(conf, (tp*tn - fn*fp)/((tp+fp)*(fn+tn))), na.rm=TRUE) # maximum TSS (any threshold)
  # Needed to convert these into numerics for tssmax as the integers get too large for R to calculate
  modeldiag$tssmax[i] <- max(with(conf, (as.numeric(tp)*as.numeric(tn) - as.numeric(fn)*as.numeric(fp))/((as.numeric(tp)+as.numeric(fp))*(as.numeric(fn)+as.numeric(tn)))), na.rm=TRUE) # maximum TSS (any threshold)
  modeldiag$accmax[i] <- max(with(conf, (tp+tn)/(tp+fp+fn+tn)), na.rm=TRUE) # maximum overall accuracy
  modeldiag$kappamax[i] <- max(e@kappa, na.rm=TRUE) # maximum kappa
  modeldiag$rpb[i] <- cor(preds1, spdata$presfit) # point biserial correlation
  # testing data
  if(length(testindsp)>0 & fittrain){ # need presences in the test dataset
    e.tt <- evaluate(p=as.vector(preds1tt[spdata$presfit[testinds]]), a=as.vector(preds1tt[!spdata$presfit[testinds]]))
    modeldiag$auc.tt[i] <- e.tt@auc
    conf.tt <- as.data.frame(e.tt@confusion) # confusion matrices (all thresholds)
    modeldiag$tssmax.tt[i] <- max(with(conf.tt, (as.numeric(tp)*as.numeric(tn) - as.numeric(fn)*as.numeric(fp))/((as.numeric(tp)+as.numeric(fp))*(as.numeric(fn)+as.numeric(tn)))), na.rm=TRUE) # maximum TSS (any threshold)
    modeldiag$accmax.tt[i] <- max(with(conf.tt, (tp+tn)/(tp+fp+fn+tn)), na.rm=TRUE) # maximum overall accuracy
    modeldiag$kappamax.tt[i] <- max(e.tt@kappa, na.rm=TRUE) # maximum kappa
  }
   
  # true skill statistic, accuracy, kappa, and other stats that require a threshold
  modeldiag$tss[i] = 	with(conf[e.ind,], (as.numeric(tp)*as.numeric(tn) - as.numeric(fn)*as.numeric(fp))/((as.numeric(tp)+as.numeric(fp))*(as.numeric(fn)+as.numeric(tn)))) # TSS for chosen threshold
  modeldiag$acc[i] = 	with(conf[e.ind,], (tp+tn)/(tp+fp+fn+tn)) # overall accuracy
  modeldiag$sens[i] = with(conf[e.ind,], (tp)/(tp+fn)) # sensitivity: fraction of correctly predicted presences
  modeldiag$spec[i] = with(conf[e.ind,], (tn)/(tn+fp)) # specificity: fraction of correctly predicted absences
  modeldiag$kappa[i] = e@kappa[e.ind] # Cohen's kappa
  # THESE are BASED ON THRESHOLDS DEVELOPED FROM THE TRAINING DATA
  if(length(testindsp)>0 & fittrain){ # need presences in the test dataset
    # Find which threshold from the training data best matches the confusion matrices from the testing data
    thresh.ind.tt <- data.frame(t=e.tt@t, thresh.tr=modeldiag$thresh.tr[i])
    thresh.ind.tt$diff <- abs(thresh.ind.tt$t - thresh.ind.tt$thresh.tr) 
    thresh.ind.tt$index <- c(1:nrow(thresh.ind.tt))
    e.ind.tt <- thresh.ind.tt$index[thresh.ind.tt$diff == min(thresh.ind.tt$diff)[1]] # index for the chosen threshold
    modeldiag$tss.tt[i] = with(conf.tt[e.ind.tt,], (as.numeric(tp)*as.numeric(tn) - as.numeric(fn)*as.numeric(fp))/((as.numeric(tp)+as.numeric(fp))*(as.numeric(fn)+as.numeric(tn))))
    modeldiag$acc.tt[i] = with(conf.tt[e.ind.tt,], (tp+tn)/(tp+fp+fn+tn)) # overall accuracy
    modeldiag$sens.tt[i] = with(conf.tt[e.ind.tt,], (tp)/(tp+fn)) # sensitivity: fraction of correctly predicted presences
    modeldiag$spec.tt[i] = with(conf.tt[e.ind.tt,], (tn)/(tn+fp)) # specificity: fraction of correctly predicted absences
    modeldiag$kappa.tt[i] = e.tt@kappa[e.ind.tt]
  }
     
  # abundance model diagnostics
  modeldiag$dev.biomass[i] = summary(mygam2)$dev.expl
  # predicted vs. observed biomass when present
  modeldiag$r2.biomass[i] = cor(log(preds2[spdata$presfit & !is.na(spdata$wtcpue)]), spdata$logwtcpue[spdata$presfit & !is.na(spdata$wtcpue)])^2 # correlation of log(biomass) where present_need to remove NAs as some species have a few
  if(length(testindsp)>0 & fittrain){
    modeldiag$r2.biomass.tt[i] = cor(preds2tt[which(testinds %in% testindsp)], spdata$logwtcpue[testindsp], use='complete.obs')^2 # only if presences exist in the test dataset
  }
    
  # full model diagnostics
  spdata$wtcpue[spdata$presfit == FALSE] <- 0 
  modeldiag$r2.all[i] = cor(preds[!is.na(spdata$wtcpue)], spdata$wtcpue[!is.na(spdata$wtcpue)])^2 # overall biomass correlation
  if(length(testindsp)>0 & fittrain) modeldiag$r2.all.tt[i] = cor(predstt, spdata$wtcpue[testinds], use='complete.obs')^2 # overall biomass correlation. only makes sense to do this if the species is present at least once in the testing dataset
  
  #Compare to models without temperature to ultimately calculation %explained by temp terms
  modeldiag$dev.pres.null[i] = summary(mygam1null)$dev.expl
  modeldiag$dev.biomass.null[i] = summary(mygam2null)$dev.expl
  
  test<-cbind(spdata,preds1,preds)
  # Plots of predicted vs. observed, with and w/o the smear_for the full model_this helps examine how realistic predictions are
  par(mfrow=c(2,2), mar=c(5, 4, 3, 1))
  if(fittrain){
    par(mfrow=c(3,2), mar=c(5, 4, 3, 1)) # A couple extra graphs are produced for species with training models
  }
  
  plot(wtcpue~preds, ylab='observed cpue', xlab='predicted cpue', main='Full Model_Obs. vs pred. (no smear)', data=test)
  mtext(paste('Median diff=', summary(test$preds - test$wtcpue)[3], sep=''), side=3, line=-1.5)
  mtext(paste('Mean diff=', summary(test$preds - test$wtcpue)[4], sep=''), side=3, line=-3)
  # mtext(paste('Smear=', round(smear, digits=2), sep=''), side=3, line=-4.5)
  # plot(wtcpue~preds.nosmear, ylab='observed cpue', xlab='predicted cpue (w/o smear)', main='Full Model_Obs. vs pred. (no smear)', data=test)
  # mtext(paste('Median diff=', summary(test$preds.nosmear - test$wtcpue)[3], sep=''), side=3, line=-1.5)
  # mtext(paste('Mean diff=', summary(test$preds.nosmear - test$wtcpue)[4], sep=''), side=3, line=-3)

  modeldiag$smear[i] = smear
  modeldiag$pred_obsMedian[i] = summary(test$preds - test$wtcpue)[3]
  modeldiag$pred_obsMean[i] = summary(test$preds - test$wtcpue)[4]
  # modeldiag$pred_obsMedianNosmear[i] = summary(test$preds.nosmear - test$wtcpue)[3]
  # modeldiag$pred_obsMeanNosmear[i] = summary(test$preds.nosmear - test$wtcpue)[4]
   
  t1<-tapply(test$preds1 > modeldiag$thresh[i], list(test$year,test$surveyfact),mean) #average predicted occurence based on full model threshold
  t1b<-tapply(test$preds1 > modeldiag$thresh.kappa[i], list(test$year,test$surveyfact),mean) 
  t2<-tapply(test$presfit, list(test$year,test$surveyfact),mean) #proportion of hauls with observed presence
  t3<-tapply(test$preds,list(test$year,test$surveyfact),mean) #average predicted abundance
  t4<-tapply(test$wtcpue,list(test$year,test$surveyfact),mean) #average observed abundance
  
  presr2<-round(cor(stack(as.data.frame(t2))[,1],stack(as.data.frame(t1))[,1],use="p")^2,2)
  presr2kap <-round(cor(stack(as.data.frame(t2))[,1],stack(as.data.frame(t1b))[,1],use="p")^2,2)
  abunr2<-round(cor(stack(as.data.frame(t4))[,1],stack(as.data.frame(t3))[,1],use="p")^2,2)
  
  plot(stack(as.data.frame(t2))[,1],stack(as.data.frame(t1))[,1],xlab="Proportion of hauls with species present (by year and survey)",ylab="Predicted proportion of cells with presences", cex=0.5,main='Full PA mod_survey/year means (black/prev; red/kappa')
  points(stack(as.data.frame(t2))[,1],stack(as.data.frame(t1b))[,1],cex=0.5,col='red')
  mtext(paste("r^2 =",presr2), side=3, line=-1.5)
  mtext(paste("r^2 =",presr2kap), side=3, line=-2.5)
  plot(stack(as.data.frame(t4))[,1],(stack(as.data.frame(t3))[,1]),xlab="Average wtcpue (by year and survey)",ylab="Average predicted wtcpue", cex=0.5,main='Full model predictions_survey/year means')
  mtext(paste("r^2 =",abunr2), side=3, line=-1.5)
   
  if(fittrain){
    plot(spdata$wtcpue[testinds]~predstt, ylab='observed cpue', xlab='predicted cpue (test data set)', main='Predictions from training model')
    mtext(paste('r2=', round(cor(predstt, spdata$wtcpue[testinds], use='complete.obs')^2, digits=3)), side=3, line=-1.5)
  
    test <- cbind(spdata[testinds,],preds1tt, predstt)
    t1<-tapply(test$preds1tt > modeldiag$thresh.tr[i], list(test$year,test$surveyfact),mean) #average predicted occurence based on training model threshold
    t1b<-tapply(test$preds1tt > modeldiag$thresh.trKap[i], list(test$year,test$surveyfact),mean) #average predicted occurence based on training model threshold
    t2<-tapply(test$presfit, list(test$year,test$surveyfact),mean) #proportion of hauls with observed presence
    t3<-tapply(test$predstt,list(test$year,test$surveyfact),mean) #average predicted abundance
    t4<-tapply(test$wtcpue,list(test$year,test$surveyfact),mean) #average observed abundance
    
    
    presrTTr2<-round(cor(stack(as.data.frame(t2))[,1],stack(as.data.frame(t1))[,1],use="p")^2,2)
    presrTTr2Kap<-round(cor(stack(as.data.frame(t2))[,1],stack(as.data.frame(t1b))[,1],use="p")^2,2)
    abundTTr2 <- round(cor(stack(as.data.frame(t4))[,1],stack(as.data.frame(t3))[,1],use="p")^2,2)
    modeldiag$prev.obs[i] = mean(stack(as.data.frame(t2))[,1], na.rm=T)
    modeldiag$prev.pred.prev[i] = mean(stack(as.data.frame(t1))[,1], na.rm=T)
    modeldiag$prev.pred.kap[i] = mean(stack(as.data.frame(t1b))[,1], na.rm=T)
    
    plot(stack(as.data.frame(t2))[,1],stack(as.data.frame(t1))[,1],xlab="Proportion of hauls with species present (by year and survey)",ylab="Predicted proportion of cells with presences", cex=0.5,main='Testing PA mod_survey/year means (black/prev; red/kappa')
    points(stack(as.data.frame(t2))[,1],stack(as.data.frame(t1b))[,1], cex=0.5,col='red') 
    mtext(paste("r^2 =",presrTTr2), side=3, line=-1.5)
    mtext(paste("r^2 =",presrTTr2Kap), side=3, line=-2.5)
    modeldiag$r2.predPATT.surv_year[i]<-presrTTr2
    modeldiag$r2.predPATTkap.surv_year[i]<-presrTTr2Kap
    
    plot(stack(as.data.frame(t4))[,1],(stack(as.data.frame(t3))[,1]),xlab="Average wtcpue (by year and survey)",ylab="Average predicted wtcpue", cex=0.5,main='Training/testing model_survey/year means')
    mtext(paste("r^2 =", abundTTr2), side=3, line=-1.5)
    modeldiag$r2.predTT.surv_year[i]<-abundTTr2
  }
       
  modeldiag$r2.pres.surv_year[i]<-presr2
  modeldiag$r2.pres.surv_year.kap[i]<-presr2kap
  modeldiag$r2.abun.surv_year[i]<-abunr2
  par(mfrow=c(1,1))
   
  ####################################################
  #### Save models for later projections
  ####################################################
  
  mods = list(mygam1 = mygam1, mygam2 = mygam2)
  
  sp <- gsub('/', '', sp) # would mess up saving the file
   
  save(mods, file=paste(modfolder, 'CEmods_Nov2017_GMEXdrop_', runname, '_', sp, '.RData', sep='')) 
  
  # write these files each time through the loop so that we can watch progress
  save(modeldiag, file=paste("output/modeldiag_Nov2017_GMEXdrop_", runname,".Rdata",sep=""))
  write.csv(modeldiag, file=paste("output/modeldiag_Nov2017_GMEXdrop_", runname,".csv",sep=""))
  
  write.csv(allwarnings, file=paste('output/warnings_Nov2017_GMEXdrop_', runname, '.csv', sep=''))
}
dev.off()
    
# BE SURE TO DROP i=260 from modeldiag and allspp
modeldiag <- modeldiag[-260,]
save(modeldiag, file=paste("output/modeldiag_Nov2017_", runname,".Rdata",sep=""))
write.csv(modeldiag, file=paste("output/modeldiag_Nov2017_", runname,".csv",sep=""))

modeldiag_gmex <- modeldiag
rm(modeldiag)
load(paste("output/modeldiag_Nov2017_", runname,".Rdata",sep=""))
modeldiag <- modeldiag[!(modeldiag$sppocean %in% allspp),]
modeldiag <- rbind(modeldiag, modeldiag_gmex)
modeldiag <- modeldiag[order(modeldiag$sppocean),]
save(modeldiag, file=paste("output/modeldiag_Nov2017_", runname,".Rdata",sep=""))
write.csv(modeldiag, file=paste("output/modeldiag_Nov2017_", runname,".csv",sep=""))

