## Set working directories
if(Sys.info()["nodename"] == "pinsky-macbookair"){
  setwd('~/Documents/Rutgers/Range projections/proj_ranges/')
  modfolder <- '../CEmodels/'
}     
if(Sys.info()["nodename"] == "amphiprion.deenr.rutgers.edu"){ # This works for Jim's amphiprion directory too.
  setwd('~/Documents/range_projections/')
  #.libPaths(new='~/R/x86_64-redhat-linux-gnu-library/3.1/') # For Jim's amphiprion this needs to be muted
  modfolder <- 'CEmodels/'
} 
if(Sys.info()["user"] == "jamesmorley"){
  setwd('/Users/jamesmorley/Documents/project_velocity')
  modfolder <- 'output/CEmodels/'
}
    
# options(repos='http://cran.rstudio.com/') # to get packages installed on amphiprion
# install.packages("dismo")

# Loop through species and fit models.
library(mgcv);library(dismo)
runname <- "fitallreg_2017" # use all regions in each fit that are from the same ocean

load('data/master_hauls_March7_2017.RData') # import master hauls file
load('data/dat_selectedspp_Feb_1_2017.Rdata')# load species catch data
dat <- dat[!(dat$wtcpue == 0 & dat$region == 'DFO_SoGulf'),] # the zeros in SoGulf are actual zeros (ie not just a scale issue) and thus are true absences
dat$wtcpue[dat$wtcpue == 0] <- 0.00002 # 'zeros' in dat are now species too light to register on scales_here a value below the lowest non-zero value is assigned_for transforming data
dat$logwtcpue <- log(dat$wtcpue)
# trim columns that are already in master hauls file, which will be merged in below with the hauls data
dat <- data.frame(haulid = dat$haulid, sppocean = dat$sppocean, Freq = dat$Freq, wtcpue = dat$wtcpue, logwtcpue = dat$logwtcpue, presfit = TRUE, stringsAsFactors = F)

# How common are 'NA' for wtcpue
# zeros <- data.frame(table(dat$sppocean, is.na(dat$wtcpue)))
# zeros <- zeros[zeros$Var2 == T,]
dat <- dat[!is.na(dat$wtcpue),] # drop NAs as it creates errors for some species. May want to go back and manually do 'oncorhynchus tshawytscha_Pac' as almost 10% of presence records have NA for wtcpue (tagging study?)
 
## small test to see which species are only marginally present in Newfoundland or WCAnn surveys
#nosurfregions <- c("DFO_NewfoundlandSpring","DFO_NewfoundlandFall","NWFSC_WCAnn")
#nsrspp <- sort(unique(dat$sppocean[dat$region %in% nosurfregions])) # spp present in at least one of these regions
#length(nsrspp)
#length(unique(dat$sppocean))
#tab <- table(dat$sppocean[dat$wtcpue>0 & dat$sppocean %in% nsrspp], dat$region[dat$wtcpue>0 & dat$sppocean %in% nsrspp]) # counts presences in all regions
#dim(tab)
#colnames(tab) <- c('AI', 'EBS', 'GOA', 'WCTri', 'Newf_F', 'Newf_S', 'Scot', 'SoGulf', 'NEUS_F', 'NEUS_S', 'WCAnn', 'GoMex') # shorter column names for easier printing to screen
#cmax <- apply(tab, 1, which.max) # region with the most catches
#nsrspp2 <- which(!(cmax %in% c(5,6,11))) # index for species that didn't have the highest presence count in Newf or WCAnn
#tab[nsrspp2,] # spp like Clupea harengus are common in many places
 
######################
# Start the big loop #
######################
     
# Create table to store model diagnostics
allspp = sort(unique(dat$sppocean))
n = rep(NA, length(allspp))
modeldiag = data.frame(sppocean=n, npres=n, fakepres=n, npres.tr=n, npres.te=n, ntot=n, thresh=n, thresh.tt=n, auc=n, auc.tt=n, tss=n, tss.tt=n, tssmax=n, tssmax.tt=n, acc=n, acc.tt=n, accmax=n, accmax.tt=n, sens=n, sens.tt=n, spec=n, spec.tt=n, kappa=n, kappa.tt=n, kappamax=n, kappamax.tt=n, rpb=n, smear=n, pred_obsMedian=n, pred_obsMean=n, r2.biomass=n, r2.biomass.tt=n, r2.all=n, r2.all.tt=n, r2.pres.surv_year=n, r2.abun.surv_year=n, r2.predTT.surv_year=n, dev.pres=n, dev.biomass=n, dev.pres.null=n, dev.biomass.null=n, stringsAsFactors=FALSE) # pred_obsMedianNosmear=n, pred_obsMeanNosmear=n, # tt is for training/testing model
 
#Open pdf to print figures 
pdf(file=paste("figures/CEmodelGAMsmooths/GAMs_PART5",runname,".pdf",sep=""),width=10,height=10)
 
options(warn=1) # print warnings as they occur
allwarnings = NULL
print(paste(length(allspp), 'models to fit'))
  
for(i in 652:length(allspp)){ #seq(from=41, to=705, by=20)){ 
  fittrain = TRUE
  mygam1tt <- mygam2tt <- mygam1 <- mygam2 <- preds <- preds1 <- preds2 <- predstt <- preds1tt <- preds2tt <- NULL 
       
  sp<-allspp[i]
  print(paste(i,sp, Sys.time()))
  mydat<-dat[dat$sppocean==sp,] 
 
  ####################################################
  # Add records for when there was a haul but no fish_note: this changed from previous version, same end product
  ####################################################
  
  spdata <- merge(hauls, mydat, by='haulid', all.x = T, sort=F) # Add empty hauls
  myocean <- head(spdata$ocean[spdata$presfit == TRUE])[1] # identify if this is a Pacific or Atlantic species
  spdata <- spdata[spdata$ocean == myocean,] # trim master hauls file to the ocean of interest 
  spdata$presfit[is.na(spdata$presfit)] <- FALSE
  spdata <- droplevels(spdata) # drop the west or east coast 'regionfact' levels
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
    testinds <- setdiff(testinds, which(as.character(spdata$regionfact) == regstodrop))
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
  # table(spdata$year[spdata$presfit]) # output number of presences by year
  
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
     
  ####################################	
  # Fit the training/testing models
  ####################################	
  
  # set 'gamma' penalty levels for gam to prevent overfitting_got this from a presentation by Simon Wood (https://people.maths.bris.ac.uk/~sw15190/mgcv/tampere/mgcv-advanced.pdf)
  gammaPA <- log(nrow(spdata[traininds,])) / 2
  gammaAbun <- log(nrow(spdata[trainindsp,])) / 2
 
    if(fittrain){
    try1 <- tryCatch({
      mygam1tt<-gam(mypresmod, family="binomial",data=spdata[traininds,], select=TRUE, gamma=gammaPA) 
      mygam2tt<-gam(myabunmod, data=spdata[trainindsp,], select=TRUE, gamma=gammaAbun) # only fit where species is present
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
    
  # fill in diagnostics
  modeldiag$sppocean[i] = sp
  modeldiag$npres[i] = freq
  modeldiag$fakepres[i] = nrow(spdata[spdata$presfit.pad,]) - nrow(spdata[spdata$presfit,])
  if(fittrain){
    modeldiag$npres.tr[i] = sum(spdata$presfit[traininds])
    modeldiag$npres.te[i] = sum(spdata$presfit[testinds])
  }
  modeldiag$ntot[i] = dim(spdata)[1]
 
  # evaluate model (use dismo package)
  # pick a threshold for pres/abs model evaluation (where needed)
  e <- evaluate(p=as.vector(preds1[spdata$presfit]), a=as.vector(preds1[!spdata$presfit]))
  modeldiag$thresh[i] <- threshold(e, stat='prevalence')
  e.ind <- which(e@t == modeldiag$thresh[i]) # index for the chosen threshold
  conf <- as.data.frame(e@confusion) # confusion matrices (all thresholds)
   
  if(length(testindsp)>0 & fittrain){ # need presences in the test dataset
    e.tt <- evaluate(p=as.vector(preds1tt[spdata$presfit[testinds]]), a=as.vector(preds1tt[!spdata$presfit[testinds]]))
    modeldiag$thresh.tt[i] <- threshold(e.tt, stat='prevalence') # for testing/training
    e.ind.tt <- which(e.tt@t == modeldiag$thresh.tt[i]) # index for the chosen threshold
    conf.tt <- as.data.frame(e.tt@confusion) # confusion matrices (all thresholds)
  }
   
  # pres/abs model diagnostics (no threshold needed)
  modeldiag$dev.pres[i] = summary(mygam1)$dev.expl
  modeldiag$auc[i] <- e@auc
  modeldiag$tssmax[i] <- max(with(conf, (tp*tn - fn*fp)/((tp+fp)*(fn+tn))), na.rm=TRUE) # maximum TSS (any threshold)
  modeldiag$accmax[i] <- max(with(conf, (tp+tn)/(tp+fp+fn+tn)), na.rm=TRUE) # maximum overall accuracy
  modeldiag$kappamax[i] <- max(e@kappa, na.rm=TRUE) # maximum kappa
  modeldiag$rpb[i] <- cor(preds1, spdata$presfit) # point biserial correlation
  
  if(length(testindsp)>0 & fittrain){ # need presences in the test dataset
    modeldiag$auc.tt[i] <- e.tt@auc
    modeldiag$tssmax.tt[i] <- max(with(conf.tt, (tp*tn - fn*fp)/((tp+fp)*(fn+tn))), na.rm=TRUE) # maximum TSS (any threshold)
    modeldiag$accmax.tt[i] <- max(with(conf.tt, (tp+tn)/(tp+fp+fn+tn)), na.rm=TRUE) # maximum overall accuracy
    modeldiag$kappamax.tt[i] <- max(e.tt@kappa, na.rm=TRUE) # maximum kappa
  }
   
  # true skill statistic, accuracy, kappa, and other stats that require a threshold
  modeldiag$tss[i] = 	with(conf[e.ind,], (tp*tn - fn*fp)/((tp+fp)*(fn+tn))) # TSS for chosen threshold
  modeldiag$acc[i] = 	with(conf[e.ind,], (tp+tn)/(tp+fp+fn+tn)) # overall accuracy
  modeldiag$sens[i] = with(conf[e.ind,], (tp)/(tp+fn)) # sensitivity: fraction of correctly predicted presences
  modeldiag$spec[i] = with(conf[e.ind,], (tn)/(tn+fp)) # specificity: fraction of correctly predicted absences
  modeldiag$kappa[i] = e@kappa[e.ind] # Cohen's kappa
  if(length(testindsp)>0 & fittrain){ # need presences in the test dataset
    modeldiag$tss.tt[i] = with(conf.tt[e.ind.tt,], (tp*tn - fn*fp)/((tp+fp)*(fn+tn)))
    modeldiag$acc.tt[i] = with(conf.tt[e.ind.tt,], (tp+tn)/(tp+fp+fn+tn)) # overall accuracy
    modeldiag$sens.tt[i] = with(conf.tt[e.ind.tt,], (tp)/(tp+fn)) # sensitivity: fraction of correctly predicted presences
    modeldiag$spec.tt[i] = with(conf.tt[e.ind.tt,], (tn)/(tn+fp)) # specificity: fraction of correctly predicted absences
    modeldiag$kappa.tt[i] = e.tt@kappa[e.ind.tt]
  }
   
  # abundance model diagnostics
  modeldiag$dev.biomass[i] = summary(mygam2)$dev.expl
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
  
  t1<-tapply(test$preds1,list(test$year,test$surveyfact),mean) #average predicted p(occur)
  t2<-tapply(test$presfit,list(test$year,test$surveyfact),mean) #proportion of hauls with presence
  t3<-tapply(test$preds,list(test$year,test$surveyfact),mean) #average predicted abundance
  t4<-tapply(test$wtcpue,list(test$year,test$surveyfact),mean) #average observed abundance
  
  presr2<-round(cor(stack(as.data.frame(t2))[,1],stack(as.data.frame(t1))[,1],use="p")^2,2)
  abunr2<-round(cor(stack(as.data.frame(t4))[,1],stack(as.data.frame(t3))[,1],use="p")^2,2)
  
  plot(stack(as.data.frame(t2))[,1],stack(as.data.frame(t1))[,1],xlab="Proportion of hauls with species present (by year and survey)",ylab="Mean predicted probability of occurrence", cex=0.5,main='Full PA model_survey/year means')
  mtext(paste("r^2 =",presr2), side=3, line=-1.5)
  plot(stack(as.data.frame(t4))[,1],(stack(as.data.frame(t3))[,1]),xlab="Average wtcpue (by year and survey)",ylab="Average predicted wtcpue", cex=0.5,main='Full model predictions_survey/year means')
  mtext(paste("r^2 =",abunr2), side=3, line=-1.5)
   
  if(fittrain){
    plot(spdata$wtcpue[testinds]~predstt, ylab='observed cpue', xlab='predicted cpue (test data set)', main='Predictions from training model')
    mtext(paste('r2=', round(cor(predstt, spdata$wtcpue[testinds], use='complete.obs')^2, digits=3)), side=3, line=-1.5)
  
    test <- cbind(spdata[testinds,],predstt)
    t3<-tapply(test$predstt,list(test$year,test$surveyfact),mean) #average predicted abundance
    t4<-tapply(test$wtcpue,list(test$year,test$surveyfact),mean) #average observed abundance
    
    abundTTr2 <- round(cor(stack(as.data.frame(t4))[,1],stack(as.data.frame(t3))[,1],use="p")^2,2)
    plot(stack(as.data.frame(t4))[,1],(stack(as.data.frame(t3))[,1]),xlab="Average wtcpue (by year and survey)",ylab="Average predicted wtcpue", cex=0.5,main='Training/testing model_survey/year means')
    mtext(paste("r^2 =", abundTTr2), side=3, line=-1.5)
    modeldiag$r2.predTT.surv_year[i]<-abundTTr2
  }
     
  modeldiag$r2.pres.surv_year[i]<-presr2
  modeldiag$r2.abun.surv_year[i]<-abunr2
  par(mfrow=c(1,1))
  
  ####################################################
  #### Save models for later projections
  ####################################################
  
  mods = list(mygam1 = mygam1, mygam2 = mygam2)
  
  sp <- gsub('/', '', sp) # would mess up saving the file
  
  save(mods, myregions, file=paste(modfolder, 'CEmods_',runname, '_', sp, '.RData', sep='')) 
  
  #think about figures to output - thermal response curves? spatial prediction by 1 deg square?
  #think about other data to save - number of pres/abs by region (?) 
  
  # write these files each time through the loop so that we can watch progress
  save(modeldiag,file=paste("output/modeldiag_PART5_",runname,".Rdata",sep=""))
  write.csv(modeldiag, file=paste("output/modeldiag_PART5_",runname,".csv",sep=""))
  
  write.csv(allwarnings, file=paste('output/warnings_PART5_', runname, '.csv', sep=''))
   
}
dev.off()
