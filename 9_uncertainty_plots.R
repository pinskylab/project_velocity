# This uses output from '6_model_fitting_UNCERTAINTY.R' and makes figures
   
setwd('/Users/jamesmorley/Documents/project_velocity')
modfolder <- 'output/CEmodels_Uncertainty_2018/'
  
library(mgcv);library(dismo)
library(latticeExtra); library(gbm)
library(brglm) # some biased reduction for binomial data
library(MASS) 
  
sp <- 'paralichthys dentatus_Atl'
   
load('data/hauls_catch_Dec2017.RData') # hauls and catch data
load(paste(modfolder, 'nicheMods_GAMs_', sp, '_uncertainty2018.RData', sep=''))
load(paste(modfolder, 'nicheMods_GLMs_', sp, '_uncertainty2018.RData', sep=''))
load(paste(modfolder, 'nicheMods_BRTs_', sp, '_uncertainty2018.RData', sep=''))

filename <- paste(modfolder, 'RawData_', sp, '_uncertainty2018.RData', sep='')  
load(filename)

# ===========================================================
# MODEL DIAGNOSTICS TO TESTING DATA =========================
# ===========================================================

# Create table to store model diagnostics
n = rep(NA, 3)
modeldiag = data.frame(sppocean=n, model=n, npres=n, fakepres=n, npres.tr=n, npres.te=n, ntot=n, thresh.kappa=n, thresh=n, thresh.tss=n, thresh.trKap=n, thresh.tr=n, thresh.trtss=n, auc=n, auc.tt=n, tss.tt_prev=n, tss.tt_kap=n, tss.tt_tss=n, tssmax=n, tssmax.tt=n, kappa.tt_prev=n, kappa.tt_kap=n, kappa.tt_tss=n, kappamax=n, kappamax.tt=n, rpb=n, prev.obs=n, prev.pred.prev=n, prev.pred.kap=n, prev.pred.tss=n, pred_obsMedian=n, pred_obsMean=n, r2.biomass=n, r2.biomass.tt=n, r2.all=n, r2.all.tt=n, r2.pres.surv_year=n, r2.pres.surv_year.kap=n, r2.pres.surv_year.tss=n, r2.predPA.surv_year_mean=n, r2.abun.surv_year=n, r2.predPATT.surv_year=n, r2.predPATTkap.surv_year=n, r2.predPATTtss.surv_year=n, r2.predPATT.surv_year_mean=n, r2.predTT.surv_year=n, dev.pres=n, dev.biomass=n, stringsAsFactors=FALSE)

pdf(height=11, width=10, file=paste('figures/Uncertainty_mods/ModDiag_GLM-GAM-BRT_', sp, '.pdf', sep=''))
mods <- c('glm', 'gam', 'brt')

for(i in 1:length(mods)){
  if(mods[i]=='glm'){
    PAmod <- myglm1
    BIOMmod <- myglm2
    PAttmod <- myglm1tt
    BIOMttmod <- myglm2tt
  }
  if(mods[i]=='gam'){
    PAmod <- mygam1
    BIOMmod <- mygam2
    PAttmod <- mygam1tt
    BIOMttmod <- mygam2tt
  }
  if(mods[i]=='brt'){
    #PAmod <- myBRT1
    #BIOMmod <- myBRT2
    PAmod <- myBRT1
    BIOMmod <- myBRT2
    PAttmod <- myBRT1tt
    BIOMttmod <- myBRT2tt
  }
  
  ####################################################
  # Compare predictions to observations to assess model performance
  ####################################################
  modeldiag$sppocean[i] <- sp
  modeldiag$model[i] <- mods[i]
  modeldiag$npres[i] = sum(spdata$presfit)
  modeldiag$fakepres[i] = nrow(spdata[spdata$presfit.pad,]) - nrow(spdata[spdata$presfit,])
  modeldiag$npres.tr[i] = sum(spdata$presfit[traininds])
  modeldiag$npres.te[i] = sum(spdata$presfit[testinds])
  modeldiag$ntot[i] = nrow(spdata) # total rows of data for full prediction model
  
  if((mods[i]=='glm')|(mods[i]=='gam')){
    # For FULL model
    preds1 <- PAmod$fitted.values
    preds2 <- exp(predict(BIOMmod, newdata = spdata, type='response', na.action='na.pass')) # abundance predictions
    preds <- preds1*preds2
    
    # And for training/testing data set
    preds1tt <- predict(PAttmod, newdata = spdata[testinds,], type="response") 
    preds2tt <- exp(predict(BIOMttmod, newdata = spdata[testinds,], type='response'))
    predstt <- preds1tt*preds2tt
    # Need the predictions for the training data PAmod to develope an appropriate threshold value for the testing data
    preds1tr <- PAttmod$fitted.values 
  }
  if(mods[i]=='brt'){
    # For FULL model
    preds1 <- PAmod$fitted
    preds2 <- exp(predict.gbm(BIOMmod, newdata = spdata, type='response', n.trees=BIOMmod$n.trees, na.action='na.pass')) # abundance predictions
    preds <- preds1*preds2
    
    # And for training/testing data set
    preds1tt <- predict.gbm(PAttmod, newdata = spdata[testinds,], type="response", n.trees=PAttmod$n.trees) 
    preds2tt <- exp(predict.gbm(BIOMttmod, newdata = spdata[testinds,], type='response', n.trees=BIOMttmod$n.trees))
    predstt <- preds1tt*preds2tt
    # Need the predictions for the training data PAmod to develope an appropriate threshold value for the testing data
    preds1tr <- PAttmod$fitted
  }
  
  # evaluate model (use dismo package) ===================================
  
  # pick a threshold for pres/abs model evaluation (where needed)
  # evaluate function acts on a column of predicted values at true presences (p); and predicted values at true absences
  # evaluate provides the AUC_some summary stats_and one type of threshold value
  e <- evaluate(p=as.vector(preds1[spdata$presfit]), a=as.vector(preds1[!spdata$presfit]))
  #plot(e, 'ROC')
  #str(e) # shows options stored in 'e'
  #e@t # how you look at the different items w/n 'e'
  modeldiag$thresh.kappa[i] <- threshold(e, stat='kappa') # threshold based on full data set
  modeldiag$thresh[i] <- threshold(e, stat='prevalence') 
  e.ind.prev <- which(e@t == modeldiag$thresh[i]) # index for the chosen threshold's confusion matrix
  e.ind.kap <- which(e@t == modeldiag$thresh.kappa[i]) 
  conf <- as.data.frame(e@confusion) # confusion matrices (all thresholds consecutively from e@t)
  conf <- data.frame(cbind(tp=as.numeric(conf$tp), fp=as.numeric(conf$fp), fn=as.numeric(conf$fn), tn=as.numeric(conf$tn)))
  conf$threshold <- e@t
  conf$n <- conf$tp + conf$fp + conf$fn + conf$tn
  conf$tss <- with(conf, ((tp/(tp+fn)) + (tn/(tn+fp)))-1)
  conf$kappa <- with(conf, (((tp+tn)/n) - (((tp+fp)*(tp+fn))+((fn+tn)*(tn+fp)))/(n^2)) / (1 - ((((tp+fp)*(tp+fn))+((fn+tn)*(tn+fp)))/(n^2))))
  modeldiag$thresh.tss[i] <-  conf$threshold[conf$tss == max(conf$tss)]# threshold based on full data set
  
  
  # THIS DEVELOPES A THRESHOLD BASED ON THE TRAINING DATA
  e.tr <- evaluate(p=as.vector(preds1tr[spdata$presfit[traininds]]), a=as.vector(preds1tr[!spdata$presfit[traininds]]))
  modeldiag$thresh.tr[i] <- threshold(e.tr, stat='prevalence') # for testing/training
  modeldiag$thresh.trKap[i] <- threshold(e.tr, stat='kappa') # for testing/training
  e.ind.prev.tr <- which(e.tr@t == modeldiag$thresh.tr[i]) # index for the chosen threshold
  e.ind.kap.tr <- which(e.tr@t == modeldiag$thresh.trKap[i]) # index for the chosen threshold
  conf.tr <- as.data.frame(e.tr@confusion) # confusion matrices (all thresholds)
  conf.tr <- data.frame(cbind(tp=as.numeric(conf.tr$tp), fp=as.numeric(conf.tr$fp), fn=as.numeric(conf.tr$fn), tn=as.numeric(conf.tr$tn)))
  conf.tr$threshold <- e.tr@t
  conf.tr$n <- conf.tr$tp + conf.tr$fp + conf.tr$fn + conf.tr$tn
  conf.tr$tss <- with(conf.tr, ((tp/(tp+fn)) + (tn/(tn+fp)))-1)
  conf.tr$kappa <- with(conf.tr, (((tp+tn)/n) - (((tp+fp)*(tp+fn))+((fn+tn)*(tn+fp)))/(n^2)) / (1 - ((((tp+fp)*(tp+fn))+((fn+tn)*(tn+fp)))/(n^2))))
  modeldiag$thresh.trtss[i] <-  conf.tr$threshold[conf.tr$tss == max(conf.tr$tss)]# threshold based on full data set
  
  # Full model diagnostics_non-threshold ========================
  if(mods[i]=='glm'){
    library(modEvA) # this could also be avoided and just calculate %deviance manually
    modeldiag$dev.pres[i] = Dsquared(PAmod)
    modeldiag$dev.biomass[i] = Dsquared(BIOMmod)
    detach("package:modEvA", unload=TRUE) # this masks some of the dismo stuff_may have to reload dismo
    library(dismo) # just in case
  }
  if(mods[i]=='gam'){
    modeldiag$dev.pres[i] = summary(PAmod)$dev.expl
    modeldiag$dev.biomass[i] = summary(BIOMmod)$dev.expl
  }
  if(mods[i]=='brt'){
    # PAmod$self.statistics$mean.null #null deviance
    # PAmod$cv.statistics$deviance.mean #residual deviance
    modeldiag$dev.pres[i] = (PAmod$self.statistics$mean.null - PAmod$cv.statistics$deviance.mean)/PAmod$self.statistics$mean.null
    modeldiag$dev.biomass[i] = (BIOMmod$self.statistics$mean.null - BIOMmod$cv.statistics$deviance.mean)/BIOMmod$self.statistics$mean.null
  }
  
  modeldiag$auc[i] <- e@auc
  modeldiag$tssmax[i] <- max(conf$tss)
  modeldiag$kappamax[i] <- max(e@kappa, na.rm=TRUE) # maximum kappa
  modeldiag$rpb[i] <- cor(preds1, spdata$presfit) # point biserial correlation
  
  # testing data
  e.tt <- evaluate(p=as.vector(preds1tt[spdata$presfit[testinds]]), a=as.vector(preds1tt[!spdata$presfit[testinds]]))
  modeldiag$auc.tt[i] <- e.tt@auc
  conf.tt <- as.data.frame(e.tt@confusion) # confusion matrices (all thresholds)
  conf.tt <- data.frame(cbind(tp=as.numeric(conf.tt$tp), fp=as.numeric(conf.tt$fp), fn=as.numeric(conf.tt$fn), tn=as.numeric(conf.tt$tn)))
  conf.tt$threshold <- e.tt@t
  conf.tt$n <- conf.tt$tp + conf.tt$fp + conf.tt$fn + conf.tt$tn
  conf.tt$tss <- with(conf.tt, ((tp/(tp+fn)) + (tn/(tn+fp)))-1)
  conf.tt$kappa <- with(conf.tt, (((tp+tn)/n) - (((tp+fp)*(tp+fn))+((fn+tn)*(tn+fp)))/(n^2)) / (1 - ((((tp+fp)*(tp+fn))+((fn+tn)*(tn+fp)))/(n^2))))
  modeldiag$tssmax.tt[i] <- max(conf.tt$tss)
  modeldiag$kappamax.tt[i] <- max(e.tt@kappa, na.rm=TRUE) # maximum kappa
  
  # BASED ON THRESHOLDS DEVELOPED FROM THE TRAINING DATA
  # Determine which threshold from the training data best matches the confusion matrices from the testing data
  # Do this for the three methods of estimating the threshold
  thresh.ind.tt <- data.frame(t=e.tt@t, thresh.tr=modeldiag$thresh.tr[i], thresh.trKap=modeldiag$thresh.trKap[i], thresh.trtss=modeldiag$thresh.trtss[i])
  thresh.ind.tt$diff.prev <- abs(thresh.ind.tt$t - thresh.ind.tt$thresh.tr) 
  thresh.ind.tt$diff.kap <- abs(thresh.ind.tt$t - thresh.ind.tt$thresh.trKap) 
  thresh.ind.tt$diff.tss <- abs(thresh.ind.tt$t - thresh.ind.tt$thresh.trtss) 
  thresh.ind.tt$index <- c(1:nrow(thresh.ind.tt))
  e.ind.tt.prev <- thresh.ind.tt$index[thresh.ind.tt$diff.prev == min(thresh.ind.tt$diff.prev)[1]] # index for the chosen threshold
  e.ind.tt.kap <- thresh.ind.tt$index[thresh.ind.tt$diff.kap == min(thresh.ind.tt$diff.kap)[1]] # index for the chosen threshold
  e.ind.tt.tss <- thresh.ind.tt$index[thresh.ind.tt$diff.tss == min(thresh.ind.tt$diff.tss)[1]] # index for the chosen threshold
  # Based on threshold determined from training data_what is the kappa and tss value in testing data_for three methods of threshold calculations
  modeldiag$tss.tt_prev[i] = with(conf.tt[e.ind.tt.prev,], ((tp/(tp+fn)) + (tn/(tn+fp)))-1)
  modeldiag$kappa.tt_prev[i] = e.tt@kappa[e.ind.tt.prev]
  modeldiag$tss.tt_kap[i] = with(conf.tt[e.ind.tt.kap,], ((tp/(tp+fn)) + (tn/(tn+fp)))-1)
  modeldiag$kappa.tt_kap[i] = e.tt@kappa[e.ind.tt.kap]
  modeldiag$tss.tt_tss[i] = with(conf.tt[e.ind.tt.tss,], ((tp/(tp+fn)) + (tn/(tn+fp)))-1)
  modeldiag$kappa.tt_tss[i] = e.tt@kappa[e.ind.tt.tss]
  
  # DIAGNOSTICS FOR THE HURDLE MODELS
  # predicted vs. observed biomass when present
  modeldiag$r2.biomass[i] = cor(log(preds2[spdata$presfit & !is.na(spdata$wtcpue)]), spdata$logwtcpue[spdata$presfit & !is.na(spdata$wtcpue)])^2 # correlation of log(biomass) where present_need to remove NAs as some species have a few
  modeldiag$r2.biomass.tt[i] = cor(log(preds2tt[which(testinds %in% testindsp)]), spdata$logwtcpue[testindsp], use='complete.obs')^2 # only if presences exist in the test dataset
  
  # full model diagnostics
  spdata$wtcpue[spdata$presfit == FALSE] <- 0 
  modeldiag$r2.all[i] = cor(preds[!is.na(spdata$wtcpue)], spdata$wtcpue[!is.na(spdata$wtcpue)])^2 # overall biomass correlation
  modeldiag$r2.all.tt[i] = cor(predstt, spdata$wtcpue[testinds], use='complete.obs')^2 # overall biomass correlation. only makes sense to do this if the species is present at least once in the testing dataset
  
  test<-cbind(spdata,preds1,preds)
  # Plots of predicted vs. observed, with and w/o the smear_for the full model_this helps examine how realistic predictions are
  par(mfrow=c(3,2), mar=c(5, 4, 3, 1)) # A couple extra graphs are produced for species with training models
  plot(wtcpue~preds, ylab='observed cpue', xlab='predicted cpue', main=paste('<< ', mods[i], ' >>    ', 'Full Mod_Obs. vs pred. (haul-level)', sep=''), data=test)
  mtext(paste('Median diff=', summary(test$preds - test$wtcpue)[3], sep=''), side=3, line=-1.5)
  mtext(paste('Mean diff=', summary(test$preds - test$wtcpue)[4], sep=''), side=3, line=-3)
  
  modeldiag$pred_obsMedian[i] = summary(test$preds - test$wtcpue)[3]
  modeldiag$pred_obsMean[i] = summary(test$preds - test$wtcpue)[4]
  
  t1<-tapply(test$preds1 > modeldiag$thresh[i], list(test$year,test$surveyfact),mean) #average predicted occurence based on full model threshold
  t1b<-tapply(test$preds1 > modeldiag$thresh.kappa[i], list(test$year,test$surveyfact),mean) 
  t1c<-tapply(test$preds1 > modeldiag$thresh.tss[i], list(test$year,test$surveyfact),mean) 
  t2<-tapply(test$presfit, list(test$year,test$surveyfact),mean) #proportion of hauls with observed presence
  t2b<-tapply(test$preds1,list(test$year,test$surveyfact),mean)
  t3<-tapply(test$preds,list(test$year,test$surveyfact),mean) #average predicted abundance
  t4<-tapply(test$wtcpue,list(test$year,test$surveyfact),mean) #average observed abundance
  
  modeldiag$r2.pres.surv_year[i] <- round(cor(stack(as.data.frame(t2))[,1],stack(as.data.frame(t1))[,1],use="p")^2,2)
  modeldiag$r2.pres.surv_year.kap[i] <- round(cor(stack(as.data.frame(t2))[,1],stack(as.data.frame(t1b))[,1],use="p")^2,2)
  modeldiag$r2.pres.surv_year.tss[i] <- round(cor(stack(as.data.frame(t2))[,1],stack(as.data.frame(t1c))[,1],use="p")^2,2)
  modeldiag$r2.abun.surv_year[i] <- round(cor(stack(as.data.frame(t4))[,1],stack(as.data.frame(t3))[,1],use="p")^2,2)
  modeldiag$r2.predPA.surv_year_mean[i] <- round(cor(stack(as.data.frame(t2))[,1],stack(as.data.frame(t2b))[,1],use="p")^2,2) # no thresholds, just mean proportion of predicted as a continuous variable against mean proportion observed
  
  plot(stack(as.data.frame(t2))[,1],stack(as.data.frame(t1))[,1], col='green3', xlab="Prop. hauls w/ species present (by year and survey)",ylab="Predicted proportion of cells with presences", cex=0.5,main='Full PA mod_surv./yr means (green/prev; red/kappa; blue/tss)')
  points(stack(as.data.frame(t2))[,1],stack(as.data.frame(t1b))[,1],cex=0.5,col='red')
  points(stack(as.data.frame(t2))[,1],stack(as.data.frame(t1c))[,1],cex=0.5,col='blue')
  abline(a=0, b=1, lwd=.5, lty=3)
  mtext(paste("prev-r^2 =",modeldiag$r2.pres.surv_year[i]), side=3, line=-1.5)
  mtext(paste("kap-r^2 =",modeldiag$r2.pres.surv_year.kap[i]), side=3, line=-2.5)
  mtext(paste("tss-r^2 =",modeldiag$r2.pres.surv_year.tss[i]), side=3, line=-3.5)
  
  plot(stack(as.data.frame(t4))[,1],(stack(as.data.frame(t3))[,1]),xlab="Average wtcpue (by year and survey)",ylab="Average predicted wtcpue", cex=0.7,main='Full model predictions_survey/year means')
  mtext(paste("r^2 =",modeldiag$r2.abun.surv_year[i]), side=3, line=-1.5)
  
  plot(spdata$wtcpue[testinds]~predstt, ylab='observed cpue', xlab='predicted cpue (test data set)', main='Predictions from training model')
  mtext(paste('r2=', round(cor(predstt, spdata$wtcpue[testinds], use='complete.obs')^2, digits=3)), side=3, line=-1.5)
  
  plot(stack(as.data.frame(t2))[,1]~stack(as.data.frame(t2b))[,1],xlab="Mean probability of presence",ylab="Prop. of hauls with presences", cex=0.5,main='Full model PA mod_surv./yr means')
  abline(a=0, b=1, lwd=.5, lty=3)
  mtext(paste("r^2 =",modeldiag$r2.predPA.surv_year_mean[i]), side=3, line=-1.5)
  
  test <- cbind(spdata[testinds,],preds1tt, predstt)
  t1<-tapply(test$preds1tt > modeldiag$thresh.tr[i], list(test$year,test$surveyfact),mean) #average predicted occurence based on training model threshold
  t1b<-tapply(test$preds1tt > modeldiag$thresh.trKap[i], list(test$year,test$surveyfact),mean) #average predicted occurence based on training model threshold
  t1c<-tapply(test$preds1tt > modeldiag$thresh.trtss[i], list(test$year,test$surveyfact),mean) #average predicted occurence based on training model threshold
  t2<-tapply(test$presfit, list(test$year,test$surveyfact),mean) #proportion of hauls with observed presence
  t2b<-tapply(test$preds1tt,list(test$year,test$surveyfact),mean)
  t3<-tapply(test$predstt,list(test$year,test$surveyfact),mean) #average predicted abundance
  t4<-tapply(test$wtcpue,list(test$year,test$surveyfact),mean) #average observed abundance
  
  modeldiag$r2.predPATT.surv_year[i] <- round(cor(stack(as.data.frame(t2))[,1],stack(as.data.frame(t1))[,1],use="p")^2,2)
  modeldiag$r2.predPATTkap.surv_year[i] <- round(cor(stack(as.data.frame(t2))[,1],stack(as.data.frame(t1b))[,1],use="p")^2,2)
  modeldiag$r2.predPATTtss.surv_year[i] <- round(cor(stack(as.data.frame(t2))[,1],stack(as.data.frame(t1c))[,1],use="p")^2,2)
  modeldiag$r2.predTT.surv_year[i] <- round(cor(stack(as.data.frame(t4))[,1],stack(as.data.frame(t3))[,1],use="p")^2,2)
  modeldiag$r2.predPATT.surv_year_mean[i] <- round(cor(stack(as.data.frame(t2))[,1],stack(as.data.frame(t2b))[,1],use="p")^2,2) # no thresholds, just mean proportion of predicted as a continuous variable against mean proportion observed
  
  plot(stack(as.data.frame(t2))[,1]~stack(as.data.frame(t2b))[,1],xlab="Mean probability of presence",ylab="Prop. of hauls with presences", cex=0.5,main='Testing PA mod_surv./yr means')
  abline(a=0, b=1, lwd=.5, lty=3)
  mtext(paste("r^2 =",modeldiag$r2.predPATT.surv_year_mean[i]), side=3, line=-1.5)
  
  modeldiag$prev.obs[i] = mean(stack(as.data.frame(t2))[,1], na.rm=T)
  modeldiag$prev.pred.prev[i] = mean(stack(as.data.frame(t1))[,1], na.rm=T)
  modeldiag$prev.pred.kap[i] = mean(stack(as.data.frame(t1b))[,1], na.rm=T)
  modeldiag$prev.pred.tss[i] = mean(stack(as.data.frame(t1c))[,1], na.rm=T)
  
  plot(stack(as.data.frame(t2))[,1],stack(as.data.frame(t1))[,1], col='green3', xlab="Prop. of hauls with species present (by year and survey)",ylab="Predicted prop. of cells with presences", cex=0.5,main='Testing PA mod_surv./yr means (green/prev; red/kappa; blue/tss)')
  points(stack(as.data.frame(t2))[,1],stack(as.data.frame(t1b))[,1], cex=0.5,col='red') 
  points(stack(as.data.frame(t2))[,1],stack(as.data.frame(t1c))[,1], cex=0.5,col='blue') 
  abline(a=0, b=1, lwd=.5, lty=3)
  mtext(paste("prev-r^2 =",modeldiag$r2.predPATT.surv_year[i]), side=3, line=-1.5)
  mtext(paste("kap-r^2 =",modeldiag$r2.predPATTkap.surv_year[i]), side=3, line=-2.5)
  mtext(paste("tss-r^2 =",modeldiag$r2.predPATTtss.surv_year[i]), side=3, line=-3.5)
  
  plot(stack(as.data.frame(t4))[,1],(stack(as.data.frame(t3))[,1]),xlab="Average wtcpue (by year and survey)",ylab="Average predicted wtcpue", cex=0.7,main='Training/testing model_survey/year means')
  mtext(paste("r^2 =", modeldiag$r2.predTT.surv_year[i]), side=3, line=-1.5)
  
}
dev.off()

# write these files each time through the loop so that we can watch progress
save(modeldiag, file=paste("output/modeldiag_Uncer_", sp, "Feb2018.Rdata",sep=""))
write.csv(modeldiag, file=paste("output/modeldiag_Uncer_", sp, "Feb2018.csv",sep=""))



#===========================================
# Plot up all niche model outputs ==========
#==========================================
region <- 'NEFSC_NEUS'
load('data/prediction_files_Nov2017/predictionEASTnoBias_rcp85_jas_GFDL-CM3.RData')
abc <- pred.bathE[pred.bathE$year==2015,]
abc$regionfact <- as.factor(region)

abc$predPAgam <- predict(mygam1, newdata = abc, type='response', na.action='na.pass')
abc$predBIOMgam <- exp(predict(mygam2, newdata = abc, type='response', na.action='na.pass'))
abc$predBIOMgam_log <- predict(mygam2, newdata = abc, type='response', na.action='na.pass')
abc$predPAglm <- predict(myglm1, newdata = abc, type='response', na.action='na.pass')
abc$predBIOMglm <- exp(predict(myglm2, newdata = abc, type='response', na.action='na.pass'))
abc$predBIOMglm_log <- predict(myglm2, newdata = abc, type='response', na.action='na.pass')
abc$predPAbrt <- predict.gbm(myBRT1, newdata = abc, type='response', n.trees=myBRT1$n.trees, na.action='na.pass')
abc$predBIOMbrt <- exp(predict.gbm(myBRT2, newdata = abc, type='response', n.trees=myBRT2$n.trees, na.action='na.pass'))
abc$predBIOMbrt_log <- predict.gbm(myBRT2, newdata = abc, type='response', n.trees=myBRT2$n.trees, na.action='na.pass')

region <- as.factor(region)

pdf(file=paste("figures/Uncertainty_mods/NicheMods_", sp, ".pdf",sep=""),width=10,height=10)
 
# ============ GLM P-A 
par(mfrow=c(4,2),mar=c(2.5,2.5,1,.1), mgp=c(1.5,.5,0))
newdat <- data.frame(SBT.seasonal=seq(range(spdata$SBT.seasonal)[1], range(spdata$SBT.seasonal)[2], by=.1), SST.seasonal.mean=mean(spdata$SST.seasonal.mean[spdata$presfit]), SBT.min=mean(spdata$SBT.min[spdata$presfit]), SBT.max=mean(spdata$SBT.max[spdata$presfit]), SST.max=mean(spdata$SST.max[spdata$presfit]), rugosity=mean(spdata$rugosity[spdata$presfit]), GRAINSIZE=mean(spdata$GRAINSIZE[spdata$presfit]), regionfact=region)
plot(predict(myglm1, newdata=newdat, type='response')~newdat$SBT.seasonal, ylab='Prob. of Presence')
mtext('GLM', side=3, line=-1, adj=0, cex=.5)

newdat <- data.frame(SBT.seasonal=mean(spdata$SBT.seasonal[spdata$presfit]), SST.seasonal.mean=seq(range(spdata$SST.seasonal.mean)[1], range(spdata$SST.seasonal.mean)[2], by=.1), SBT.min=mean(spdata$SBT.min[spdata$presfit]), SBT.max=mean(spdata$SBT.max[spdata$presfit]), SST.max=mean(spdata$SST.max[spdata$presfit]), rugosity=mean(spdata$rugosity[spdata$presfit]), GRAINSIZE=mean(spdata$GRAINSIZE[spdata$presfit]), regionfact=region)
plot(predict(myglm1, newdata=newdat, type='response')~newdat$SST.seasonal.mean, ylab='Prob. of Presence')

newdat <- data.frame(SBT.seasonal=mean(spdata$SBT.seasonal[spdata$presfit]), SST.seasonal.mean=mean(spdata$SST.seasonal.mean[spdata$presfit]), SBT.min=seq(range(spdata$SBT.min)[1], range(spdata$SBT.min)[2], by=.1), SBT.max=mean(spdata$SBT.max[spdata$presfit]), SST.max=mean(spdata$SST.max[spdata$presfit]), rugosity=mean(spdata$rugosity[spdata$presfit]), GRAINSIZE=mean(spdata$GRAINSIZE[spdata$presfit]), regionfact=region)
plot(predict(myglm1, newdata=newdat, type='response')~newdat$SBT.min, ylab='Prob. of Presence')

newdat <- data.frame(SBT.seasonal=mean(spdata$SBT.seasonal[spdata$presfit]), SST.seasonal.mean=mean(spdata$SST.seasonal.mean[spdata$presfit]), SBT.min=mean(spdata$SBT.min[spdata$presfit]), SBT.max=seq(range(spdata$SBT.max)[1], range(spdata$SBT.max)[2], by=.1), SST.max=mean(spdata$SST.max[spdata$presfit]), rugosity=mean(spdata$rugosity[spdata$presfit]), GRAINSIZE=mean(spdata$GRAINSIZE[spdata$presfit]), regionfact=region)
plot(predict(myglm1, newdata=newdat, type='response')~newdat$SBT.max, ylab='Prob. of Presence')

newdat <- data.frame(SBT.seasonal=mean(spdata$SBT.seasonal[spdata$presfit]), SST.seasonal.mean=mean(spdata$SST.seasonal.mean[spdata$presfit]), SBT.min=mean(spdata$SBT.min[spdata$presfit]), SBT.max=mean(spdata$SBT.max[spdata$presfit]), SST.max=seq(range(spdata$SST.max)[1], range(spdata$SST.max)[2], by=.1), rugosity=mean(spdata$rugosity[spdata$presfit]), GRAINSIZE=mean(spdata$GRAINSIZE[spdata$presfit]), regionfact=region)
plot(predict(myglm1, newdata=newdat, type='response')~newdat$SST.max, ylab='Prob. of Presence')

newdat <- data.frame(SBT.seasonal=mean(spdata$SBT.seasonal[spdata$presfit]), SST.seasonal.mean=mean(spdata$SST.seasonal.mean[spdata$presfit]), SBT.min=mean(spdata$SBT.min[spdata$presfit]), SBT.max=mean(spdata$SBT.max[spdata$presfit]), SST.max=mean(spdata$SST.max[spdata$presfit]), rugosity=seq(range(spdata$rugosity)[1], range(spdata$rugosity)[2], by=.1), GRAINSIZE=mean(spdata$GRAINSIZE[spdata$presfit]), regionfact=region)
plot(predict(myglm1, newdata=newdat, type='response')~newdat$rugosity, ylab='Prob. of Presence')

newdat <- data.frame(SBT.seasonal=mean(spdata$SBT.seasonal[spdata$presfit]), SST.seasonal.mean=mean(spdata$SST.seasonal.mean[spdata$presfit]), SBT.min=mean(spdata$SBT.min[spdata$presfit]), SBT.max=mean(spdata$SBT.max[spdata$presfit]), SST.max=mean(spdata$SST.max[spdata$presfit]), rugosity=mean(spdata$rugosity[spdata$presfit]), GRAINSIZE=seq(range(spdata$GRAINSIZE)[1], range(spdata$GRAINSIZE)[2], by=.1), regionfact=region)
plot(predict(myglm1, newdata=newdat, type='response')~newdat$GRAINSIZE, ylab='Prob. of Presence')

# ============ GLM CPUE
par(mfrow=c(4,2),mar=c(2.5,2.5,1,.1), mgp=c(1.5,.5,0))
newdat <- data.frame(SBT.seasonal=seq(range(spdata$SBT.seasonal)[1], range(spdata$SBT.seasonal)[2], by=.1), SST.seasonal.mean=mean(spdata$SST.seasonal.mean[spdata$presfit]), SBT.min=mean(spdata$SBT.min[spdata$presfit]), SBT.max=mean(spdata$SBT.max[spdata$presfit]), SST.max=mean(spdata$SST.max[spdata$presfit]), rugosity=mean(spdata$rugosity[spdata$presfit]), GRAINSIZE=mean(spdata$GRAINSIZE[spdata$presfit]), regionfact=region)
plot(exp(predict(myglm2, newdata=newdat, type='response'))~newdat$SBT.seasonal, ylab='CPUE')
mtext('GLM', side=3, line=-1, adj=0, cex=.5)

newdat <- data.frame(SBT.seasonal=mean(spdata$SBT.seasonal[spdata$presfit]), SST.seasonal.mean=seq(range(spdata$SST.seasonal.mean)[1], range(spdata$SST.seasonal.mean)[2], by=.1), SBT.min=mean(spdata$SBT.min[spdata$presfit]), SBT.max=mean(spdata$SBT.max[spdata$presfit]), SST.max=mean(spdata$SST.max[spdata$presfit]), rugosity=mean(spdata$rugosity[spdata$presfit]), GRAINSIZE=mean(spdata$GRAINSIZE[spdata$presfit]), regionfact=region)
plot(exp(predict(myglm2, newdata=newdat, type='response'))~newdat$SST.seasonal.mean, ylab='CPUE')

newdat <- data.frame(SBT.seasonal=mean(spdata$SBT.seasonal[spdata$presfit]), SST.seasonal.mean=mean(spdata$SST.seasonal.mean[spdata$presfit]), SBT.min=seq(range(spdata$SBT.min)[1], range(spdata$SBT.min)[2], by=.1), SBT.max=mean(spdata$SBT.max[spdata$presfit]), SST.max=mean(spdata$SST.max[spdata$presfit]), rugosity=mean(spdata$rugosity[spdata$presfit]), GRAINSIZE=mean(spdata$GRAINSIZE[spdata$presfit]), regionfact=region)
plot(exp(predict(myglm2, newdata=newdat, type='response'))~newdat$SBT.min, ylab='CPUE')

newdat <- data.frame(SBT.seasonal=mean(spdata$SBT.seasonal[spdata$presfit]), SST.seasonal.mean=mean(spdata$SST.seasonal.mean[spdata$presfit]), SBT.min=mean(spdata$SBT.min[spdata$presfit]), SBT.max=seq(range(spdata$SBT.max)[1], range(spdata$SBT.max)[2], by=.1), SST.max=mean(spdata$SST.max[spdata$presfit]), rugosity=mean(spdata$rugosity[spdata$presfit]), GRAINSIZE=mean(spdata$GRAINSIZE[spdata$presfit]), regionfact=region)
plot(exp(predict(myglm2, newdata=newdat, type='response'))~newdat$SBT.max, ylab='CPUE')

newdat <- data.frame(SBT.seasonal=mean(spdata$SBT.seasonal[spdata$presfit]), SST.seasonal.mean=mean(spdata$SST.seasonal.mean[spdata$presfit]), SBT.min=mean(spdata$SBT.min[spdata$presfit]), SBT.max=mean(spdata$SBT.max[spdata$presfit]), SST.max=seq(range(spdata$SST.max)[1], range(spdata$SST.max)[2], by=.1), rugosity=mean(spdata$rugosity[spdata$presfit]), GRAINSIZE=mean(spdata$GRAINSIZE[spdata$presfit]), regionfact=region)
plot(exp(predict(myglm2, newdata=newdat, type='response'))~newdat$SST.max, ylab='CPUE')

newdat <- data.frame(SBT.seasonal=mean(spdata$SBT.seasonal[spdata$presfit]), SST.seasonal.mean=mean(spdata$SST.seasonal.mean[spdata$presfit]), SBT.min=mean(spdata$SBT.min[spdata$presfit]), SBT.max=mean(spdata$SBT.max[spdata$presfit]), SST.max=mean(spdata$SST.max[spdata$presfit]), rugosity=seq(range(spdata$rugosity)[1], range(spdata$rugosity)[2], by=.1), GRAINSIZE=mean(spdata$GRAINSIZE[spdata$presfit]), regionfact=region)
plot(exp(predict(myglm2, newdata=newdat, type='response'))~newdat$rugosity, ylab='CPUE')

newdat <- data.frame(SBT.seasonal=mean(spdata$SBT.seasonal[spdata$presfit]), SST.seasonal.mean=mean(spdata$SST.seasonal.mean[spdata$presfit]), SBT.min=mean(spdata$SBT.min[spdata$presfit]), SBT.max=mean(spdata$SBT.max[spdata$presfit]), SST.max=mean(spdata$SST.max[spdata$presfit]), rugosity=mean(spdata$rugosity[spdata$presfit]), GRAINSIZE=seq(range(spdata$GRAINSIZE)[1], range(spdata$GRAINSIZE)[2], by=.1), regionfact=region)
plot(exp(predict(myglm2, newdata=newdat, type='response'))~newdat$GRAINSIZE, ylab='CPUE')

# ============ GAMs
par(mfrow=c(1,1))
plot(mygam1,pages=1,scale=0, all.terms=TRUE, shade=T); mtext(paste(sp,"presence"),outer=T,line=-2)
plot(mygam2,pages=1,scale=0, all.terms=TRUE, shade=T); mtext(paste(sp,"abundance"),outer=T,line=-2)

# ============ BRTs

par(mfrow=c(1,1))
gbm.plot(myBRT1) # accounts for the average effects of all other variables, for each variable plot_I think effects average to zero for each subplot
gbm.plot(myBRT2) # accounts for the average effects of all other variables, for each variable plot_I think effects average to zero for each subplot

BRT1.int <- gbm.interactions(myBRT1) # returns a list
BRT2.int <- gbm.interactions(myBRT2) # returns a list
par(mfrow=c(3,2))
gbm.perspec(myBRT1,x=BRT1.int$rank.list[1,1], y=BRT1.int$rank.list[1,3], main='PA mods')
gbm.perspec(myBRT1,x=BRT1.int$rank.list[2,1], y=BRT1.int$rank.list[2,3], main='PA mods')
gbm.perspec(myBRT1,x=BRT1.int$rank.list[3,1], y=BRT1.int$rank.list[3,3], main='PA mods')
gbm.perspec(myBRT2,x=BRT2.int$rank.list[1,1], y=BRT2.int$rank.list[1,3], main='CPUE mods')
gbm.perspec(myBRT2,x=BRT2.int$rank.list[2,1], y=BRT2.int$rank.list[2,3], main='CPUE mods')
gbm.perspec(myBRT2,x=BRT2.int$rank.list[3,1], y=BRT2.int$rank.list[3,3], main='CPUE mods')

a <- levelplot(predPAgam~lonBathgrid*latBathgrid, main='PA gam', abc)
b <- levelplot(predPAglm~lonBathgrid*latBathgrid, main='PA glm', abc)
c <- levelplot(predPAbrt~lonBathgrid*latBathgrid, main='PA brt', abc)
d <- levelplot(predBIOMgam~lonBathgrid*latBathgrid, main='Biom gam', abc)
e <- levelplot(predBIOMglm~lonBathgrid*latBathgrid, main='Biom glm', abc)
f <- levelplot(predBIOMbrt~lonBathgrid*latBathgrid, main='Biom brt', abc)
g <- levelplot((predBIOMgam*predPAgam)~lonBathgrid*latBathgrid, main='Hurd. gam', abc)
h <- levelplot((predBIOMglm*predPAglm)~lonBathgrid*latBathgrid, main='Hurd. glm', abc)
i <- levelplot((predBIOMbrt*predPAbrt)~lonBathgrid*latBathgrid, main='Hurd. brt', abc)

# distribution of wtcpue to inform how/why biomass models are effective or not
par(mfrow=c(2,1))
hist(spdata$logwtcpue[spdata$presfit])
plot(lat~lon, pch=20, cex=.2, spdata[spdata$presfit==T,])

par(mfrow=c(2,2))
print(a, split=c(1,1,2,2), more=TRUE)
print(b, split=c(1,2,2,2), more=TRUE)
print(c, split=c(2,1,2,2))
print(d, split=c(1,1,2,2), more=TRUE)
print(e, split=c(1,2,2,2), more=TRUE)
print(f, split=c(2,1,2,2))
print(g, split=c(1,1,2,2), more=TRUE)
print(h, split=c(1,2,2,2), more=TRUE)
print(i, split=c(2,1,2,2))

# Presence absence binary maps_each with threshold
par(mfrow=c(2,2))
plot(latBathgrid~lonBathgrid, cex=.05, col='gray80', main='P/A_GAM_prev_thr.', abc)
points(latBathgrid~lonBathgrid, cex=.05, col='red', abc[abc$predPAgam > modeldiag$thresh[2],])
plot(latBathgrid~lonBathgrid, cex=.05, col='gray80', main='P/A_BRT_prev_thr.', abc)
points(latBathgrid~lonBathgrid, cex=.05, col='red', abc[abc$predPAbrt > modeldiag$thresh[3],])
plot(latBathgrid~lonBathgrid, cex=.05, col='gray80', main='P/A_GLM_prev_thr.', abc)
points(latBathgrid~lonBathgrid, cex=.05, col='red', abc[abc$predPAglm > modeldiag$thresh[1],])
par(mfrow=c(2,2))
plot(latBathgrid~lonBathgrid, cex=.05, col='gray80', main='P/A_GAM_kappa_thr.', abc)
points(latBathgrid~lonBathgrid, cex=.05, col='red', abc[abc$predPAgam > modeldiag$thresh.kappa[2],])
plot(latBathgrid~lonBathgrid, cex=.05, col='gray80', main='P/A_BRT_kappa_thr.', abc)
points(latBathgrid~lonBathgrid, cex=.05, col='red', abc[abc$predPAbrt > modeldiag$thresh.kappa[3],])
plot(latBathgrid~lonBathgrid, cex=.05, col='gray80', main='P/A_GLM_kappa_thr.', abc)
points(latBathgrid~lonBathgrid, cex=.05, col='red', abc[abc$predPAglm > modeldiag$thresh.kappa[1],])
par(mfrow=c(2,2))
plot(latBathgrid~lonBathgrid, cex=.05, col='gray80', main='P/A_GAM_tss_thr.', abc)
points(latBathgrid~lonBathgrid, cex=.05, col='red', abc[abc$predPAgam > modeldiag$thresh.tss[2],])
plot(latBathgrid~lonBathgrid, cex=.05, col='gray80', main='P/A_BRT_tss_thr.', abc)
points(latBathgrid~lonBathgrid, cex=.05, col='red', abc[abc$predPAbrt > modeldiag$thresh.tss[3],])
plot(latBathgrid~lonBathgrid, cex=.05, col='gray80', main='P/A_GLM_tss_thr.', abc)
points(latBathgrid~lonBathgrid, cex=.05, col='red', abc[abc$predPAglm > modeldiag$thresh.tss[1],])

par(mfrow=c(1,1))
pairs(abc[,c('predPAgam','predPAglm','predPAbrt','predBIOMgam_log','predBIOMglm_log','predBIOMbrt_log')])
dev.off()
  
 