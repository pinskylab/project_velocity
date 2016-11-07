
# A bunch of fragmented code used to develop and test models... in progress.


library(mgcv);library(ROCR)
source("CSquareCode.r") #taken from VMStools in googlecode
#start with dat file from fit_climate_envelopes.R

spp = sort(unique(dat$sppocean))
for(sp in spp){

#sp<-"gadus morhua_Atl"
#sp<- "merluccius albidus_Atl"
mydat<-dat[dat$sppocean==sp,] 
myocean<-mydat[1,"ocean"]
#Need to add records for when there was a haul but no fish
zs<-dat[!dat$haulid %in% mydat$haulid & dat$ocean==myocean,] #extract haulids where missing
# Extract one record of each haulid
matchpos<-match(unique(zs$haulid),zs$haulid)
zeros<-zs[matchpos,]

#Add/change relevant columns --> zero catch for target spp.
zeros$spp<-mydat$spp[1]
zeros$sppl<-mydat$sppl[1]
zeros$sppnew<-mydat$sppnew[1]
zeros$sppregion<-paste(zeros$sppnew,zeros$region,sep="_")
zeros$sppocean<-mydat$sppocean[1] #may need to add "ocean"
zeros$wtcpue<-0
zeros$presfit<-F
zeros$wtcpuena<-1e-4 #for consistency, but should be zero
zeros$wtcpuenal<-log(zeros$wtcpuena)

mydatw0<-rbind(mydat,zeros) #combine positive hauls with zero hauls

mydatw0$cs1<-as.factor(CSquare(mydatw0$lon,mydatw0$lat,1)) #1 degree squares

#Calculate mean catch per haul (by region??). For true biomass estimate, would need to stratify
ave.catch.wt<-tapply(mydatw0$wtcpue,list(mydatw0$year,mydatw0$region),mean,na.rm=T)
#avecatchyr<-cbind("year"=as.numeric(names(ave.catch.wt)),"biomassmean"=ave.catch.wt)
avecatchyr<-cbind("year"=as.numeric(rownames(ave.catch.wt)),ave.catch.wt)
avecatchyrreg<-cbind(stack(as.data.frame(ave.catch.wt)),rep(rownames(ave.catch.wt),6))
colnames(avecatchyrreg)<-c("biomassmean","region","year")
spdata<-merge(mydatw0,avecatchyrreg)
spdata<-spdata[complete.cases(spdata[,c("surftemp","bottemp","rugosity")]),]

#subset training and testing data by stratum: should be roughly 80/20

allstrata<-unique(spdata$stratum)
nstrata<-length(allstrata)
trainstrata<-sample(allstrata,round(nstrata*0.8))
teststrata<- allstrata[! allstrata %in% trainstrata]
dim(spdata[spdata$stratum %in% trainstrata,])[1]/dim(spdata)[1] #generally btw 75-85%

#subset training and testing data by 5-deg square: should be roughly 80/20

# allstrata<-unique(spdata$cs1)
# nstrata<-length(allstrata)
# trainstrata<-sample(allstrata,round(nstrata*0.8))
# teststrata<- allstrata[! allstrata %in% trainstrata]
# dim(spdata[spdata$cs1 %in% trainstrata,])[1]/dim(spdata)[1] 

#subset training and testing data by region (use 1-4 to predict 5th)

allstrata<-unique(spdata$region)
nstrata<-length(allstrata)
trainstrata<-sample(allstrata,round(nstrata*0.8))
teststrata<- allstrata[! allstrata %in% trainstrata]
dim(spdata[spdata$region %in% trainstrata,])[1]/dim(spdata)[1] 

#subset training and testing data by year (use 1963-2002 to predict 2003-2014)

allstrata<-unique(spdata$year)
trainstrata<-1963:2002
teststrata<- 2003:2014
dim(spdata[spdata$year %in% trainstrata,])[1]/dim(spdata)[1] 

#subset training and testing data by year (use first 80% to predict last 20%)

spdata1<-spdata[order(spdata$year,spdata$month),]
ninds<-dim(spdata)[1]
traininds<-1:round(ninds*0.8)
testinds<- (round(ninds*0.8)+1):ninds

#Fit model with training dataset
mygamm<-gam(presfit~s(bottemp,k=4)+s(surftemp,k=4)+s(logrugosity,k=4) + region ,family="binomial",data=spdata[traininds,]) #+ s(biomassmean,k=4)

#mygamm2<-gam(presfit~s(lon,lat),family="binomial",data=spdata)
#mygamm3<-gam(presfit~s(rugosity,k=4)+s(biomassmean),family="binomial",data=spdata) 

#Test model with testing data set
preds1<- predict(mygamm,spdata[testinds,],type="response")
preds1.rocr = prediction(predictions=as.numeric(preds1), labels=spdata$presfit[testinds])
perf = performance(preds1.rocr, 'tpr', 'fpr')
plot(perf)
auc = performance(preds1.rocr, 'auc')@y.values[[1]]

#Using all data to train and validate:
mygammA<-gam(presfit~s(bottemp,k=4)+s(surftemp,k=4)+s(logrugosity,k=4)+regionfact+s(biomassmean,k=4) ,family="binomial",data=spdata) #+ s(biomassmean,k=4)
preds1A<- predict(mygammA,spdata,type="response")
preds1.rocrA = prediction(predictions=as.numeric(preds1A), labels=spdata$presfit)
perfA = performance(preds1.rocrA, 'tpr', 'fpr')
plot(perfA,add=T,col=2)
aucA = performance(preds1.rocrA, 'auc')@y.values[[1]] #worse than with stratified??

#table(preds1>0.3,spdata$presfit[spdata$stratum %in% teststrata]) #0.878

#test abundance models
gama2 = gam(wtcpuenal ~ s(surftemp,k=4) + s(bottemp,k=4) +s(logrugosity,k=4) +s(biomassmean,k=4)+regionfact, family=gaussian, data=spdata[spdata$year %in% trainstrata & spdata$presfit,], select=TRUE, control=list(mgcv.half=30)) # 
gama2B = gam(wtcpuenal ~ s(surftemp,k=4) + s(bottemp,k=4) +s(logrugosity,k=4), family=gaussian, data=spdata[spdata$year %in% trainstrata & spdata$presfit,], select=TRUE, control=list(mgcv.half=30)) # 
preds2<- predict(gama2B,spdata[spdata$year %in% teststrata & spdata$presfit,],type="response")

plot(spdata[spdata$year %in% teststrata & spdata$presfit,"wtcpuenal"],preds2)
cor(spdata[spdata$year %in% teststrata & spdata$presfit,"wtcpuenal"],preds2)^2

#with full data
gama2A = gam(wtcpuenal ~ s(surftemp,k=4) + s(bottemp,k=4) +s(logrugosity,k=4) +s(biomassmean,k=4) +regionfact, family=gaussian, data=spdata[spdata$presfit,], select=TRUE, control=list(mgcv.half=30)) # 
#Without temperature covariates
gama2A2 = gam(wtcpuenal ~ s(logrugosity,k=4) +s(biomassmean,k=4)+regionfact, family=gaussian, data=spdata[spdata$presfit,], select=TRUE, control=list(mgcv.half=30)) # 
gama2A3 = gam(wtcpuenal ~ s(surftemp,k=4) + s(bottemp,k=4) , family=gaussian, data=spdata[spdata$presfit,], select=TRUE, control=list(mgcv.half=30)) # 

preds2A<- predict(gama2A3,spdata[spdata$presfit,],type="response")

plot(spdata[spdata$presfit,"wtcpuenal"],preds2A)
cor(spdata[spdata$presfit,"wtcpuenal"],preds2A)^2

predsA<-predict(gama2A,spdata,type="response")
preds<-preds1A*exp(predsA)
plot(spdata$bottemp,preds,col=spdata$regionfact)

plot(spdata[spdata$stratum %in% teststrata,"wtcpuenal"],preds) 


test<-cbind(spdata[spdata$year %in% teststrata,],preds1)
t1<-tapply(test$preds1,list(test$year,test$cs1),mean)
t2<-tapply(test$presfit,list(test$year,test$cs1),mean)
plot(stack(as.data.frame(t2))[,1],stack(as.data.frame(t1))[,1],xlab="Proportion of hauls with species present (by 1 deg square)",ylab="Mean predicted probability of occurrence", cex=0.5)
cor(stack(as.data.frame(t2))[,1],stack(as.data.frame(t1))[,1],use="p")

#What are reasonable performance measures for these models? For CEM, we want to know that temperature does a reasonably good job of predicting presence/absence before we use it to make projections. 
#ROC curves are overly optimistic if training and testing data are a) the same, or b) randomly selected without regard to spatial autocorrelation. If species data and predictors are spatially autocorrelated, then model may be prone to overfitting and predictive ability over estimated. This is especially important for models which are used to extrapolate in space or time. We want to do both.
#We should test model predictive ability (AUC for pres/abs, R^2 for abundance) by using a spatial split in the data. For instance, strata could be left out, or latitudinal bands such that approx 80% of data are used for fitting and 20% used for testing. This would necessitate not including a stratum effect, which would eventually ease projections.
#Another option is to use the C-square codes to divide region into blocks. Larger blocks may be more useful for model testing than the sometimes small strata.
#An extreme option is to leave out one region and predict with other regions (for species that occur in multiple regions). Performance is pretty terrible in this case.
#The tendency to overfit under spatial autocorrelation is a good reason to constrain GAM splines to only a few knots.
#Need to figure out how to use abundance data to predict pres/absence. In general, prediction of abundance (wtcpue) is poor, and probably largely driven by the biomassmean term. I'm inclined to leave biomass term out if we can't use it for projections, even though it obviously improves the model fit. Need to test how much temperature improves the fit of the abundance models.
#AUC curves generally show better (or similar) performance when models are trained and tested on different data rather than on the same data. This is contradictory.
#If probability of occurrance in a given haul is predicted to be 30%, what if we summarized these over a broader area (1/4 degree squares, stratum) and compared to the observed (mean) probability of occurrance? 