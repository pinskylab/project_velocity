## Set working directories
if(Sys.info()["nodename"] == "pinsky-macbookair"){
	setwd('~/Documents/Rutgers/Range projections/proj_ranges/')
	modfolder <- '../CEmodels/'
	}
if(Sys.info()["nodename"] == "amphiprion.deenr.rutgers.edu"){
	setwd('~/Documents/range_projections/')
	.libPaths(new='~/R/x86_64-redhat-linux-gnu-library/3.1/') # so that it can find my old packages
	modfolder <- 'CEmodels/'
	}
if(Sys.info()["user"] == "lauren"){
	setwd('~/backup/NatCap/proj_ranges/')
	modfolder <- 'output/CEmodels/'
}




# Loop through species and fit models.
library(mgcv);library(ROCR)

runname <- "fitallreg" # use all regions in each fit that are from the same ocean


load("data/dat_selectedspp.Rdata")
	dat$logwtcpue <- log(dat$wtcpue)
	dat$survey <- dat$region #keeps NEFSC as two separate surveys in "survey"
	dat$surveyfact <-as.factor(dat$survey)
	dat$region[dat$region %in% c("NEFSC_NEUSFall","NEFSC_NEUSSpring")] <- "NEFSC_NEUS" #NEFSC as one region in "region"
	dat$regionfact <- as.factor(dat$region) # the version ot use for model fitting
	# if using season
	dat$season <- as.factor(c(rep('wi', 3), rep('sp', 3), rep('su', 3), rep('fa', 3))[dat$month])
	dat$regseas <- as.factor(paste(dat$region,dat$season,sep="_"))

allspp = sort(unique(dat$sppocean))
n = rep(NA, length(allspp))
modeldiag = data.frame(sppocean=n, npres=n, npres.tr=n, npres.te=n, ntot=n, auc=n, auc.tt=n, tss=n, tss.tt=n, r2.biomass=n, r2.biomass.tt=n, r2.all=n, r2.all.tt=n, r2.pres.1deg=n, r2.abun.1deg=n, dev.pres=n, dev.biomass=n,dev.pres.null=n, dev.biomass.null=n, stringsAsFactors=FALSE) # tt is for training/testing model


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

#Open pdf to print figures 
pdf(file=paste("figures/CEmodelGAMsmooths/GAMs_",runname,".pdf",sep=""),width=8,height=6)

options(warn=1) # print warnings as they occur
allwarnings = NULL
print(paste(length(allspp), 'models to fit'))

for(i in 1:length(allspp)){ 
#for(i in c(18,22,25,27,28)){ # for testing on trouble species
	fittrain = TRUE
	mygam1tt <- mygam2tt <- mygam1 <- mygam2 <- preds <- preds1 <- preds2 <- predstt <- preds1tt <- preds2tt <- NULL

	sp<-allspp[i]
	print(paste(i,sp, Sys.time()))

	mydat<-dat[dat$sppocean==sp,] 
	mydat$logwtcpue[is.infinite(mydat$logwtcpue)] <- NA
	myregions<-unique(mydat$region)
	myhauls<-unique(mydat$haulid)
	myseasons<-unique(mydat$season)
	myocean <- unique(mydat$ocean)




	####################################################
	# Add records for when there was a haul but no fish
	####################################################
	# without season. expand to all regions in the same ocean.
	zs<-dat[!(dat$haulid %in% myhauls) & dat$ocean %in% myocean,] #extract haulids where this species is missing

	matchpos<-match(unique(zs$haulid),zs$haulid) # Extract one record of each haulid
	zeros<-zs[matchpos,]
	#Add/change relevant columns --> zero catch for target spp.
	zeros$spp<-mydat$spp[1]
	zeros$sppl<-mydat$sppl[1]
	zeros$sppnew<-mydat$sppnew[1]
	zeros$sppocean<-mydat$sppocean[1] #may need to add "ocean"
	zeros$wtcpue<-0
	zeros$logwtcpue <- NA
	zeros$presfit<-FALSE

	mydatw0<-rbind(mydat,zeros) #combine positive hauls with zero hauls

	##########################################################
	# Calculate mean catch per haul by region for this taxon #
	##########################################################

	# (For true biomass estimate, would need to stratify the mean)
	ave.catch.wt<-tapply(mydatw0$wtcpue,list(mydatw0$year,mydatw0$region),mean,na.rm=T)
	# transform into a data.frame
	if(dim(ave.catch.wt)[2]<2) { # only one region
		avecatchyrreg<-cbind(as.data.frame(ave.catch.wt), rep(colnames(ave.catch.wt), dim(ave.catch.wt)[2]), rownames(ave.catch.wt))
		colnames(avecatchyrreg)<-c("biomassmean","region","year")
	} else { # multiple regions
		avecatchyrreg<-cbind(stack(as.data.frame(ave.catch.wt)), rep(rownames(ave.catch.wt), dim(ave.catch.wt)[2]))
		colnames(avecatchyrreg)<-c("biomassmean","region","year")
	}
	spdata<-merge(mydatw0,avecatchyrreg)
	
	#Also save average biomassmean for each region (across all years) to use in later predictions.
	avemeanbiomass<-apply(ave.catch.wt,2,mean,na.rm=T)

	####################################################
	# Trim data to complete cases
	####################################################

	spdata<-spdata[complete.cases(spdata[,c("surftemp","bottemp","rugosity","presfit")]),]
	spdata <- droplevels(spdata)


	##############################################
	# Set up data in regions with no presences
	##############################################
	spdata$logwtcpue.pad <- spdata$logwtcpue # has some zeros changed to -18 to allow abundance model fitting
	spdata$presfit.pad <- spdata$presfit # has some FALSE changed to TRUE to allow abundance model fitting


	npres <- table(spdata$regionfact[spdata$presfit])
	if(any(npres < 1)){
		mywarn <- paste('Zero presences for', i, sp, 'in', paste(names(npres)[npres<1], collapse=', '), 'so adding some as -23 (1e-10) to allow abundance model fitting')
		allwarnings <- c(allwarnings, mywarn)
		warning(mywarn)
		regstofill <- names(npres)[as.numeric(npres) == 0]

		spdata$regionfact[spdata$region %in% regstofill] <- names(npres)[which.max(npres)] # in regions with no observations, replace the region ID with that from a region with observations. this prevents a low region coefficient from explaining the zero observations.

		# if a region has no presences
		# pick some absences and force them to low abundance presences for abundance model fitting
		for(j in 1:length(regstofill)){
			theseinds <- spdata$region == regstofill[j]
			fake0s <- sample(which(theseinds), size = round(0.1 * sum(theseinds)))
			spdata$logwtcpue.pad[fake0s] <- -23
			spdata$presfit.pad[fake0s] <- TRUE
			print(paste(regstofill[j], ': Added', length(fake0s), 'fake zeros'))
		}
	}
	
	spdata <- droplevels(spdata)


	####################################################
	#Set up data for training and testing to evaluate performance
	####################################################

	#Subset training and testing data by year (use first 80% to predict last 20%)
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

	# test if we have enough presences in testing and training sets (at least one per region as set up for model fitting)
	nprestrain <- table(spdata$regionfact[trainindsp])
	nprestest <- table(spdata$regionfact[testindsp])
	if(any(nprestrain < 1)){
		mywarn <- paste('Zero training presences for', i, sp, 'in', paste(names(nprestrain)[nprestrain<1], collapse=', '))
		allwarnings <- c(allwarnings, mywarn)
		warning(mywarn)
		regstofill <- names(nprestrain)[as.numeric(nprestrain) == 0]
	}
	if(any(nprestest < 1)){
		mywarn <- paste('Zero testing presences for', i, sp, 'in', paste(names(nprestest)[nprestest<1], collapse=', '))
		allwarnings <- c(allwarnings, mywarn)
		warning(mywarn)
	}
	
	# make sure we have at least 6 unique levels for each variable (necessary to fit gam with 4 knots)
	# look at training presence indices, since the most constraining (for mygam2tt)
	levs <- apply(spdata[trainindsp,c('bottemp', 'surftemp', 'logrugosity', 'biomassmean')], 2, FUN=function(x) length(unique(x)))
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
	# note that biomassmean is now linear (not smoothed)
	if(length(levels(spdata$regionfact))==1){
			mypresmod<-formula(presfit ~ s(bottemp,k=6)+s(surftemp,k=6)+s(logrugosity,k=4)+biomassmean)
			myabunmod<-formula(logwtcpue.pad ~ s(bottemp,k=6)+s(surftemp,k=6)+s(logrugosity,k=4)+biomassmean)
			mynullpresmod<-formula(presfit ~ s(logrugosity,k=4)+biomassmean) #Null model w/o temp
			mynullabunmod<-formula(logwtcpue.pad ~ s(logrugosity,k=4)+biomassmean) #Null model w/o temp
	} else {
			mypresmod<-formula(presfit ~ s(bottemp,k=6)+s(surftemp,k=6)+s(logrugosity,k=4)+regionfact+biomassmean-1)
			myabunmod<-formula(logwtcpue.pad ~ s(bottemp,k=6)+s(surftemp,k=6)+s(logrugosity,k=4)+regionfact+biomassmean-1)
			mynullpresmod<-formula(presfit ~ s(logrugosity,k=4)+regionfact+biomassmean-1) #Null model w/o temp
			mynullabunmod<-formula(logwtcpue.pad ~ s(logrugosity,k=4)+regionfact+biomassmean-1) #Null model w/o temp
	}


	####################################	
	# Fit the training/testing models
	####################################	

	# We could use select=TRUE so that terms can be smoothed out of the model (a model selection algorithm),
	if(fittrain){
		try1 <- tryCatch({
			mygam1tt<-gam(mypresmod, family="binomial",data=spdata[traininds,]) 
			mygam2tt<-gam(myabunmod, data=spdata[trainindsp,]) # only fit where species is present
		}, error = function(e) { # ignore warnings, since no function to catch them
			mywarn <- paste('Error in training gam fitting for', i, sp, ':', e)
			assign('allwarnings', c(get('allwarnings', envir=.GlobalEnv), mywarn), envir=.GlobalEnv) # these assigns outside the local scope are poor form in R. But not sure how else to do it here...
			assign('fittrain', FALSE, envir=.GlobalEnv) # if we hit an error in predictions, we can't calculate performance stats
			warning(mywarn)
		})
	}

			

	####################################################
	#Fit models to All data (no test/training split)
	####################################################

	try2 <- tryCatch({
		mygam1<-gam(mypresmod,family="binomial",data=spdata)
		mygam2<-gam(myabunmod,data=spdata[spdata$presfit.pad,]) # only fit where spp is present
		mygam1null<-gam(mynullpresmod,family="binomial",data=spdata)
		mygam2null<-gam(mynullabunmod,data=spdata[spdata$presfit.pad,]) # only fit where spp is present
	
	}, error = function(e) {
		mywarn <- paste('Error in gam fitting for', i, sp, ':', e)
		assign('allwarnings', c(get('allwarnings', envir=.GlobalEnv), mywarn), envir=.GlobalEnv) # these assigns outside the local scope are poor form in R. But not sure how else to do it here...
		warning(mywarn)
	})


	####################################################
	# Plot gam smooths to check for unrealistic out-of-range responses 
	####################################################
	
	#Should write out to PDF
	plot(mygam1,pages=1,scale=0,all.terms=TRUE);mtext(paste(sp,"presence"),outer=T,line=-2)
	plot(mygam2,pages=1,scale=0,all.terms=TRUE);mtext(paste(sp,"abundance"),outer=T,line=-2)


	####################################################
	# Compare predictions to observations to assess model performance
	####################################################

	# For FULL model
	preds1 <- predict(mygam1,spdata,type="response") #can also use mygam1$fitted.values
	preds2 <- exp(predict(mygam2, newdata = spdata, type='response')) # abundance predictions
	smear = mean(exp(mygam2$residuals)) # smearing estimator for re-transformation bias (see Duan 1983, http://www.herc.research.va.gov/include/page.asp?ID=cost-regression)
	preds <- preds1*preds2*smear # adds the bias correction as well
	preds[preds<0] = 0

	# And for training/testing data set
	if(fittrain){
		try3 <- tryCatch({
			preds1tt <- predict(mygam1tt,spdata[testinds,],type="response") 
			preds2tt <- exp(predict(mygam2tt, newdata = spdata[testinds,], type='response'))
			smear = mean(exp(mygam2tt$residuals)) # smearing estimator for re-transformation bias (see Duan 1983, http://www.herc.research.va.gov/include/page.asp?ID=cost-regression)
			predstt <- preds1tt*preds2tt*smear
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
	modeldiag$npres[i] = sum(spdata$presfit)
	if(fittrain){
		modeldiag$npres.tr[i] = sum(spdata$presfit[traininds])
		modeldiag$npres.te[i] = sum(spdata$presfit[testinds])
	}
	modeldiag$ntot[i] = dim(spdata)[1]
	# fill in myregions and myseasons too? would be useful for projections

	# calculate performance using AUC (in part using ROCR)
	preds1.rocr = prediction(predictions=as.numeric(preds1), labels=spdata$presfit)
	modeldiag$auc[i] = performance(preds1.rocr, 'auc')@y.values[[1]] # area under the ROC curve
	if(length(testindsp)>0 & fittrain){ # need presences in the test dataset
		preds1tt.rocr = prediction(predictions=as.numeric(preds1tt), labels=spdata$presfit[testinds])
		modeldiag$auc.tt[i] = performance(preds1tt.rocr, 'auc')@y.values[[1]] #
	}

	# true skill statistic
	a = performance(preds1.rocr, 'tpr')@y.values[[1]] # true positive
	b = performance(preds1.rocr, 'fnr')@y.values[[1]] # false negative
	c = performance(preds1.rocr, 'fpr')@y.values[[1]] # false pos
	d = performance(preds1.rocr, 'tnr')@y.values[[1]] # true neg
	modeldiag$tss[i] = 	max((a*d - b*c)/((a+c)*(b+d)), na.rm=TRUE)
	if(length(testindsp)>0 & fittrain){ # need presences in the test dataset
		preds1tt.rocr = prediction(predictions=as.numeric(preds1tt), labels=spdata$presfit[testinds])
		a = performance(preds1tt.rocr, 'tpr')@y.values[[1]] # true positive
		b = performance(preds1tt.rocr, 'fnr')@y.values[[1]] # false negative
		c = performance(preds1tt.rocr, 'fpr')@y.values[[1]] # false pos
		d = performance(preds1tt.rocr, 'tnr')@y.values[[1]] # true neg
		modeldiag$tss.tt[i] = max((a*d - b*c)/((a+c)*(b+d)), na.rm=TRUE)
	}

	modeldiag$r2.biomass[i] = cor(log(preds2[spdata$presfit]), spdata$logwtcpue[spdata$presfit])^2 # correlation of log(biomass) where present
	if(length(testindsp)>0 & fittrain) modeldiag$r2.biomass.tt[i] = cor(preds2tt[which(testinds %in% testindsp)], spdata$logwtcpue[testindsp])^2 # only if presences exist in the test dataset
	modeldiag$r2.all[i] = cor(preds, spdata$wtcpue)^2 # overall biomass correlation
	if(length(testindsp)>0 & fittrain) modeldiag$r2.all.tt[i] = cor(predstt, spdata$wtcpue[testinds])^2 # overall biomass correlation. only makes sense to do this if the species is present at least once in the testing dataset
	modeldiag$dev.pres[i] = summary(mygam1)$dev.expl
	modeldiag$dev.biomass[i] = summary(mygam2)$dev.expl
	
	#Compare to models without temperature to ultimately calculation %explained by temp terms
	modeldiag$dev.pres.null[i] = summary(mygam1null)$dev.expl
	modeldiag$dev.biomass.null[i] = summary(mygam2null)$dev.expl

	# Some metrics at a spatially aggregated level (1x1deg square) (by year) may be more informative:
	test<-cbind(spdata,preds1,preds)
	t1<-tapply(test$preds1,list(test$year,test$cs1),mean) #average predicted p(occur)
	t2<-tapply(test$presfit,list(test$year,test$cs1),mean) #proportion of hauls with presence
	t3<-tapply(test$preds,list(test$year,test$cs1),mean) #average predicted abundance
	t4<-tapply(test$wtcpue,list(test$year,test$cs1),mean) #average observed abundance

	presr2<-round(cor(stack(as.data.frame(t2))[,1],stack(as.data.frame(t1))[,1],use="p")^2,2)
	abunr2<-round(cor(stack(as.data.frame(t4))[,1],(stack(as.data.frame(t3))[,1]),use="p")^2,2)
	#par(mfrow=c(1,2))
	#plot(stack(as.data.frame(t2))[,1],stack(as.data.frame(t1))[,1],xlab="Proportion of hauls with species present (by 1 deg square)",ylab="Mean predicted probability of occurrence", cex=0.5,main=sp)
	#mtext(paste("r^2 =",presr2))
	#plot(stack(as.data.frame(t4))[,1],(stack(as.data.frame(t3))[,1]),xlab="Average log(WTCPUE) (by 1 deg square)",ylab="Average predicted log(WTCPUE)", cex=0.5,main=sp)
	#mtext(paste("r^2 =",abunr2))

	modeldiag$r2.pres.1deg[i]<-presr2
	modeldiag$r2.abun.1deg[i]<-abunr2

	####################################################
	#### Save models for later projections
	####################################################

	mods = list(mygam1=mygam1, mygam2 = mygam2)
	
	sp <- gsub('/', '', sp) # would mess up saving the file
	
	save(mods, avemeanbiomass, myregions, file=paste(modfolder, 'CEmods_',runname, '_', sp, '.RData', sep='')) # ~4mb file

	#think about figures to output - thermal response curves? spatial prediction by 1 deg square?
	#think about other data to save - number of pres/abs by region (?) 

	# write these files each time through the loop so that we can watch progress
	save(modeldiag,file=paste("output/modeldiag_",runname,".Rdata",sep=""))
	write.csv(modeldiag, file=paste("output/modeldiag_",runname,".csv",sep=""))

	write.csv(allwarnings, file=paste('output/warnings_', runname, '.csv', sep=''))

}
dev.off()

