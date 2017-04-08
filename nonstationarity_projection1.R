## A script to test whether thermal envelopes change through time (non-stationarity)
## This fits models to first decade of dataset and projects them later to see if the fit is still good
## Only for pres/abs model
## Works with Jim Morley's Feb 2017 data and model format

## Set working directories depending on computer
if(Sys.info()["nodename"] == "pinsky-macbookair"){
	setwd('~/Documents/Rutgers/Range projections/project_velocity/')
	modfolder <- '../CEmodels_nonstationaritytest/'
	nthreads=2
	}
if(Sys.info()["nodename"] == "amphiprion.deenr.rutgers.edu"){
	setwd('~/Documents/range_projections/')
	.libPaths(new='~/R/library/') # so that it can find my old packages
	modfolder <- 'CEmodels_nonstationaritytest/'
	nthreads=10
	}




# Load packages
library(mgcv)
#library(MASS)
#library(ggplot2)
#library(reshape2)
#library(lattice)  
library(ROCR)
library(data.table)
library(plyr)
library(dismo) # used for model evaluation
#require(Hmisc)

# Load data
load('data/master_hauls_March7_2017.RData') # import master hauls file
load('data/dat_selectedspp_Feb_1_2017.Rdata')# load species catch data
#load('data/ProjectionBathGrid_Feb27_2017.RData')# load projection grid (temporary until we receive climate projection data)

# Trim only to NEFSC for the non-stationarity test
dat <- dat[dat$region %in% c('NEFSC_NEUSSpring', 'NEFSC_NEUSFall'),]
hauls <- hauls[hauls$region == 'NEFSC_NEUS',]


# Create survey and region cols
# NOT NEEDED?
#	dat$survey <- dat$region #keeps NEFSC as two separate surveys in "survey"
#	dat$surveyfact <-as.factor(dat$survey)
#	dat$region[dat$region %in% c("NEFSC_NEUSFall","NEFSC_NEUSSpring")] <- "NEFSC_NEUS" #NEFSC as one region in "region"
#	dat$regionfact <- as.factor(dat$region) # the version ot use for model fitting

# Fix rounding-error zeros issue
dat <- dat[!(dat$wtcpue == 0 & dat$region == 'DFO_SoGulf'),] # the zeros in SoGulf are actual zeros (ie not just a scale issue) and thus are true absences. remove here so not confused with rounding-error zeros. true zeros are added back in later.
# NOTE: All chinook salmon ('oncorhynchus tshawytscha_Pac') caught in WC_ANN 2005-2014 (n = 100, out of 1195 total catches in Pacific) were not weighed_probably a tagging study or something....
# Weight could potentially be estimated from #caught? Would need to redo that from script 1.....maybe if a major revision is required (as it is an important species)
dat$wtcpue[dat$wtcpue == 0] <- 0.00001 # 'zeros' in dat are now species too light to register on scales_here a value below the lowest non-zero value is assigned_for transforming data

# Add columns
dat$logwtcpue <- log(dat$wtcpue)

# Trim columns 
# Removed cols are already in master hauls file, which will be merged in below with the hauls data
dat <- data.frame(haulid = dat$haulid, sppocean = dat$sppocean, Freq = dat$Freq, wtcpue = dat$wtcpue, logwtcpue = dat$logwtcpue, presfit = TRUE, stringsAsFactors = F)

# Classify decades
hauls$decade <- floor(hauls$year/10)*10
hauls$decade[hauls$decade==1970] <- 1960 # merge 1960s and 1970s (not much data in 1960s)

# Extract species
allspp = sort(unique(dat$sppocean))

# set up dataframe to hold model diagnostics
decs <- sort(unique(hauls$decade))
n = NA
modeldiag = expand.grid(decade=decs, sppocean=allspp)
modeldiag <- merge(modeldiag, data.frame(ntot=n, npres=n, thresh=n, dev.pres=n, tss=n, acc=n, sens=n, spec=n, kappa=n, auc=n, tssmax=n, accmax=n, kappamax=n, rpb=n, thresh.gr=n, dev.pres.gr=n, tss.gr=n, acc.gr=n, sens.gr=n, spec.gr=n, kappa.gr=n, auc.gr=n, tssmax.gr=n, accmax.gr=n, kappamax.gr=n, rpb.gr=n, stringsAsFactors=FALSE))


######################
# Start the big loop #
######################

#Open pdf to print figures 
pdf(file=paste("figures/GAMs_nonstationarity_projection_1960+1970.pdf",sep=""),width=8,height=6)

options(warn=1) # print warnings as they occur
allwarnings = NULL
print(paste(length(allspp), 'models to fit'))

for(i in 1:length(allspp)){ 
	fittrain = TRUE
	mygam1 <- NULL

	sp<-allspp[i]
	print(paste(i,sp, Sys.time()))

	mydat<-dat[dat$sppocean==sp,] 

#	mydat$logwtcpue[is.infinite(mydat$logwtcpue)] <- NA
#	myregions<-unique(mydat$region)
#	myhauls<-unique(mydat$haulid)
	ocean <- strsplit(sp, split='_')[[1]][2]




	####################################################
	# Add records for when there was a haul but no fish
	####################################################
	# OLD APPROACH FROM LAUREN
	# expand to all regions in the same ocean.
#	zs<-dat[!(dat$haulid %in% myhauls) & dat$ocean %in% myocean,] #extract haulids where this species is missing
#
#	matchpos<-match(unique(zs$haulid),zs$haulid) # Extract one record of each haulid
#	zeros<-zs[matchpos,]
#	#Add/change relevant columns --> zero catch for target spp.
#	zeros$spp<-mydat$spp[1]
#	zeros$sppl<-mydat$sppl[1]
#	zeros$sppnew<-mydat$sppnew[1]
#	zeros$sppocean<-mydat$sppocean[1] #may need to add "ocean"
#	zeros$wtcpue<-0
#	zeros$logwtcpue <- NA
#	zeros$presfit<-FALSE
#
#	mydatw0<-rbind(mydat,zeros) #combine positive hauls with zero hauls

	# NEW APPROACH FROM JIM
	haulsMod <- hauls[hauls$ocean == ocean,] # trim master hauls file to the ocean of interest 
	haulsMod <- merge(haulsMod, mydat, by='haulid', all.x = T, sort=F)   # Add empty hauls from the relavant ocean
	haulsMod$presfit[is.na(haulsMod$presfit)] <- FALSE
	spdata <- droplevels(haulsMod) # drop the west coast 'regionfact' levels

	#Order by time
	spdata<-spdata[order(spdata$year,spdata$month),]


	##############################################################
	# Check that we have enough data for fitting GAMs
	##############################################################

	# check the decades
	if(sum(is.na(spdata$decade))>0){ # should assign every row to a decade
		mywarn <- paste('Some decades undefined for', i, sp)
		allwarnings <- c(allwarnings, mywarn)
		warning(mywarn)
	}
	

	# make sure at least 10 presences in fitting decade
	ninds <- table(spdata$presfit, spdata$decade)
	if(!('TRUE' %in% rownames(ninds))){ # means no presences ever!
		mywarn <- paste('Not enough presences at all for', i, sp)
		allwarnings <- c(allwarnings, mywarn)
		warning(mywarn)
		fittrain = FALSE
	} else {
		if(ninds['TRUE','1960'] < 10){ # if some presences, makes sure enough for fitting in 1960
			mywarn <- paste('Not >=10 presences in 1960s for', i, sp)
			allwarnings <- c(allwarnings, mywarn)
			warning(mywarn)
			fittrain = FALSE		
		}
		if(any(ninds['TRUE',] < 4)){ # if some presences, makes sure enough in any decade for evaluation
			mywarn <- paste('Not >=4 presence in all decades for', i, sp)
			allwarnings <- c(allwarnings, mywarn)
			warning(mywarn)
			fittrain = FALSE		
		}
	}
	
	# make sure we have at least 6 unique levels for each variable (necessary to fit gam with 4 knots) in the fitting decade
	levs1 <- apply(spdata[spdata$decade==1960,c('SBT.seasonal', 'SST.seasonal.mean', 'SBT.min', 'SBT.max', 'SST.min', 'SBT.max', 'rugosity')], 2, FUN=function(x) length(unique(x)))
	if(any(levs1 < 6)){
		mywarn <- paste("Not enough (>=6) unique levels in first decade presence set for", i, sp, ". Won't fit models")
		allwarnings <- c(allwarnings, mywarn)
		warning(mywarn)
		fittrain = FALSE
	}	
	# table(spdata$year[spdata$presfit]) # output number of presences by year
	
	####################################################
	# Set up model formula
	####################################################

	#Default models. Leave out region factor since only NEFSC
	mypresmod<-formula(presfit ~ s(SBT.seasonal) + s(SST.seasonal.mean) + s(SBT.min) + s(SBT.max) + s(SST.max) + s(rugosity) + s(GRAINSIZE)) # no regionfact in these models

	# set 'gamma' penalty levels for gam to prevent overfitting_got this from a presentation by Simon Wood (https://people.maths.bris.ac.uk/~sw15190/mgcv/tampere/mgcv-advanced.pdf)
	gammaPA <- log(nrow(spdata)) / 2


	####################################################
	#Fit models to first decade
	####################################################

	if(fittrain){
		try2 <- tryCatch({
			mygam1<-gam(mypresmod, family="binomial", select=TRUE, gamma=gammaPA, data=spdata[spdata$decade==1960,], control=list(maxit=500, nthreads=nthreads)) # add more iterations to help convergence, plus use more threads
		}, error = function(e) {
			mywarn <- paste('Error in gam fitting for', i, sp, ':', e)
			assign('allwarnings', c(get('allwarnings', envir=.GlobalEnv), mywarn), envir=.GlobalEnv) # this assigns outside the local scope are poor form in R. But not sure how else to do it here...
			warning(mywarn)
		})


		####################################################
		# Plot gam smooths to check visually
		####################################################
	
		#Should write out to PDF opened before the loop
		plot(mygam1,pages=1,scale=0,all.terms=TRUE);mtext(paste(sp,"presence"),outer=T,line=-2)


		#################
		# Model evaluation
		#################
		k <- which(modeldiag$sppocean==sp)

		# fill in basic model data
		for(j in 1:length(decs)){
			k2 <- intersect(k, which(modeldiag$decade==decs[j]))
			modeldiag$ntot[k2] = sum(spdata$decade==decs[j])
			modeldiag$npres[k2] <- sum(spdata$presfit[spdata$decade==decs[j]])
		}
		
		# make predictions for each observation point
		spdata <- as.data.table(spdata)
		spdata[,predpres := predict(mygam1, newdata=spdata, type='response')]

		# evaluate each decade and set up threshold for pres/abs calls
		e <- dlply(.data=spdata, .fun=function(x) evaluate(p=as.vector(x$predpres[x$presfit]), a=as.vector(x$predpres[!x$presfit])), .variables=~decade) # evaluate each decade separately
		modeldiag$thresh[k] <- threshold(e[[1]], stat='prevalence') # set threshold from fitting decade
		e.ind <- sapply(e, function(x) which.min(abs(x@t - modeldiag$thresh[k][1]))) # index for the chosen threshold, for each decade
		conf <- lapply(e, function(x) as.data.frame(x@confusion)) # confusion matrices (all thresholds)
		conftrim <- conf
		for(j in 1:length(conftrim)) conftrim[[j]] <- conf[[j]][e.ind[j],] # trim conf to only the matrix for the chosen threshold, in each decade
		conftrim <- as.data.frame(do.call(rbind, conftrim)) # collapse to a df
		conftrim$decade <- rownames(conftrim)

		# true skill statistic, accuracy, kappa, and other stats that require a threshold
		for(j in 1:length(decs)){
			k2 <- intersect(k, which(modeldiag$decade==decs[j]))
			myconf <- conftrim[as.character(decs[j]),]
			modeldiag$tss[k2] = with(myconf, (tp*tn - fn*fp)/((tp+fp)*(fn+tn))) # TSS for chosen threshold
			modeldiag$acc[k2] = with(myconf, (tp+tn)/(tp+fp+fn+tn)) # overall accuracy
			modeldiag$sens[k2] = with(myconf, (tp)/(tp+fn)) # sensitivity: fraction of correctly predicted presences
			modeldiag$spec[k2] = with(myconf, (tn)/(tn+fp)) # specificity: fraction of correctly predicted absences
			modeldiag$kappa[k2] = e[[j]]@kappa[e.ind[j]] # Cohen's kappa
		}

		# pres/abs model diagnostics (no threshold needed)
		modeldiag$dev.pres[k] = summary(mygam1)$dev.expl

		for(j in 1:length(decs)){
			k2 <- intersect(k, which(modeldiag$decade==decs[j]))
			modeldiag$auc[k2] <- e[[j]]@auc
			modeldiag$tssmax[k2] <- max(with(conf[[j]], (tp*tn - fn*fp)/((tp+fp)*(fn+tn))), na.rm=TRUE) # maximum TSS (any threshold)
			modeldiag$accmax[k2] <- max(with(conf[[j]], (tp+tn)/(tp+fp+fn+tn)), na.rm=TRUE) # maximum overall accuracy
			modeldiag$kappamax[k2] <- max(e[[j]]@kappa, na.rm=TRUE) # maximum kappa
			modeldiag$rpb[k2] <- e[[j]]@cor # point biserial correlation
		}
		
		# test against future decades by grid cell
		gridsz <- 1/20
		spdata[, latgrid := ceiling(lat/gridsz)*gridsz + gridsz/2]
		spdata[, longrid := ceiling(lon/gridsz)*gridsz + gridsz/2]
		spdatbydec <- spdata[,.(presfit=mean(presfit), n=length(presfit), predpres=mean(predpres)), by=.(latgrid, longrid, decade)]

		# evaluate each decade and set up threshold for pres/abs calls (by grid cell)
		e2 <- dlply(.data=spdatbydec, .fun=function(x) evaluate(p=as.vector(x$predpres[round(x$presfit)==1]), a=as.vector(x$predpres[round(x$presfit)==0])), .variables=~decade) # evaluate each decade separately
		modeldiag$thresh.gr[k] <- threshold(e2[[1]], stat='prevalence') # set threshold from fitting decade
		e2.ind <- sapply(e2, function(x) which.min(abs(x@t - modeldiag$thresh.gr[k][1]))) # index for the chosen threshold, for each decade
		conf2 <- lapply(e2, function(x) as.data.frame(x@confusion)) # confusion matrices (all thresholds)
		conftrim2 <- conf2
		for(j in 1:length(conftrim2)) conftrim2[[j]] <- conf2[[j]][e2.ind[j],] # trim conf to only the matrix for the chosen threshold, in each decade
		conftrim2 <- as.data.frame(do.call(rbind, conftrim2)) # collapse to a df
		conftrim2$decade <- rownames(conftrim2)

		# true skill statistic, accuracy, kappa, and other stats that require a threshold (by grid cell)
		for(j in 1:length(decs)){
			k2 <- intersect(k, which(modeldiag$decade==decs[j]))
			myconf <- conftrim2[as.character(decs[j]),]
			modeldiag$tss.gr[k2] = with(myconf, (tp*tn - fn*fp)/((tp+fp)*(fn+tn))) # TSS for chosen threshold
			modeldiag$acc.gr[k2] = with(myconf, (tp+tn)/(tp+fp+fn+tn)) # overall accuracy
			modeldiag$sens.gr[k2] = with(myconf, (tp)/(tp+fn)) # sensitivity: fraction of correctly predicted presences
			modeldiag$spec.gr[k2] = with(myconf, (tn)/(tn+fp)) # specificity: fraction of correctly predicted absences
			modeldiag$kappa.gr[k2] = e2[[j]]@kappa[e2.ind[j]] # Cohen's kappa
		}

		# pres/abs model diagnostics (no threshold needed) (by grid cell)
		modeldiag$dev.pres.gr[k] = summary(mygam1)$dev.expl

		for(j in 1:length(decs)){
			k2 <- intersect(k, which(modeldiag$decade==decs[j]))
			modeldiag$auc.gr[k2] <- e2[[j]]@auc
			modeldiag$tssmax.gr[k2] <- max(with(conf2[[j]], (tp*tn - fn*fp)/((tp+fp)*(fn+tn))), na.rm=TRUE) # maximum TSS (any threshold)
			modeldiag$accmax.gr[k2] <- max(with(conf2[[j]], (tp+tn)/(tp+fp+fn+tn)), na.rm=TRUE) # maximum overall accuracy
			modeldiag$kappamax.gr[k2] <- max(e2[[j]]@kappa, na.rm=TRUE) # maximum kappa
			modeldiag$rpb.gr[k2] <- e2[[j]]@cor # point biserial correlation
		}


		####################################################
		#### Save models for later projections
		####################################################

		mods = list(mygam1=mygam1)
	
		sp <- gsub('/', '', sp) # would mess up saving the file
	
		save(mods, file=paste(modfolder, 'CEmods_nonstationarity_projection2_1960+1970', sp, '.RData', sep='')) # ~4mb file

		# write these files each time through the loop so that we can watch progress
		write.csv(modeldiag, file="output/modeldiag_nonstationarity_projection2_1960+1970.csv", row.names=FALSE)

		write.csv(allwarnings, file='output/warnings_nonstationarity_projection2_1960+1970.csv')
	}
}


dev.off() # close the figure

