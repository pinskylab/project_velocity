## Set working directories
if(Sys.info()["nodename"] == "pinsky-macbookair"){
	setwd('~/Documents/Rutgers/Range projections/proj_ranges/')
	}
if(Sys.info()["nodename"] == "amphiprion.deenr.rutgers.edu"){
	setwd('~/Documents/range_projections/')
	.libPaths(new='~/R/x86_64-redhat-linux-gnu-library/3.1/') # so that it can find my old packages
	ncores=15
	}
# could add code for Lauren's working directory here

############
## Flags
############

## choose which runs to use
## runtype refers to how the Species Distribution Models (SDMs) were fit
## projtype refers to how the SDM projections were done
#runtype <- 'test'; projtype=''
#runtype <- ''; projtype=''`
#runtype <- 'testK6noSeas'; projtype='_xreg'
runtype <- 'fitallreg'; projtype='_xreg'

# choose the rcp
rcp <- 85
# rcp <- 45


####################
## helper functions
####################
require(Hmisc) # for summarize()
require(RColorBrewer)
require(parallel) # for multi-core calculations

lu <- function(x) return(length(unique(x)))

# expects x to have columns 'period' and 'rich' (in that order)
# each as used below (period has a dash, and we use the midpoint here)
calcrichtrend <- function(x){
	pds <- as.numeric(unlist(strsplit(as.character(x[,1]), split='-')))
	dim(pds) <- c(2,nrow(x))
	mids <- colMeans(pds)
	mod <- lm(x[,2] ~ mids)
	return(mod$coefficients[2])
}

# normalize to 0-1
norm01 <- function(x){
	mn <- min(x)
	mx <- max(x)
	return((x-mn)/(mx-mn))
}

# turnover metrics
# for use with Hmisc::summarize() but doesn't work for some reason
# expects region in column 1, lat (2), lon (3), time period in column 4, spp in col 5, pres in col 6
# only present species included
turnover <- function(x){
	pers <- sort(unique(x[,4]))
	if(length(pers) != 2) stop(paste('expecting 2 time periods:', unique(x[,1]), unique(x[,2]), unique(x[,3]), paste(unique(x[,4], collapse=','))))
	initcomm <- x[x[,4] == pers[1] & x[,6]==TRUE,5]
	finalcomm <- x[x[,4] == pers[2] & x[,6]==TRUE,5]
	
	a <- length(intersect(initcomm, finalcomm)) # number of shared species
	lost <- length(setdiff(initcomm, finalcomm)) # number of species lost
	gained <- length(setdiff(finalcomm, initcomm))
	return(c(nstart=length(initcomm), nend=length(finalcomm), nlost=lost, ngained=gained, flost=lost/length(initcomm), fgained=gained/length(initcomm), beta_sor=2*a/(2*a+gained+lost)))
}

findthresh <- function(counter, model, sppocean, region, bm, ncount){
	print(paste(counter, 'of', ncount, Sys.time()))
	wts <- sort(bm$wtcpue.proj[bm$model == model & bm$sppocean == sppocean & bm$region == region])
	ind<-max(which(cumsum(wts)/sum(wts)<0.05))
	return(wts[ind])
}

# version that relies on bmave2020 being in the environment
findthresh2 <- function(counter, model, sppocean, region, ncount){ 
	print(paste(counter, 'of', ncount, Sys.time()))
	wts <- sort(bmave2020$wtcpue.proj[bmave2020$model == model & bmave2020$sppocean == sppocean & bmave2020$region == region])
	ind<-max(which(cumsum(wts)/sum(wts)<0.05))
	return(wts[ind])
}
	

############################################
## Summarize community change by grid cell
## for each model separately
############################################

load(paste('data/biomassavemapbymod_', runtype, projtype, '_rcp', rcp, '.RData', sep='')) # loads biomassavemapbymod data.frame. SLOW (a 1GB file).
	dim(biomassavemapbymod) # 101,196,030 x 7
#	summary(biomassavemapbymod) # some NAs
#		table(biomassavemapbymod$region[is.na(biomassavemapbymod$wtcpue.proj)]) # WCTri, ScotianShelf, NEUSFall, NEUSSpring, GOMex
#		table(biomassavemapbymod$period[is.na(biomassavemapbymod$wtcpue.proj)]) # evenly across all periods
#		table(biomassavemapbymod$sppocean[is.na(biomassavemapbymod$wtcpue.proj)], biomassavemapbymod$region[is.na(biomassavemapbymod$wtcpue.proj)]) # same number missing for all species in a region: a purely spatial issue

	# trim out NAs
	biomassavemapbymod <- biomassavemapbymod[!is.na(biomassavemapbymod$wtcpue.proj),]
	dim(biomassavemapbymod) # 94,935,750 x 7

# find abundance threshold to count as present for each taxon
# use cumulative 5% of wtcpue from earliest timeperiod, by region
# this would likely be faster with data.table or dplyr, but I'm having trouble installing either package on Amphiprion
i <- !duplicated(biomassavemapbymod[,c('model', 'region', 'sppocean')])
taxthreshbymod <- biomassavemapbymod[i,c('model', 'region', 'sppocean')]
	taxthreshbymod$thresh<-NA

#taxthreshbymod$thresh <- mcmapply(findthresh, counter=1:nrow(taxthreshbymod), model=taxthreshbymod$model, sppocean=taxthreshbymod$sppocean, region=taxthreshbymod$region, MoreArgs=list(bm=biomassavemapbymod[biomassavemapbymod$period == '2006-2020',], ncount=nrow(taxthreshbymod)), SIMPLIFY = TRUE, mc.cores=ncores) # initiating each run is slow: perhaps because it creates many copies of biomassavemapbymod? each instance uses 17GB memory, yikes. takes about 5 hours?

bmave2020 <- biomassavemapbymod[biomassavemapbymod$period == '2006-2020',] # smaller version with just the first period
	nrow(bmave2020) # 18,987,150
	
i <- 1:nrow(taxthreshbymod) # a simple way to run just a part of the data as a test
taxthreshbymod$thresh[i] <- mcmapply(findthresh2, counter=(1:nrow(taxthreshbymod))[i], model=taxthreshbymod$model[i], sppocean=taxthreshbymod$sppocean[i], region=taxthreshbymod$region[i], MoreArgs=list(ncount=nrow(taxthreshbymod[i,])), SIMPLIFY = TRUE, mc.cores=ncores) # each instance uses 17GB memory, yikes, but some fo this is shared. takes about 6 hours?


# save threshold
save(taxthreshbymod, file=paste('data/taxthreshbymod_', runtype, projtype, '_rcp', rcp, '.RData', sep=''))


# apply threshold to projections to determine presence/absence
presmapbymod <- merge(biomassavemapbymod, taxthreshbymod) # surprisingly fast (10 minutes?)
presmapbymod$pres <- presmapbymod$wtcpue.proj > presmapbymod$thresh
	
	# examine the results
	table(presmapbymod$sppocean, presmapbymod$pres) # how many marked present vs. absent

#	i = presmapbymod$sppocean == 'gadus morhua_Atl' & presmapbymod$period == '2006-2020'
#	i = presmapbymod$sppocean == 'gadus morhua_Atl' & presmapbymod$period == '2081-2100'
#	i = presmapbymod$sppocean == 'theragra chalcogramma_Pac' & presmapbymod$period == '2006-2020'
#	i = presmapbymod$sppocean == 'theragra chalcogramma_Pac' & presmapbymod$period == '2081-2100'
#	plot(presmapbymod$lon[i], presmapbymod$lat[i], col=c('black', 'red')[presmapbymod$pres[i]+1])

# write out pres/abs
save(presmapbymod, file=paste('data/presmapbymod_', runtype, projtype, '_rcp', rcp, '.RData', sep=''))
	
# richness by grid cell by time period
i <- presmapbymod$pres # only count where a spp is present
richbymod <- aggregate(list(rich = presmapbymod$sppocean[i]), by=list(model=presmapbymod$model[i], region=presmapbymod$region[i], period=presmapbymod$period[i], lat=presmapbymod$lat[i], lon=presmapbymod$lon[i]), FUN=lu) # calculate richness (# taxa) by grid cell in each time period for each model

	# examine
	nrich <- aggregate(list(nrich = richbymod$rich), by=list(model=richbymod$model, region=richbymod$region, lat=richbymod$lat, lon=richbymod$lon), FUN=lu)
		summary(nrich) # up to five values: good (most 4-5)
	
# write out richness
save(richbymod, file=paste('data/richbymod_', runtype, projtype, '_rcp', rcp, '.RData', sep=''))
	
# load in richness calcs
#load(paste('data/richbymod_', runtype, projtype, '_rcp', rcp, '.RData', sep='')) # loads rich data.frame

	
# calculate trend in richness by grid cell
richtrendbymod <- Hmisc::summarize(X=richbymod[,c('period', 'rich')], by=list(model=richbymod$model, region=richbymod$region, lat=richbymod$lat, lon=richbymod$lon, dummy=richbymod$lon), FUN=calcrichtrend, stat.name='trend') # calculate richness (# taxa) by grid cell in each time period for each model. for an unknown resion, the last column in the by list gets NAs, so padded with a dummy column here
	richtrendbymod <- richtrendbymod[,-grep('dummy', names(richtrendbymod))]
	
	# examine
	
# write out richness trends
save(richtrendbymod, file=paste('data/richtrendbymod_', runtype, projtype, '_rcp', rcp, '.RData', sep=''))
