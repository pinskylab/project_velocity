## Set working directories
if(Sys.info()["nodename"] == "pinsky-macbookair"){
	setwd('~/Documents/Rutgers/Range projections/proj_ranges/')
	projfolder <- '../CEmodels_proj' # holds model projections (outside Git)
	modfolder <- '../CEModels' # holds the models (outside Git)
	}
if(Sys.info()["nodename"] == "amphiprion.deenr.rutgers.edu"){
	setwd('~/Documents/range_projections/')
	projfolder <- 'CEmodels_proj'
	modfolder <- 'CEmodels'
	.libPaths(new='~/R/x86_64-redhat-linux-gnu-library/3.1/') # so that it can find my old packages (chron and ncdf4)
	}
# could add code for Lauren's working directory here

###################
## Flags
###################

## choose which runs to use
## runtype refers to how the Species Distribution Models (SDMs) were fit
## projtype refers to how the SDM projections were done
#runtype <- 'test'; projtype=''
#runtype <- 'testK6noSeas'; projtype=''
#runtype <- 'testK6noSeas'; projtype='_xreg'
runtype <- 'fitallreg'; projtype='_xreg' # for cross-region models and projections

# choose the RCP
rcp <- 85
# rcp <- 45

## choose the time periods
timeperiods <- data.frame(year = 2006:2100, period = c(rep('2006-2020', 15), rep('2021-2040', 20), rep('2041-2060', 20), rep('2061-2080', 20), rep('2081-2100', 20)))
periods <- sort(unique(timeperiods$period))
nt <- length(unique(periods))

# number of climate models
nm <- 13


#############
## Functions
#############
#require(mgcv)
require(Hmisc)

# weighted mean function to use with summarize()
wmean <- function(x){ # values in col 1, weights in col 2
	inds <- !is.na(x[,1]) & !is.na(x[,2])
	return(weighted.mean(x[inds,1], x[inds,2])) 
}



############################################################
## Summarize distributions by two-decade periods
## Keep each model distinct (don't average across models)
############################################################

# list all projections from this run
files <- list.files(path = projfolder, pattern=paste('summproj_', runtype, projtype, '_rcp', rcp, sep=''))


# set up dataframes
# Don't know how many regions for each taxon, so can't pre-allocate the rows (but could guess high... might speed this up)
# could also calc mean depth, but would need to pull from rugosity file
biomassavemapbymod <- data.frame(model=numeric(0), sppocean = character(0), region = character(0), period = character(0), lat = numeric(0), lon = numeric(0), wtcpue.proj = numeric(0))

options(warn=1) # print warnings as they occur

# loop through all files
for(i in 1:length(files)){ # takes a while (a couple hours ?)
	# load data for this species
	load(paste(projfolder, '/', files[i], sep='')) # load summproj for this taxon
	myregions <- sort(unique(summproj$region))
	mysppocean <- gsub('.Rdata', '', gsub(paste('summproj_', runtype, projtype, '_rcp', rcp, '_', sep=''), '', files[i]))

	print(paste(i, 'of', length(files), mysppocean, paste(myregions, collapse=', '), Sys.time()))

	summproj <- merge(summproj, timeperiods)

	# set up dataframe for this taxon
	# need a block of regions/locations for each time period and each model
	inds <- !duplicated(summproj[,c('region', 'lat', 'lon')]) # unique locations to record
	mymap <- data.frame(model=rep(1, sum(inds)*nt), sppocean = rep(mysppocean, sum(inds)*nt), region = rep(summproj$region[inds], nt), period = rep(sort(unique(timeperiods$period)), rep(sum(inds), nt)), lat = rep(summproj$lat[inds], nt), lon = rep(summproj$lon[inds], nt), wtcpue.proj = NA)
		row.names(mymap) <- 1:nrow(mymap)
	mymap <- mymap[rep(1:nrow(mymap),nm),] # expand mymap so there is a row for each model
	mymap$model <- rep(1:nm, rep(sum(inds)*nt, nm)) # assign model IDs to each row

	mymap <- mymap[order(mymap$model, mymap$region, mymap$period, mymap$lat, mymap$lon),]
	summproj <- summproj[order(summproj$region, summproj$year, summproj$lat, summproj$lon),]

	# Summarize projections across each model within time periods
	cols <- grep('wtcpue.proj', names(summproj), value=TRUE)
	for(k in 1:nm){ # for each model
		thiscol <- cols[k]
		for(j in 1:nt){ # for each time period
			inds <- summproj$period == periods[j]
			inds2 <- mymap$period == periods[j] & mymap$model == k
			temp <- aggregate(list(wtcpue.proj = summproj[inds,thiscol]), by=list(region = summproj$region[inds], period = summproj$period[inds], lat = summproj$lat[inds], lon = summproj$lon[inds]), FUN=mean, na.rm=TRUE)
			temp <- temp[order(temp$region, temp$period, temp$lat, temp$lon),]

			if(all(temp$region == mymap$region[inds2] & temp$period == mymap$period[inds2] & temp$lat == mymap$lat[inds2] & temp$lon == mymap$lon[inds2])){
				mymap$wtcpue.proj[inds2] <- temp$wtcpue.proj
			} else {
				warning(paste('rows do not match on j=', j))		
			}
		}
	}
	biomassavemapbymod <- rbind(biomassavemapbymod, mymap) # an inefficient way to do this: better to pre-allocate

}

summary(biomassavemapbymod)
dim(biomassavemapbymod)


### Save the biomass maps by model
save(biomassavemapbymod, file = paste('data/biomassavemapbymod_', runtype, projtype, '_rcp', rcp, '.RData', sep=''))
