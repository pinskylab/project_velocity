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


#######################
## Script-wide flags ##
#######################

## choose which run and time periods to use
#runtype <- 'test'
#runtype <- 'testK6noSeas'; projtype <- ''
#runtype <- 'testK6noSeas'; projtype <- '_xreg' # with cross-region projections
runtype <- 'fitallreg'; projtype <- '_xreg' # with cross-region projections

# choose the RCP
rcp <- 85
# rcp <- 45


###########################
## Script-wide functions ##
###########################
#require(mgcv)
require(Hmisc)
require(lattice)
require(gridExtra)

# weighted mean function for summarize()
wmean <- function(x){
	i <- !is.na(x[,1]) & !is.na(x[,2])
	if(sum(i) == 0) return(NA)
	else return(weighted.mean(x[i,1], x[i,2]))
}


##########################################################
## Calculate mean lat/depth and sum of biomass by year  ##
## within each region                                   ##
## and across regions (if xreg)                         ##
##########################################################

# list all projections from this run
files <- list.files(path = projfolder, pattern=paste('summproj_', runtype, projtype, '_rcp', rcp, sep=''))

# list to save predicted positions for lat, lon, and depth
meanpos <- list(0)

# summarize within regions (whether xreg or not)
# set up dataframes
# Don't know how many regions for each taxon, so can't pre-allocate the rows (but could guess high... might speed this up)
# could also calc mean depth, but would need to pull from rugosity file
biomasssum <- data.frame(sppocean = character(0), region = character(0), year = numeric(0)) # sum of average wtcpue across the region (2006-2100) for each survey in each model (columns)
	for(i in 1:13) biomasssum[[paste('summwtcpue', i, sep='_')]] <- numeric(0)
meanlat <- data.frame(sppocean = character(0), region = character(0), year = numeric(0)) # biomass-weighted mean latitude across the region (2006-2100) for each survey in each model (columns)
	for(i in 1:13) meanlat[[paste('lat', i, sep='_')]] <- numeric(0)
meanlon <- data.frame(sppocean = character(0), region = character(0), year = numeric(0))
	for(i in 1:13) meanlon[[paste('lon', i, sep='_')]] <- numeric(0)


options(warn=1) # print warnings as they occur

# loop through all files
for(i in 1:length(files)){ # takes a while (a couple hours ?)
	# load data for this species
	load(paste(projfolder, '/', files[i], sep='')) # load summproj for this taxon
	myregions <- sort(unique(summproj$region))
	mysppocean <- gsub('.Rdata', '', gsub(paste('summproj_', runtype, projtype, '_rcp', rcp, '_', sep=''), '', files[i])) # strip sppocean out of file name

	print(paste(i, 'of', length(files), mysppocean, paste(myregions, collapse=', '), Sys.time()))

	# set up dataframes for this taxon
	mybiomasssum <- data.frame(sppocean = mysppocean, region = rep(myregions, rep(95,length(myregions))), year = rep(2006:2100, length(myregions)))
	mymeanlat <- data.frame(sppocean = mysppocean, region = rep(myregions, rep(95,length(myregions))), year = rep(2006:2100, length(myregions)))
	mymeanlon <- data.frame(sppocean = mysppocean, region = rep(myregions, rep(95,length(myregions))), year = rep(2006:2100, length(myregions)))

	# Summarize projections
	for(j in 1:13){
		snm <- paste('wtcpue.proj', j, sep='_')
		temp <- aggregate(summproj[,snm], by=list(year = summproj$year, region = summproj$region), FUN=sum, na.rm=TRUE) # remove NAs: assumes that NA gridcells are constant through time, which they should be, since determined by the GCM grid
			names(temp)[3] <- paste('summwtcpue', j, sep='_')
		mybiomasssum <- merge(mybiomasssum, temp)

		temp <- summarize(summproj[,c('lat', snm)], by=list(region = summproj$region, year = summproj$year), FUN=wmean)
			names(temp)[3] <- paste('lat', j, sep='_')
		mymeanlat <- merge(mymeanlat, temp)

		temp <- summarize(summproj[,c('lon', snm)], by=list(region = summproj$region, year = summproj$year), FUN=wmean)
			names(temp)[3] <- paste('lon', j, sep='_')
		mymeanlon <- merge(mymeanlon, temp)

	}

	biomasssum <- rbind(biomasssum, mybiomasssum) # an inefficient way to do this: better to pre-allocate
	meanlat <- rbind(meanlat, mymeanlat)
	meanlon <- rbind(meanlon, mymeanlon)

}

### Save meanlat, meanlon and biomasssum
save(biomasssum, meanlat, meanlon, file = paste('data/meanlat,lon,biomass_', runtype, projtype, '_rcp', rcp, '.RData', sep=''))


# if projections across regions, summarize by spp overall as well
if(projtype == '_xreg'){
	# set up dataframes (can preallocate since only one ocean per species)
	n <- rep(NA, length(files) * 80)
	biomasssumbyspp <- data.frame(sppocean = n, year = n, stringsAsFactors=FALSE) # sum of average wtcpue across the region (2020-2099) for each survey in each model (columns)
		for(i in 1:13) biomasssumbyspp[[paste('summwtcpue', i, sep='_')]] <- n
	meanlatbyspp <- data.frame(sppocean = n, year = n, stringsAsFactors=FALSE) # biomass-weighted mean latitude across the region (2020-2099) for each survey in each model (columns)
		for(i in 1:13) meanlatbyspp[[paste('lat', i, sep='_')]] <- n
	meanlonbyspp <- data.frame(sppocean = n, year = n, stringsAsFactors=FALSE)
		for(i in 1:13) meanlonbyspp[[paste('lon', i, sep='_')]] <- n


	options(warn=1) # print warnings as they occur

	# loop through all files
	for(i in 1:length(files)){ # takes a while (a couple hours ?)
		# load data for this species
		load(paste(projfolder, '/', files[i], sep='')) # load summproj for this taxon
		mysppocean <- gsub('.Rdata', '', gsub(paste('summproj_', runtype, projtype, '_rcp', rcp, '_', sep=''), '', files[i])) # strip sppocean out of file name

		print(paste(i, 'of', length(files), mysppocean, Sys.time()))

		# set up dataframes for this taxon
		mybiomasssumbyspp <- data.frame(sppocean = mysppocean, year = 2020:2099, stringsAsFactors=FALSE)
		mymeanlatbyspp <- data.frame(sppocean = mysppocean, year = 2020:2099, stringsAsFactors=FALSE)
		mymeanlonbyspp <- data.frame(sppocean = mysppocean, year = 2020:2099, stringsAsFactors=FALSE)

		# Summarize projections (loop through each model)
		for(j in 1:13){
			snm <- paste('wtcpue.proj', j, sep='_')
			temp <- aggregate(summproj[,snm], by=list(year = summproj$year), FUN=sum, na.rm=TRUE) # remove NAs: assumes that NA gridcells are constant through time, which they should be, since determined by the GCM grid
				names(temp)[2] <- paste('summwtcpue', j, sep='_')
			mybiomasssumbyspp <- merge(mybiomasssumbyspp, temp)

			temp <- summarize(summproj[,c('lat', snm)], by=list(year = summproj$year), FUN=wmean)
				names(temp)[2] <- paste('lat', j, sep='_')
			mymeanlatbyspp <- merge(mymeanlatbyspp, temp)

			temp <- summarize(summproj[,c('lon', snm)], by=list(year = summproj$year), FUN=wmean)
				names(temp)[2] <- paste('lon', j, sep='_')
			mymeanlonbyspp <- merge(mymeanlonbyspp, temp)

		}

		mybiomasssumbyspp <- mybiomasssumbyspp[,names(biomasssumbyspp)]
		mymeanlatbyspp <- mymeanlatbyspp[,names(meanlatbyspp)]
		mymeanlonbyspp <- mymeanlonbyspp[,names(meanlonbyspp)]

		startind <- (i-1)*80+1
		biomasssumbyspp[startind:(startind+79),] <- mybiomasssumbyspp
		meanlatbyspp[startind:(startind+79),] <- mymeanlatbyspp
		meanlonbyspp[startind:(startind+79),] <- mymeanlonbyspp

	}
}

### Save meanlatbyspp, meanlonbyspp and biomasssumbyspp
save(biomasssumbyspp, meanlatbyspp, meanlonbyspp, file = paste('data/meanlat,lon,biomassbyspp_', runtype, projtype, '_rcp', rcp, '.RData', sep=''))



###################################################
## Summarize distributions by two-decade periods
###################################################

timeperiods <- data.frame(year = 2006:2100, period = c(rep('2006-2020', 15), rep('2021-2040', 20), rep('2041-2060', 20), rep('2061-2080', 20), rep('2081-2100', 20)))
periods <- sort(unique(timeperiods$period))
nt <- length(unique(periods))

# list all projections from this run
files <- list.files(path = projfolder, pattern=paste('summproj_', runtype, projtype, '_rcp', rcp, sep=''))


# set up dataframes
# Don't know how many regions for each taxon, so can't pre-allocate the rows (but could guess high... might speed this up)
# could also calc mean depth, but would need to pull from rugosity file
biomassavemap <- data.frame(sppocean = character(0), region = character(0), period = character(0), lat = numeric(0), lon = numeric(0), wtcpue.proj = numeric(0))

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
	inds <- !duplicated(summproj[,c('region', 'lat', 'lon')]) # unique locations to record
	mymap <- data.frame(sppocean = rep(mysppocean, sum(inds)*nt), region = rep(summproj$region[inds], nt), period = rep(sort(unique(timeperiods$period)), rep(sum(inds), nt)), lat = rep(summproj$lat[inds], nt), lon = rep(summproj$lon[inds], nt), wtcpue.proj = NA)

	mymap <- mymap[order(mymap$region, mymap$period, mymap$lat, mymap$lon),]
	summproj <- summproj[order(summproj$region, summproj$year, summproj$lat, summproj$lon),]

	# Summarize projections across all models within time periods
	cols <- grep('wtcpue.proj', names(summproj))
	for(j in 1:nt){
		inds <- summproj$period == periods[j]
		inds2 <- mymap$period == periods[j]
		temp <- apply(summproj[inds,cols], MARGIN=1, FUN=mean) # average across models
		temp2 <- aggregate(list(wtcpue.proj = temp), by=list(region = summproj$region[inds], period = summproj$period[inds], lat = summproj$lat[inds], lon = summproj$lon[inds]), FUN=mean, na.rm=TRUE) # average within time periods
		temp2 <- temp2[order(temp2$region, temp2$period, temp2$lat, temp2$lon),]

		if(all(temp2$region == mymap$region & temp2$period == mymap$period & temp2$lat == mymap$lat & temp2$lon == mymap$long)){
			mymap$wtcpue.proj[inds2] <- temp2$wtcpue.proj
		} else {
			warning(paste('rows do not match on j=', j))		
		}
	}

	biomassavemap <- rbind(biomassavemap, mymap) # an inefficient way to do this: better to pre-allocate

}

summary(biomassavemap)
dim(biomassavemap)


### Save the file
save(biomassavemap, file = paste('data/biomassavemap_', runtype, projtype, '_rcp', rcp, '.RData', sep=''))

