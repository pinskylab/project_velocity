## Set working directories
if(Sys.info()["nodename"] == "pinsky-macbookair"){
	setwd('~/Documents/Rutgers/Range projections/proj_ranges/')
	projfolder <- '../CEmodels_proj' # holds model projections (outside Git)
	modfolder <- '../CEModels' # holds the models (outside Git)
	climgridfolder <- '../data/'
	}
if(Sys.info()["nodename"] == "amphiprion.deenr.rutgers.edu"){
	setwd('~/Documents/range_projections/')
	projfolder <- 'CEmodels_proj'
	modfolder <- 'CEmodels'
	climgridfolder <- 'data/'
	}
# could add code for Lauren's working directory here


##############################
## Process NatCap Wind data ##
##############################
require(rgdal)

# read in data
wind_AKw <- readGDAL('cmsp_data/Wind/output/npv_US_millions_Alaska_west1050.tif') # read in SpatialGridDataFrame
	# image(wind_AKw)
wind_AKe <- readGDAL('./cmsp_data/Wind/output/npv_US_millions_Alaska_east1050.tif') # read in SpatialGridDataFrame
wind_WC <- readGDAL('./cmsp_data/Wind/output/npv_US_millions_WestCoast1050.tif')
wind_GM <- readGDAL('./cmsp_data/Wind/output/npv_US_millions_GoMex1050.tif')
wind_NE <- readGDAL('./cmsp_data/Wind/output/npv_US_millions_Northeast1050.tif')
	# image(wind_NE)

load(paste(climgridfolder, 'climGrid_rcp85.proj2.RData', sep='')) # loads clim. has the grid cells we want to summarize to

# combine wind output into a list
wind <- list(wind_AKw, wind_AKe, wind_WC, wind_GM, wind_NE)

# project to LL. also converts to SpatialPointsDataFrame
wind.t <- wind
for(i in 1:length(wind)){
	wind.t[[i]] <- spTransform(wind[[i]], CRS('+proj=longlat +data=WGS84')) # produces warnings about coercing to points. this is ok
#		colfun <- colorRamp(c('white', 'blue'))
#		sc <- (wind.t[[i]]$band1- min(wind.t[[i]]$band1))/(max(wind.t[[i]]$band1) - min(wind.t[[i]]$band1))
#		plot(wind.t[[i]], col=rgb(colfun(sc), maxColorValue=255), pch=16, cex=0.05) # plots
#		plot(wind.t[[i]], col=c('blue', 'red')[1+(wind.t[[i]]$band1> -600)], pch=16, cex=0.05) # plots of <> a threshold
}
	
# Add a grid indicator
gridsize=0.25 # size of grid of the climate data, in degrees
for(i in 1:length(wind.t)){
	wind.t[[i]]$lat <- floor(coordinates(wind.t[[i]])[,2]/gridsize)*gridsize + gridsize/2 # round to nearest grid center
	wind.t[[i]]$lon <- floor(coordinates(wind.t[[i]])[,1]/gridsize)*gridsize + gridsize/2
}

# Summarize by grid cell
wind.sum <- wind.t
for(i in 1:length(wind.t)){
	wind.sum[[i]] <- aggregate(list(wind_npv = wind.t[[i]]$band1), by=list(lat = wind.t[[i]]$lat, lon = wind.t[[i]]$lon), FUN=mean)
}
	# plot to make sure it worked
	i <- 5 # pick the region to plot
	par(mfrow=c(1,2))
	plot(wind.t[[i]], col=c('blue', 'red')[1+(wind.t[[i]]$band1> 0)], pch=16, cex=0.05) # plots of <> a threshold
	plot(wind.sum[[i]]$lon, wind.sum[[i]]$lat, col=c('blue', 'red')[1+(wind.sum[[i]]$wind_npv > 0)], pch=16, cex=0.2)	

# concatenate the regions together
wind.out <- wind.sum[[1]]
if(length(wind.sum)>1){
	for(i in 2:length(wind.sum)){
		wind.out <- rbind(wind.out, wind.sum[[i]])
	}
}

# convert to positive longitude, to match climgrid
wind.out$lon[wind.out$lon<0] <- wind.out$lon[wind.out$lon<0] + 360
	range(wind.out$lon)
	
	# whole plot
	plot(wind.out$lon, wind.out$lat, col=c('blue', 'red')[1+(wind.out$wind_npv > 0)], pch=16, cex=0.3)	

# mark which grids are in clim
wind.out$keep <- FALSE
wind.out$keep[paste(wind.out$lat, wind.out$lon) %in% paste(clim$lat, clim$lon)] <- TRUE

	# compared untrimmed and trimmed
	par(mfrow=c(1,2))
	plot(wind.out$lon, wind.out$lat, col=c('blue', 'red')[1+(wind.out$wind_npv > 0)], pch=16, cex=0.3)	
	plot(wind.out$lon[wind.out$keep], wind.out$lat[wind.out$keep], col=c('blue', 'red')[1+(wind.out$wind_npv[wind.out$keep] > 0)], pch=16, cex=0.3)	

# remove grids not in clim
wind.out <- wind.out[wind.out$keep, c('lat', 'lon', 'wind_npv')]

# write out
write.csv(wind.out, 'cmsp_data/wind_npv.csv')

##############################
## Process NatCap Wave data ##
##############################
require(rgdal)

# read in data
wave_AKw <- readGDAL('cmsp_data/Wave/output/npv_usd_Alaska_west.tif') # read in SpatialGridDataFrame
wave_AKe <- readGDAL('./cmsp_data/Wave/output/npv_usd_Alaska_east.tif') # read in SpatialGridDataFrame
wave_WC <- readGDAL('./cmsp_data/Wave/output/npv_usd_WestCoast.tif')
wave_GM <- readGDAL('./cmsp_data/Wave/output/npv_usd_GoMex.tif')
wave_GMg <- readGDAL('./cmsp_data/Wave/output/npv_usd_GoMex_global.tif') # with global wave data
wave_NE <- readGDAL('./cmsp_data/Wave/output/npv_usd_Northeast.tif')
wave_NEg <- readGDAL('./cmsp_data/Wave/output/npv_usd_Northeast_global.tif') # with global wave data
	# image(wave_AKw)
	# image(wave_AKe)
	# image(wave_WC)
	# image(wave_GM)
	# image(wave_GMg)
	# image(wave_NE)
	# image(wave_NEg)

load(paste(climgridfolder, 'climGrid_rcp45.proj2.RData', sep='')) # loads clim. has the grid cells we want to summarize to

# trim clim to only non-duplicated rows (don't need every year) and those with climate data
dim(clim)
clim <- clim[!duplicated(clim[,c('lat', 'lon', 'region')]) & !is.na(clim$bottemp.proj_1) & !is.na(clim$surftemp.proj_1),]
dim(clim) # 6156

# combine wave output into a list
wave <- list(wave_AKw, wave_AKe, wave_WC, wave_GM, wave_GMg, wave_NE, wave_NEg)
names(wave) <- c('wave_AKw', 'wave_AKe', 'wave_WC', 'wave_GM', 'wave_GMg', 'wave_NE', 'wave_NEg')

# project to LL. also converts to SpatialPointsDataFrame
wave.t <- wave
for(i in 1:length(wave)){
	wave.t[[i]] <- spTransform(wave[[i]], CRS('+proj=longlat +data=WGS84')) # produces warnings about coercing to points. this is ok
#		colfun <- colorRamp(c('white', 'blue'))
#		sc <- (wave.t[[i]]$band1- min(wave.t[[i]]$band1))/(max(wave.t[[i]]$band1) - min(wave.t[[i]]$band1))
#		plot(wave.t[[i]], col=rgb(colfun(sc), maxColorValue=255), pch=16, cex=0.05) # plots
#		plot(wave.t[[i]], col=c('blue', 'red')[1+(wave.t[[i]]$band1> -600)], pch=16, cex=0.05) # plots of <> a threshold
} 
	
# Add a grid indicator
gridsize=0.25 # size of grid of the climate data, in degrees
for(i in 1:length(wave.t)){
	wave.t[[i]]$lat <- floor(coordinates(wave.t[[i]])[,2]/gridsize)*gridsize + gridsize/2 # round to nearest grid center
	wave.t[[i]]$lon <- floor(coordinates(wave.t[[i]])[,1]/gridsize)*gridsize + gridsize/2
}

# Summarize by grid cell
wave.sum <- wave.t
for(i in 1:length(wave.t)){
	wave.sum[[i]] <- aggregate(list(wave_npv = wave.t[[i]]$band1), by=list(lat = wave.t[[i]]$lat, lon = wave.t[[i]]$lon), FUN=mean)
}

	# plot to make sure it worked
	i <- 4 # pick the region to plot
	thresh <- median(wave.t[[i]]$band1)
	print(thresh)
	par(mfrow=c(1,2))
	plot(wave.t[[i]], col=c('blue', 'red')[1+(wave.t[[i]]$band1> thresh)], pch=16, cex=0.05) # plots of <> a threshold
	plot(wave.sum[[i]]$lon, wave.sum[[i]]$lat, col=c('blue', 'red')[1+(wave.sum[[i]]$wave_npv > thresh)], pch=16, cex=0.2)	

# concatenate the regions together
wave.out <- wave.sum[[1]]
wave.out$wave_region <- names(wave.sum)[1]
if(length(wave.sum)>1){
	for(i in 2:length(wave.sum)){
		new <- wave.sum[[i]]
		new$wave_region <- names(wave.sum)[i]
		wave.out <- rbind(wave.out, new)
	}
}

# convert to positive longitude, to match climgrid
wave.out$lon[wave.out$lon<0] <- wave.out$lon[wave.out$lon<0] + 360
	range(wave.out$lon)
	
	# whole plot
	plot(wave.out$lon, wave.out$lat, col=c('blue', 'red')[1+(wave.out$wave_npv > 0)], pch=16, cex=0.3)	# only positive NPV on west coast, SE Alaska, and off Newfoundland

# compare regional to global data in Gulf of Mexico and Northeast
#	gomcomp <- merge(wave.out[wave.out$wave_region=='wave_GM',], wave.out[wave.out$wave_region=='wave_GMg',], by=c('lat', 'lon'))
#		plot(gomcomp$wave_npv.x, gomcomp$wave_npv.y) # slight indications that global dataset is a little low on the low side, and a little high on the high side
#		abline(a=0,b=1)
#
#	necomp <- merge(wave.out[wave.out$wave_region=='wave_NE',], wave.out[wave.out$wave_region=='wave_NEg',], by=c('lat', 'lon'))
#		plot(necomp$wave_npv.x, necomp$wave_npv.y) # more scatter than for Gulf of Mexico
#		abline(a=0,b=1)


# mark which grids are in which clim region
dim(wave.out)
wave.out <- merge(wave.out, clim[,c('region', 'lat', 'lon')], all.x=TRUE)
dim(wave.out) # more rows. some lat/lon grids are in multiple regions

sum(is.na(wave.out$region)) # 8063 rows that don't match clim
table(wave.out$region, wave.out$wave_region) # which wave regions match which climate projection regions

	# plots of clim grid and wave region overlaps
#	plot(clim$lon, clim$lat, pch=16, cex=0.5, col=as.numeric(as.factor(clim$region)))
#	regwavedat <- wave.out$wave_region %in% c('wave_WC', 'wave_NE', 'wave_GM')
#	points(wave.out$lon[regwavedat], wave.out$lat[regwavedat], pch=16, cex=0.3, col='brown')
#	points(wave.out$lon, wave.out$lat, pch=16, cex=0.1) # only the west coast is fully covered by a high resolution wave dataset. Northeast regional data miss the offshore part of Georges Bank in NEUS. Gulf of Mexico regional data miss a few inshore and offshore points.

	# plot of zoomed in region to verify grid spacing
#	plot(clim$lon, clim$lat, pch=16, cex=1, xlim=c(285, 290), ylim=c(39, 42))
#	points(wave.out$lon, wave.out$lat, pch=16, cex=0.5, col='red') # yes, wave points are 0.25 deg apart, like the climate points

# remove global wave data from wave.out for Northeast and Gulf of Mexico where regional wave data are present
nrow(wave.out)
reginds <- wave.out$wave_region %in% c('wave_GM', 'wave_NE') # rows with regional data for NE or GoM
globdups <- (paste(wave.out$lat, wave.out$lon) %in% paste(wave.out$lat[reginds], wave.out$lon[reginds])) & (wave.out$wave_region %in% c('wave_GMg', 'wave_NEg'))
	# make sure I selected what I expected to select
#	sort(unique(wave.out$wave_region[globdups]))
#	plot(wave.out$lon[reginds], wave.out$lat[reginds]) # where I have regional data
#	points(wave.out$lon[globdups], wave.out$lat[globdups], col='red', pch=16, cex=0.5) # global data to be removed
#	globnotdups <- !globdups & wave.out$wave_region %in% c('wave_GMg', 'wave_NEg')
#	points(wave.out$lon[globnotdups], wave.out$lat[globnotdups], col='blue', pch=16, cex=0.5) # global data to not be removed
wave.out <- wave.out[!globdups,]
nrow(wave.out)
sort(unique(wave.out$wave_region))

tokeep <- !is.na(wave.out$region)

	# compared untrimmed and trimmed wave.out
	par(mfrow=c(1,2))
	plot(clim$lon, clim$lat, pch=16, cex=0.3)
	points(wave.out$lon, wave.out$lat, col=c('blue', 'red')[1+(wave.out$wave_npv > 0)], pch=16, cex=0.2)	
	plot(clim$lon, clim$lat, pch=16, cex=0.3)
	points(wave.out$lon[tokeep], wave.out$lat[tokeep], col=c('blue', 'red')[1+(wave.out$wave_npv[tokeep] > 0)], pch=16, cex=0.2)	

# remove grids not in clim
nrow(wave.out)
wave.out <- wave.out[tokeep, c('lat', 'lon', 'wave_npv')]
nrow(wave.out)

# order wave.out
wave.out <- wave.out[order(wave.out$lon, wave.out$lat),]

# remove rows duplicated by merge with clim regions
dups <- duplicated(wave.out)
#	sum(dups)
#	head(wave.out[paste(wave.out$lat, wave.out$lon) %in% paste(wave.out$lat[dups], wave.out$lon[dups]),])

	nrow(wave.out)
wave.out <- wave.out[!dups,]
	nrow(wave.out) # 4596

# write out
write.csv(wave.out, 'cmsp_data/wave_npv.csv')
