###############################################################################
## Create climatologies for each region: gridded averages with interpolation ##
###############################################################################
setwd('~/Documents/Rutgers/Range projections/proj_ranges')
load("data/trawl_allregionsforprojections_wSST_2015-06-02.RData") # loads dat

# remove NA lat/lon
	dat = dat[!is.na(dat$lat) & !is.na(dat$lon),]

# output lat/lons from trawl tows for Lauren Roger to match against benthic habitat data
#	outlatlons = dat[!duplicated(dat[,c('lat', 'lon', 'depth')]),c('lat', 'lon', 'depth')]
#	outlatlons = outlatlons[order(outlatlons$lat, outlatlons$lon, outlatlons$depth),]
#	write.csv(outlatlons, file=paste('data/trawl_latlons_forLauren_', Sys.Date(), '.csv', sep=''), row.names=FALSE)

# Fix lon to only positive to match CMIP5
	dat$lon[dat$lon < 0] = dat$lon[dat$lon < 0] + 360 # fix lons to only positive to match CMIP5 data
	
# Prep temperature data
	temp = dat[complete.cases(dat[,c('lat', 'lon', 'depth')]),c('region', 'stratum', 'lat', 'lon', 'depth', 'year', 'yearsurv', 'month', 'bottemp', 'surftemp')]
		temp = droplevels(temp)
		rm(dat)

	# Add grid indicator for the climatology
	gridsize=0.25 # size of grid of the climate data, in degrees
	temp$latgrid = floor(temp$lat/gridsize)*gridsize + gridsize/2 # round to nearest grid center
	temp$longrid = floor(temp$lon/gridsize)*gridsize + gridsize/2
	
	# Trim to no later than December 16, 2005 (to match historical period in CMIP5)
	temp = temp[temp$year<2006,]
	
	# Label by decade within each region (using start date within region as start of the decades)
	regs = sort(unique(temp$region))
	temp$decade = NA
	for(i in 1:length(regs)){
		inds = temp$region == regs[i]
		yrs = sort(unique(temp$year[inds]))
		decs = seq(min(yrs), max(yrs), by=10) + 5 # center of each decade
		temp$decade[inds] = floor((temp$yearsurv[inds] - min(decs)+5)/10)*10+min(decs)
		decdf <- data.frame(decindex=1:length(unique(temp$decade[inds])), decade=sort(unique(temp$decade[inds])))
			row.names(decdf) <- decdf$decade
			
		if(nrow(decdf)>1 & (max(yrs) - max(decdf$decade) < 0)){ # if the survey ends less than half-way through the last decade, add its survey points to the previous decade (if more than 1 decade)
			decdf$decindex[nrow(decdf)] = decdf$decindex[nrow(decdf)-1]
		}
		temp$decindex[inds] <- decdf[as.character(temp$decade[inds]), 'decindex'] # index by rownames

		# look at the choices made
		print(regs[i])
		print(range(yrs))
		print(decdf)
	}
	
		table(temp$region, temp$decindex)

	# simplistic averaging by month within regions to look at cold vs. warm seasons
#	tempsbyreg <- aggregate(list(bottemp=temp$bottemp, surftemp=temp$surftemp), by=list(region=temp$region, month=temp$month), FUN=mean, na.rm=TRUE)
#		require(lattice)
#		xyplot(bottemp~month|region, data=tempsbyreg, type=c('o', 'g'), scales=list(x=list(relation='same'), y=list(relation='free')))
#
	
# Average by grid cell and by decade within seasons. lat, lon are now grid centers
	climdec = aggregate(list(bottemp.clim = temp$bottemp, surftemp.clim = temp$surftemp, depth=temp$depth), by=list(lat = temp$latgrid, lon = temp$longrid, decade = temp$decindex, region=temp$region), FUN=mean, na.rm=TRUE)
		dim(climdec) # 13,130
		names(climdec)
		
# Make a list of all grid cells in all regions in all decades in all season
	full <- data.frame(lat=numeric(0), lon=numeric(0), decade=numeric(0), region=character(0))
	regs <- sort(unique(climdec$region))
	for(i in 1:length(regs)){
		inds <- climdec$region == regs[i]
		decs <- data.frame(decade = sort(unique(climdec$decade[inds])))
		inds2 <- climdec$region == regs[i] & !duplicated(climdec[,c('lat', 'lon', 'region')]) # unique lat/lon
		newregdf <- merge(climdec[inds2,c('lat', 'lon', 'region')], decs)
		full <- rbind(full, newregdf)
	}
	dim(full) # 15,058
	
# Outer join to add NAs for temperatures not measured
	climdec2 <- merge(climdec, full, all=TRUE)
		dim(climdec2) # 15,058

# Average across decades (done as a second step to avoid one decade with lots of data swamping the average)
	clim = aggregate(list(bottemp.clim = climdec2$bottemp.clim, surftemp.clim = climdec2$surftemp.clim), by=list(lat = climdec2$lat, lon = climdec2$lon,  region=climdec2$region), FUN=mean) # don't remove NAs because it would mean dropping a decade
		dim(clim) # 6,535
		names(clim)

	climdepth = aggregate(list(depth = climdec2$depth), by=list(lat = climdec2$lat, lon = climdec2$lon,  region=climdec2$region), FUN=mean, na.rm=TRUE) # remove NAs since not worried about decadal variability in depth
		dim(climdepth) # 6,535
		summary(climdepth)
	
	clim = merge(clim, climdepth) # merge in depth data
		dim(clim)
	
		# plot each region in a separate figure. color axis is relative within each region
		require(lattice)
		cols = colorRampPalette(colors = c('blue', 'white', 'red'))
		pdf(width=10, height=9, file=paste('figures/climBT_grid.pdf', sep=''))
		levelplot(bottemp.clim ~ lon*lat|region, scales = list(x='free', y='free'), data=clim, at = seq(0,1, length.out=40), col.regions=cols, panel=function(...,z, subscripts, at){
				panel.fill(col='light grey')
				panel.levelplot(..., z=z, subscripts=subscripts, at=seq(min(z[subscripts], na.rm=TRUE), max(z[subscripts], na.rm=TRUE), length.out=20))})

		dev.off()

		pdf(width=10, height=9, file=paste('figures/climSST_grid.pdf', sep=''))
		levelplot(surftemp.clim ~ lon*lat|region, scales = list(x='free', y='free'), data=clim, at = seq(0,1, length.out=40), col.regions=cols, panel=function(...,z, subscripts, at){
				panel.fill(col='light grey')
				panel.levelplot(..., z=z, subscripts=subscripts, at=seq(min(z[subscripts], na.rm=TRUE), max(z[subscripts], na.rm=TRUE), length.out=20))})

		dev.off()


# Interpolate temperatures
	# useful functions from http://www.r-bloggers.com/great-circle-distance-calculations-in-r/
	deg2rad <- function(deg) return(deg*pi/180)
	# Haversine formula (hf)
	gcd.hf <- function(long1, lat1, long2, lat2) {
		long1 = deg2rad(long1); long2 = deg2rad(long2); lat1 = deg2rad(lat1); lat2=deg2rad(lat2)
		R <- 6371 # Earth mean radius [km]
		delta.long <- (long2 - long1)
		delta.lat <- (lat2 - lat1)
		a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
		c <- 2 * asin(pmin(1,sqrt(a)))
		return(d = R * c) # Distance in km
	}

	sum(is.na(clim$bottemp.clim)) # = 1371: some missing values
	inds = which(is.na(clim$bottemp.clim)) # points to fill for BT
	clim$bottemp.clim.int = clim$bottemp.clim # interpolated BT
	for(i in 1:length(inds)){
		lt = clim$lat[inds[i]]
		ln = clim$lon[inds[i]]
		nn = clim$lat %in% c(lt, lt-gridsize, lt+gridsize, lt-gridsize*2, lt+gridsize*2, lt-gridsize*3, lt+gridsize*3) & clim$lon %in% c(ln, ln-gridsize, ln+gridsize, ln-gridsize*2, ln+gridsize*2, ln-gridsize*3, ln+gridsize*3) & clim$region == unique(clim$region[inds[i]]) # nearest neighbors from same region (including the focal grid cell)
		dist = pmax(1, gcd.hf(ln, lt, clim$lon[nn], clim$lat[nn])) # distances in km, but set 1km as a floor (so no zeros)
		clim$bottemp.clim.int[inds[i]] = weighted.mean(x=clim$bottemp.clim[nn], w=1/dist, na.rm=TRUE) # weight values from 3 grid cells in every direction by 1/distance
	}	
		sum(is.na(clim$bottemp.clim.int)) # = 379: some missing values
		unique(clim$region[is.na(clim$bottemp.clim.int)]) # WCTri, ScotianShelf, NEUSSpring and Fall, GOMex
		table(clim$region, is.na(clim$bottemp.clim.int))

	sum(is.na(clim$surftemp.clim)) # = 1173 missing SST values
	inds = which(is.na(clim$surftemp.clim)) # points to fill for SST
	clim$surftemp.clim.int = clim$surftemp.clim # interpolated ST
	for(i in 1:length(inds)){
		lt = clim$lat[inds[i]]
		ln = clim$lon[inds[i]]
		nn = clim$lat %in% c(lt, lt-gridsize, lt+gridsize, lt-gridsize*2, lt+gridsize*2, lt-gridsize*3, lt+gridsize*3) & clim$lon %in% c(ln, ln-gridsize, ln+gridsize, ln-gridsize*2, ln+gridsize*2, ln-gridsize*3, ln+gridsize*3) & clim$region == unique(clim$region[inds[i]])# nearest neighbors (including the focal grid cell)
		dist = pmax(1,gcd.hf(ln, lt, clim$lon[nn], clim$lat[nn])) # distances in km (>=1)
		clim$surftemp.clim.int[inds[i]] = weighted.mean(x=clim$surftemp.clim[nn], w=1/dist, na.rm=TRUE)
	}	
		sum(is.na(clim$surftemp.clim.int)) # = 378 missing SST values
		unique(clim$region[is.na(clim$surftemp.clim.int)]) # WCTri, ScotianShelf, NEUSSpring and Fall, GOMex
		table(clim$region, is.na(clim$surftemp.clim.int))

	inds = which(is.na(clim$depth)) # points to fill for depth (0)
		length(inds)

	# plot each region in a separate figure (interpolated temp)
		#BT
	require(lattice)
	cols = colorRampPalette(colors = c('blue', 'white', 'red'))
	pdf(width=10, height=9, file=paste('figures/climBT_grid_interp.pdf', sep=''))
	levelplot(bottemp.clim.int ~ lon*lat|region, scales = list(x='free', y='free'), data=clim, at = seq(0,1, length.out=40), col.regions = cols, panel=function(...,z, subscripts, at){
			panel.fill(col='light grey')
			panel.levelplot(..., z=z, subscripts=subscripts, at=seq(min(z[subscripts], na.rm=TRUE), max(z[subscripts], na.rm=TRUE), length.out=20))})
	dev.off()

		#SST
	pdf(width=10, height=9, file=paste('figures/climSST_grid_interp.pdf', sep=''))
	levelplot(surftemp.clim.int ~ lon*lat|region, scales = list(x='free', y='free'), data=clim, at = seq(0,1, length.out=40), col.regions = cols, panel=function(...,z, subscripts, at){
			panel.fill(col='light grey')
			panel.levelplot(..., z=z, subscripts=subscripts, at=seq(min(z[subscripts], na.rm=TRUE), max(z[subscripts], na.rm=TRUE), length.out=20))})

	dev.off()

		#Depth
	pdf(width=10, height=6, file=paste('figures/climDepth_grid_interp.pdf', sep=''))
	levelplot(depth ~ lon*lat|region, scales = list(x='free', y='free'), data=clim, at = seq(0,1, length.out=40), col.regions = terrain.colors, panel=function(...,z, subscripts, at){
			panel.fill(col='light grey')
			panel.levelplot(..., z=z, subscripts=subscripts, at=quantile(z[subscripts], probs = seq(0,1, length.out=40), na.rm=TRUE))})

	dev.off()

# Add grid indicator for merging with climate deltas (1deg instead of 1/4deg for climatology)
	climgridsize=1 # size of grid of the climate data, in degrees
	clim$latgrid = floor(clim$lat/climgridsize)*climgridsize + climgridsize/2 # round to nearest grid center
	clim$longrid = floor(clim$lon/climgridsize)*climgridsize + climgridsize/2

## Add stratum based on majority within the grid (a slow process since line-by-line, 15+ min)
#	clim$stratum = character(nrow(clim))
#	ties = numeric(0)
#	for(i in 1:nrow(clim)){
#		if(i %% 100 == 0) print(i)
#		strats = temp$stratum[temp$latgrid == clim$lat[i] & temp$longrid == clim$lon[i] & temp$region == clim$region[i]]
#		t = sort(table(strats), decreasing=TRUE) # sort by frequency
#		if(length(t)>1){ # check for ties
#			if(t[1] == t[2]){
#				print(paste('there was a tie for i=', i))
#				ties = c(ties, i)
#			}
#		}
#		clim$stratum[i] = names(t[1])
#		
#		# this would be for nearest neighbor instead
#		#dist = gcd.hf(clim$lon[i], clim$lat[i], temp$lon, temp$lat) # distances in km to all hauls
#		#clim$stratum[i] = as.character(temp$stratum[which.min(dist)])
#	}
#	length(ties) # number of grid cells with ties between which stratum had the most points (34)
#	sum(is.na(clim$stratum)) # 0

# write out climatology
## WOULD BE GOOD TO TRIM OUT LINES WITH NA VALUES HERE (E.G., BOTTEMP.CLIM.INT OR SURFTEMP.CLIM.INT)
	write.csv(clim, file=paste('data/climGrid.csv', sep=''))
	
# write out all grid cell lat/lons for Lauren Rogers
#	clim = read.csv('data/climGrid_2015-02-02.csv')
#	outlatlons2 = clim[!duplicated(clim[,c('lat', 'lon')]),c('lat', 'lon', 'depth')]
#	outlatlons2 = outlatlons2[order(outlatlons2$lat, outlatlons2$lon),]
#	write.csv(outlatlons2, file=paste('data/projectiongrid_latlons_forLauren_', Sys.Date(), '.csv', sep=''), row.names=FALSE)
#
#	# split each grid cell into 4 sub-cells (from 1/4 to 1/8 degree)
#	a = numeric(nrow(outlatlons2)*4)
#	outlatlons3 = data.frame(lat=a, lon=a, depth=a)
#	for(i in 1:nrow(outlatlons2)){
#		outlatlons3$lat[((i-1)*4+1):(i*4)] = outlatlons2$lat[i] + rep(c(-0.0625, 0.0625),2)
#		outlatlons3$lon[((i-1)*4+1):(i*4)] = outlatlons2$lon[i] + rep(c(-0.0625, 0.0625),c(2,2))
#	}
#	
#	write.csv(outlatlons3, file=paste('data/projectiongrid_latlons1.8th_forLauren_', Sys.Date(), '.csv', sep=''), row.names=FALSE)
#	
#	# split each 1/8 grid cell into 4 sub-cells (from 1/8 to 1/16 degree)
#	a = numeric(nrow(outlatlons3)*4)
#	outlatlons4 = data.frame(lat=a, lon=a, depth=a)
#	nrow(outlatlons3)
#	for(i in 1:nrow(outlatlons3)){ # takes a few minutes
#		if(i %% 1000 == 0) print(i) # a little progress meter
#		outlatlons4$lat[((i-1)*4+1):(i*4)] = outlatlons3$lat[i] + rep(c(-0.03125, 0.03125),2)
#		outlatlons4$lon[((i-1)*4+1):(i*4)] = outlatlons3$lon[i] + rep(c(-0.03125, 0.03125),c(2,2))
#	}
#	
#	write.csv(outlatlons4, file=paste('data/projectiongrid_latlons1.16th_forLauren_', Sys.Date(), '.csv', sep=''), row.names=FALSE)
	