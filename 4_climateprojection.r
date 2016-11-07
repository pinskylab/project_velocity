## Set working directories
if(Sys.info()["nodename"] == "pinsky-macbookair"){
	setwd('~/Documents/Rutgers/Range projections/proj_ranges/')
	deltafolder <- '../data/' # not stored in Git
	}
if(Sys.info()["nodename"] == "amphiprion.deenr.rutgers.edu"){
	setwd('~/Documents/range_projections/')
	deltafolder <- 'data/'
	.libPaths(new='~/R/x86_64-redhat-linux-gnu-library/3.1/') # so that it can find my old packages (chron and ncdf4)
	}
# could add code for Lauren's working directory here

######################
# useful functions 
######################
require(Hmisc)

deg2rad <- function(deg) return(deg*pi/180)

# Haversine formula (hf) for distance between 2 points on a sphere
#from http://www.r-bloggers.com/great-circle-distance-calculations-in-r/
gcd.hf <- function(long1, lat1, long2, lat2) {
	long1 = deg2rad(long1); long2 = deg2rad(long2); lat1 = deg2rad(lat1); lat2=deg2rad(lat2)
	R <- 6371 # Earth mean radius [km]
	delta.long <- (long2 - long1)
	delta.lat <- (lat2 - lat1)
	a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
	c <- 2 * asin(pmin(1,sqrt(a)))
	return(d = R * c) # Distance in km
}

# weighted mean for use with summarize(). values in col 1, weights in col 2
wgtmean = function(x, na.rm=FALSE) wtd.mean(x=x[,1], weights=x[,2], na.rm=na.rm)

# rounds x to nearest number in y
roundto = function(x,y){r = which.min(abs(x-y)); return(y[r])}

	
####################################################
## Add deltas to climatology (run this overnight) ##
####################################################

# Load data needed to run
clim = read.csv('data/climGrid.csv', row.names=1, stringsAsFactors=FALSE); type='Grid' # the gridded climatology

# Choose which RCP to process
#rcp <- 85
rcp <- 45

if(rcp==85){
	load(paste(deltafolder, 'delta2100rcp85long.RData', sep='')) # loads delta2100rcp85long. slow (378M) 
	delta2100long <- delta2100rcp85long # script is written to operate on delta2100long
	rm(delta2100rcp85long)
}
if(rcp==45){
	load(paste(deltafolder, 'delta2100rcp45long.RData', sep='')) # loads delta2100rcp45long. slow (378M) 
	delta2100long <- delta2100rcp45long # script is written to operate on delta2100long
	rm(delta2100rcp45long)
}

regs = sort(unique(clim$region))

#####
# Start calculations
#####
# add depthgrid
for(i in 1:length(delta2100long)){ # for each climate model
	print(i)
	# make depth grid
	dgnm = paste('depthgrid', i, sep='') # different grid for each model
	clim[[dgnm]] = NA # to hold depth rounded to the GCM grid
	for(j in 1:length(regs)){ # for each region
		inds = clim$region == regs[j]
		#add depth grid specific to each region and each model
		inds2 = delta2100long[[i]]$region == regs[j]
		deps = sort(unique(delta2100long[[i]]$depth[inds2]))
		clim[[dgnm]][inds] = sapply(clim$depth[inds], FUN=roundto, y=deps)
	}
}


# expand clim to include each year in delta2100long
clim$year = min(delta2100long[[1]]$year)
clim_old = clim
nrow(clim)*80

for(i in (min(clim$year)+1):max(delta2100long[[1]]$year)){
	cat(i); cat(' ') # watch the loop slow down as clim gets big...
	clim_old$year = i
	clim = rbind(clim, clim_old)
		#dim(clim)
}
dim(clim)

# save a temporary copy in case the next step goes awry
clim_old2 = clim 
save(clim, file=paste(deltafolder, 'climGrid.proj_step1.RData', sep=''))
# load(paste(deltafolder, 'climGrid.proj_step1.RData', sep=''))

# merge in the bottemp deltas from delta2100long (~10 min)
for(i in 1:length(delta2100long)){ # for each climate model
	print(i)
	dgnm = paste('depthgrid', i, sep='') # different grid for each model

	# make sure all lats, lons, and depths match: should be numeric(0). Otherwise indicates places where clim spatial footprint exceeds our delta footprint
	lats1 = sort(unique(clim$latgrid))
	lats2 = sort(unique(delta2100long[[i]]$lat))
		print(paste('lat', paste(c(setdiff(lats1,lats2), setdiff(lats2, lats1)), collapse=', ')))
	lons1 = sort(unique(clim$longrid))
	lons2 = sort(unique(delta2100long[[i]]$lon))
		print(paste('lon', paste(c(setdiff(lons1,lons2), setdiff(lons2, lons1)), collapse=', ')))
	depths1 = sort(unique(clim[[dgnm]]))
	depths2 = sort(unique(delta2100long[[i]]$depth))
		print(paste('depth', paste(c(setdiff(depths1,depths2), setdiff(depths2, depths1)), collapse=', ')))

	# do the merge
		print(dim(clim))
	clim = merge(clim, delta2100long[[i]], by.x = c('region', 'latgrid', 'longrid', 'year', dgnm), by.y=c('region', 'lat', 'lon', 'year', 'depth'), all.x=TRUE) # this step takes many minutes each iteration
		print(dim(clim))
	dnm = paste('bottemp.delta_', i, sep='')
	names(clim)[names(clim)=='delta'] = dnm

}	
names(clim)

# save a temporary copy in case the next step goes awry
clim_old3 = clim 
save(clim, file=paste(deltafolder, 'climGrid.proj_step2.RData', sep=''))
# load(paste(deltafolder, 'climGrid.proj_step2.RData', sep=''))

# calculate future bottemp by adding delta to climatology
# look for NAs in bottemp.delta2100 and fix to nearest available delta (same lat/lon, shallower depth)
# occurs because GCM topography is coarse
# NOTE: This takes many hours (run overnight)
options(warn=2) # turn warnings to errors
clim = clim[order(clim$region, clim$lat, clim$lon, clim$depth, clim$year),]
for(i in 1:length(delta2100long)){ # for each climate model
	print(paste('Model', i))

	dnm = paste('bottemp.delta_', i, sep='')
	dgnm = paste('depthgrid', i, sep='') # different grid for each model
	missinds = which(is.na(clim[[dnm]]) & clim$year == 2020) # ignore different years: we'll do them all together. Also do all points in the same latgrid/longrid together
	missinds = missinds[!duplicated(clim[missinds, c('region', 'latgrid', 'longrid', dgnm)])] # only take unique latgrid/longrid/depthgrid/region, since all points that match on these criteria will be filled together
	
	print(paste(length(missinds), 'missing for 2020'))
	if(length(missinds)>0){
		full = delta2100long[[i]] # copy of this model output so that we can add missing years
		if(!all(2020:2100 %in% unique(full$year))){ # make sure we have every year (we don't for some climate models). we need all years (even if NA) in order to merge seemlessly into clim
			missingyears = setdiff(2020:2100, full$year)
			temp = full[!duplicated(full[,c('lon', 'lat', 'depth', 'region')]),]
			temp$delta = NA
			for(j in 1:length(missingyears)){ # add in the datapoints for each year if missing a year
				print(paste('adding', missingyears[j]))
				temp$year = missingyears[j]
				full = rbind(full,temp)
			}
			print(nrow(delta2100long[[i]]))
			print(nrow(full))
			print(nrow(full)/80) # should be an even number if equal numbers of lines per year
		}
		#full = full[order(full$region, full$lat, full$lon, full$depth, full$year),] # sort in same order as clim. slow and not needed?

		level1 = 0 # keep track of how often we had to go out to adjacent lat/lon to get a delta value

		for(k in 1:length(missinds)){ # for each missing delta value in clim
			if(k %% 100 == 0) print(k)
			fill = FALSE # indicator for whether to fill this missing delta value
			thisreg = clim$region[missinds[k]]
			thislat = clim$latgrid[missinds[k]]
			thislon = clim$longrid[missinds[k]]
			thisdepth = clim[[dgnm]][missinds[k]]
				#print(paste('k=', k, 'of', length(missinds), thisreg, thislat, thislon, thisdepth))
			useclim = which(clim$region == thisreg & clim$latgrid == thislat & clim$longrid == thislon & clim[[dgnm]] == thisdepth) # get all rows in clim at this lat/lon/depth (across all years). These are the rows in clim to fill with delta values.
			availinds2 = full$region == thisreg & full$lat == thislat & full$lon == thislon  # trim full (the climate model output) to focal region and lat/lon
				# full[availinds2,]; thisdepth
			if(!all(is.na(clim[[dnm]][useclim]))) stop('not all deltas are NA') # check that useclim's really are NAs to be filled
			if(any(!is.na(full$delta[availinds2]))){ # first, try to use the closest depth at the same place (if some values are not NA
				fulldeps = sort(unique(full$depth[availinds2 & !is.na(full$delta)])) # find the depths at this location that have non-NA deltas
				usedepth = fulldeps[abs(fulldeps - thisdepth) == min(abs(fulldeps - thisdepth))] # find the closest depth
				
				if(length(usedepth)>1) warning(paste('too many usedepths for i=', i, 'k=', k))
				usefull = which(full$depth[availinds2] == usedepth) # an index to the depth we want to use
					#print(paste(length(usefull), length(useclim))) # can be different if clim has multiple lat/lons that correspond to one model grid cell. but useclim must be a multiple of usefull
					#print(thisdepth)
				fill = TRUE # indicator for whether to fill this missing delta value
			}
			if(!any(!is.na(full$delta[availinds2]))){ # Level 1: if all deltas at this lat/lon are NA, use average of +/-10deg surrounding deltas at the same or closest depth and weight by 1/distance
				print(paste('going to level 1: i=', i, 'k=', k))
				level1 = level1 + 1
				availinds2 = full$region == thisreg & full$lat %in% c(thislat + -10:10) & full$lon %in% c(thislon + -10:10) # any depth
					#print(full[availinds2,])
				if(any(!is.na(full$delta[availinds2]))){ # try to use the closest depth in the 400-cell region
					fulldeps = sort(unique(full$depth[availinds2 & !is.na(full$delta)])) # find the depths at this location that have non-NA deltas
					usedepth = fulldeps[abs(fulldeps - thisdepth) == min(abs(fulldeps - thisdepth))] # find the closest depth

					if(length(usedepth)>1) warning(paste('too many usedepths for i=', i, 'k=', k))
					usefull = which(full$depth[availinds2] == usedepth) # an index to the depth we want to use. do it this way so that years with NAs are still included
						#print(paste(length(usefull), length(useclim))) # can be different if clim has multiple lat/lons that correspond to one model grid cell. but useclim must be a multiple of usefull
						#print(full[availinds2,][usefull,])
						#print(thisdepth)
					fill = TRUE # indicator for whether to fill this missing delta value
				} else {
					print(paste("still couldn't find a non-NA value +/- 10deg for k=", k))
				}

			} 
			if(fill){
				dists = pmax(1, gcd.hf(thislon, thislat, full$lon[availinds2][usefull], full$lat[availinds2][usefull])) # distances in km, but set 1km as a floor (so no zeros)

#				filler = aggregate(list(delta = full$delta[availinds2][usefull]), by = list(year = full$year[availinds2][usefull]), FUN=weighted.mean, na.rm=TRUE, w=1/dists) # summarize full deltas by year (average across cells available in each year), weighted by inverse distance
				filler = summarize(cbind(full$delta[availinds2][usefull], 1/dists), by = list(year = full$year[availinds2][usefull]), FUN=wgtmean, na.rm=TRUE, stat.name='delta')
				rownames(filler) = filler$year # make rows of filler easy to reference
				filler2 = filler[as.character(clim$year[useclim]), ] # re-order and expand filler to match years of clim rows to fill
				if(all(clim$year[useclim] == filler2$year)){
					clim[[dnm]][useclim] = filler2$delta
				} else {
					stop(paste('clim and filler year not in the same order, k=', k)) # make sure that years line up perfectly before we do the assign	
				}		
			}			
			if(k %% 100 == 0) cat(paste(k, ''))
		}
		cat('\n')	
		# make sure it worked
		print(paste('still missing', sum(is.na(clim[[dnm]])), 'on', dnm)) # should be 0
		print(paste('years are', paste(sort(unique(clim$year[is.na(clim[[dnm]])])), collapse=', ')))
		print(paste('Level1:', level1))
	}
	# calculate bottom temp projection for this model
	colnm = paste('bottemp.proj_', i, sep='')
	clim[[colnm]] = clim[[dnm]] + clim$bottemp.clim.int
	print(paste('still missing', sum(is.na(clim[[colnm]])), 'on', colnm))
	print(paste('years are', paste(sort(unique(clim$year[is.na(clim[[colnm]])])), collapse=', ')))
}

unique(clim$region)

# save a temporary copy in case the next step goes awry
clim_old4 = clim 
save(clim, file=paste(deltafolder, 'climGrid.proj_step3.RData', sep=''))
# load(paste(deltafolder, 'climGrid.proj_step3.RData', sep=''))

# merge in the surftemp deltas from delta2100long
for(i in 1:length(delta2100long)){ # for each climate model
	print(i)

	# choose only the surface
	depths2b = sort(unique(delta2100long[[i]]$depth))
	surfindsb = delta2100long[[i]]$depth == depths2b[1]

	# make sure all lats, lons, and depths match: should be numeric(0)
	lats1 = sort(unique(clim$latgrid))
	lats2 = sort(unique(delta2100long[[i]]$lat))
		print(c(setdiff(lats1,lats2), setdiff(lats2, lats1)))
	lons1 = sort(unique(clim$longrid))
	lons2 = sort(unique(delta2100long[[i]]$lon))
		print(c(setdiff(lons1,lons2), setdiff(lons2, lons1)))


	# do the merge
	print(dim(clim))
	cols = c("lon", "lat", "delta", "region", "year")
	clim = merge(clim, delta2100long[[i]][surfindsb,cols], by.x = c('region', 'latgrid', 'longrid', 'year'), by.y=c('region', 'lat', 'lon', 'year'), all.x=TRUE)
		print(dim(clim))
	dnm = paste('surftemp.delta_', i, sep='')
	names(clim)[names(clim)=='delta'] = dnm

}	
names(clim)

# save a temporary copy in case the next step goes awry
clim_old4 = clim 
save(clim, file=paste(deltafolder, 'climGrid.proj_step4.RData', sep=''))
# load(paste(deltafolder, 'climGrid.proj_step4.RData', sep=''))


# calculate future surftemp
# look for NAs in surftemp.delta and fix to nearest available delta (same lat/lon, shallower depth)
# occurs because GCM topography is coarse
# NOTE: This takes 30 min or so
options(warn=2) # turn warnings to errors
clim = clim[order(clim$region, clim$lat, clim$lon, clim$depth, clim$year),]
for(i in 1:length(delta2100long)){ # for each climate model
	print(paste('Model', i))

	dnm = paste('surftemp.delta_', i, sep='')
	missinds = which(is.na(clim[[dnm]]) & clim$year == 2020) # ignore different years: we'll do them all together. Also do all points in the same latgrid/longrid together
	missinds = missinds[!duplicated(clim[missinds, c('region', 'latgrid', 'longrid')])] # only take unique latgrid/longrid/region, since all points that match on these criteria will be filled together

	print(paste(length(missinds), 'missing for 2020'))
	if(length(missinds)>0){
		surfindsb = delta2100long[[i]]$depth == min(delta2100long[[i]]$depth) # get only the surface values
		cols = c("lon", "lat", "delta", "region", "year")
		full = delta2100long[[i]][surfindsb,cols] # only surface temperatures
		if(!all(2020:2100 %in% unique(full$year))){ # make sure we have every year (we don't for some climate models), and add them if not
			missingyears = (2020:2100)[!(2020:2100 %in% unique(full$year))]
			temp = full[!duplicated(full[,c('lon', 'lat', 'region')]),]
			temp$delta = NA
			for(j in 1:length(missingyears)){ # add in the datapoints for each year if missing a year
				print(paste('adding', missingyears[j]))
				temp$year = missingyears[j]
				full = rbind(full,temp)
			}
			print(nrow(delta2100long[[i]]))
			print(nrow(full))
			print(nrow(full)/80) # should be an even number if equal numbers of lines per year
		}
		full = full[order(full$region, full$lat, full$lon, full$year),] # sort in same order as clim

		level1 = 0 # keep track of how far I had to go to fill in the missing values

		for(k in 1:length(missinds)){
			fill = FALSE # indicator for whether to fill this missing delta value
			thisreg = clim$region[missinds[k]]
			thislat = clim$latgrid[missinds[k]]
			thislon = clim$longrid[missinds[k]]
				#print(paste('k=', k, thisreg, thislat, thislon))
			useclim = which(clim$region == thisreg & clim$latgrid == thislat & clim$longrid == thislon) # get all rows in clim at this lat/lon/region (across all years). These are the rows in clim to fill with delta values.
			availinds2 = full$region == thisreg & full$lat == thislat & full$lon == thislon  # trim full (the climate model output) to focal region and lat/lon
				#full[availinds2,]
			if(!all(is.na(clim[[dnm]][useclim]))) stop('not all deltas are NA') # check that useclim's really are NAs to be filled
			if(any(!is.na(full$delta[availinds2]))){ # first, try to use values at the same place (if some values are not NA
					print(paste('k=', k, 'VERY strange that a surface value exists here')) # we shouldn't end up here
					sum(availinds2)
					length(useclim) # can be different if clim has multiple lat/lons that correspond to one model grid cell. but useclim length must be a multiple of usefull length
					full[availinds2,]
					fill = TRUE
			}
			if(!any(!is.na(full$delta[availinds2]))){ # Level 1: if all deltas at this lat/lon are NA, use average of +/-10deg surrounding deltas
				#print(paste('going to level 1: i=', i, 'k=', k))
				level1 = level1 + 1
				#availinds2 = full$region == thisreg & full$lat %in% c(thislat, thislat-1, thislat+1) & full$lon %in% c(thislon, thislon-1, thislon+1) # 9 grid cells surroudning
				availinds2 = full$region == thisreg & full$lat %in% c(thislat + -10:10) & full$lon %in% c(thislon + -10:10) # +/-10deg

					#full[availinds2,]
				if(any(!is.na(full$delta[availinds2]))){ # first, try to use the closest depth at the same place
					fill = TRUE
				} else {
					print(paste("still couldn't find a non-NA value +/- 10deg for k=", k))
				}
			}
			if(fill){
				dists = pmax(1, gcd.hf(thislon, thislat, full$lon[availinds2], full$lat[availinds2])) # distances in km, but set 1km as a floor (so no zeros)

				filler = summarize(cbind(full$delta[availinds2], 1/dists), by = list(year = full$year[availinds2]), FUN=wgtmean, na.rm=TRUE, stat.name='delta') # summarize full deltas by year (average across cells available in each year and weight by inverse distance)
				rownames(filler) = filler$year # make rows of filler easy to reference
				filler2 = filler[as.character(clim$year[useclim]), ] # re-order and expand filler to match years of clim rows to fill
				if(all(clim$year[useclim] == filler2$year)){
					clim[[dnm]][useclim] = filler2$delta		
				} else {
					stop(paste('clim and filler year not in the same order, k=', k)) # make sure that years line up perfectly before we do the assign
				}
			}
			
			if(k %% 100 == 0) cat(paste(k, ''))
		}
		cat('\n')	
		# make sure it worked
		print(paste('still missing', sum(is.na(clim[[dnm]])), 'on', dnm))
		print(paste('years are', paste(sort(unique(clim$year[is.na(clim[[dnm]])])), collapse=', '))) # should be 0
		print(paste('Level1:', level1))
	}
	# calculate surface temp projection
	colnm = paste('surftemp.proj_', i, sep='')
	clim[[colnm]] = clim[[dnm]] + clim$surftemp.clim.int
	print(paste('still missing', sum(is.na(clim[[colnm]])), 'on', colnm))
	print(paste('years are', paste(sort(unique(clim$year[is.na(clim[[colnm]])])), collapse=', '))) # should only be 2100
}

unique(clim$region)

# write out clim dataframe with projections
	save(clim, file=paste(deltafolder, 'clim', type, '_rcp', rcp, '.proj2.RData', sep=''))


