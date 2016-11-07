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
	}
# could add code for Lauren's working directory here

#######################
## Process WDPA data ##
#######################
# Read in a shapefile of protecte area polygons
# Calculate area in each grid cell
# NEED TO TRIM WDPA TO RELEVANT MPAs of various types

require(maptools)
require(RColorBrewer)
require(rgdal)
#require(PBSmapping)
require(rgeos)

# read in the climate grid
clim = read.csv('data/climGrid.csv', row.names=1, stringsAsFactors=FALSE) # the climatology grid
	clim$lon[clim$lon>180] = clim$lon[clim$lon>180] - 360 # convert lon to shp format (-180 to 180)
	gridsz = 0.25 # grid size

# read in the LMEs
wdpa = readShapePoly('data/WDPA/NA_marine_MPA/mpinsky-search-1382225374362.shp', proj4string = CRS('+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs'))
	# plot(wdpa) # very slow
	wdpa = wdpa[wdpa$marine == 1,] # trim to only marine (remove 20 of 3023)
	
	# some data exploration
	#table(wdpa$marine)
	#sort(unique(wdpa$desig))
	#sort(unique(wdpa$desig_eng))
	#as.matrix(table(wdpa$iucn_cat))
	#table(wdpa$is_point==1) # 448 are points, so turn these into circles where possible
	#	table(wdpa$rep_m_area[wdpa$is_point]) # all points report an area of 20 km2. Suspicious.
	#	plot(wdpa[wdpa$is_point==1 & !is.na(wdpa$is_point),][448,]) # plot one. they seem to be approximated as parallelograms
	#table(is.na(wdpa$is_point)) # 668 report NA for is_point.. what to make of these?
	#	plot(wdpa[is.na(wdpa$is_point),][660,]) # plot one. it seems to have a shape
	#as.matrix(table(wdpa$no_take))
		
	# subset wdpa? e.g., to no-take areas
	
	# merge wdpa into one polygon?
	#wdpa_merge = unionSpatialPolygons(wdpa, rep("wdpa", nrow(wdpa)), threshold=0.0001) # threshold avoids slivers of area<0.0001 # extremely slow (>3 min)
		
# Make the grid based on corner points
	inds <- !duplicated(clim[,c('lat', 'lon')]) # clim has multiple entries per location for different seasons and regions
	lats = clim$lat[inds]; length(lats); y = numeric(5*length(lats)) # all the Y coords, in order
	# lower right, lower left, upper left, upper right, lower right
	for(i in 1:length(lats)){ y[(5*(i-1)+1):(5*(i-1)+5)] = c(lats[i]-gridsz/2, lats[i]-gridsz/2, lats[i]+gridsz/2, lats[i]+gridsz/2, lats[i]-gridsz/2)}

	lons = clim$lon[inds]; x = numeric(5*length(lons))
	# lower right, lower left, upper left, upper right, lower right: clockwise so that sp sees it as an island, not a hole
	for(i in 1:length(lons)){ x[(5*(i-1)+1):(5*(i-1)+5)] = c(lons[i]+gridsz/2, lons[i]-gridsz/2, lons[i]-gridsz/2, lons[i]+gridsz/2, lons[i]+gridsz/2) }

	# Create a SpatialPolygonsDataFrame (a list of "Polygons", each of which is a list of "Polygon")
	pgns = vector('list', length(lats))
	for(i in 1:length(pgns)){
		inds2 = (5*(i-1)+1):(5*(i-1)+5)
		pgns[[i]] = Polygons(list(Polygon(cbind(x[inds2],y[inds2]))), i)
	}
	SP <- SpatialPolygons(pgns, proj4string=CRS(proj4string(wdpa)))
	SPdata <- data.frame(gridpolyID = as.numeric(sapply(slot(SP, 'polygons'), slot, 'ID')), lat = clim$lat[inds], lon = clim$lon[inds], region=clim$region[inds]) # would be better to put this in with SP as a SpatialPolygonsDataFrame
		length(SP)
		#plot(SP[1:10])
		#plot(SP[1:1000]) # slow for so many
		#plot(SP) # very slow

# Intersect the grid and the PAs (slow: can skip if done earlier)
	# see https://stat.ethz.ch/pipermail/r-sig-geo/2012-June/015340.html and https://stat.ethz.ch/pipermail/r-sig-geo/2011-June/012099.html
	gI <- gIntersects(SP, wdpa, byid=TRUE) # slowish (~1 min). rownames are wdpa IDs, colnames are SP ID names
		sum(gI); dim(gI)
		ng = sum(colSums(gI)>0); ng # number of SP polygons that intersect wdpa
		cols = which(colSums(gI)>0) # ids of SP polygons that intersect wdpa
	out <- vector(mode="list", length=ng) # one entry for each SP grid cell
	for(i in 1:length(out)){ # 20 min? steps through each grid cell. quickly for some, slowly for others.
		print(paste(i, 'of', ng))
		out[[i]] <- gIntersection(SP[cols[i],], wdpa[gI[,cols[i]],], byid=TRUE) # only calc intersections for those polygons that actually intersect (according to gIntersects)
		if(class(out[[i]]) == 'SpatialCollections'){
			out[[i]] = out[[i]]@polyobj # extract only the polygon part if the intersection has points, lines, or rings as well
			print(paste('converted at', i))
		}
	}
	table(sapply(out, class)) # should only be SpatialPolygons
	out1 <- do.call("rbind", out) # make a big SpatialPolygons out of the intersections
	save(out1, file=paste('data/wdpa_by_grid', gridsz, '_intersect.RData', sep='')) # save
	
# Read in intersection (if desired)
load('data/wdpa_by_grid0.25_intersect.RData'); gridsz <- 0.25
	
# Calc fraction of each grid covered by a PA
	# Merge together wdpa pieces in the same grid cell
	rn = row.names(out1) # the rownames of rn have the rownames of SP and wdpa (which we now use as IDs) that make up the intersection
	nrn = do.call('rbind', strsplit(rn, " ")) # split apart the ids
	out2 <- gUnaryUnion(out1, id = nrn[,1]) # do the merge by grid cell ID (=rowname in SP)

	save(out2, file=paste('data/grid', gridsz, '_int_by_wdpa.RData', sep='')) # save the SpatialPolygons object

	# Calculate area of the grid intersected by PAs
	rn = as.numeric(row.names(out2)) # the rownames have the ids of SP (=rownames of SP)
	df = data.frame(grid = rn, area_wdpa=sapply(slot(out2, 'polygons'), slot, 'area')) # area of all intersected wdpa pieces	
	df$prop_wdpa = df$area_wdpa/gridsz^2 # turn to a proportion of the grid cell (e.g., 0.25x0.25°). would be even better to project to an equal-area projection...
	df <- df[order(df$grid),]
		head(df)
		summary(df$prop_wdpa)
		head(names(SP)) # match against these grid cells
		
		all(as.numeric(names(SP)) == 1:length(SP)) # TRUE: row names (and df$grid) == index
	
		# use SP row names (== df$grid) as index into lats and lons (originally used for grid centers)
		df$lat = clim$lat[df$grid]
		df$lon = clim$lon[df$grid]
		df$region = clim$region[df$grid]
	
	# add grid cells with 0 coverage
	newinds = setdiff(1:nrow(clim), df$grid) # all the rows in clim that didn't find a matching PA
	new = data.frame(grid=newinds, area_wdpa = rep(0, length(newinds)), prop_wdpa = rep(0, length(newinds)), lat = clim$lat[newinds], lon = clim$lon[newinds], region = clim$region[newinds])
	wdpacov = rbind(df, new)		
		dim(wdpacov)
		nrow(clim)
		hist(wdpacov$prop_wdpa, xlab='Proportion of grid cell covered by WDPA', main='')

	# re-order and add rownames
	wdpacov = wdpacov[order(wdpacov$grid),]
	row.names(wdpacov) = 1:nrow(wdpacov)
		head(wdpacov)

	# Write out
	write.csv(wdpacov, file=paste('data/grid', gridsz, '_cov_by_wdpa.csv', sep=''))

# Calculate the fraction of each PA in each grid cell
	rn <- row.names(out1) # the rownames have the ids of SP and wdpa that make up the intersection
	nrn <- do.call('rbind', strsplit(rn, " ")) # split apart the ids
	wdpa.by.grid <- data.frame(wdpapolyID = as.numeric(nrn[,2]), gridpolyID = as.numeric(nrn[,1]), area_grid=sapply(slot(out1, 'polygons'), slot, 'area')) # pull out area. would be even better to project to an equal-area projection...
		dim(wdpa.by.grid)
	wdpaareas <- data.frame(wdpapolyID=sapply(slot(wdpa, 'polygons'), slot, 'ID'), area_wdpa = sapply(slot(wdpa, 'polygons'), slot, 'area'))
	wdpa.by.grid <- merge(wdpa.by.grid, wdpaareas)
		dim(wdpa.by.grid)
		head(wdpa.by.grid)
	wdpa.by.grid$prop_grid <- wdpa.by.grid$area_grid/wdpa.by.grid$area_wdpa
		summary(wdpa.by.grid$prop_grid)

	# add in other wdpa metadata
	wdpatomerge = wdpa
	wdpatomerge@data <- cbind(wdpa@data, wdpapolyID = sapply(slot(wdpa, 'polygons'), slot, 'ID')) # add a column of polygon ID, to allow merging
		intersect(names(wdpa.by.grid), names(wdpatomerge@data)) # wdpapolyID
	wdpa.by.grid2 <- merge(wdpa.by.grid, wdpatomerge@data) # merges on wdpapolyID
		dim(wdpa.by.grid2)
		head(wdpa.by.grid2)

	# add in grid cell metadata
		intersect(names(wdpa.by.grid2), names(SPdata)) # gridpolyID
	wdpa.by.grid2 <- merge(wdpa.by.grid2, SPdata)
			
	# re-order
	wdpa.by.grid2 <- wdpa.by.grid2[order(wdpa.by.grid2$wdpapolyID, wdpa.by.grid2$gridpolyID),]
		head(wdpa.by.grid2)
	
	# write out
	write.csv(wdpa.by.grid2, file=paste('data/wdpa_cov_by_grid', gridsz, '.csv', sep=''))