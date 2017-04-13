library(data.table)
library(raster)
library(ncdf4)
library(rasterVis)
library(chron)
library(ggplot2)
setwd('/Users/abigail/Desktop/Jim_Imac/Documents/proj_ranges/Jim')

# First plot up some info on haul locations so I get the projection grid right
load('haulInfo_Dec26_2016.RData')
ggplot(hauls, aes(x=depth)) + geom_histogram(binwidth = 50) + facet_wrap(~region)
ggplot(hauls[hauls$region == 'DFO_NewfoundlandFall',], aes(x=lat)) + geom_histogram(binwidth = 1)

# ==============================================================================================================================================
# Develope the projection grid==============================================================================================================================================
# ==============================================================================================================================================

# Load one SODA data file
soda <- nc_open('SODA temp climatology/soda3.3.1_5dy_ocean_or_2008_01_01.nc')
print(soda)
names(soda$var) # Gives the names of all the variables provided (like temperature)
sodaTemp <- ncvar_get(soda, 'temp') # Extracts one of the 26 variables in this nc file_an array_temperature
dim(sodaTemp)
sodaTemp[300,400,10] # can get cell values within the array
sodaSurf <- sodaTemp[,,1] # creates a matrix of just the surface layer
dim(sodaSurf) # 1440 by 1070

#extract some of the 'dimensions' within the soda set to set up a projection grid
lon <- as.vector(ncvar_get(soda, 'xt_ocean')) 
lat <- as.vector(ncvar_get(soda, 'yt_ocean'))
# round lat so that it has a reasonable number of significant figures
lat <- signif(lat, digits = 6) # one more sigfig than lon, for no real reason other than it is slightly higher resolution

# image(lon, lat, sodaSurf, col = rev(brewer.pal(11, "RdBu")))

grid <- expand.grid(lon = lon, lat = lat)

# Different type of plot
# cutpts <- c(-5, 0, 5, 10, 15, 20, 25, 30, 35)
# levelplot(sodaSurf ~ lon * lat, data = grid, at = cutpts, cuts = 9, pretty = T, col.regions = (rev(brewer.pal(10, "RdBu"))))

# turn into data.frame
sodaSurfVec <- as.vector(sodaSurf)
SODA <- as.data.frame(cbind(grid, sst = sodaSurfVec))
# plot(sst~lat, SODA) # looks good, dome shaped and warm at equator
rm(grid, sodaSurf, lat, lon, soda, sodaSurfVec, sodaTemp)

# ==============================================================================================================================================
# Now working with the bathymetry EAST data, just within the .nc format b/c we don't need to calculate rugosity for this
# ==============================================================================================================================================

bathy <- nc_open('Rugosity/Rugosity_data/GEBCO_2014_2D_-100.1942_23.6311_-42.3301_62.0777.nc')
print(bathy)
names(bathy$var) # Elevation is only variable
bathyArray <- ncvar_get(bathy, 'elevation') 
dim(bathyArray)
bathyArray[100,240]# can get cell values within the matrix

lon <- ncvar_get(bathy, 'lon') 
lat <- ncvar_get(bathy, 'lat')
# image(lon, lat, bathyArray, col = rev(brewer.pal(11, "RdBu")))

grid <- expand.grid(lon = lon, lat = lat)

# turn into data.frame
bathyVec <- as.vector(bathyArray)
bathy.grid <- as.data.frame(cbind(grid, elev = bathyVec))

# First drop obvious rows of non interest, i.e. positive values or too deep
bathy.grid <- bathy.grid[bathy.grid$elev < 0,]
bathy.grid <- bathy.grid[bathy.grid$elev > -401,]
# Create longitude bins with at .25 degree resolution to match SODA
gridsize <- 0.25 
bathy.grid$longrid <- floor(bathy.grid$lon/gridsize)*gridsize + gridsize/2

# Bin to soda lat grid_tougher b/c resolution changes with latitude
# First cut down the SODA grid to the right dimensions in bathygrid to save time
SODAeast <- SODA[SODA$lon < (max(bathy.grid$lon)+.5) & SODA$lon > (min(bathy.grid$lon)-.5),] # added/subtracted the .5 to buffer a bit
SODAeast <- SODAeast[SODAeast$lat < (max(bathy.grid$lat)+.5) & SODAeast$lat > (min(bathy.grid$lat)-.5),] # added/subtracted the .5 to buffer a bit
SODAeast <- SODAeast[!is.na(SODAeast$sst),] # Drop rows with no temp data (ie on land)

#For loop to assign bathymetry values to SODA lat grid bins
SODAlat <- unique(as.vector(SODAeast$lat)) # Vector of all the lats in the SODA grid_to be used in a loop
bathyLat <- unique(as.vector(bathy.grid$lat)) # all the lat values that need to be associated with the correct SODA bins with a forloop

latBins <- data.frame(lat=numeric(), latgrid=numeric()) 
  
for(i in 1:length(bathyLat)){
  lat = bathyLat[i]
  bins = data.frame(lat = lat, latBin = SODAlat)
  bins$diffs = abs(bins$lat - bins$latBin) # column showing how far away a lat value is from each SODA bin
  latBin = bins$latBin[bins$diffs == min(bins$diffs)] # Choose the SODA lat bin that the bathymetry lat value is closest to
  latBins[i,] = data.frame(lat=lat, latgrid=latBin)
}

# Merge the lat bins over to the master bathymetry file and then Aggregate depths to the right projection grid
proj.grid <- merge(bathy.grid, latBins, all.x=T, by='lat', sort=F)
proj.grid <- data.table(proj.grid) # easier to do aggregations
proj.grid.ag <- proj.grid[, list(depth = mean(elev), depthrange = (max(elev) - min(elev))), by = list(longrid, latgrid)]

#plot up depth data 
plot(latgrid~longrid, data=proj.grid.ag) # can see I'll need to 'prune' out some lat/lon boxes
levelplot(depth ~ longrid * latgrid, data = proj.grid.ag) 
levelplot(depthrange ~ longrid * latgrid, data = proj.grid.ag) 

# Near the coastlines some grid cells occur that have no SODA data
# change names of SODA lat/lon columns and merge, then look for cells with NA for sst
setnames(SODAeast, c('longrid', 'latgrid', 'sst'))
proj.grid.east <- merge(proj.grid.ag, SODAeast, all.x=T, by=c('longrid', 'latgrid'), sort=F)
summary(proj.grid.east) # there are 1310 NAs for sst
# Where are the NAs
plot(latgrid~longrid, data=proj.grid.east[is.na(proj.grid.east$sst),]) # makes sense, most nas in embayments and freshwater, try dropping them and seeing if coastline looks better
plot(latgrid~longrid, data=proj.grid.east[!is.na(proj.grid.east$sst),]) 
levelplot(depth ~ longrid * latgrid, data = proj.grid.east[!is.na(proj.grid.east$sst),])# much improved coastlines
# better resolution than with the original trawl data method (based on looking at google maps and old climatology), although some of the cells right on the coast may need to be removed later (e.g. pamlico sound, chesapeake bay)

proj.grid.east <- proj.grid.east[!is.na(proj.grid.east$sst),]
levelplot(depth ~ longrid * latgrid, data = proj.grid.east)# much improved coastlines
# Do some pruning
proj.grid.east <- proj.grid.east[longrid < -43.5 & latgrid > 24 & latgrid < 61]
proj.grid.east <- proj.grid.east[!(longrid < -72 & latgrid > 50)]
proj.grid.east <- proj.grid.east[!(longrid > -50 & latgrid > 57)]
proj.grid.east <- proj.grid.east[!(longrid > -67 & latgrid < 37)]
proj.grid.east <- proj.grid.east[!(longrid > -90 & longrid < -87 & latgrid < 25)]
# Kept the bahamas
levelplot(depth ~ longrid * latgrid, data = proj.grid.east)# much improved coastlines

save(proj.grid.east, file='projectionGridEast_Jan9_2017.RData')
rm(SODAeast)

# ==============================================================================================================================================
# Now do the identical with west coast and merge the two
# ==============================================================================================================================================

bathy <- nc_open('Rugosity/Rugosity_data/GEBCO_2014_2D_-179.4175_13.6893_-113.3981_70.3883.nc')
print(bathy)
names(bathy$var) # Elevation is only variable
bathyArray <- ncvar_get(bathy, 'elevation') 
dim(bathyArray)

lon <- ncvar_get(bathy, 'lon') 
lat <- ncvar_get(bathy, 'lat')
image(lon, lat, bathyArray, col = rev(brewer.pal(11, "RdBu")))

grid <- expand.grid(lon = lon, lat = lat)

# turn into data.frame
bathyVec <- as.vector(bathyArray)
bathy.grid <- as.data.frame(cbind(grid, elev = bathyVec))

# First drop obvious rows of non interest, i.e. positive values or too deep
bathy.grid <- bathy.grid[bathy.grid$elev < 0,]
bathy.grid <- bathy.grid[bathy.grid$elev > -401,]
# Create longitude bins with at .25 degree resolution to match SODA
gridsize <- 0.25 
bathy.grid$longrid <- floor(bathy.grid$lon/gridsize)*gridsize + gridsize/2

# Bin to soda lat grid_tougher b/c resolution changes with latitude
# First cut down the SODA grid to the right dimensions in bathygrid to save time
SODAwest <- SODA[SODA$lon < (max(bathy.grid$lon)+.5) & SODA$lon > (min(bathy.grid$lon)-.5),] # added/subtracted the .5 to buffer a bit
SODAwest <- SODAwest[SODAwest$lat < (max(bathy.grid$lat)+.5) & SODAwest$lat > (min(bathy.grid$lat)-.5),] # added/subtracted the .5 to buffer a bit
SODAwest <- SODAwest[!is.na(SODAwest$sst),] # Drop rows with no temp data (ie on land)

#For loop to assign bathymetry values to SODA lat grid bins
SODAlat <- unique(as.vector(SODAwest$lat)) # Vector of all the lats in the SODA grid_to be used in a loop
bathyLat <- unique(as.vector(bathy.grid$lat)) # all the lat values that need to be associated with the correct SODA bin with a forloop

latBins <- data.frame(lat=numeric(), latgrid=numeric()) 

for(i in 1:length(bathyLat)){
  lat = bathyLat[i]
  bins = data.frame(lat = lat, latBin = SODAlat)
  bins$diffs = abs(bins$lat - bins$latBin) # column showing how far away a lat value is from each SODA bin
  latBin = bins$latBin[bins$diffs == min(bins$diffs)] # Choose the SODA lat bin that the bathymetry lat value is closest to
  latBins[i,] = data.frame(lat=lat, latgrid=latBin)
}

# Merge the lat bins over to the master bathymetry file and then Aggregate depths to the right projection grid
proj.grid <- merge(bathy.grid, latBins, all.x=T, by='lat', sort=F)
proj.grid <- data.table(proj.grid) # easier to do aggregations
proj.grid.ag <- proj.grid[, list(depth = mean(elev), depthrange = (max(elev) - min(elev))), by = list(longrid, latgrid)]

#plot up depth data 
plot(latgrid~longrid, data=proj.grid.ag) # can see I'll need to 'prune' out some lat/lon boxes
levelplot(depth ~ longrid * latgrid, data = proj.grid.ag) 
levelplot(depthrange ~ longrid * latgrid, data = proj.grid.ag) 

# I think there are some lat/lon combinations that may not be in SODA, near the coastlines
# change names of SODA lat/lon columns and merge, then look for cells with NA for sst
proj.grid.west <- merge(proj.grid.ag, SODAwest, all.x=T, by=c('longrid', 'latgrid'), sort=F)
summary(proj.grid.west) # there are 2587 NAs for sst
# Where are the NAs
plot(latgrid~longrid, data=proj.grid.west[is.na(proj.grid.west$sst),]) # makes sense, most nas in embayments and freshwater, try dropping them and seeing if coastline looks better
plot(latgrid~longrid, data=proj.grid.west[!is.na(proj.grid.west$sst),]) 
levelplot(depth ~ longrid * latgrid, data = proj.grid.west[!is.na(proj.grid.west$sst),])# much improved coastlines
# better resolution than with the original trawl data method (based on looking at google maps and original climatology), although some of the cells right on the coast may need to be removed later

proj.grid.west <- proj.grid.west[!is.na(proj.grid.west$sst),]
levelplot(depth ~ longrid * latgrid, data = proj.grid.west) # Need to prune some areas out
proj.grid.west <- proj.grid.west[longrid < -116 & latgrid > 31.5 & latgrid < 65]
proj.grid.west <- proj.grid.west[!(longrid < -140 & longrid > -153 & latgrid > 50 & latgrid < 55)]
proj.grid.west <- proj.grid.west[!(longrid < -130 & latgrid < 47)]
proj.grid.west <- proj.grid.west[!(longrid < -135 & longrid > -145 & latgrid < 54)]
proj.grid.west <- proj.grid.west[!(longrid < -130.5 & latgrid < 50.8)]
proj.grid.west <- proj.grid.west[!(longrid < -179 & latgrid < 55 & latgrid > 52)]

save(proj.grid.west, file='projectionGridWest_Jan9_2017.RData')
rm(bathy.grid, bathyArray, bins, grid, latBins, proj.grid, proj.grid.ag, SODA, SODAwest, bathy, bathyLat, bathyVec, gridsize, lon, SODAlat, i, lat, latBin)

# ==============================================================================================================================================
# now to get the rest of the Aleutian islands on there
# ==============================================================================================================================================

bathy.aleut <- nc_open('Rugosity/Rugosity_data/GEBCO_2014_2D_162.0_47.0_180.0_59.0.nc')
bathy.aleut2 <- nc_open('Rugosity/Rugosity_data/GEBCO_2014_2D_-180.0_46.0_-173.0_57.0.nc')

bathyArray <- ncvar_get(bathy.aleut, 'elevation') 
bathyArray2 <- ncvar_get(bathy.aleut2, 'elevation') 

lon <- ncvar_get(bathy.aleut, 'lon') 
lon <- -180 - (180-lon) # to convert to the SODA longitude_not necessary with bathyArray2
lat <- ncvar_get(bathy.aleut, 'lat')
# image(lon, lat, bathyArray, col = rev(brewer.pal(11, "RdBu")))

grid <- expand.grid(lon = lon, lat = lat)

# turn into data.frame
bathyVec <- as.vector(bathyArray)
bathy.grid <- as.data.frame(cbind(grid, elev = bathyVec))

# First drop obvious rows of non interest, i.e. positive values or too deep
bathy.grid <- bathy.grid[bathy.grid$elev < 0,]
bathy.grid <- bathy.grid[bathy.grid$elev > -401,]
# Create longitude bins with at .25 degree resolution to match SODA
gridsize <- 0.25 
bathy.grid$longrid <- floor(bathy.grid$lon/gridsize)*gridsize + gridsize/2

# Bin to soda lat grid_tougher b/c resolution changes with latitude
# First cut down the SODA grid to the right dimensions in bathygrid to save time
SODAaleut <- SODA[SODA$lon < (max(bathy.grid$lon)+.5) & SODA$lon > (min(bathy.grid$lon)-.5),] # added/subtracted the .5 to buffer a bit
SODAaleut <- SODAaleut[SODAaleut$lat < (max(bathy.grid$lat)+.5) & SODAaleut$lat > (min(bathy.grid$lat)-.5),] # added/subtracted the .5 to buffer a bit
levelplot(sst ~ lon * lat, data = SODAaleut) # looks like I converted the longitude correctly
SODAaleut <- SODAaleut[!is.na(SODAaleut$sst),] # Drop rows with no temp data (ie on land)

#For loop to assign bathymetry values to SODA lat grid bins
SODAlat <- unique(as.vector(SODAaleut$lat)) # Vector of all the lats in the SODA grid_to be used in a loop
bathyLat <- unique(as.vector(bathy.grid$lat)) # all the lat values that need to be associated with the correct SODA bins with a forloop

latBins <- data.frame(lat=numeric(), latgrid=numeric()) 

for(i in 1:length(bathyLat)){
  lat = bathyLat[i]
  bins = data.frame(lat = lat, latBin = SODAlat)
  bins$diffs = abs(bins$lat - bins$latBin) # column showing how far away a lat value is from each SODA bin
  latBin = bins$latBin[bins$diffs == min(bins$diffs)] # Choose the SODA lat bin that the bathymetry lat value is closest to
  latBins[i,] = data.frame(lat=lat, latgrid=latBin)
}

# Merge the lat bins over to the master bathymetry file and then Aggregate depths to the right projection grid
proj.grid <- merge(bathy.grid, latBins, all.x=T, by='lat', sort=F)
proj.grid <- data.table(proj.grid) # easier to do aggregations
proj.grid.ag <- proj.grid[, list(depth = mean(elev), depthrange = (max(elev) - min(elev))), by = list(longrid, latgrid)]

#plot up depth data 
plot(latgrid~longrid, data=proj.grid.ag) # can see I'll need to 'prune' out some lat/lon boxes
levelplot(depth ~ longrid * latgrid, data = proj.grid.ag) 

# Near the coastlines some grid cells occur that have no SODA data
# change names of SODA lat/lon columns and merge, then look for cells with NA for sst
setnames(SODAaleut, c('longrid', 'latgrid', 'sst'))
proj.grid.aleut <- merge(proj.grid.ag, SODAaleut, all.x=T, by=c('longrid', 'latgrid'), sort=F)
summary(proj.grid.aleut) # there are 25 NAs for sst

proj.grid.aleut <- proj.grid.aleut[!is.na(proj.grid.aleut$sst),]
levelplot(depth ~ longrid * latgrid, data = proj.grid.aleut)# much improved coastlines
# Do some pruning
proj.grid.aleut <- proj.grid.aleut[!(longrid < -196 | latgrid > 56)]
plot(latgrid~longrid, data=proj.grid.aleut) # can see I'll need to 'prune' out some lat/lon boxes

# Now fill in that last little data gap in the aleutians==============================================================================================================================================

bathy.aleut2 <- nc_open('Rugosity/Rugosity_data/GEBCO_2014_2D_-180.0_46.0_-173.0_57.0.nc')
bathyArray2 <- ncvar_get(bathy.aleut2, 'elevation') 

lon <- ncvar_get(bathy.aleut2, 'lon') 
lat <- ncvar_get(bathy.aleut2, 'lat')
# image(lon, lat, bathyArray, col = rev(brewer.pal(11, "RdBu")))

grid <- expand.grid(lon = lon, lat = lat)

# turn into data.frame
bathyVec <- as.vector(bathyArray2)
bathy.grid <- as.data.frame(cbind(grid, elev = bathyVec))

# First drop obvious rows of non interest, i.e. positive values or too deep
bathy.grid <- bathy.grid[bathy.grid$elev < 0,]
bathy.grid <- bathy.grid[bathy.grid$elev > -401,]
# Create longitude bins with at .25 degree resolution to match SODA
bathy.grid$longrid <- floor(bathy.grid$lon/gridsize)*gridsize + gridsize/2

# Bin to soda lat grid_tougher b/c resolution changes with latitude
# First cut down the SODA grid to the right dimensions in bathygrid to save time
SODAaleut2 <- SODA[SODA$lon < (max(bathy.grid$lon)+.5) & SODA$lon > (min(bathy.grid$lon)-.5),] # added/subtracted the .5 to buffer a bit
SODAaleut2 <- SODAaleut2[SODAaleut2$lat < (max(bathy.grid$lat)+.5) & SODAaleut2$lat > (min(bathy.grid$lat)-.5),] # added/subtracted the .5 to buffer a bit
levelplot(sst ~ longrid * latgrid, data = SODAaleut2) # looks like I converted the longitude correctly
SODAaleut2 <- SODAaleut2[!is.na(SODAaleut2$sst),] # Drop rows with no temp data (ie on land)

#For loop to assign bathymetry values to SODA lat grid bins
SODAlat <- unique(as.vector(SODAaleut2$lat)) # Vector of all the lats in the SODA grid_to be used in a loop
bathyLat <- unique(as.vector(bathy.grid$lat)) # all the lat values that need to be associated with the correct SODA bins with a forloop

latBins <- data.frame(lat=numeric(), latgrid=numeric()) 

for(i in 1:length(bathyLat)){
  lat = bathyLat[i]
  bins = data.frame(lat = lat, latBin = SODAlat)
  bins$diffs = abs(bins$lat - bins$latBin) # column showing how far away a lat value is from each SODA bin
  latBin = bins$latBin[bins$diffs == min(bins$diffs)] # Choose the SODA lat bin that the bathymetry lat value is closest to
  latBins[i,] = data.frame(lat=lat, latgrid=latBin)
}

# Merge the lat bins over to the master bathymetry file and then Aggregate depths to the right projection grid
proj.grid <- merge(bathy.grid, latBins, all.x=T, by='lat', sort=F)
proj.grid <- data.table(proj.grid) # easier to do aggregations
proj.grid.ag <- proj.grid[, list(depth = mean(elev), depthrange = (max(elev) - min(elev))), by = list(longrid, latgrid)]

#plot up depth data 
plot(latgrid~longrid, data=proj.grid.ag) # can see I'll need to 'prune' out some lat/lon boxes
levelplot(depth ~ longrid * latgrid, data = proj.grid.ag) 

# Near the coastlines some grid cells occur that have no SODA data
# change names of SODA lat/lon columns and merge, then look for cells with NA for sst
setnames(SODAaleut2, c('longrid', 'latgrid', 'sst'))
proj.grid.aleut2 <- merge(proj.grid.ag, SODAaleut2, all.x=T, by=c('longrid', 'latgrid'), sort=F)
summary(proj.grid.aleut2) # there are 1 NAs for sst

proj.grid.aleut2 <- proj.grid.aleut2[!is.na(proj.grid.aleut2$sst),]
levelplot(depth ~ longrid * latgrid, data = proj.grid.aleut2)# much improved coastlines
# Do some pruning
proj.grid.aleut2 <- proj.grid.aleut2[!(longrid > -176 | latgrid > 54)]
plot(latgrid~longrid, data=proj.grid.aleut2) # can see I'll need to 'prune' out some lat/lon boxes

rm(bathy.grid, bathyArray, bins, grid, latBins, proj.grid, proj.grid.ag, SODA, SODAaleut, SODAaleut2, bathy, bathyLat, bathyVec, gridsize, lon, SODAlat, i, lat, latBin)

# ==============================================================================================================================================
# Now merge all regions and expand the grid_will need to drop a handful of Aleutian duplicates
# ==============================================================================================================================================
#load('projectionGridWest_Jan9_2017.RData')
#load('projectionGridEast_Jan9_2017.RData')

proj.grid <- data.table(rbind(proj.grid.east, proj.grid.west, proj.grid.aleut2, proj.grid.aleut))
plot(latgrid~longrid, cex=.05, data=proj.grid)
levelplot(depth ~ longrid * latgrid, data = proj.grid) 
# Find duplicates in the Aleutian region
setkey(proj.grid, longrid, latgrid)
proj.grid <- unique(proj.grid)

proj.grid$rowIndex <- 1:nrow(proj.grid)
#drop 'depthrange' and 'sst'
proj.grid$depthrange <- proj.grid$sst <- NULL
rowIndex <- as.vector(1:nrow(proj.grid))
year <- as.vector(2006:2100)
grid <- expand.grid(rowIndex = rowIndex, year = year)
proj.grid <- merge(proj.grid, grid, by="rowIndex", all=T, sort=F)
rm(grid, rowIndex, year)
proj.grid$rowIndex <- NULL
 
save(proj.grid, file='projectionGrid_Jan9_2017.RData')
write.csv(proj.grid, file='projectionGrid_Jan9_2017.csv', row.names=F)

# Add a season to the projection grid
load('projectionGrid_Jan9_2017.RData')
proj.grid$id <- 1:nrow(proj.grid) # add a row ID

# Create grid of row IDs and seasons to add back to proj.grid
newproj <- expand.grid(id=proj.grid$id, season=1:4)
nrow(proj.grid) *4 # same as newproj so seemed to work
newproj <- merge(newproj, proj.grid)
newproj <- newproj[order(newproj$longrid, newproj$latgrid, newproj$year, newproj$season), c('longrid', 'latgrid', 'depth', 'year', 'season')]
save(newproj, file='projectionGrid_Jan14_2017.RData')
write.csv(newproj, file='projectionGrid_Jan14_2017.csv', row.names=F)

# Make some pdfs of the grids etc.
plot(latgrid~longrid, cex=.05, data=proj.grid[year==2006])
levelplot(depth ~ longrid * latgrid, data = proj.grid[year==2006]) 
 
png("figures/projection_gridDepthDist.png")
ggplot(proj.grid[year==2006], aes(x=depth)) + geom_histogram(binwidth = 25)
dev.off()
 