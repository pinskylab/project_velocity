# THIS SCRIPT STILL NEEDS TO BE CLEANED UP

library(data.table)
library(raster)
library(ncdf4)
library(rasterVis)
library(chron)
library(ggplot2)
require(maps)
require(RColorBrewer)
library(geosphere)
library(plyr)
setwd('/Users/jamesmorley/Documents/project_velocity')

# Bring in as raster
rugos.east <- raster('Rugosity/Rugosity_data/GEBCO_2014_2D_-100.1942_23.6311_-42.3301_62.0777.nc')
rugos.east # get raster attributes_shows that it is .00833 resolution, which is 30-second resolution
print(rugos.east) # characteristics of source netCDF file_including useful meta data for writing about etc.
names(rugos.east)
values(rugos.east)[1:10] # Can see values in first 10 cells
res(rugos.east) # get the resolution
plot(rugos.east, main = 'Elevation/bathymetry') # plots data_doesn't take long

# Create new layers for rugosity_do as a RasterStack
tri.east <- terrain(rugos.east, opt='TRI') #TRI (Terrain Ruggedness Index) is the mean of the absolute differences between the value of a cell and the value of its 8 surrounding cells.
rough.east <- terrain(rugos.east, opt='roughness') #Roughness is the difference between the maximum and the minimum value of a cell and its 8 surrounding cells.
bathyEast <- stack(rugos.east, tri.east, rough.east)
nlayers(bathyEast)
# elev <- raster(bathyEast, layer=1) # Would extract one of the layers in a stack if needed

# Convert to a datatable (rasterToPoints creates a matrix)
bathy.grid <- data.table(rasterToPoints(bathyEast)) # all the lat/lon expanded and all three values in the rasterStack included

# Need to quality control to see if the depth values match up when imported as an ncdf file...maybe to see if the lat lon vectors end up identical as they should and depth values correspond with what I saw with the raster import
# bathy <- nc_open('Rugosity/Rugosity_data/GEBCO_2014_2D_-100.1942_23.6311_-42.3301_62.0777.nc')
# print(bathy)
# names(bathy$var) # Elevation is only variable
# bathyArray <- ncvar_get(bathy, 'elevation') 
# dim(bathyArray); dim(bathyEast) # they match up, but are flipped 90 degreesy
# bathyArray[100,240]; bathyArray[240,100]; bathyEast[240,100]; bathyEast[100,240]# the two don't match up, doesn't necessarily mean something is wrong...

# lon <- ncvar_get(bathy, 'lon') 
# lat <- ncvar_get(bathy, 'lat')
# rm(bathy)
# grid <- expand.grid(lon = lon, lat = lat)
# identical(nrow(bathy.grid), nrow(grid)) # True, so dimensions match up_how about specific lat/lon values for depth?

# turn ncdf file into data.frame
# bathyVec <- as.vector(bathyArray)
# bathy.grid.ncdf <- as.data.frame(cbind(grid, elev = bathyVec))

# Compare unique values of lat and lon, to see if everything matches up
setnames(bathy.grid, c('lon', 'lat', 'elev', 'tri', 'roughness'))
# identical(unique(bathy.grid.ncdf$lat), unique(sort(bathy.grid$lat))) # Need to sort bathy.grid as is high to low, still is False_maybe due to sigfigs? Seems to match up
# unique(bathy.grid.ncdf$lat)
# unique(sort(bathy.grid$lat))
# don't include the code here, but found same thing with lon
# Now check some specific lat/lon values, the elev column should match up for both files
# rm(bathyArray, grid, bathyVec, lat, lon)
# Drop positive elevations to make these files more manageable
bathy.grid <- bathy.grid[elev < 0]
# bathy.grid.ncdf <- bathy.grid.ncdf[bathy.grid.ncdf$elev < 0,]

# the lines below I recycled a bunch to point check some specic lat lons
# The convoluted search is needed (< and >) b/c of the number of hidden significant figures
# bathy.grid.ncdf[200000,]
# bathy.grid[16000000]# 
# bathy.grid[lon < -96.40416 & lon > -96.40418 & lat > 26.60416 & lat < 26.60418]
# bathy.grid.ncdf[bathy.grid.ncdf$lon < -96.40416 & bathy.grid.ncdf$lon > -96.40418 & bathy.grid.ncdf$lat > 26.60416 & bathy.grid.ncdf$lat < 26.60418,]
# The above code shows that the datatable with depths and rugosity measures matches up with the basic ncdf file extraction for depth_so everything checks out so far
# rm(bathy.grid.ncdf)

# make some levelplots to check on rugosity see if it makes sense_First drop some of the deepest areas to save room
# Don't want to drop any cells that may appear in climate projection grid which was set to 400m, so I'll be conservative for now
bathy.grid <- bathy.grid[elev > -3000]
rm(bathyEast, rough.east, rugos.east, tri.east)
# Can trim out some areas (hudson bay, great lakes, greenland shelf, etc.)
bathy.grid <- bathy.grid[lon < -43.5 & lat > 24 & lat < 61]
bathy.grid <- bathy.grid[!(lon < -72 & lat > 50)]
bathy.grid <- bathy.grid[!(lon > -50 & lat > 57)]
bathy.grid <- bathy.grid[!(lon > -67 & lat < 37)]
bathy.grid <- bathy.grid[!(lon > -90 & lon < -87 & lat < 25)]
bathy.grid <- bathy.grid[!(lon < -75.5 & lat > 43)]
 
levelplot(elev ~ lon * lat, data = bathy.grid[elev > -1000]) 
nrow(bathy.grid[tri > 300 & elev > -1000]) # 460 of 3.9 million have these high TRI values_in shallower waters anyway
plot(lat ~ lon, data = bathy.grid[tri > 300 & elev > -1000]) # See where these huge rugosity-TRI values are_level plot doesn't work for this sparse data
plot(elev ~ lat, data = bathy.grid[tri > 300 & elev > -1000]) # most (~80%) are in depths below 400m
# Most of these values are clustered around bahamas, with a handful more dispersed along the shelf break all over.
# This plot gives the best look at TRI rugosity
cutpts <- c(0, 0.25, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 5.7)
levelplot(log(tri+1) ~ lon * lat, data = bathy.grid[elev > -1000], at = cutpts, cuts = 11, pretty = T, col.regions = (rev(brewer.pal(11, "RdBu"))))
# Clearly the rugosity data match up well with what I expected, so the above data manipulations seem to work
# MAY NEED TO TRANSFORM RUGOSITY AS I DID IN FIGURE (LOG)_I THINK I NEED TO SEE WHAT TRI VALUES ACTUALL SHOW UP IN THE HAULS FIRST_SO TRANSFORMATION WILL TAKE PLACE IN HABITAT MODEL FITTING STAGE
save(bathy.grid, file='bathy_rugos_grid_Jan19_2017.RData')

# Make high quality pdf
bathy.shall <- bathy.grid[elev > -1000] # Make a separate datatable for the plot
bathy.shall$logTRI <- log(bathy.shall$tri + 1)
bathy.shall$TRIplot <- bathy.shall$logTRI
bathy.shall <- bathy.shall[logTRI > 5.7, TRIplot := 6]# Only 482 points are above a value of 5.7, so this pools them all together 
# Remove a handful of points removed from the shelf or on land
bathy.shall <- bathy.shall[!(lon > -64 & lat < 40)]
bathy.shall <- bathy.shall[!(lat < 34 & lat > 32 & lon > -75)]
bathy.shall <- bathy.shall[!(lat < 46 & lat > 45 & lon < -71)]

# normalize to 0-1
norm01 <- function(x){
  mn <- min(x)
  mx <- max(x)
  return((x-mn)/(mx-mn))
}
colfun <- colorRamp(rev(brewer.pal(11, 'Spectral')))

pdf(width=7, height=7, file='figures/rugosity_map_east.pdf')
par(mai=c(0.7,0.7,0.3, 0.1), las=1, mgp=c(2,1,0))
map(database='world', fill=TRUE, col='dark gray', border=FALSE, xlim=range(bathy.grid$lon), ylim=range(bathy.grid$lat), xlab='Longitude', ylab='Latitude')
map.axes()
title('Rugosity')
points(bathy.shall$lon, bathy.shall$lat, col=rgb(colfun(norm01(bathy.shall$TRIplot)), maxColorValue=255), pch=16, cex=0.05)
legend('left', legend=round(seq(min(bathy.shall$TRIplot), max(bathy.shall$TRIplot), length.out=10),2), col=rgb(colfun(norm01(seq(min(bathy.shall$TRIplot), max(bathy.shall$TRIplot), length.out=10))), maxColorValue=255), pch=16, cex=0.9, title='Log(rugosity)', bty='n')
dev.off()
rm(bathy.shall, cutpts, xlims, ylims)

# ==============================================================================
# Now bring in west coast data, which consists of 1 big, and 2 small bathymetry files
# ==============================================================================

# First the main west coast file
rugos.west <- raster('Rugosity/Rugosity_data/GEBCO_2014_2D_-179.4175_13.6893_-113.3981_70.3883.nc')
tri.west <- terrain(rugos.west, opt='TRI')
rough.west <- terrain(rugos.west, opt='roughness') 
bathyWest <- stack(rugos.west, tri.west, rough.west)
bathy.gridW <- data.table(rasterToPoints(bathyWest)) # all the lat/lon expanded and all three values in the rasterStack included
setnames(bathy.gridW, c('lon', 'lat', 'elev', 'tri', 'roughness'))
bathy.gridW <- bathy.gridW[elev < 0]
bathy.gridW <- bathy.gridW[elev > -3000]
rm(bathyWest, rough.west, rugos.west, tri.west)

# Prune some uneeded areas
levelplot(elev ~ lon * lat, data = bathy.gridW[elev > -1000]) 
bathy.gridW <- bathy.gridW[lon < -116 & lat > 31.5 & lat < 65]
bathy.gridW <- bathy.gridW[!(lon < -140 & lon > -153 & lat > 50 & lat < 55)]
bathy.gridW <- bathy.gridW[!(lon < -130 & lat < 47)]
bathy.gridW <- bathy.gridW[!(lon < -135 & lon > -145 & lat < 54)]
bathy.gridW <- bathy.gridW[!(lon < -130.5 & lat < 50.8)]
bathy.gridW <- bathy.gridW[!(lon < -179 & lat < 55 & lat > 52)]
# still a couple more to get before making a plot

# Second is the western Aleutians file
rugos.aleut <- raster('Rugosity/Rugosity_data/GEBCO_2014_2D_162.0_47.0_180.0_59.0.nc')
tri.aleut <- terrain(rugos.aleut, opt='TRI')
rough.aleut <- terrain(rugos.aleut, opt='roughness') 
bathyAleut <- stack(rugos.aleut, tri.aleut, rough.aleut)
bathy.gridA <- data.table(rasterToPoints(bathyAleut)) # all the lat/lon expanded and all three values in the rasterStack included
setnames(bathy.gridA, c('lon', 'lat', 'elev', 'tri', 'roughness'))
bathy.gridA <- bathy.gridA[elev < 0]
bathy.gridA <- bathy.gridA[elev > -3000]
rm(bathyAleut, rough.aleut, rugos.aleut, tri.aleut)
bathy.gridA$lon <- -180 - (180-bathy.gridA$lon) # Need to convert longitude to different scale here at the edge of the map

# Third is a slot of data in between the easter and western Aleutians files_there is some overlap here
rugos.slot <- raster('Rugosity/Rugosity_data/GEBCO_2014_2D_-180.0_46.0_-173.0_57.0.nc')
tri.slot <- terrain(rugos.slot, opt='TRI')
rough.slot <- terrain(rugos.slot, opt='roughness') 
bathySlot <- stack(rugos.slot, tri.slot, rough.slot)
bathy.gridS <- data.table(rasterToPoints(bathySlot)) # all the lat/lon expanded and all three values in the rasterStack included
setnames(bathy.gridS, c('lon', 'lat', 'elev', 'tri', 'roughness'))
bathy.gridS <- bathy.gridS[elev < 0]
bathy.gridS <- bathy.gridS[elev > -3000]
rm(bathySlot, rough.slot, rugos.slot, tri.slot)

# Combine all of west coast, do some pruning, and make a nice graph
bathy.grid.W <- rbind(bathy.gridW, bathy.gridA, bathy.gridS)
# drop repeat values from some overlap b/w the data regions_for some lat/lon will need to make sure the rows with rugosity values are kept, as margin cells from the parent rugosity rasters have NA
bathy.grid.W <- bathy.grid.W[!is.na(tri)]
setkey(bathy.grid.W, lon, lat)
bathy.grid.W <- unique(bathy.grid.W)
rm(bathy.gridA, bathy.gridS, bathy.gridW)
levelplot(elev ~ lon * lat, data = bathy.grid.W[elev > -1000]) 
bathy.grid.W <- bathy.grid.W[!lon < -196]
bathy.grid.W <- bathy.grid.W[!(lat > 57 & lon < -180)]
bathy.grid.W <- bathy.grid.W[!(lon > -148 & lon < -138 & lat < 57)]
bathy.grid.W <- bathy.grid.W[!(lon > -140 & lon < -137.2 & lat < 55.5)] # this didn't seem to get anything...

save(bathy.grid, bathy.grid.W, file='bathy_rugos_grid_Jan20_2017.RData')

# Plot up the west coast data
bathy.shall <- bathy.grid.W[elev > -1000] # Make a separate datatable for the plot
bathy.shall$logTRI <- log(bathy.shall$tri + 1)
bathy.shall$TRIplot <- bathy.shall$logTRI
bathy.shall <- bathy.shall[logTRI > 5.7, TRIplot := 6]
# Remove a handful of points removed from the shelf or on land
bathy.shall <- bathy.shall[!(lon < -122.5 & lat < 35)]
bathy.shall <- bathy.shall[!(lat < 55 & lon < -136 & lon > -140)]
bathy.shall <- bathy.shall[!(lat < 50 & lon < -130)]
bathy.shall <- bathy.shall[!(lat < 52.5 & lon < -134 & lon > -140)]
bathy.shall <- bathy.shall[!(lon > -118 & lat > 35.7)]
#bathy.shall <- bathy.shall[!(lon > -116 & lat > 34)] # this was an accident, but I don't think it removed anything on the shelf

pdf(width=10, height=7, file='figures/rugosity_map_west.pdf')
par(mai=c(0.7,0.7,0.3, 0.1), las=1, mgp=c(2,1,0))
map(database='world', fill=TRUE, col='dark gray', border=FALSE, xlim=range(bathy.grid.W$lon), ylim=range(bathy.grid.W$lat), xlab='Longitude', ylab='Latitude')
map.axes()
title('Rugosity')
points(bathy.shall$lon, bathy.shall$lat, col=rgb(colfun(norm01(bathy.shall$TRIplot)), maxColorValue=255), pch=16, cex=0.05)
legend('bottomleft', legend=round(seq(min(bathy.shall$TRIplot), max(bathy.shall$TRIplot), length.out=10),2), col=rgb(colfun(norm01(seq(min(bathy.shall$TRIplot), max(bathy.shall$TRIplot), length.out=10))), maxColorValue=255), pch=16, cex=0.9, title='Log(rugosity)', bty='n')
dev.off()
rm(bathy.shall, colfun, norm01)
  
# Combine east and west coast into one file
bath.grid <- rbind(bathy.grid, bathy.grid.W)
rm(bathy.grid, bathy.grid.W)

# ==============================================================================

# correlate the rugosity levels to compare the two metrics, and plot the regression_ultimately only really interested in <400m
bath.sub <- bath.grid[elev > -401] # N = 6,477,921 data points
cor(bath.sub$tri, bath.sub$roughness) # correlation of 0.972 b/w these two metrics
summary(lm(tri~roughness, data=bath.sub)) # r2 of 0.9453_p<.000...0001
cor(log(bath.sub$tri + 1), log(bath.sub$roughness + 1)) # correlation of 0.981 b/w these two metrics
summary(lm(log(tri + 1)~log(roughness + 1), data=bath.sub)) # r2 of 0.962_p<.000...0001
# seems more than reasonable just to go with TRI, which seems more appealing anyway
rm(bath.sub)
bath.grid$roughness <- NULL
save(bath.grid, file='bathy_grid_Jan22_2017.RData')

# ==============================================================================
# Repeat the above code that organizes bathymetry data, but do so at a courser resolution
# ==============================================================================
  
# Bring in as raster
rugos.east <- raster('Rugosity/Rugosity_data/GEBCO_2014_2D_-100.1942_23.6311_-42.3301_62.0777.nc')
res(rugos.east) # 0.008333333 0.008333333
dim(rugos.east) # 4615 6945 1
range(rugos.east) # see the start and end values for lat and lon, want these to be even values to make nice grid
ymax(rugos.east) # 62.08333_can't redifine ymax as it changes the resolution while keeping same number of rows
xmin(rugos.east) # -100.2
rugos.east <- crop(rugos.east, extent(rugos.east, 5, 4615, 19, 6945)) # rerunning the above code shows this worked to get lat and lon at same decimal to make even grid with the aggregate below

# Create new layers for rugosity_do as a RasterStack
tri.east <- terrain(rugos.east, opt='TRI') #TRI (Terrain Ruggedness Index) is the mean of the absolute differences between the value of a cell and the value of its 8 surrounding cells.
rugos.east[rugos.east > 0] <- NA # remove the positive (terrestrial) values
tri.trimmed <- mask(tri.east, mask = rugos.east) # Remove rugosity values that correspond to terrestrial values
bathyEast <- stack(rugos.east, tri.trimmed)

# Aggregate to a courser resolution_1/20 of a degree lat and lon
bathyAgg <- aggregate(bathyEast, fact=6, fun=mean, expand=FALSE, na.rm=TRUE)
dim(bathyEast); res(bathyEast) # 4611 6927    2; 0.008333333 0.008333333
dim(bathyAgg); res(bathyAgg) #768 1154    2; 0.05 0.05
plot(bathyAgg[[1]])
plot(bathyAgg[[2]])
  
# Convert to a datatable (rasterToPoints creates a matrix)
bathy.grid <- data.table(rasterToPoints(bathyAgg)) # all the lat/lon expanded and all three values in the rasterStack included_when all values are NA they are dropped?
setnames(bathy.grid, c('lon', 'lat', 'elev', 'tri'))
# I won't drop any cells here b/c the file is not too big_I'll let the climate projection grid decide what stays
sort(unique(bathy.grid$lon)) # nice even grid
sort(unique(bathy.grid$lat)) # nice even grid_matches with lon
rm(bathyAgg, bathyEast, rugos.east, tri.east, tri.trimmed)

cutpts <- c(0, 0.25, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6)
levelplot(log(tri+1) ~ lon * lat, data = bathy.grid, at = cutpts, cuts = 11, pretty = T, col.regions = (rev(brewer.pal(11, "RdBu"))))

# ==============================================================================
# Now bring in west coast data, which consists of 1 big, and 2 small bathymetry files
# ==============================================================================

# First the main west coast file
rugos.west <- raster('Rugosity/Rugosity_data/GEBCO_2014_2D_-179.4175_13.6893_-113.3981_70.3883.nc')
res(rugos.west) # 0.008333333 0.008333333
dim(rugos.west) # 6805 7924    1
range(rugos.west) # see the start and end values for lat and lon, want these to be even values to make nice grid
ymax(rugos.west) # 70.39167
xmin(rugos.west) # -179.425 # need to trim these so they are both with a 5 in second decimal
rugos.west <- crop(rugos.west, extent(rugos.west, 6, 6805, 10, 7924)) # rerunning the above code shows this worked to get lat and lon at same decimal to make even grid with the aggregate below

tri.west <- terrain(rugos.west, opt='TRI')
rugos.west[rugos.west > 0] <- NA # remove the positive (terrestrial) values
tri.trimmed <- mask(tri.west, mask = rugos.west) # Remove rugosity values that correspond to terrestrial values
bathyWest <- stack(rugos.west, tri.trimmed)

bathyAgg <- aggregate(bathyWest, fact=6, fun=mean, expand=FALSE, na.rm=TRUE)
dim(bathyWest); res(bathyWest) # 6800 7915    2; 0.008333333 0.008333333
dim(bathyAgg); res(bathyAgg) #1133 1319    2; 0.05 0.05
plot(bathyAgg[[1]])
plot(bathyAgg[[2]])

bathy.gridW <- data.table(rasterToPoints(bathyAgg)) # all the lat/lon expanded and all three values in the rasterStack included
setnames(bathy.gridW, c('lon', 'lat', 'elev', 'tri'))
sort(unique(bathy.gridW$lon))
sort(unique(bathy.gridW$lat))
rm(bathyWest, rugos.west, tri.west, tri.trimmed, bathyAgg)

# Second is the western Aleutians file
rugos.aleut <- raster('Rugosity/Rugosity_data/GEBCO_2014_2D_162.0_47.0_180.0_59.0.nc')
res(rugos.aleut) # 0.008333333 0.008333333
dim(rugos.aleut) # 1441 2160    1
range(rugos.aleut) # see the start and end values for lat and lon, want these to be even values to make nice grid
ymax(rugos.aleut) #  59.00833
xmin(rugos.aleut) # 162 need to trim these so they are both with a 5 or a 0 in second decimal
rugos.aleut <- crop(rugos.aleut, extent(rugos.aleut, 2, 1441, 1, 2160)) # rerunning the above code shows this worked to get lat and lon at same decimal to make even grid with the aggregate below

tri.aleut <- terrain(rugos.aleut, opt='TRI')
rugos.aleut[rugos.aleut > 0] <- NA # remove the positive (terrestrial) values
tri.trimmed <- mask(tri.aleut, mask = rugos.aleut) # Remove rugosity values that correspond to terrestrial values
bathyAleut <- stack(rugos.aleut, tri.trimmed)

bathyAgg <- aggregate(bathyAleut, fact=6, fun=mean, expand=FALSE, na.rm=TRUE)
dim(bathyAleut); res(bathyAleut) # 1440 2160    2; 0.008333333 0.008333333
dim(bathyAgg); res(bathyAgg) #240 360   2; 0.05 0.05
plot(bathyAgg[[1]])
plot(bathyAgg[[2]])

bathy.gridA <- data.table(rasterToPoints(bathyAgg)) # all the lat/lon expanded and all three values in the rasterStack included
setnames(bathy.gridA, c('lon', 'lat', 'elev', 'tri'))
sort(unique(bathy.gridA$lon))
sort(unique(bathy.gridA$lat))

rm(bathyAleut, rugos.aleut, tri.aleut, tri.trimmed, bathyAgg)
bathy.gridA$lon <- -180 - (180-bathy.gridA$lon) # Need to convert longitude to different scale here at the edge of the map

# Third is a slot of data in between the easter and western Aleutians files_there is some overlap here
rugos.slot <- raster('Rugosity/Rugosity_data/GEBCO_2014_2D_-180.0_46.0_-173.0_57.0.nc')
res(rugos.slot) # 0.008333333 0.008333333
dim(rugos.slot) # 1321  841    1
range(rugos.slot) # see the start and end values for lat and lon, want these to be even values to make nice grid
ymax(rugos.slot) #  57.00833
xmin(rugos.slot) # -180 need to trim these so they are both with a 5 or a 0 in second decimal
rugos.slot <- crop(rugos.slot, extent(rugos.slot, 2, 1321, 1, 841)) # rerunning the above code shows this worked to get lat and lon at same decimal to make even grid with the aggregate below

tri.slot <- terrain(rugos.slot, opt='TRI')
rugos.slot[rugos.slot > 0] <- NA # remove the positive (terrestrial) values
tri.trimmed <- mask(tri.slot, mask = rugos.slot) # Remove rugosity values that correspond to terrestrial values
bathySlot <- stack(rugos.slot, tri.trimmed)

bathyAgg <- aggregate(bathySlot, fact=6, fun=mean, expand=FALSE, na.rm=TRUE)
dim(bathySlot); res(bathySlot) # 1320  841    2; 0.008333333 0.008333333
dim(bathyAgg); res(bathyAgg) #220 140   2; 0.05 0.05
plot(bathyAgg[[1]])
plot(bathyAgg[[2]])

bathy.gridS <- data.table(rasterToPoints(bathyAgg)) # all the lat/lon expanded and all three values in the rasterStack included
setnames(bathy.gridS, c('lon', 'lat', 'elev', 'tri'))
sort(unique(bathy.gridS$lon))
sort(unique(bathy.gridS$lat))
rm(bathySlot, rugos.slot, tri.slot, tri.trimmed, bathyAgg)

# Combine all of west coast
bathy.grid.W <- rbind(bathy.gridW, bathy.gridA, bathy.gridS)
# drop repeat values from some overlap b/w the data regions_for some lat/lon will need to make sure the rows with rugosity values are kept, as margin cells from the parent rugosity rasters have NA
bathy.grid.W <- bathy.grid.W[!is.na(tri)]
setkey(bathy.grid.W, lon, lat)
bathy.grid.W <- unique(bathy.grid.W)
rm(bathy.gridA, bathy.gridS, bathy.gridW)
levelplot(elev ~ lon * lat, data = bathy.grid.W) 

# Take a closer look at the 'slot' area of the Aleutians to make sure there is full coverage
levelplot(elev ~ lon * lat, data = bathy.grid.W[lon < -179 & lon > -180.5 & lat > 45 & lat < 60]) 
plot(lat~lon, cex=.1, data = bathy.grid.W[lon < -179 & lon > -180.5 & lat > 45 & lat < 60]) # No evidence of a gap in the grid from binding these together

# Combine east and west coast into one file
bath.grid <- rbind(bathy.grid, bathy.grid.W)
rm(bathy.grid, bathy.grid.W)
save(bath.grid, file='bathy_course_untrimmed_Feb8_2017.RData') # Note that this was originally 'bathy_grid_course_Feb7_2017.RData', but I had to redo as I accidentally replaced it below with the file used for projections

bath.grid$logRugos <- log(bath.grid$tri + 1)
cutpts <- c(0, 0.25, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 5.7)
pdf(width=15, height=8, file='/Users/abigail/Desktop/rugosity_course_untrimmed.pdf')
par(mai=c(0.7,0.7,0.3, 0.1), las=1, mgp=c(2,1,0))
levelplot(logRugos ~ lon * lat, data = bath.grid, at = cutpts, cuts = 11, pretty = T, col.regions = (rev(brewer.pal(11, "RdBu"))))
#map(database='world', fill=TRUE, col='dark gray', border=FALSE, xlim=range(west$lonClimgrid), ylim=range(west$latClimgrid), xlab='Longitude', ylab='Latitude')
#points(latBathgrid~lonBathgrid, cex = .05, data=west) 
dev.off()
  
# ==============================================================================
# # Assign lat/lon bins from bathymetry to hauls
# # Note that used this code for course and fine bathymetry resolution, currently code specific to fine res is muted
# ==============================================================================

#load('bathy_grid_Jan22_2017.RData')
load('bathy_course_untrimmed_Feb8_2017.RData')
# Are there obvious-repeated bins within each degree of lat or lon to match up with?
sort(unique(bath.grid$lat)) # 10 values w/n each latitude, and they repeat_an even grid
#abc <- unique(data.frame(lon = bath.grid$lon)) # almost same deal with lon, but sigfigs of the decimal bins differ here depending if it's a two or three digit lon (e.g. -78 vs. -120)
#rm(abc)
# I'll try binning hauls using the gridsize method
gridsize <- 0.05 
# specify the lat/lon grid, b/c there are some very small-hidden decimal places
bath.grid$latgrid <- signif(floor(bath.grid$lat/gridsize)*gridsize + gridsize/2, digits = 6) 
bath.grid$longrid <- signif(floor(bath.grid$lon/gridsize)*gridsize + gridsize/2, digits = 7) 
summary(bath.grid$latgrid - bath.grid$lat); summary(bath.grid$longrid - bath.grid$lon) # very small differences_definitely no grid cell shifting

load('haulInfo_Dec26_2016.RData') # bring in haul info file
# originally I used a for loop below to assign grid values to hauls, but now I don't think that is necessary
hauls$latgrid <- signif(floor(hauls$lat/gridsize)*gridsize + gridsize/2, digits = 6)
hauls$longrid <- signif(floor(hauls$lon/gridsize)*gridsize + gridsize/2, digits = 7) 
summary(hauls$latgrid - hauls$lat); summary(hauls$longrid - hauls$lon) # actual haul values never more than 0.025 away from grid, makes sense for a 0.05 res
# Now merge in the bathymetry data
abc <- copy(bath.grid)
abc$lat <- abc$lon <- NULL # don't need these for the merge, only the grid values without all the sigfigs
haulsTest <- merge(hauls, abc, all.x = TRUE, by=c('latgrid', 'longrid'), sort = FALSE)
haulsTest[(is.na(haulsTest$tri) & !is.na(haulsTest$lat)),] # single haul that wasn't present in bathymetry grid
bath.grid[latgrid == 57.525 & longrid == -152.675] # no value for this in the bathymetry grid_based on figure below, is reasonable to assign rugosity value one click to the south
# where is this haul?
pdf(width=13, height=7, file='/Users/abigail/Desktop/trawl_coords.pdf')
par(mai=c(0.7,0.7,0.3, 0.1), las=1, mgp=c(2,1,0))
map(database='world', fill=TRUE, col='dark gray', border=FALSE, xlim=range(haulsTest$longrid, na.rm=T), ylim=range(haulsTest$latgrid, na.rm=T), xlab='Longitude', ylab='Latitude')
points(lat~lon, cex=.1, data=haulsTest[!is.na(haulsTest$tri),])
points(lat~lon, cex=.3, col='red', data=haulsTest[(is.na(haulsTest$tri) & !is.na(haulsTest$lat)),])
dev.off()

haulsTest$tri[haulsTest$haulid == '024-198701-036'] <- 40.2 # assign rugosity value from grid cell to the south
plot(depth~abs(elev), cex=.1, data=haulsTest)
summary(lm(depth~abs(elev), data=haulsTest)) # slope = 0.9037606; intercept = 6.7194116; p < .00001 (~15 zeros?); r2 = 0.92
# seems reasonable to replace missing depth values with elev if we use depth
levelplot(tri ~ lon * lat, data = haulsTest) 

save(haulsTest, file='haulInfo_course_Feb8_2017.RData') # save the hauls file
# LOAD UP THE FINE SCALE BATHYMETRY HAULS, AND CORRELATE THE TWO SCALES OF RUGOSITY
haulsCourse <- haulsTest # rename to bring in the fine resolution file
rm(abc, hauls, gridsize, haulsTest)
 
load('haulInfo_Jan27_2017.RData') # load the hauls file with the finer resolution
course <- data.frame(haulid = haulsCourse$haulid, rugosCourse = haulsCourse$tri)
fine <- data.frame(haulid = haulsTest$haulid, rugosFine = haulsTest$rugosity)
rugos <- merge(course, fine, all=T, by='haulid')
plot(rugosFine~rugosCourse, cex=.1, data=rugos)
summary(lm(rugosFine~rugosCourse, data=rugos)) # slope=0.872(0.0012); intercep=0.595; p(slope) < 2.2e-16; r2= 0.774; df=148533 (n=148535)
rm(course, fine, haulsTest, rugos)

# THE CODE BELOW WAS TO GRID HAULS TO THE FINE (1/120) RESOLUTION_SO NO LONGER NEEDED_PLUS I THINK IT COULD HAVE BEEN DONE MORE SIMPLY (LIKE ABOVE)
# make longitude bins
haulLon <- as.vector(hauls$lon) # Vector of all the lons in hauls_don't use 'unique' as I can just cbind the result of the forloop to hauls
bathyLon <- unique(as.vector(bath.grid$lon)) 
lonBins <- data.frame(lon=numeric(), longrid=numeric()) 

for(i in 1:length(haulLon)){
  lon = haulLon[i]
  bins = data.frame(lon = lon, lonBin = bathyLon)
  bins$diffs = abs(bins$lon - bins$lonBin) # column showing how far away a lon value is from each proj.grid bin
  min = min(bins$diffs) # What is the minimum value for differences between the haul lon and the bathy lon grid
  minBins = bins[bins$diffs == min,] # Make a small dataframe with all the bathygrids that equal the minimum difference, there should only be 1 or 2 values for this_this is necessary as there were a handful of ties occuring b/c of the course res of some of the haul lat/lon values
  lonBin = ifelse(nrow(minBins) == 1, minBins$lonBin, max(minBins$lonBin))# If there is a tie, this chooses the higher longitude value, somewhat arbitrary_I did this b/c some really shallow effort on east coast, so less values will be pushed to land bath grid values, maybe
  lonBins[i,] = data.frame(lon=lon, longrid=lonBin)
}

# check to see if order is maintained with hauls and then cbind it, in case some rounding difficulties prevent a good merge
identical(haulLon, lonBins$lon) # TRUE, so order is conserved as expected so can cbind rather than merge back to hauls
setnames(lonBins, c('lonVec', 'longrid'))
hauls <- cbind(hauls, lonBins)
identical(hauls$lon, hauls$lonVec) # TRUE, again order is preserved
hauls$lonVec <- NULL
hauls$diff <- hauls$lon - hauls$longrid # Everything looks good, all 'diff' values are within 0.00417, which is half of the bathymetry resolution of 0.00833
rm(bathyLon, haulLon, lon, lonBin)

# latitude
haulLat <- as.vector(hauls$lat) # Vector of all the lats in hauls_don't use 'unique' as I can just cbind the result of the forloop to hauls
bathyLat <- unique(as.vector(bath.grid$lat)) 
latBins <- data.frame(lat=numeric(), latgrid=numeric()) 

for(i in 1:length(haulLat)){
  lat = haulLat[i]
  bins = data.frame(lat = lat, latBin = bathyLat)
  bins$diffs = abs(bins$lat - bins$latBin) # column showing how far away a lat value is from each proj.grid bin
  min = min(bins$diffs) # What is the minimum value for differences between the haul lat and the bathy lat grid
  minBins = bins[bins$diffs == min,] # Make a small dataframe with all the bathygrids that equal the minimum difference, there should only be 1 or 2 values for this
  latBin = ifelse(nrow(minBins) == 1, minBins$latBin, min(minBins$latBin))# If there is a tie, this chooses the lower latitude value, somewhat arbitrary, but on average on both coasts this will push further from shore reducing values forced to land bins
  latBins[i,] = data.frame(lat=lat, latgrid=latBin)
  print(i)
}
 
# check to see if order is maintained with hauls and then cbind it, in case some rounding difficulties prevent a good merge
identical(haulLat, latBins$lat) # TRUE, so order is conserved as expected so can cbind rather than merge back to hauls
setnames(latBins, c('latVec', 'latgrid'))
hauls <- cbind(hauls, latBins)
identical(hauls$lat, hauls$latVec) # TRUE, again order is preserved
hauls$latVec <- NULL
hauls$diff <- hauls$lat - hauls$latgrid # also looks good, difference between actual and grid all range b/w half the bathy resolution above or below
hauls$diff <- NULL
rm(bathyLat, haulLat, i, lat, latBin, min, minBins, bins, latBins, lonBins)

# ============================================================================================================================

# Need to merge in bath.grid based on latgrid longrid need to adjust some names first
setnames(bath.grid, c('longrid', 'latgrid', 'bathgriddepth', 'rugosity'))
haulsTest <- data.table(merge(hauls, bath.grid, by=c('latgrid', 'longrid'), all.x=T, sort=F))

# Do some plots to see how things look, and make some temporary columns to show that things seem to merge well, plus compare with original hauls to make sure things are right.
# Do some spot checking_I recycled (changed the number) this code a bunch to see if the merge generally worked as expected
hauls[7000,]
haulsTest[haulsTest$haulid == hauls$haulid[7000],] # this gives row x on original file, and corresponding row on haulsTest based on haulid
bath.grid[latgrid == hauls$latgrid[7000] & longrid == hauls$longrid[7000]] # everything checks out here
identical(length(unique(hauls$haulid)), length(unique(haulsTest$haulid))) # TRUE, so all hauls still present and uniquely represented
summary(hauls); summary(haulsTest) # It seems that only 59 hauls didn't match up with a bathy.grid value
haulsTest[is.na(rugosity) & !is.na(lat)]

# ~59 NAs for rugosity (+ 69 more that were hauls with no lat/lon), which suggests some lat/lon combinations were not in bath grid
cde <- haulsTest[is.na(rugosity) & !is.na(lat)]# These rows have all NAs for bathymetry stuff, but they did have lat/longrids...
# maybe these are unique combos of lat/lon that are not actually on the bathy grid???
bath.grid[latgrid == cde$latgrid[40] & longrid == cde$longrid[40]] # none seem to be on bathy.grid as expected
# Find out where these points are
pdf(width=11, height=7, file='/Users/jamesmorley/Desktop/test_graph.pdf')
par(mai=c(0.7,0.7,0.3, 0.1), las=1, mgp=c(2,1,0))
map(database='world', fill=TRUE, col='dark gray', border=FALSE, xlim=range(bath.grid$longrid), ylim=range(bath.grid$latgrid), xlab='Longitude', ylab='Latitude')
points(lat~lon, cex=.1, data=haulsTest[!is.na(rugosity)])
points(lat~lon, cex=.1, col='red', data=cde)
dev.off()
# All these hauls are right next to shore, so must have snapped to a bathy grid that happened to be on land and dropped earlier
# For those hauls only I'll separately identify the bathy grid cell based on a loop using lat and lon coords

haulCoords <- data.frame(lat = cde$lat, lon = cde$lon) # dataframe of all lat lon coordinates from the hauls to be assigned_don't use 'unique' as I can just cbind the result of the forloop to hauls
bathCoords <- unique(data.frame(latgrid = bath.grid$latgrid, longrid = bath.grid$longrid)) 
# Splitting the east from west I think will save some time in the loop
bathCoordsE <- bathCoords[bathCoords$longrid > -105,]
bathCoordsW <- bathCoords[bathCoords$longrid < -105,]
grid <- data.frame(lat=numeric(), lon=numeric(), latgrid=numeric(), longrid=numeric()) 

diff <- function(data) {
  dist = with(data, distHaversine(p1=c(lon, lat), p2=c(lonBin, latBin)))
}

for(i in 1:nrow(haulCoords)){
  coord = haulCoords[i,]
  if(coord$lon > -105)
  {
    bins = data.frame(lat=coord$lat, lon=coord$lon, latBin=bathCoordsE$latgrid, lonBin=bathCoordsE$longrid)
  }  
  if(coord$lon < -105)
  {
    bins = data.frame(lat=coord$lat, lon=coord$lon, latBin=bathCoordsW$latgrid, lonBin=bathCoordsW$longrid)
  }
  bins = bins[bins$latBin < coord$lat + 0.2 & bins$latBin > coord$lat - 0.2,] # to restrict the data that the following function works on
  bins = adply(bins, .margins=1, .fun=diff) # This still takes awhile, even with only <90k of data
  min = min(bins$V1) # What is the minimum value for differences between the haul lat and the bathy lat grid
  minBins = bins[bins$V1 == min,] # Make a small dataframe with all the bathygrids that equal the minimum difference, there should only be 1 or 2 values for this
  coordLat = ifelse(nrow(minBins) == 1, minBins$latBin, max(minBins$latBin))# If there is a tie, this chooses the lower latitude value, somewhat arbitrary, but on average on both coasts this will push further from shore reducing values forced to land bins
  coordLon = ifelse(nrow(minBins) == 1, minBins$lonBin, max(minBins$lonBin))
  grid[i,] = data.frame(lat=coord$lat, lon=coord$lon, latgrid=coordLat, longrid=coordLon)
  print(i)
}
  
# cbind grid to cde_I checked and order was preserved
cde$latgrid <- cde$longrid <- cde$rugosity <- cde$bathgriddepth <- NULL # drop old values that were off grid
grid$lat <- grid$lon <- NULL
cde <- cbind(cde, grid)
cde <- merge(cde, bath.grid, by=c('latgrid', 'longrid'), all.x=T, sort=F)# merge bath.grid to cde to get rugosity data
# delete cde from haulsTest and rbind cde
haulsTest <- haulsTest[!(is.na(rugosity) & !is.na(lat))] # remove the hauls that I'm about to readd
setcolorder(cde, c("latgrid", "longrid","haulid","region","year", "yearsurv","month", "stratum","lat","lon","depth","bottemp","surftemp_orig", "adjSST", "surftemp","bathgriddepth", "rugosity"))
haulsTest <- rbind(haulsTest, cde)
length(unique(haulsTest$haulid))# make sure all haulid's still unique_yes
summary(haulsTest$lat - haulsTest$latgrid); summary(haulsTest$lon - haulsTest$longrid) # looks good
save(haulsTest, file='haulInfo_Jan27_2017.RData') # save the hauls file
rm(bathCoords, bathCoordsE, bathCoordsW, bins, cde, coord, grid, haulCoords, minBins, coordLat, coordLon, i, min, diff)

pdf(width=11, height=7, file='/Users/jamesmorley/Desktop/trawl_coords.pdf')
par(mai=c(0.7,0.7,0.3, 0.1), las=1, mgp=c(2,1,0))
map(database='world', fill=TRUE, col='dark gray', border=FALSE, xlim=range(bath.grid$longrid), ylim=range(bath.grid$latgrid), xlab='Longitude', ylab='Latitude')
points(lat~lon, cex=.1, data=haulsTest)
dev.off()
 
# Test with actual trawl data depths vs. elevation column
summary(lm(depth~abs(bathgriddepth), data=haulsTest)) # r2 = 0.95
pdf(width=7, height=7, file='/Users/jamesmorley/Desktop/obs_v_bathy_depth.pdf')
plot(depth~abs(bathgriddepth), data=haulsTest)
dev.off()
# How does rugosity relate to depth 
summary(lm(rugosity~abs(bathgriddepth), data=bath.grid)) # p < 0.00...01_r2 = 0.11_too many points for a figure

# ==============================================================================
# Creating a bathymetry grid for use with projections
# Note, I reused this section both for the fine and course bathymetry grid. So currently, lines for the fine resulotion are muted
# ==============================================================================

#load('bathy_grid_Jan22_2017.RData') # fine resolution
load('bathy_grid_course_Feb7_2017.RData') # aggregated to 1/20 of a degree, courser res
# For the bathymetry grid, need to limit the extent of the file to match the climate projection grid anything that doesn't fit precisely into a projection.grid cell is dropped
load('projectionGrid_Jan9_2017.RData')
proj.grid <- proj.grid[year == 2006] # Only need one 'copy' of the lat lons
proj.grid$year <- NULL
# Need to assign grid cell values to bath.grid, based on the proj.grid format, which is SODA
# This is easy for longitude, as the resolution is 0.25 degrees throughout
gridsize <- 0.25 
bath.grid$longrid <- floor(bath.grid$lon/gridsize)*gridsize + gridsize/2
rm(gridsize)
  
#For loop to assign bathymetry values to projection grid lat grid bins_resolution changes with latitude
# bath.grid$latAdj <- signif(bath.grid$lat, digits = 7) # Only necessary for fine res_specify decimal places otherwise things get confused later_this yields 5 decimal places, which is more than adequate for the 0.0083 resolution
projLat <- unique(as.vector(proj.grid$latgrid)) # Vector of all the lats in the SODA grid_to be used in a loop
#bathyLat <- unique(as.vector(bath.grid$latAdj)) # all the lat values that need to be associated with the correct SODA bins with a forloop
bathyLat <- unique(as.vector(bath.grid$lat))
#latBins <- data.frame(latAdj=numeric(), latgrid=numeric()) 
latBins <- data.frame(lat=numeric(), latgrid=numeric())

for(i in 1:length(bathyLat)){
  lat = bathyLat[i]
  bins = data.frame(lat = lat, latBin = projLat)
  bins$diffs = abs(bins$lat - bins$latBin) # column showing how far away a lat value is from each proj.grid bin
  latBin = bins$latBin[bins$diffs == min(bins$diffs)] # Choose the proj.grid lat bin that the bathymetry lat value is closest to
  #latBins[i,] = data.frame(latAdj=lat, latgrid=latBin)
  latBins[i,] = data.frame(lat=lat, latgrid=latBin)
}

# bath.grid <- merge(bath.grid, latBins, all.x=T, by='latAdj', sort=F) # merge in latBins to bath.grid
bath.gridB <- as.data.frame(bath.grid) # had issues with this merge as a data.table, not clear why
bath.gridB <- merge(bath.gridB, latBins, all.x=T, by='lat', sort=F) # merge in latBins to bath.grid
summary(bath.gridB$lon - bath.gridB$longrid)
summary(bath.gridB$lat - bath.gridB$latgrid) # some are off, as the bathymetry grid extended beyond the projection grid for latitude_correction for this:

#summary(latBins$latAdj - latBins$latgrid) # makes sense, all points of latitude w/n 0.12 degrees
summary(latBins$lat - latBins$latgrid) # some are way off, this is expected as bath.grid is greater geographic coverage
# Need to figure out what the northern/southern limits of proj.grid should be, based on the top and bottom rows
sort(projLat) # did some simple math here to determine the min and max values of latitude that should be included
# min value in bath.grid should be 24.0145_max should be 64.9734
bath.gridB <- bath.gridB[bath.gridB$lat < 64.9734 & bath.gridB$lat > 24.0145,]
summary(bath.gridB$lat - bath.gridB$latgrid) # all values are close to grid values now

# Need to limit bath.grid to latgrid/longrid that overlap with proj.grid
proj.grid.bath <- merge(proj.grid, bath.gridB, all.x=T, by = c('longrid', 'latgrid'))
nrow(unique(proj.grid.bath, by=c('longrid', 'latgrid')))# 13,637; this makes sense as is same value as proj.grid
# also extracting the oringal bathymetry lon and lat should show all rows are unique
nrow(unique(data.frame(lon = proj.grid.bath$lon, lat = proj.grid.bath$lat))) # = 199,246, so all rows are unique
# proj.grid.bath$lat <- NULL # Stick with 'latAdj', the rounded version, there is a negligible (.000001 degree) difference
setnames(proj.grid.bath, c('lonClimgrid', 'latClimgrid', 'meandepth', 'latBathgrid', 'lonBathgrid', 'depth', 'rugosity')) # Note that 'meandepth' is not fully accurate, as for the projection grid I culled out deep values; it could easily be redone here accurately as each grid cell has all bathymetry values
proj.grid.bath[is.na(lonBathgrid)] # there is 5 climate grid cells that have no bathymetry....
# ...must've been when aggregating bathymetry, in some nearshore areas this may have led to previously sparse areas_plot up to see_they are in very shallow areas and ok to drop these cells
proj.grid.bath <- proj.grid.bath[!is.na(lonBathgrid)]
proj.grid.bath$meandepth <- NULL # drop this as it's a bit misleading (a carryover column from developing projection grid)
 
#save(proj.grid.bath, file='bathy_grid_Jan30_2017.RData') # final bathymetry-rugosity grid for projections!
save(proj.grid.bath, file='bathy_grid_course_Feb7_2017.RData')
rm(bins, latBins, bathyLat, i, lat, latBin, projLat)

proj.grid.bath$logRugos <- log(proj.grid.bath$rugosity + 1)
east <- proj.grid.bath[lonClimgrid > -105]
west <- proj.grid.bath[lonClimgrid < -105]

cutpts <- c(0, 0.25, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 5.5)
# Plot up all bathymetry grid points east coast then west
pdf(width=7, height=7, file='/Users/abigail/Desktop/rugosity_course_east.pdf')
par(mai=c(0.7,0.7,0.3, 0.1), las=1, mgp=c(2,1,0))
#map(database='world', fill=TRUE, col='dark gray', border=FALSE, xlim=range(east$lonClimgrid), ylim=range(east$latClimgrid), xlab='Longitude', ylab='Latitude')
levelplot(logRugos ~ lonBathgrid * latBathgrid, data = east, at = cutpts, cuts = 11, pretty = T, col.regions = (rev(brewer.pal(11, "RdBu"))))
#points(latBathgrid~lonBathgrid, cex = .05, data=east) 
dev.off()

pdf(width=9, height=7, file='/Users/abigail/Desktop/rugosity_course_west.pdf')
par(mai=c(0.7,0.7,0.3, 0.1), las=1, mgp=c(2,1,0))
levelplot(logRugos ~ lonBathgrid * latBathgrid, data = west, at = cutpts, cuts = 11, pretty = T, col.regions = (rev(brewer.pal(11, "RdBu"))))
#map(database='world', fill=TRUE, col='dark gray', border=FALSE, xlim=range(west$lonClimgrid), ylim=range(west$latClimgrid), xlab='Longitude', ylab='Latitude')
#points(latBathgrid~lonBathgrid, cex = .05, data=west) 
dev.off()

pdf(width=15, height=7, file='/Users/abigail/Desktop/rugosity_course_NorthAmer.pdf')
par(mai=c(0.7,0.7,0.3, 0.1), las=1, mgp=c(2,1,0))
levelplot(logRugos ~ lonBathgrid * latBathgrid, data = proj.grid.bath, at = cutpts, cuts = 11, pretty = T, col.regions = (rev(brewer.pal(11, "RdBu"))))
#map(database='world', fill=TRUE, col='dark gray', border=FALSE, xlim=range(west$lonClimgrid), ylim=range(west$latClimgrid), xlab='Longitude', ylab='Latitude')
#points(latBathgrid~lonBathgrid, cex = .05, data=west) 
dev.off()
 
# ==============================================================================
# Bring sediment data into the projection grid_combine with rugosity==================================
# ==============================================================================

load('data/bathy_grid_course_Feb7_2017.RData') # loads 'proj.grid.bath'
proj.sed <- readRDS('data/benthic_hab_proj.rds')
setkey(proj.grid.bath, latBathgrid, lonBathgrid); setkey(proj.sed, latBathgrid, lonBathgrid)
proj.grid <- merge(proj.grid.bath, proj.sed[,c('latBathgrid', 'lonBathgrid', 'GRAINSIZE')], all.x=TRUE)

proj.grid <- proj.grid[complete.cases(proj.grid),]
# proj.grid$habitatFact <- as.character(proj.grid$habitat)
# proj.grid$habitatFact[proj.grid$habitatFact %in% c("rocky_reef", "shelf_hard", "slope_hard", "deep_hard")] <- "hard"
# proj.grid$habitatFact[proj.grid$habitatFact %in% c("subt_soft", "slope_soft", "shelf_soft", "deep_soft")] <- "soft"
# proj.grid$habitatFact <- factor(proj.grid$habitatFact, levels=c('soft', 'hard'), ordered=FALSE) 

# Drop a handful of grid cells that are apart from the main prediction areas
proj.grid <- proj.grid[!(latBathgrid < 25 & lonBathgrid > -100 & lonBathgrid < -95)]
proj.grid <- proj.grid[!(latBathgrid < 55 & latBathgrid > 53.5 & lonBathgrid < -175)]
# Drop deep grid cells
proj.grid <- proj.grid[depth > -401]
proj.grid <- proj.grid[latBathgrid > 24.27]
# Drop the bahamas
proj.grid <- proj.grid[!(latBathgrid < 27.7 & lonBathgrid > -79.7)]

proj.grid$rugosity <- log(proj.grid$rugosity + 1) #Log rugosity to match habitat models

# Plot up sediment and rugosity data
cutpts <- c(-5, -3, -1, 0, 1, 2, 3, 4, 5, 7, 8)
pdf(width=15, height=9, file='figures/projection_grid_sediment.pdf')
par(mai=c(0.7,0.7,0.3, 0.1), las=1, mgp=c(2,1,0))
levelplot(GRAINSIZE ~ lonBathgrid * latBathgrid, data = proj.grid, at = cutpts, cuts = 11, pretty = T, col.regions = (rev(brewer.pal(11, "RdBu"))))
dev.off()
 
cutpts <- c(0, .5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 7)
pdf(width=15, height=9, file='figures/projection_grid_rugosity.pdf')
par(mai=c(0.7,0.7,0.3, 0.1), las=1, mgp=c(2,1,0))
levelplot(rugosity ~ lonBathgrid * latBathgrid, data = proj.grid, at = cutpts, cuts = 11, pretty = T, col.regions = (rev(brewer.pal(11, "RdBu"))))
dev.off()

rm(cutpts, proj.sed, proj.grid.bath)
save(proj.grid, file='data/ProjectionBathGrid_Feb27_2017.RData')

# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================

# map the depth data with a color ramp
mapTheme <- colorRampPalette(rev(brewer.pal(10, "Spectral")))(100)
plot(rugos.east, col=mapTheme) 

cols <- colorRampPalette(rev(brewer.pal(10, "Spectral")))(50)
cols = colorRampPalette(colors = c('white', 'orange', 'dark orange', 'black'))(40)
plot(rough.east, col=cols)
map(database="world", fill=T, col="gray", xlim=c(-100.2, -42.325), ylim=c(23.625, 62.083), add=T) # mar=c(2,3,0,0))
map(database="state", fill=T, col="gray", xlim=c(-100.2, -42.325), ylim=c(23.625, 62.083), mar=c(2,3,0,0), add=TRUE)

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
plot(latgrid~longrid, data=proj.grid.aleut)

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
