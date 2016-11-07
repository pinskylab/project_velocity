# Make a shapefile of our Areas of Interest, for use in the NatCap's InVest software
# Also make a csv of lat/lon for landing points and grid connection points

## Set working directories
if(Sys.info()["nodename"] == "pinsky-macbookair"){
	setwd('~/Documents/Rutgers/Range projections/proj_ranges/')
	natcapfolder <- '../NatCap'
	}
if(Sys.info()["nodename"] == "amphiprion.deenr.rutgers.edu"){
	setwd('~/Documents/range_projections/')
	}
# could add code for Lauren's working directory here

#######################
## Load libraries
#######################
require(maptools)
require(RColorBrewer)
require(rgdal)
require(PBSmapping)
require(rgeos) # for gBuffer


##############################################################
## Read in and process data to create Areas of Interest (AOI)
##############################################################
climgrid <- read.csv('data/climGrid.csv', row.names=1)
#land <- readOGR(dsn='cmsp_data/', layer='global_polygon') # global land layer from NatCap
#	plot(land) # takes too long

# label by area of interest regions
#climgrid$AOI <- 'northamerica'
climgrid$AOI <- NA
climgrid$AOI[climgrid$region %in% c('AFSC_Aleutians', 'AFSC_EBS', 'AFSC_GOA')] <- 'Alaska'
#climgrid$AOI[climgrid$region %in% c('AFSC_Aleutians')] <- 'Alaska_Aleutians'
#climgrid$AOI[climgrid$region %in% c('AFSC_Aleutians') & climgrid$lon > 183] <- 'Alaska_Aleutians_east'
#climgrid$AOI[climgrid$region %in% c('AFSC_Aleutians') & climgrid$lon >= 177 & climgrid$lon <= 183] <- 'Alaska_Aleutians_middle'
#climgrid$AOI[climgrid$region %in% c('AFSC_Aleutians') & climgrid$lon < 177] <- 'Alaska_Aleutians_west'
#climgrid$AOI[climgrid$region %in% c('AFSC_EBS')] <- 'Alaska_EBS'
#climgrid$AOI[climgrid$region %in% c('AFSC_GOA')] <- 'Alaska_GOA'
#climgrid$AOI[climgrid$region %in% c('AFSC_WCTri', 'NWFSC_WCAnn', 'SEFSC_GOMex', 'DFO_NewfoundlandFall', 'DFO_NewfoundlandSpring', 'DFO_ScotianShelf', 'DFO_SoGulf', 'NEFSC_NEUSFall', 'NEFSC_NEUSSpring')] <- 'notAlaska'
climgrid$AOI[climgrid$region %in% c('AFSC_WCTri', 'NWFSC_WCAnn')] <- 'WestCoast'
climgrid$AOI[climgrid$region %in% c('SEFSC_GOMex')] <- 'GoMex'
climgrid$AOI[climgrid$region %in% c('DFO_NewfoundlandFall', 'DFO_NewfoundlandSpring', 'DFO_ScotianShelf', 'DFO_SoGulf', 'NEFSC_NEUSFall', 'NEFSC_NEUSSpring')] <- 'Northeast'
sum(is.na(climgrid$AOI))
sort(unique(climgrid$AOI))

# draw 1/4deg boxes around all grid cell centers
PIDs <- rep(1:nrow(climgrid), rep(4,nrow(climgrid))) # polygon ids
POSs <- rep(1:4, nrow(climgrid)) # vertex ID within a polygon
Xs <- rep(climgrid$lon, rep(4, nrow(climgrid))) + c(-0.125, 0.125, 0.125, -0.125)
Ys <- rep(climgrid$lat, rep(4, nrow(climgrid))) + c(-0.125, -0.125, 0.125, 0.125) # lower left, upper left, upper right, lower right

gridPS <- as.PolySet(data.frame(PID = PIDs, POS=POSs, X=Xs, Y=Ys), projection='LL')
	plotPolys(gridPS) # appears to cover a somewhat wider area than we project to. that's OK
gridSP <- PolySet2SpatialPolygons(gridPS) # convert to spatial polygon for union

# merge all into multipart polygons based on AOI
gridSP.m <- unionSpatialPolygons(gridSP, climgrid$AOI)	
	plot(gridSP.m, lwd=0.2, axes=TRUE, col=1:length(gridSP.m))
#	plot(gridSP.m, axes=TRUE, ylim=c(40,70), xlim=c(-170,-60)) # doesn't work
	
# project into meters using Lambert Conformal Conic
gridSP.p <- spTransform(gridSP.m, CRS('+proj=lcc +lat_1=32 +lat_2=44 +lat_0=40 +lon_0=-96 +datum=WGS84'))
	plot(gridSP.p, lwd=0.2, axes=TRUE, col=1:4)

# buffer
gridSP.b <- gBuffer(gridSP.p, width=100000, byid=TRUE) # 100km
	plot(gridSP.b, lwd=0.2, axes=TRUE, col=1:4)
	
# convert to spatialpolygonsdataframe for writing
pid <- sapply(slot(gridSP.b, "polygons"), function(x) slot(x, "ID"))
gridSPD <- SpatialPolygonsDataFrame(gridSP.b, data.frame(region=pid, row.names=pid))

# split Alaska along 179 and 181deg longitude (3 pieces)
maskPS <- as.PolySet(data.frame(PID=rep(c(1:3),rep(4,3)), POS=rep(1:4, 3), X=rep(c(150,179.99,180.01,230), c(2,4,4,2)), Y=rep(c(45,65,45,65,45,65,45), c(1,2,2,2,2,2,1))), projection='LL')
mask <- PolySet2SpatialPolygons(maskPS)
mask.t <- spTransform(mask, proj4string(gridSPD))
	#plot(mask.t, col=1:length(mask.t))
	plot(gridSP.b, lwd=0.2, axes=TRUE, col=1:4)
	plot(mask.t, add=TRUE)
gridSPAK <- gIntersection(gridSPD[grep('Alaska', gridSPD$region),], mask.t, byid=TRUE)
	plot(gridSPAK, col=1:length(gridSPAK))

	# convert Alaska to spatialpolygonsdataframe for writing
	pid <- sapply(slot(gridSPAK, "polygons"), function(x) slot(x, "ID"))
	regs <- c('Alaska_west', 'Alaska_middle', 'Alaska_east')
	gridSPDAK <- SpatialPolygonsDataFrame(gridSPAK, data.frame(region=regs, row.names=pid))

# write out each region as a separate shapefile
#writePolyShape(merged3, fn='cmsp_data/AOI_northamerica') # won't write .prj
notAK <- grep('Alaska', gridSPD$region, invert=TRUE)
for(i in notAK){ # for all except Alaska
	writeOGR(gridSPD[i,], "./cmsp_data/", paste('AOI', gridSPD$region[i], sep='_'), driver="ESRI Shapefile", overwrite_layer=TRUE) # write the .prj file as part of the shapefile
}

# use this if using 'notAlaska' region
#i <- which(gridSPD$region=='notAlaska')
#writeOGR(gridSPD[i,], "./cmsp_data/", paste('AOI', gridSPD$region[i], sep='_'), driver="ESRI Shapefile", overwrite_layer=TRUE) # write the .prj file as part of the shapefile

for(i in 1:length(gridSPDAK)){ # for Alaska
	writeOGR(gridSPDAK[i,], "./cmsp_data/", paste('AOI', gridSPDAK$region[i], sep='_'), driver="ESRI Shapefile", overwrite_layer=TRUE) # write the .prj file as part of the shapefile
}


#######################################
## Make land and grid connection points
#######################################
# set parameters
crs <- CRS('+proj=lcc +lat_1=32 +lat_2=44 +lat_0=40 +lon_0=-96 +datum=WGS84')
crslatlong = CRS("+init=epsg:4326")

# read in files
coastlines <- readOGR(dsn=natcapfolder, layer='NAmainland_lines') # global coast layer from NatCap
#	plot(coastlines)
towns <- readOGR(dsn=paste(natcapfolder, 'ne_10m_populated_places', sep='/'), layer='ne_10m_populated_places') # populated places layer from Natural Earth data
#	plot(towns, pch=16, cex=0.2)
	# plot(towns, pch=16, cex=0.2, add=TRUE) # to add to coastlines plot

# calculate coastline length
ln <- SpatialLinesLengths(coastlines, longlat=TRUE) # returns answer in km for each line
	ln
	
# trim towns to those >1000 people
towns1000 <- towns[which(towns@data$POP_MAX >= 1000),]
	length(towns) # 7343
	length(towns1000) # 6933
	
# convert shps to planar coordinates for spatial sampling
coastlines.p <- spTransform(coastlines, CRSobj=crs)
	plot(coastlines.p, lwd=0.2, axes=TRUE)
towns.p <- spTransform(towns1000, CRSobj=crs)
#	plot(towns.p, pch=16) # odd plot: some towns seem to have extreme coordinats
	plot(towns.p, pch=16, cex=0.5, add=TRUE, col='blue') # but plots on top of NA well

# trim town to those <20km from the NA coast
nearcoast <- gWithinDistance(coastlines.p, towns.p, byid=TRUE, dist=50*1000) # slow (a few min). units in meters (50km). returns a matrix with columns corresponding to each coastline (includes a few major islands)
nearcoastany <- rowSums(nearcoast) > 0 # sum across rows: we don't care which coast a town is close to
	sum(nearcoast) # 370
	sum(nearcoastany) # 347. shows that some towns were close to multiple coastlines
	points(towns.p[which(nearcoastany),], col='red', pch=16, cex=1) # adds to the plot before

towns.nearcoast <- towns.p[which(nearcoastany),]


# add points along the coastline. 1 every km
landpts <- spsample(coastlines.p, type='regular', n=round(ln))
	plot(coastlines.p, lwd=0.2, axes=TRUE)
	points(landpts, col='red', pch=16, cex=0.2)	
	
	plot(coastlines.p, lwd=0.2, axes=TRUE, ylim=c(0, 6e5), xlim=c(-2.4e6, -2e6)) # zoom in on CA coast
	points(landpts, col='red', pch=16, cex=0.2)	

# find nearest landpt to each town
landdist <- gDistance(landpts, towns.nearcoast, byid=TRUE) # takes a couple minutes. rows are towns. cols are landpts

# project back to latlong in order to make a table for NatCap
landpts.ll <- spTransform(landpts, crslatlong)
towns.nearcoast.ll <- spTransform(towns.nearcoast, crslatlong)

# make output table of town and nearest landing point locations
towns.coords <- coordinates(towns.nearcoast.ll)
land.coords <- coordinates(landpts.ll)

n <- numeric(nrow(landdist))
outgrid <- data.frame(ID=1:nrow(landdist), LAT=towns.coords[,2], LONG=towns.coords[,1], TYPE='GRID', LOCATION=towns.nearcoast.ll@data$NAME)
outland <- data.frame(ID=1:nrow(landdist), LAT=n, LONG=n, TYPE='LAND', LOCATION=towns.nearcoast.ll@data$NAME)

for(i in 1:nrow(landdist)){ # fill in nearest landpt for each town
	j <- which.min(landdist[i,]) # index for nearest landpt
	outland$LAT[i] <- land.coords[j,2]
	outland$LONG[i] <- land.coords[j,1]
}

	# make sure it looks OK
	plot(outland$LONG, outland$LAT, pch=16, cex=0.5)
	points(outgrid$LONG, outgrid$LAT, col='red', pch=16, cex=0.5) # the towns
	
# combine town and landings points
out <- rbind(outland, outgrid)

	head(out)
	tail(out)
	
# write out
write.csv(out, file='cmsp_data/LandGridPts_NorthAmerica.csv', row.names=FALSE)