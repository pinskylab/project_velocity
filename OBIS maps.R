library(maps)
    
setwd('/Users/jamesmorley/Documents/project_velocity')
SPECIES <- c('homarus americanus_Atl','paralichthys dentatus_Atl','centropristis striata_Atl','anoplopoma fimbria_Pac','doryteuthis opalescens_Pac','hippoglossus stenolepis_Pac','sebastes alutus_Pac')
z=3
 
OBIS <- data.frame(read.csv(paste('OBIS_files/', SPECIES[z], '_OBIS.csv', sep='')), stringsAsFactors = FALSE)
OBIS <- OBIS[,c(2,5,6,7,8,4,142,35,60)]
#OBIS <- OBIS[,c(3,6,7,1,9,5,145,38,10)] # for ocean perch only, I think b/c I may have manipulated the csv file in excel

# Need to remove unwanted data sources
OBIS <- OBIS[!is.na(OBIS$year),]
OBIS <- OBIS[!OBIS$year < 1970,]
OBIS <- droplevels(OBIS)

nrow(OBIS[is.na(OBIS$institutionCode),]) # make sure everything has a source code
table(OBIS$institutionCode)
table(OBIS$collectionCode)
drops <- c('SeamountsOnline','EGG COLLECTION','LARVAL COLLECTION',"Neuston net (16', 0.947)",'Bongo Net','EAISSNA','MOC10GB','noaa_cpr_zoo','FALL NMFS NEFSC BOTTOM TRAWL SURVEY', 'NEAMAP', 'SEAMAP South Atlantic','SPRING NMFS NEFSC BOTTOM TRAWL SURVEY') 
drops2 <- c('USNM','Bedford Institute of Oceanography (BIO)', 'DFO Gulf Fisheries Centre') #For 'institutionCode' drops
#plot(decimalLatitude~decimalLongitude, data=OBIS[OBIS$institutionCode=='Pacific Biological Station',])
table(OBIS$scientificName)

#OBIS <- OBIS[OBIS$scientificName=='Homarus americanus',]
OBIS <- OBIS[!(OBIS$collectionCode %in% drops),]
OBIS <- OBIS[!(OBIS$institutionCode %in% drops2),]
OBIS <- droplevels(OBIS)
hist(OBIS$depth) # need to drop deep areas?

# Add back trawl data that we have_need to format columns
modfolder <- 'output/CEmodels_Uncertainty_2018/'
filename <- paste(modfolder, 'RawData_', SPECIES[z], '_uncertainty2018.RData', sep='')  
load(filename)
rm(spdataB, spdataC, testinds, testindsp, traininds, trainindsp)
spdata <- spdata[spdata$presfit==TRUE,]

# Combine the lat and lon records from survey data and OBIS, then snap to a lat and lon grid
colnames(OBIS) <- c("scientificName","lon","lat","institutionCode","collectionCode","eventDate","year","depth","individualCount") 
obsRecs <- data.frame(lat=spdata$lat, lon=spdata$lon)
obsRecs_OBIS <- data.frame(lat=OBIS$lat, lon=OBIS$lon)
obsRecs <- data.frame(rbind(obsRecs, obsRecs_OBIS))

gridsize=0.25 # size of grid of the climate data, in degrees
obsRecs$latgrid = floor(obsRecs$lat/gridsize)*gridsize + gridsize/2 # round to nearest grid center
obsRecs$longrid = floor(obsRecs$lon/gridsize)*gridsize + gridsize/2

obsRecsAgg <- aggregate(obsRecs$lat, by=list(obsRecs$latgrid, obsRecs$longrid), FUN=length)
obsRecsAgg$logAbun <- log(obsRecsAgg$x)

# Create a scale for a levelplot type deal
cols = colorRampPalette(colors = c('gray80', 'blue', 'dark blue', 'black'))(n=20)
scaleOBIS <- seq(0, max(obsRecsAgg$logAbun), length.out=20)

colKey <- character()
colVal <- numeric()
for(i in 1:nrow(obsRecsAgg)){
  value = obsRecsAgg$logAbun[i]
  scaleVals <- value-scaleOBIS
  indexVals <- which(scaleVals>=0)
  colKey[i] <- cols[max(indexVals)]
  colVal[i] <- max(indexVals)
}
 
obsRecsAgg$colKey <- colKey
obsRecsAgg$colVal <- colVal

# Need to load pred.agg_PA and use that to create identical map boundaries as the main figure
load(paste('figures/Uncertainty_MS_2018/Final_proj_figure_', SPECIES[z], '.RData', sep=''))
rm(plot1,plot2)

# Grab the same map limits from the other figure
if(SPECIES[z]  == 'homarus americanus_Atl'){ 
  pred.agg_PA <- pred.agg_PA[pred.agg_PA$latitude < 52,]
  pred.agg_PA <- pred.agg_PA[pred.agg_PA$latitude > 34,]
  xlimit = c(min(pred.agg_PA$longitude), max(pred.agg_PA$longitude))
  ylimit = c(min(pred.agg_PA$latitude), max(pred.agg_PA$latitude))
}
if(SPECIES[z] == 'paralichthys dentatus_Atl'){
  pred.agg_PA <- pred.agg_PA[pred.agg_PA$latitude < 51,]
  pred.agg_PA <- pred.agg_PA[!pred.agg_PA$longitude < -82,]
  pred.agg_PA <- pred.agg_PA[!(pred.agg_PA$longitude < -80.75 & pred.agg_PA$latitude < 27),]
  xlimit = c(min(pred.agg_PA$longitude)-.5, max(pred.agg_PA$longitude))
  ylimit = c(min(pred.agg_PA$latitude), max(pred.agg_PA$latitude))
  obsRecsAgg <- obsRecsAgg[-c(6,57,191),]
}
if(SPECIES[z] == 'centropristis striata_Atl'){
  pred.agg_PA <- pred.agg_PA[pred.agg_PA$latitude < 51,]
  xlimit = c(min(pred.agg_PA$longitude), max(pred.agg_PA$longitude))
  ylimit = c(min(pred.agg_PA$latitude), max(pred.agg_PA$latitude))
}
if(SPECIES[z] == 'anoplopoma fimbria_Pac'){
  xlimit = c(min(pred.agg_PA$longitude), max(pred.agg_PA$longitude) + 0.5)
  ylimit = c(min(pred.agg_PA$latitude), max(pred.agg_PA$latitude))
}
if(SPECIES[z] == 'doryteuthis opalescens_Pac'){
  xlimit = c(min(pred.agg_PA$longitude), max(pred.agg_PA$longitude) + 0.5)
  ylimit = c(min(pred.agg_PA$latitude), max(pred.agg_PA$latitude))
}
if(SPECIES[z] == 'hippoglossus stenolepis_Pac'){
  xlimit = c(min(pred.agg_biom$longitude), max(pred.agg_biom$longitude) + 0.5)
  ylimit = c(min(pred.agg_biom$latitude), max(pred.agg_biom$latitude))
}
if(SPECIES[z] == 'sebastes alutus_Pac'){
  xlimit = c(min(pred.agg_biom$longitude), max(pred.agg_biom$longitude) + 0.5)
  ylimit = c(min(pred.agg_biom$latitude), max(pred.agg_biom$latitude))
}
 
save(obsRecsAgg, xlimit, ylimit, file=paste('OBIS_files/', SPECIES[z], '_OBIS_plot.RData', sep='')) 

# Now put the figure together with all the save files

pdf(width=6.5, height=8, file='figures/Uncertainty_MS_2018/OBISmaps_Jan2019.pdf')
# Halibut
par(mfrow=c(4,2), mar=c(0,0,0,0), oma=c(0.1,0.1,0.1,0.1))
load(paste('OBIS_files/', SPECIES[6], '_OBIS_plot.RData', sep=''))
plot(Group.1~Group.2, xlim=xlimit, ylim=ylimit, xaxt='n', yaxt='n', xlab='', ylab='', cex=.2, pch=15, col=obsRecsAgg$colKey, obsRecsAgg)
map('world', xlim=xlimit, ylim=ylimit, fill=TRUE, lwd=.4, col='dark gray', add=TRUE)
box()
points(Group.1~Group.2, cex=.2, pch=15, col=obsRecsAgg$colKey, obsRecsAgg)
mtext('(a)',line=-1.25,adj=0.01,font=2)
# Perch
load(paste('OBIS_files/', SPECIES[7], '_OBIS_plot.RData', sep=''))
plot(Group.1~Group.2, xlim=xlimit, ylim=ylimit, xaxt='n', yaxt='n', xlab='', ylab='', cex=.2, pch=15, col=obsRecsAgg$colKey, obsRecsAgg)
map('world', xlim=xlimit, ylim=ylimit, fill=TRUE, lwd=.4, col='dark gray', add=TRUE)
box()
points(Group.1~Group.2, cex=.2, pch=15, col=obsRecsAgg$colKey, obsRecsAgg)
mtext('(b)',line=-1.25,adj=0.02,font=2)
# Summer flounder
load(paste('OBIS_files/', SPECIES[2], '_OBIS_plot.RData', sep=''))
plot(Group.1~Group.2, xlim=xlimit, ylim=ylimit, xaxt='n', yaxt='n', xlab='', ylab='', cex=.4, pch=15, col=obsRecsAgg$colKey, obsRecsAgg)
map('world', xlim=xlimit, ylim=ylimit, fill=TRUE, lwd=.4, col='dark gray', add=TRUE)
box()
points(Group.1~Group.2, cex=.4, pch=15, col=obsRecsAgg$colKey, obsRecsAgg)
mtext('(c)',line=-1.25,adj=0.01,font=2)
# Lobster
load(paste('OBIS_files/', SPECIES[1], '_OBIS_plot.RData', sep=''))
plot(Group.1~Group.2, xlim=xlimit, ylim=ylimit, xaxt='n', yaxt='n', xlab='', ylab='', cex=.5, pch=15, col=obsRecsAgg$colKey, obsRecsAgg)
map('world', xlim=xlimit, ylim=ylimit, fill=TRUE, lwd=.4, col='dark gray', add=TRUE)
box()
points(Group.1~Group.2, cex=.5, pch=15, col=obsRecsAgg$colKey, obsRecsAgg)
mtext('(d)',line=-1.25,adj=0.02,font=2)
# Sablefish
load(paste('OBIS_files/', SPECIES[4], '_OBIS_plot.RData', sep=''))
plot(Group.1~Group.2, xlim=xlimit, ylim=ylimit, xaxt='n', yaxt='n', xlab='', ylab='', cex=.2, pch=15, col=obsRecsAgg$colKey, obsRecsAgg)
map('world', xlim=xlimit, ylim=ylimit, fill=TRUE, lwd=.4, col='dark gray', add=TRUE)
box()
points(Group.1~Group.2, cex=.3, pch=15, col=obsRecsAgg$colKey, obsRecsAgg)
mtext('(e)',line=-1.25,adj=0.01,font=2)
# Squid
load(paste('OBIS_files/', SPECIES[5], '_OBIS_plot.RData', sep=''))
plot(Group.1~Group.2, xlim=xlimit, ylim=ylimit, xaxt='n', yaxt='n', xlab='', ylab='', cex=.2, pch=15, col=obsRecsAgg$colKey, obsRecsAgg)
map('world', xlim=xlimit, ylim=ylimit, fill=TRUE, lwd=.4, col='dark gray', add=TRUE)
box()
points(Group.1~Group.2, cex=.3, pch=15, col=obsRecsAgg$colKey, obsRecsAgg)
mtext('(f)',line=-1.25,adj=0.02,font=2)
# BSB
load(paste('OBIS_files/', SPECIES[3], '_OBIS_plot.RData', sep=''))
plot(Group.1~Group.2, xlim=xlimit, ylim=ylimit, xaxt='n', yaxt='n', xlab='', ylab='', cex=.3, pch=15, col=obsRecsAgg$colKey, obsRecsAgg)
map('world', xlim=xlimit, ylim=ylimit, fill=TRUE, lwd=.4, col='dark gray', add=TRUE)
box()
points(Group.1~Group.2, xlim=xlimit, ylim=ylimit, cex=.3, pch=15, col=obsRecsAgg$colKey, obsRecsAgg)
mtext('(g)',line=-1.25,adj=0.01,font=2)

dev.off()

