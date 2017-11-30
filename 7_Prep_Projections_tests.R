# Script to try and figure out the blocky climate projection data

  setwd('/Users/jamesmorley/Documents/project_velocity')
  climgridfolder <- 'data/'

library(latticeExtra)
library(maps)
library(geosphere)

###############################################
# Choose projection options_rerun this code for different rcp and seasons
###############################################

# rcp <- 26
rcp <- 85

# season <- 'jfm'
# season <- 'amj'
season <- 'jas'
# season <- 'ond'

#################################
# Prep environmental data
#################################
      
load('data/projectionGrid_Feb24_2017.RData')# load projection grid to get lat/lon values
clim.grid <- proj.grid # rename as a different 'proj.grid' imported below with bathymetry
clim.grid$depth <- NULL # depth not needed
rm(proj.grid)
   
modelrun <- c('bcc-csm1-1-m','bcc-csm1-1','CanESM2','CCSM4','CESM1-CAM5','CNRM-CM5','GFDL-CM3','GFDL-ESM2M','GFDL-ESM2G','GISS-E2-R','GISS-E2-H','IPSL-CM5A-LR','IPSL-CM5A-MR','MIROC-ESM','MPI-ESM-LR','NorESM1-ME')
#latest projections were missing one of the models for rcp85 SBT: 'IPSL-CM5A-MR'
#modelrun <- c('bcc-csm1-1', 'bcc-csm1-1-m','CanESM2','CCSM4','CESM1-CAM5','CNRM-CM5','GFDL-CM3', 'GFDL-ESM2M', 'GFDL-ESM2G','GISS-E2-R','GISS-E2-H','IPSL-CM5A-LR','MIROC-ESM','MPI-ESM-LR','NorESM1-ME')
#modelrun <- c('bcc-csm1-1-m','CanESM2','CCSM4','CESM1-CAM5','CNRM-CM5','GFDL-CM3', 'GFDL-ESM2G','GISS-E2-R','IPSL-CM5A-LR','MPI-ESM-LR','NorESM1-ME')

pred.metric <- c('max', 'min', 'mean')

# Functions for reading in the temperature data_need separate function for the min/max files as they have an extra column
processlinesMean <- function(x){
  if(grepl('.....', x, fixed=TRUE)){
    return(rep(as.numeric(NA), 94))
  } else {
    x <- sub('^ +', '', x) # remove leading whitespace
    return(as.numeric(unlist(strsplit(x, split=' +'))))
  }
}

processlinesMinMax <- function(x){
  if(grepl('.....', x, fixed=TRUE)){
    return(rep(as.numeric(NA), 95))
  } else {
    x <- sub('^ +', '', x) # remove leading whitespace
    return(as.numeric(unlist(strsplit(x, split=' +'))))
  }
}

pred.folder <- 'tests/test7/'

# For reading in all models for a particular type of data_minimum surface temps was not in habitat models
tempsmeansbt <- matrix(as.numeric(NA), nrow=1281878, ncol=length(modelrun))
#tempsmeansst <- matrix(as.numeric(NA), nrow=1281878, ncol=length(modelrun))
tempsminsbt <- matrix(as.numeric(NA), nrow=1281878, ncol=length(modelrun))
tempsmaxsbt <- matrix(as.numeric(NA), nrow=1281878, ncol=length(modelrun))
#tempsmaxsst <- matrix(as.numeric(NA), nrow=1281878, ncol=length(modelrun))

i = 7
# The loops below produce matrix where each column is a different climate model and the years are stacked in one vector
# seasonal mean sbt
for(i in 1:length(modelrun)){
  #print(i) 
  #filename = paste('data/Temp proj files/', pred.folder[2], modelrun[i], '_rcp', rcp, '_fill.nc_2006_2100_', season,'_', pred.metric[3], '.txt', sep="")
  filename = paste('data/Temp proj files/', pred.folder, modelrun[i], '_rcp', rcp, '_fill.nc_2006_2100_', season,'_', pred.metric[3], '.txt', sep="")
  filein <- readLines(filename)
  filein <- t(sapply(filein, processlinesMean))
  rownames(filein) <- NULL
  print(paste(i, modelrun[i], dim(filein)[1], dim(filein)[2]))
  tempsmeansbt[,i] <- as.vector(filein)
} 
# seasonal mean sst
for(i in 1:length(modelrun)){
  #print(i) 
  filename = paste('data/Temp proj files/', pred.folder[1], modelrun[i], '_rcp', rcp, '_r1i1p1_1950_2100.nc_regrid.nc_fill.nc_2006_2100_', season,'_', pred.metric[3], '_no_scaling.txt', sep="")
  filein <- readLines(filename)
  filein <- t(sapply(filein, processlinesMean))
  rownames(filein) <- NULL
  print(paste(i, modelrun[i], dim(filein)[1], dim(filein)[2]))
  tempsmeansst[,i] <- as.vector(filein)
} 
# min annual sbt
for(i in 1:length(modelrun)){
  #print(i)
  #filename = paste('data/Temp proj files/', pred.folder[2], modelrun[i], '_rcp', rcp, '_fill.nc_2006_2100_', season,'_', pred.metric[2], '.txt', sep="")
  filename = paste('data/Temp proj files/', pred.folder, modelrun[i], '_rcp', rcp, '_fill.nc_2006_2100_', season,'_', pred.metric[2], '.txt', sep="")
  filein <- readLines(filename)
  filein <- t(sapply(filein, processlinesMinMax))
  rownames(filein) <- NULL
  fileinVec <- as.vector(filein)
  fileinVec <- fileinVec[1:1281878]
  print(paste(i, modelrun[i], dim(filein)[1], dim(filein)[2]))
  tempsminsbt[,i] <- fileinVec
}
# max annual sbt
for(i in 1:length(modelrun)){
  #print(i)
  #filename = paste('data/Temp proj files/', pred.folder[2], modelrun[i], '_rcp', rcp, '_fill.nc_2006_2100_', season,'_', pred.metric[1], '.txt', sep="")
  filename = paste('data/Temp proj files/', pred.folder, modelrun[i], '_rcp', rcp, '_fill.nc_2006_2100_', season,'_', pred.metric[1], '.txt', sep="")
  filein <- readLines(filename)
  filein <- t(sapply(filein, processlinesMinMax))
  rownames(filein) <- NULL
  fileinVec <- as.vector(filein)
  fileinVec <- fileinVec[1:1281878]
  print(paste(i, modelrun[i], dim(filein)[1], dim(filein)[2]))
  tempsmaxsbt[,i] <- fileinVec
}
# max annual sst
for(i in 1:length(modelrun)){
  #print(i)
  filename = paste('data/Temp proj files/', pred.folder[1], modelrun[i], '_rcp', rcp, '_r1i1p1_1950_2100.nc_regrid.nc_fill.nc_2006_2100_', season,'_', pred.metric[1], '_no_scaling.txt', sep="")
  filein <- readLines(filename)
  filein <- t(sapply(filein, processlinesMinMax))
  rownames(filein) <- NULL
  fileinVec <- as.vector(filein)
  fileinVec <- fileinVec[1:1281878]
  print(paste(i, modelrun[i], dim(filein)[1], dim(filein)[2]))
  tempsmaxsst[,i] <- fileinVec
}
    
colnames(clim.grid) <- c('lonClimgrid', 'latClimgrid') 
clim.grid$index <- c(1:nrow(clim.grid))
load('data/ProjectionBathGrid_Feb27_2017.RData')# load bathymetry projection grid
proj.grid$depth <- NULL # To reduce file size
year.bin <- data.frame(year=2007:2100, bin=c(rep('2007-2020', 14), rep('2021-2040', 20), rep('2041-2060', 20), rep('2061-2080', 20), rep('2081-2100', 20))) # for making summary graphs in the loop
cols = colorRampPalette(colors = c('dark blue', 'blue', 'white', 'red', 'dark red'))
 
# Build and save all 16 prediction data sets_make some summary figures in 20 year bins
for(i in 1:length(modelrun)){ 
  print(paste('Beginning model number ', i, ': ', modelrun[i], sep=''))
  pred <- data.frame(cbind(clim.grid, year=rep(2007:2100, rep=94, each=13637), SBT.seasonal = tempsmeansbt[,i], SBT.min = tempsminsbt[,i], SBT.max = tempsmaxsbt[,i]))
  #pred <- data.frame(cbind(clim.grid, year=rep(2007:2100, rep=94, each=13637), SST.seasonal.mean = tempsmeansst[,i], SST.max = tempsmaxsst[,i]))
  pred <- merge(pred, year.bin, by='year', all.x=T)  # Bin into 20 year periods to make some plots
  pred.bathE <- merge(pred, proj.grid[proj.grid$lonClimgrid > -105,], by=c('lonClimgrid', 'latClimgrid'), all.y = T, sort=F)
  pred.bathW <- merge(pred, proj.grid[proj.grid$lonClimgrid < -105,], by=c('lonClimgrid', 'latClimgrid'), all.y = T, sort=F)
  pred.bathE <- pred.bathE[order(pred.bathE$year, pred.bathE$index),] # reorder by year and index
  pred.bathW <- pred.bathW[order(pred.bathW$year, pred.bathW$index),] # reorder by year and index
  
  #aggE <- aggregate(list(sbt.seasonal = pred.bathE$SBT.seasonal, SST.seasonal.mean = pred.bathE$SST.seasonal.mean, SBT.min = pred.bathE$SBT.min,
   #         SBT.max = pred.bathE$SBT.max, SST.max = pred.bathE$SST.max), by=list(year_range=pred.bathE$bin, latitude=pred.bathE$latBathgrid, longitude=pred.bathE$lonBathgrid), FUN=mean)
  #aggW <- aggregate(list(sbt.seasonal = pred.bathW$SBT.seasonal, SST.seasonal.mean = pred.bathW$SST.seasonal.mean, SBT.min = pred.bathW$SBT.min,
   #         SBT.max = pred.bathW$SBT.max, SST.max = pred.bathW$SST.max), by=list(year_range=pred.bathW$bin, latitude=pred.bathW$latBathgrid, longitude=pred.bathW$lonBathgrid), FUN=mean)
  aggE <- aggregate(list(SBT.seasonal = pred.bathE$SBT.seasonal, SBT.min = pred.bathE$SBT.min,
                         SBT.max = pred.bathE$SBT.max), by=list(year_range=pred.bathE$bin, latitude=pred.bathE$latBathgrid, longitude=pred.bathE$lonBathgrid), FUN=mean)
  aggW <- aggregate(list(SBT.seasonal = pred.bathW$SBT.seasonal, SBT.min = pred.bathW$SBT.min,
                         SBT.max = pred.bathW$SBT.max), by=list(year_range=pred.bathW$bin, latitude=pred.bathW$latBathgrid, longitude=pred.bathW$lonBathgrid), FUN=mean)
  
  # 20yr. summary figures for each model
  #pdf(width=14, height=5, file=paste('figures/Temp_projections_noBiasCorr2/', modelrun[i], '_rcp', rcp, '_20yr_tempProj.pdf', sep=''))
  pdf(width=14, height=5, file=paste('figures/TEST_7_deltas_SODAdepths_climatology', modelrun[i], '_rcp', rcp, '_20yr_tempProj.pdf', sep=''))
   
  print(levelplot(SBT.seasonal~longitude*latitude|year_range, data=aggE, xlab=NULL, ylab=NULL, main=paste('sbt.seasonal_rcp', rcp, '  season:',season, '  Model:', modelrun[i]), at = seq(min(aggE$SBT.seasonal-.01, na.rm=T),max(aggE$SBT.seasonal+.01, na.rm=T), length.out=40), 
            col.regions=cols, layout = c(5, 1), panel = function(...) {
              panel.fill(col = "light gray")
              panel.levelplot(...)
            }))  
  #print(levelplot(SST.seasonal.mean~longitude*latitude|year_range, data=aggE, xlab=NULL, ylab=NULL, main=paste('SST.seasonal.mean_rcp', rcp, '  season:',season, '  Model:', modelrun[i]), at = seq(min(aggE$SST.seasonal.mean-.1),max(aggE$SST.seasonal.mean+.1), length.out=40), 
   #         col.regions=cols, layout = c(5, 1), panel = function(...) {
    #          panel.fill(col = "light gray")
     #         panel.levelplot(...)
      #      }))  
  print(levelplot(SBT.min~longitude*latitude|year_range, data=aggE, xlab=NULL, ylab=NULL, main=paste('SBT.min_rcp', rcp, '  season:',season, '  Model:', modelrun[i]), at = seq(min(aggE$SBT.min-.01, na.rm=T),max(aggE$SBT.min+.01, na.rm=T), length.out=40), 
            col.regions=cols, layout = c(5, 1), panel = function(...) {
              panel.fill(col = "light gray")
              panel.levelplot(...)
            }))  
  print(levelplot(SBT.max~longitude*latitude|year_range, data=aggE, xlab=NULL, ylab=NULL, main=paste('SBT.max_rcp', rcp, '  season:',season, '  Model:', modelrun[i]), at = seq(min(aggE$SBT.max-.01, na.rm=T),max(aggE$SBT.max+.01, na.rm=T), length.out=40), 
            col.regions=cols, layout = c(5, 1), panel = function(...) {
              panel.fill(col = "light gray")
              panel.levelplot(...)
            }))  
  #print(levelplot(SST.max~longitude*latitude|year_range, data=aggE, xlab=NULL, ylab=NULL, main=paste('SST.max_rcp', rcp, '  season:',season, '  Model:', modelrun[i]), at = seq(min(aggE$SST.max-.1),max(aggE$SST.max+.1), length.out=40), 
   #         col.regions=cols, layout = c(5, 1), panel = function(...) {
    #          panel.fill(col = "light gray")
     #         panel.levelplot(...)
      #      }))  
  print(levelplot(SBT.seasonal~longitude*latitude|year_range, data=aggW, xlab=NULL, ylab=NULL, main=paste('sbt.seasonal_rcp', rcp, '  season:',season, '  Model:', modelrun[i]), at = seq(min(aggW$SBT.seasonal-.01, na.rm=T),max(aggW$SBT.seasonal+.01, na.rm=T), length.out=40), 
            col.regions=cols, layout = c(5, 1), panel = function(...) {
              panel.fill(col = "light gray")
              panel.levelplot(...)
            }))  
  #print(levelplot(SST.seasonal.mean~longitude*latitude|year_range, data=aggW, xlab=NULL, ylab=NULL, main=paste('SST.seasonal.mean_rcp', rcp, '  season:',season, '  Model:', modelrun[i]), at = seq(min(aggW$SST.seasonal.mean-.1),max(aggW$SST.seasonal.mean+.1), length.out=40), 
   #         col.regions=cols, layout = c(5, 1), panel = function(...) {
    #          panel.fill(col = "light gray")
     #         panel.levelplot(...)
      #      }))  
  print(levelplot(SBT.min~longitude*latitude|year_range, data=aggW, xlab=NULL, ylab=NULL, main=paste('SBT.min_rcp', rcp, '  season:',season, '  Model:', modelrun[i]), at = seq(min(aggW$SBT.min-.01, na.rm=T),max(aggW$SBT.min+.01, na.rm=T), length.out=40), 
            col.regions=cols, layout = c(5, 1), panel = function(...) {
              panel.fill(col = "light gray")
              panel.levelplot(...)
            }))  
  print(levelplot(SBT.max~longitude*latitude|year_range, data=aggW, xlab=NULL, ylab=NULL, main=paste('SBT.max_rcp', rcp, '  season:',season, '  Model:', modelrun[i]), at = seq(min(aggW$SBT.max-.01, na.rm=T),max(aggW$SBT.max+.01, na.rm=T), length.out=40), 
            col.regions=cols, layout = c(5, 1), panel = function(...) {
              panel.fill(col = "light gray")
              panel.levelplot(...)
            }))  
  #print(levelplot(SST.max~longitude*latitude|year_range, data=aggW, xlab=NULL, ylab=NULL, main=paste('SST.max_rcp', rcp, '  season:',season, '  Model:', modelrun[i]), at = seq(min(aggW$SST.max-.1),max(aggW$SST.max+.1), length.out=40), 
   #         col.regions=cols, layout = c(5, 1), panel = function(...) {
    #          panel.fill(col = "light gray")
     #         panel.levelplot(...)
      #      }))  
  dev.off()

  #filenameE <- paste('data/prediction_files_noBiasCorr/predictionEASTnoBias_rcp',rcp, '_', season, '_', modelrun[i], '.RData', sep='')
  #filenameW <- paste('data/prediction_files_noBiasCorr/predictionWESTnoBias_rcp',rcp, '_', season, '_', modelrun[i], '.RData', sep='')
  #save(pred.bathE, file=filenameE)
  #save(pred.bathW, file=filenameW)
}
  
# Plot the SODA SBT climatology 
filename = paste('data/Temp proj files/', pred.folder, 'soda_btm.txt', sep="")
clim <- read.table(filename)
clim <- data.frame(cbind(clim.grid, clim))

pdf(width=14, height=9, file='figures/TEST_3_SODA_climatology')
print(levelplot(V1~lonClimgrid*latClimgrid, data=clim, xlab=NULL, ylab=NULL, main='SODA climatology_SBT', at = seq(min(clim$V1-.01, na.rm=T),max(clim$V1+.01, na.rm=T), length.out=40), 
                col.regions=cols))
dev.off()

# To plot a specific area for troubleshooting_can delete this later 
pred2 <- pred[pred$year > 2080,]
print(levelplot(SBT.seasonal~lonClimgrid*latClimgrid|year, data=pred2, ylim=c(55,61), xlim=c(-172,-142), xlab=NULL, ylab=NULL, main=paste('sbt.seasonal_rcp', rcp, '  season:',season, '  Model:', modelrun[i]), at = seq(0,18, length.out=40), col.regions=cols))
 
cols = colorRampPalette(colors = c('dark blue', 'blue', 'white', 'red', 'dark red'))
depths <- proj.grid[proj.grid$longrid > -172 & proj.grid$longrid < -142 & proj.grid$latgrid > 55 & proj.grid$latgrid < 61,]
levelplot(depth~longrid*latgrid, at = seq(min(depths$depth-.1),max(depths$depth+.1), length.out=40), col.regions=cols, data=depths)
 

# =======================================================================================
# CODE BELOW TO LOOK AT PROJECTION DATA A DIFFERENT WAY THAN THE ABOVE 20yr. SUMMARIES_this was used mainly for QAQC
# =======================================================================================
# Below are four blocks that can be recycled to bring in any of the data_adjust 'i' for model runs_this could be condensed if necessary
proj.clim <- unique(data.frame(lonClimgrid=proj.grid$lonClimgrid, latClimgrid=proj.grid$latClimgrid))
# seasonal mean sbt
  filename = paste('data/', pred.folder[2], modelrun[3], '_rcp', rcp, '_fill.nc_2006_2100_', season,'_', pred.metric[3], '.txt', sep="")
  filein <- readLines(filename)
  filein <- t(sapply(filein, processlinesMean))
  rownames(filein) <- NULL
   
# seasonal mean sst
  filename = paste('data/', pred.folder[1], modelrun[i], '_rcp', rcp, '_r1i1p1_1950_2100.nc_regrid.nc_fill.nc_2006_2100_', season,'_', pred.metric[3], '.txt', sep="")
  filein <- readLines(filename)
  filein <- t(sapply(filein, processlinesMean))
  rownames(filein) <- NULL

# min or max sbt_for these two need to adjust 'i' but also 'pred.metric' can be 1 or 2 (max or min)
  filename = paste('data/', pred.folder[2], modelrun[1], '_rcp', rcp, '_fill.nc_2006_2100_', season,'_', pred.metric[1], '.txt', sep="")
  filein <- readLines(filename)
  filein <- t(sapply(filein, processlinesMinMax))
  rownames(filein) <- NULL

# min or max sst
  filename = paste('data/', pred.folder[1], modelrun[i], '_rcp', rcp, '_r1i1p1_1950_2100.nc_regrid.nc_fill.nc_2006_2100_', season,'_', pred.metric[2], '.txt', sep="")
  filein <- readLines(filename)
  filein <- t(sapply(filein, processlinesMinMax))
  rownames(filein) <- NULL

# Take a look at some temperature values_I recycle this code to look at different combinations of model, year, etc.
abc <- data.frame(cbind(filein, clim.grid)) # attach lat and lon of grid cells
abc <- merge(proj.clim, abc, by=c('lonClimgrid', 'latClimgrid'), all.x=T) # trim out areas we will not project to
summary(abc)

plot(latClimgrid~lonClimgrid, cex=.1, data=abc[!is.na(abc$V90),])
points(latClimgrid~lonClimgrid, col='red', cex=.5, data=abc[abc$V90 > 38, ])
points(latClimgrid~lonClimgrid, col='blue', data=abc[abc$V90 < -5, ])
levelplot(V92 ~ lonClimgrid * latClimgrid, data = abc)
# another plot option for showing just the outier areas well
forPlot <- abc$X89# choose a given year to work with the plots_this is only necessary for establishing the cutpoints
cutpts <- c(min(forPlot, na.rm=T)-1, -10, -5, 40, 80, max(forPlot, na.rm=T)+1) # add some cutpoints to show outliers_gray is relatively normal, with two levels of both negative and positive outliers
levelplot(X89 ~ lonClimgrid * latClimgrid, data = abc[!is.na(abc$X89),], at = cutpts, cuts = 6, pretty = T, col.regions = (rev(brewer.pal(9, "RdBu")))) 

# =======================================================================================
# Making file that stores the region for each species' prediction
# =======================================================================================

# Make a file that lists, for each species, the dominant 'regionfact' to use with predictions
load('data/dat_selectedspp_Feb_1_2017.Rdata')# load species catch data
dat$region[dat$region %in% c("NEFSC_NEUSFall","NEFSC_NEUSSpring")] <- "NEFSC_NEUS" #NEFSC as one region in "region"
dat$region[dat$region %in% c("DFO_NewfoundlandFall","DFO_NewfoundlandSpring")] <- "DFO_Newfoundland"
dat$region[dat$region %in% c("DFO_ScotianShelfFall","DFO_ScotianShelfSpring", "DFO_ScotianShelfSummer")] <- "DFO_ScotianShelf"
dat$region[dat$region %in% c("SCDNR_SEUSFall","SCDNR_SEUSSpring", "SCDNR_SEUSSummer")] <- "SCDNR_SEUS"
dat$region[dat$region %in% c("SEFSC_GOMexFall","SEFSC_GOMexSummer")] <- "SEFSC_GOMex"
dat$region[dat$region %in% c("VIMS_NEAMAPFall","VIMS_NEAMAPSpring")] <- "VIMS_NEAMAP"

regionFreq <- as.character()
for(i in 1:length(projspp)){
  species <- projspp[i]
  datTrim <- dat[dat$sppocean == species,]
  regionFreqs <- data.frame(table(datTrim$region))
  regionfact <- as.character(regionFreqs$Var1[regionFreqs$Freq == max(regionFreqs$Freq)])[1] # the '[1]' at the end is in case of a tie
  regionFreq[i] <- regionfact
}
rm(species, datTrim, regionFreqs, regionfact)
save(projspp, regionFreq, file='data/speciesProjectionList.RData') 

# ========================================================================
# Sums up projections of all species from script 7b
# ========================================================================

# Make figure of projections in 20 year blocks for each species
pdf(width=14, height=5, file=paste('figures/speciesProjections/', 'species_projections3.pdf', sep=''))
for(i in 401:658){
  load(paste(projfolder, '_', projspp[i], '_', modelrun[6], '_', season, '_', rcp, '_prediction_AGG.RData', sep=''))
  Projmap <- map('world', xlim=c(min(pred.agg$longitude),max(pred.agg$longitude)), ylim=c(min(pred.agg$latitude),max(pred.agg$latitude)),plot=FALSE) 
  Projmap <- data.frame(lon=Projmap$x, lat=Projmap$y)
  cols = colorRampPalette(colors = c('gray90', 'blue', 'dark blue', 'black'))
  #pred.agg$logPred <- log(pred.agg$pred.mean + 1)
  print(levelplot(pred.mean ~ longitude*latitude|year_range, data=pred.agg, xlab=NULL, ylab=NULL, main=paste('rcp', rcp, 'season' ,season, 'Model', modelrun[6], projspp[i], sep='_'), 
      at = seq(0,max(pred.agg$pred.mean+.00001), length.out=20), col.regions=cols, layout = c(5, 1)) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=.05, col='dark gray')) 
  #levelplot(logPred ~ longitude*latitude|year_range, data=pred.agg, xlab=NULL, ylab=NULL, main=paste('rcp', rcp, 'season' ,season, 'Model', modelrun[i], projspp[i], 'LOGGED', sep='_'), 
  #          at = seq(0,max(pred.agg$logPred+.0001), length.out=20), col.regions=cols, layout = c(5, 1)) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=.05, col='dark gray') 
}
dev.off()
 
# Figs for changes in thermal habitat and also centroid vectors
centPreds2 <- data.frame(centLatBase=numeric(), centLonBase=numeric(), centLatProj=numeric(), centLonProj=numeric(), habChange=numeric())
pdf(width=10, height=11, file=paste('figures/speciesProjections/', 'species_thermalHabitatChange2.pdf', sep=''))
for(i in 1:length(projspp)){
  load(paste(projfolder, '_', projspp[i], '_', modelrun[6], '_', season, '_', rcp, '_prediction_AGG.RData', sep=''))
  predBase <- pred.agg[pred.agg$year_range == '2007-2020',]
  predProj <- pred.agg[pred.agg$year_range == '2081-2100',]
  thermHab <- data.frame(cbind(latitude = predBase$latitude, longitude = predBase$longitude, predBase = predBase$pred.mean, predProj = predProj$pred.mean))
  thermHab <- thermHab[!(thermHab$latitude > 60.5 & thermHab$longitude > -64.62603),]
  thermHab <- thermHab[!(thermHab$latitude > 55 & thermHab$longitude < -64.6 & thermHab$longitude > -75),]
  thermHab$diff <- thermHab$predProj - thermHab$predBase
  thermHab$diffProp <- thermHab$diff/max(thermHab$predProj, thermHab$predBase)
  # calculate centroid for initial and end periods_also calculate latitudinal centroid
  # Do compass plots in all the major regions, with a species falling into a region based on freq. of capture
  # Do separate for GOM?_Would need to be based on if species was caught there more than x times in survey I think.
  centLatBase <- wtd.mean(thermHab$latitude, weights=thermHab$predBase) 
  centLonBase <- wtd.mean(thermHab$longitude, weights=thermHab$predBase) 
  centLatProj <- wtd.mean(thermHab$latitude, weights=thermHab$predProj) 
  centLonProj <- wtd.mean(thermHab$longitude, weights=thermHab$predProj) 
    
  Projmap <- map('world', xlim=c(min(pred.agg$longitude),max(pred.agg$longitude)), ylim=c(min(pred.agg$latitude),max(pred.agg$latitude)),plot=FALSE) 
  Projmap <- data.frame(lon=Projmap$x, lat=Projmap$y)
  cols = colorRampPalette(colors = c('dark red', 'red', 'gray90', 'blue', 'dark blue'))
  print(levelplot(diffProp ~ longitude*latitude, data=thermHab, xlab=NULL, ylab=NULL, main=paste('rcp', rcp, 'season' ,season, 'Model', modelrun[6], projspp[i], sep='_'), 
                  at = seq(-1, 1, length.out=20), col.regions=cols, panel = function(...) {
              panel.levelplot(...)
              grid.points(x=c(centLonBase), y=c(centLatBase), pch = 16)
              larrows(x0=centLonBase, y0=centLatBase, x1=centLonProj, y1=centLatProj, length=.1)
            }) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=.5, col='dark gray'))
  habChange = 100*(sum(thermHab$predProj)/sum(thermHab$predBase) - 1) # this is predicted change in biomass_as a percentage of the baseline period
  centPreds2[i,] <- data.frame(centLatBase=centLatBase, centLonBase=centLonBase, centLatProj=centLatProj, centLonProj=centLonProj, habChange=habChange)
}
dev.off()
   
centPreds2$species <- projspp # add species names
# Some other values that may be eventually useful for compass plots_still haven't figured out how to do them
centPreds$vectDist <- distHaversine(p1=matrix(cbind(centPreds$centLonBase, centPreds$centLatBase), ncol=2), p2=matrix(cbind(centPreds$centLonProj, centPreds$centLatProj), ncol=2)) / 1000# distance of vectors in km (I think in km anyway)
centPreds$vectAng <- bearingRhumb(p1=matrix(cbind(centPreds$centLonBase, centPreds$centLatBase), ncol=2), p2=matrix(cbind(centPreds$centLonProj, centPreds$centLatProj), ncol=2))
centPreds$radians <- NISTdegTOradian(centPreds$vectAng) # convert angle to radians_I think zero is the right x-axis
centPreds$slope <- (centPreds$centLatProj - centPreds$centLatBase)/(centPreds$centLonProj - centPreds$centLonBase)# This probably is not appropriate

save(centPreds, file='output/centPreds.RData')

# Plotting all the species vectors together
load('data/ProjectionBathGrid_Feb27_2017.RData')# load bathymetry projection grid
Projmap <- map('world', xlim=c(min(proj.grid$lonBathgrid),max(proj.grid$lonBathgrid)), ylim=c(min(proj.grid$latBathgrid),max(proj.grid$latBathgrid)),plot=FALSE) 
Projmap <- data.frame(lon=Projmap$x, lat=Projmap$y)
pdf(width=14, height=8, file=paste('figures/speciesProjections/', 'species_centroid_shifts2.pdf', sep=''))
proj.grid2 <- proj.grid[!(proj.grid$latClimgrid > 60.5 & proj.grid$lonClimgrid > -64.62603),]
proj.grid2 <- proj.grid2[!(proj.grid2$latClimgrid > 55 & proj.grid2$lonClimgrid < -64.6 & proj.grid2$lonClimgrid > -75),]
plot(latClimgrid~lonClimgrid, col='gray90', cex=.7, pch=16, xlab='Longitude', ylab='Latitude', main='Predicted centroid shifts for 658 species', data=proj.grid2) 
points(lat ~ lon, Projmap, type='l', lty=1, lwd=.7, col='gray50')
arrows(x0=centPreds2$centLonBase, y0=centPreds2$centLatBase, x1=centPreds2$centLonProj, y1=centPreds2$centLatProj, length=.05, lwd=.8, col='dark blue')
dev.off()

# ==============================================================================
# Figures of thermal habitat change for MAFMC_June 2017
# Figs for changes in thermal habitat and also centroid vectors
# ==============================================================================
for(i in 1:length(mafmc)){
  species = mafmc[i] 
  print(paste("Starting figures for ", species))
  pdf(width=10, height=11, file=paste('figures/MAFMC/', species, '_thermalHabitatChange.pdf', sep=''))
  for(k in 4:19){
    load(paste(projfolder, species, '_rcp', rcp[1], '_jas_prediction_AGG.RData', sep=''))
    predBase <- pred.agg[pred.agg$year_range == '2007-2020',]
    predProj <- pred.agg[pred.agg$year_range == '2081-2100',]
    thermHab <- data.frame(cbind(latitude = predBase$latitude, longitude = predBase$longitude, predBase = predBase[,k], predProj = predProj[,k]))
    thermHab <- thermHab[!(thermHab$latitude > 60.5 & thermHab$longitude > -64.62603),]
    thermHab <- thermHab[!(thermHab$latitude > 55 & thermHab$longitude < -64.6 & thermHab$longitude > -75),]
    thermHab <- thermHab[!thermHab$longitude < -82,]
    thermHab <- thermHab[!thermHab$latitude < 26.6,]
    thermHab <- thermHab[!(thermHab$latitude > 55 & thermHab$longitude < -64),]
    
    thermHab$diff <- thermHab$predProj - thermHab$predBase
    thermHab$diffProp <- thermHab$diff/max(thermHab$predProj, thermHab$predBase)
    # calculate centroid for initial and end periods_also calculate latitudinal centroid
    # Do compass plots in all the major regions, with a species falling into a region based on freq. of capture
    # Do separate for GOM?_Would need to be based on if species was caught there more than x times in survey I think.
    centLatBase <- wtd.mean(thermHab$latitude, weights=thermHab$predBase) 
    centLonBase <- wtd.mean(thermHab$longitude, weights=thermHab$predBase) 
    centLatProj <- wtd.mean(thermHab$latitude, weights=thermHab$predProj) 
    centLonProj <- wtd.mean(thermHab$longitude, weights=thermHab$predProj) 
    
    Projmap <- map('world', xlim=c(min(pred.agg$longitude),max(pred.agg$longitude)), ylim=c(min(pred.agg$latitude),max(pred.agg$latitude)),plot=FALSE) 
    Projmap <- data.frame(lon=Projmap$x, lat=Projmap$y)
    cols = colorRampPalette(colors = c('dark red', 'red', 'gray90', 'blue', 'dark blue'))
    print(levelplot(diffProp ~ longitude*latitude, data=thermHab, xlab=NULL, ylab=NULL, main=paste('rcp', rcp[1], 'season_jas_Model', modelrun[k-3], species, sep='_'), 
                    at = seq(-1, 1, length.out=20), col.regions=cols, panel = function(...) {
                      panel.levelplot(...)
                      grid.points(x=c(centLonBase), y=c(centLatBase), pch = 16)
                      larrows(x0=centLonBase, y0=centLatBase, x1=centLonProj, y1=centLatProj, length=.1)
                    }) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=.5, col='dark gray'))
    #habChange = 100*(sum(thermHab$predProj)/sum(thermHab$predBase) - 1) # this is predicted change in biomass_as a percentage of the baseline period
    #centPreds2[i,] <- data.frame(centLatBase=centLatBase, centLonBase=centLonBase, centLatProj=centLatProj, centLonProj=centLonProj, habChange=habChange)
     
    load(paste(projfolder, species, '_rcp', rcp[2], '_jas_prediction_AGG.RData', sep=''))
    predBase <- pred.agg[pred.agg$year_range == '2007-2020',]
    predProj <- pred.agg[pred.agg$year_range == '2081-2100',]
    thermHab <- data.frame(cbind(latitude = predBase$latitude, longitude = predBase$longitude, predBase = predBase[,k], predProj = predProj[,k]))
    thermHab <- thermHab[!(thermHab$latitude > 60.5 & thermHab$longitude > -64.62603),]
    thermHab <- thermHab[!(thermHab$latitude > 55 & thermHab$longitude < -64.6 & thermHab$longitude > -75),]
    thermHab <- thermHab[!thermHab$longitude < -82,]
    thermHab <- thermHab[!thermHab$latitude < 26.6,]
    thermHab <- thermHab[!(thermHab$latitude > 55 & thermHab$longitude < -64),]
    
    thermHab$diff <- thermHab$predProj - thermHab$predBase
    thermHab$diffProp <- thermHab$diff/max(thermHab$predProj, thermHab$predBase)
    # calculate centroid for initial and end periods_also calculate latitudinal centroid
    # Do compass plots in all the major regions, with a species falling into a region based on freq. of capture
    # Do separate for GOM?_Would need to be based on if species was caught there more than x times in survey I think.
    centLatBase <- wtd.mean(thermHab$latitude, weights=thermHab$predBase) 
    centLonBase <- wtd.mean(thermHab$longitude, weights=thermHab$predBase) 
    centLatProj <- wtd.mean(thermHab$latitude, weights=thermHab$predProj) 
    centLonProj <- wtd.mean(thermHab$longitude, weights=thermHab$predProj) 
    
    Projmap <- map('world', xlim=c(min(pred.agg$longitude),max(pred.agg$longitude)), ylim=c(min(pred.agg$latitude),max(pred.agg$latitude)),plot=FALSE) 
    Projmap <- data.frame(lon=Projmap$x, lat=Projmap$y)
    cols = colorRampPalette(colors = c('dark red', 'red', 'gray90', 'blue', 'dark blue'))
    print(levelplot(diffProp ~ longitude*latitude, data=thermHab, xlab=NULL, ylab=NULL, main=paste('rcp', rcp[2], 'season_jas_Model', modelrun[k-3], species, sep='_'), 
                    at = seq(-1, 1, length.out=20), col.regions=cols, panel = function(...) {
                      panel.levelplot(...)
                      grid.points(x=c(centLonBase), y=c(centLatBase), pch = 16)
                      larrows(x0=centLonBase, y0=centLatBase, x1=centLonProj, y1=centLatProj, length=.1)
                    }) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=.5, col='dark gray'))
    
  }
  dev.off()
}

# =================================================================================
# Below here is Malin's original code
# =================================================================================

load('data/ProjectionBathGrid_Feb27_2017.RData')# load projection grid 

if(!file.exists(paste(climgridfolder, 'climGrid_rcp', rcp, '.proj2_wrugos.RData', sep=''))){
	print(paste('climGrid with rugosity does not exist for RCP', rcp, '. Making it.', sep=''))
	
	load(paste(climgridfolder, 'climGrid_rcp', rcp, '.proj2.RData', sep='')) # projected temperature for each year ("clim")

	# drop unneeded columns
	clim <- clim[,!grepl('depthgrid', names(clim))] #  refer to GCM depth grids
	clim <- clim[,!grepl('bottemp.clim|surftemp.clim|delta|latgrid|longrid', names(clim))] #  the temp climatologies, deltas, and GCM lat/lon grids (1 degree)

	# add regionfact
	clim$region<- as.factor(clim$region)
	names(clim)[names(clim)=='region'] <- 'regionfact'

	# add logrugosity
	rugos <- read.csv('data/projectiongrid_latlons.1.16th_withRugosity_2015-05-06.csv')
		names(rugos)[names(rugos) == 'lon'] <- 'lon16th'
		names(rugos)[names(rugos) == 'lat'] <- 'lat16th'
		names(rugos)[names(rugos) == 'depth'] <- 'depth16th'

		gridsize=0.25 # size of grid of the climate data, in degrees
		rugos$lat <- floor(rugos$lat16th/gridsize)*gridsize + gridsize/2 # round to nearest grid center
		rugos$lon <- floor(rugos$lon16th/gridsize)*gridsize + gridsize/2

	clim <- merge(clim, rugos) # slow
	dim(clim) # 9623120 rows
	
	save(clim, file=paste(climgridfolder, 'climGrid_rcp', rcp, '.proj2_wrugos.RData', sep=''))
	
} else {
	print('climGrid with rugosity exists. Loading it')
	load(paste(climgridfolder, 'climGrid_rcp', rcp, '.proj2_wrugos.RData', sep=''))
}

############################################
## Project GAMS onto annual climate data  ##
############################################
	
options(warn=1) # print warnings as they occur

# thisprojspp <- projspp[1]

# VERY slow
doprojection <- function(thisprojspp, files, clim, projfolder, modfolder, runtype, stayinregion=TRUE){ 
	# the stayinregion flag determines whether or not we project species outside of the regions in which they were observed historically (and to which models were fit)

	# clear variables
	mods <- avemeanbiomass <- NULL
	
	# load model fits (mod and avemeanbiomass)
	fileindex <- which(grepl(gsub('/|\\(|\\)', '', thisprojspp), gsub('/|\\(|\\)', '', files)))
	print(paste(fileindex, thisprojspp, Sys.time()))

	load(paste(modfolder, '/', files[fileindex], sep='')) # loads mods and avemeanbiomass

	fitregions <- gsub('regionfact', '', grep('regionfact', names(coef(mods$mygam1)), value=TRUE)) # regions included in the model fit (if more than one region fit)
	if(length(fitregions)==0) fitregions <- names(avemeanbiomass)[which.max(avemeanbiomass)] # else, if only one region fit (and no regionfact term included), then use avemeanbiomass to pull out the name of that one region

	if(stayinregion){ # if we don't allow species to move into new regions, then we can use the average observed biomass for each region
		# add mean biomass by region
		clim$biomassmean <- 0
		clim$biomassmean[clim$regionfact %in% names(avemeanbiomass)] <- avemeanbiomass[as.character(clim$regionfact[clim$regionfact %in% names(avemeanbiomass)])] # use region to pull the correct mean biomass values
	}
	if(!stayinregion){ #else, add a standard mean biomass across all regions (use region with highest mean biomass that was in the model fit)
		clim$biomassmean <- max(avemeanbiomass[fitregions])
	}

	# smearing estimator for re-transformation bias (see Duan 1983, http://www.herc.research.va.gov/resources/faq_e02.asp)
	smear <- mean(exp(mods[['mygam2']]$residuals))

	# rows of clim to project to
	regstoproj <- names(avemeanbiomass)
	if(any(regstoproj == 'NEFSC_NEUS')){
		regstoproj <- c(regstoproj, 'NEFSC_NEUSSpring', 'NEFSC_NEUSFall') # add spring and fall surveys (how they are marked in the climatology) to the list. The models are fit to all NEUS data and don't distinguish seasons.
	}
	if(stayinregion){ # only project to regions in which this species had a climate envelope fit
		inds <- clim$regionfact %in% regstoproj
	}
	if(!stayinregion){ # project to all regions in the same ocean (Atlantic or Pacific)
		if(any(regstoproj %in% c("DFO_NewfoundlandFall", "DFO_NewfoundlandSpring", "DFO_ScotianShelf", "DFO_SoGulf", "NEFSC_NEUSFall", "NEFSC_NEUSSpring", "SEFSC_GOMex"))){ # Atlantic
			regstoproj <- c("DFO_NewfoundlandFall", "DFO_NewfoundlandSpring", "DFO_ScotianShelf", "DFO_SoGulf", "NEFSC_NEUSFall", "NEFSC_NEUSSpring", "SEFSC_GOMex")
		}
		if(any(regstoproj %in% c("AFSC_Aleutians", "AFSC_EBS", "AFSC_GOA", "AFSC_WCTri", "NWFSC_WCAnn"))){ # Pacific
			regstoproj <- c("AFSC_Aleutians", "AFSC_EBS", "AFSC_GOA", "AFSC_WCTri", "NWFSC_WCAnn")
		}
		inds <- clim$regionfact %in% regstoproj
	}
	

	# Dataframe for this species' projections
	thisproj <- clim[inds,c('regionfact', 'lat', 'lon', 'depth16th', 'year')] # dataframe to hold projections for this taxon
	for(i in 1:13) thisproj[[paste('wtcpue.proj_', i, sep='')]] <- NA 	# Add projected biomass density columns for each climate model

	# Add region factor to use in model during projection
	if(stayinregion){
		clim$regiontoproj <- as.character(clim$regionfact)
		# Adjust region names for NEUS
		if(any(regstoproj %in% 'NEFSC_NEUS')){
			clim$regiontoproj[clim$regiontoproj %in% c('NEFSC_NEUSSpring', 'NEFSC_NEUSFall')] <- 'NEFSC_NEUS'
		}
	}
	if(!stayinregion){ 
		clim$regiontoproj <- fitregions[which.max(avemeanbiomass[fitregions])] # pick the region with the highest biomass. this ensures that it had at least some presences, and so would have been in the model fit
	}

	# Calculate predictions for 2020-2100 for each model
	for(i in 1:13){ # this alone takes a long time
		print(paste(fileindex, 'model', i))
		nd <- data.frame(regionfact = clim$regiontoproj[inds], surftemp = clim[[paste('surftemp.proj_', i, sep='')]][inds], bottemp = clim[[paste('bottemp.proj_', i, sep='')]][inds], logrugosity = log(clim$rugosity[inds]+0.01), biomassmean = clim$biomassmean[inds], row.names=1:sum(inds))
		preds1 <- predict.gam(mods$mygam1, newdata = nd, type='response')
		preds2 <- exp(predict(mods$mygam2, newdata = nd, type='response'))
		preds <- preds1*preds2*smear
		preds[preds<0] <- 0
		
		pnm = paste('wtcpue.proj_', i, sep='')
		thisproj[[pnm]] = preds # can do this because clim[inds,] and thisproj are in the same order

		# set biomass projections on land to zero
		thisproj[[pnm]][thisproj$depth16th > 0] <- 0
	}
	
	# summarize by 1/4 degree grid
	summproj <- aggregate(thisproj[,grepl('wtcpue.proj', names(thisproj))], by=list(region=thisproj$regionfact, lat=thisproj$lat, lon=thisproj$lon, year=thisproj$year), FUN=mean) # also slow
	
#	print(summary(summproj))
#	print(dim(summproj))	

	thisprojspp <- gsub('/', '', thisprojspp) # would mess up saving the file if the species name had /
	outfile <- paste(projfolder, '/summproj_', runtype, '_rcp', rcp, '_', thisprojspp, '.Rdata', sep='')
	if(!stayinregion) outfile <- paste(projfolder, '/summproj_', runtype, '_xreg_rcp', rcp, '_', thisprojspp, '.Rdata', sep='') # append xreg to projections if species could be projected out of their observed regions
	save(summproj, file=outfile) # write out the projections (15MB file)
}

result <- mclapply(X= projspp, FUN=doprojection, files=files, clim=clim, projfolder=projfolder, modfolder=modfolder, runtype=runtype, stayinregion=stayinregion, mc.cores=numcorestouse) # spawn out to multiple cores. errors will be stored in results

print(result)
