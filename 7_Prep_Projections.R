# First, examine habitat model fits and decide what species to include for projections
# Second, bring in climate projection data and prep for predictions with GAMs
    
## Set working directory
if(Sys.info()["nodename"] == "amphiprion.deenr.rutgers.edu"){
	setwd('~/Documents/range_projections/')
	climgridfolder <- 'data/'
}
if(Sys.info()["user"] == "jamesmorley"){
  setwd('/Users/jamesmorley/Documents/project_velocity')
  climgridfolder <- 'data/'
}

library(latticeExtra)
library(maps)
library(geosphere)

#############################
# Choose species to project #
#############################

runtype <- 'fitallreg_2017'
# 703 species had habitat models successfully fit for them_none failed to converge (although a couple species took 2 tries)
load(paste('output/modeldiag_Nov2017_', runtype, '.Rdata', sep='')) # model diagnostics
# 3 species had zero observations in either the training or test data sets so could not evaluate with this framework
modeldiag <- modeldiag[!is.na(modeldiag$auc.tt),] # 700 species

#With additional criteria for inclusion
modeldiag <- modeldiag[!modeldiag$auc.tt < 0.75,] # 686 species

projspp <- modeldiag$sppocean 
# save(modeldiag, file='temp_prefs/sppList.RData') # For estimating temp prefs for all spp (e.g. Patrick's paper)

# Some graphing options for looking at model diagnostics
#boxplot <- data.frame(cbind(Presence_models = modeldiag$dev.pres, Biomass_models = modeldiag$dev.biomass)) 
#boxplot <- melt(boxplot)
#ggplot(boxplot, aes(x=factor(variable), y=value)) + geom_boxplot() + labs(x='', y='%deviance explained')

#	hist(modeldiag$auc.tt)
#	hist(modeldiag$auc)
#	hist(modeldiag$dev.pres - modeldiag$dev.pres.null); abline(v=0.05, col='red')
#	hist(modeldiag$dev.biomass - modeldiag$dev.biomass.null); abline(v=0.05, col='red')

# =======================================================================================
# Making file that stores the region for each species' prediction
# =======================================================================================

# Make a file that lists, for each species, the dominant 'regionfact' to use with predictions
load('data/dat_selectedspp_Feb_1_2017.Rdata')# load species catch data
dat <- dat[!(dat$wtcpue == 0 & dat$region == 'DFO_SoGulf'),] # the zeros in SoGulf are actual zeros (ie not just a scale issue) and thus are true absences
dat$sppocean[dat$sppocean=='heppasteria phygiana_Atl'] <- 'hippasteria phrygiana_Atl' # no need to aggregate as they never overlapped w/n hauls (only found in Newfoundland and they called it one way or the other)

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
 
###############################################
# Choose projection options_rerun this code for different rcp and seasons
###############################################
 
#rcp <- 26
#rcp <- 45
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
   
modelrun <- c('bcc-csm1-1-m','bcc-csm1-1','CanESM2','CCSM4','CESM1-CAM5','CNRM-CM5','GFDL-CM3','GFDL-ESM2M','GFDL-ESM2G','GISS-E2-R','GISS-E2-H','HadGEM2-ES','IPSL-CM5A-LR','IPSL-CM5A-MR','MIROC-ESM','MIROC5','MPI-ESM-LR','NorESM1-ME')
modelrun <- c('HadGEM2-ES','MIROC5')

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

if(rcp==85){
  pred.folder <- c('sst_rcp85_final/tos_Omon_','sbt_rcp85_final/')
} 
if(rcp==45){
  pred.folder <- c('sst_rcp45_final/tos_Omon_','sbt_rcp45_final/')
}
if(rcp==26){
  pred.folder <- c('sst_rcp26_final/tos_Omon_','sbt_rcp26_final/')
}

# For reading in all models for a particular type of data_minimum surface temps was not in habitat models
tempsmeansbt <- matrix(as.numeric(NA), nrow=1281878, ncol=length(modelrun))
tempsmeansst <- matrix(as.numeric(NA), nrow=1281878, ncol=length(modelrun))
tempsminsbt <- matrix(as.numeric(NA), nrow=1281878, ncol=length(modelrun))
tempsmaxsbt <- matrix(as.numeric(NA), nrow=1281878, ncol=length(modelrun))
tempsmaxsst <- matrix(as.numeric(NA), nrow=1281878, ncol=length(modelrun))

# The loops below produce matrix where each column is a different climate model and the years are stacked in one vector
# seasonal mean sbt 
for(i in 1:length(modelrun)){
  print(i) 
  filename = paste('data/Temp proj files/', pred.folder[2], modelrun[i], '_rcp', rcp, '_fill.nc_2006_2100_', season,'_', pred.metric[3], '.txt', sep="")
  filein <- readLines(filename)
  filein <- t(sapply(filein, processlinesMean))
  rownames(filein) <- NULL
  print(paste(i, modelrun[i], dim(filein)[1], dim(filein)[2]))
  tempsmeansbt[,i] <- as.vector(filein)
} 
# seasonal mean sst
for(i in 1:length(modelrun)){
  print(i) 
  filename = paste('data/Temp proj files/', pred.folder[1], modelrun[i], '_rcp', rcp, '_r1i1p1_1950_2100.nc_regrid.nc_fill.nc_2006_2100_', season,'_', pred.metric[3], '_no_scaling.txt', sep="")
  filein <- readLines(filename)
  filein <- t(sapply(filein, processlinesMean))
  rownames(filein) <- NULL
  print(paste(i, modelrun[i], dim(filein)[1], dim(filein)[2]))
  tempsmeansst[,i] <- as.vector(filein)
} 
# min annual sbt
for(i in 1:length(modelrun)){
  print(i)
  filename = paste('data/Temp proj files/', pred.folder[2], modelrun[i], '_rcp', rcp, '_fill.nc_2006_2100_', season,'_', pred.metric[2], '.txt', sep="")
  #filenameB = paste('data/Temp proj files/rcp45_nofill.txt')
  filein <- readLines(filename)
  #fileinB <- readLines(filenameB)
  #fileinB <- fileinB[!fileinB==""]
  #abc <- fileinB[218193:231828]
  #length(abc) 
  filein <- t(sapply(filein, processlinesMinMax))
  #abc <- t(sapply(abc, processlinesMinMax))
  rownames(filein) <- NULL
  #rownames(abc) <- NULL
  #summary(abc)
  #cde <- data.frame(cbind(abc, clim.grid[-13637,]))
  #levelplot(X2~longrid*latgrid, data=cde)
  #plot(latgrid~longrid, cex=.1, cde)
  #points(latgrid~longrid, col='red', cex=.1, cde[is.na(cde$X2),])
  fileinVec <- as.vector(filein)
  fileinVec <- fileinVec[1:1281878]
  print(paste(i, modelrun[i], dim(filein)[1], dim(filein)[2]))
  tempsminsbt[,i] <- fileinVec
}
# max annual sbt
for(i in 1:length(modelrun)){#length(modelrun)){
  print(i)
  filename = paste('data/Temp proj files/', pred.folder[2], modelrun[i], '_rcp', rcp, '_fill.nc_2006_2100_', season,'_', pred.metric[1], '.txt', sep="")
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
  print(i)
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
 
# Build and save all 18 prediction data sets_make some summary figures in 20 year bins
for(i in 1:length(modelrun)){ 
  print(paste('Beginning model number ', i, ': ', modelrun[i], sep=''))
  pred <- data.frame(cbind(clim.grid, year=rep(2007:2100, rep=94, each=13637), SBT.seasonal = tempsmeansbt[,i], SBT.min = tempsminsbt[,i], SBT.max = tempsmaxsbt[,i], SST.seasonal.mean = tempsmeansst[,i], SST.max = tempsmaxsst[,i]))
  pred <- merge(pred, year.bin, by='year', all.x=T)  # Bin into 20 year periods to make some plots
  pred.bathE <- merge(pred, proj.grid[proj.grid$lonClimgrid > -105,], by=c('lonClimgrid', 'latClimgrid'), all.y = T, sort=F)
  pred.bathW <- merge(pred, proj.grid[proj.grid$lonClimgrid < -105,], by=c('lonClimgrid', 'latClimgrid'), all.y = T, sort=F)
  pred.bathE <- pred.bathE[order(pred.bathE$year, pred.bathE$index),] # reorder by year and index
  pred.bathW <- pred.bathW[order(pred.bathW$year, pred.bathW$index),] # reorder by year and index
  
  aggE <- aggregate(list(SBT.seasonal = pred.bathE$SBT.seasonal, SST.seasonal.mean = pred.bathE$SST.seasonal.mean, SBT.min = pred.bathE$SBT.min,
            SBT.max = pred.bathE$SBT.max, SST.max = pred.bathE$SST.max), by=list(year_range=pred.bathE$bin, latitude=pred.bathE$latBathgrid, longitude=pred.bathE$lonBathgrid), FUN=mean)
  aggW <- aggregate(list(SBT.seasonal = pred.bathW$SBT.seasonal, SST.seasonal.mean = pred.bathW$SST.seasonal.mean, SBT.min = pred.bathW$SBT.min,
            SBT.max = pred.bathW$SBT.max, SST.max = pred.bathW$SST.max), by=list(year_range=pred.bathW$bin, latitude=pred.bathW$latBathgrid, longitude=pred.bathW$lonBathgrid), FUN=mean)
 
  # 20yr. summary figures for each model
  pdf(width=14, height=5, file=paste('figures/Temp_projections_Nov6_2017/', modelrun[i], '_rcp', rcp, '_20yr_tempProj.pdf', sep=''))

  print(levelplot(SBT.seasonal~longitude*latitude|year_range, data=aggE, xlab=NULL, ylab=NULL, main=paste('sbt.seasonal_rcp', rcp, '  season:',season, '  Model:', modelrun[i]), at = seq(min(aggE$SBT.seasonal-.1, na.rm=T),max(aggE$SBT.seasonal+.1, na.rm=T), length.out=40), 
            col.regions=cols, layout = c(5, 1), panel = function(...) {
              panel.fill(col = "light gray")
              panel.levelplot(...)
            }))  
  print(levelplot(SST.seasonal.mean~longitude*latitude|year_range, data=aggE, xlab=NULL, ylab=NULL, main=paste('SST.seasonal.mean_rcp', rcp, '  season:',season, '  Model:', modelrun[i]), at = seq(min(aggE$SST.seasonal.mean-.1),max(aggE$SST.seasonal.mean+.1), length.out=40), 
            col.regions=cols, layout = c(5, 1), panel = function(...) {
              panel.fill(col = "light gray")
              panel.levelplot(...)
            }))  
  print(levelplot(SBT.min~longitude*latitude|year_range, data=aggE, xlab=NULL, ylab=NULL, main=paste('SBT.min_rcp', rcp, '  season:',season, '  Model:', modelrun[i]), at = seq(min(aggE$SBT.min-.1, na.rm=T),max(aggE$SBT.min+.1, na.rm=T), length.out=40), 
            col.regions=cols, layout = c(5, 1), panel = function(...) {
              panel.fill(col = "light gray")
              panel.levelplot(...)
            }))  
  print(levelplot(SBT.max~longitude*latitude|year_range, data=aggE, xlab=NULL, ylab=NULL, main=paste('SBT.max_rcp', rcp, '  season:',season, '  Model:', modelrun[i]), at = seq(min(aggE$SBT.max-.1, na.rm=T),max(aggE$SBT.max+.1, na.rm=T), length.out=40), 
            col.regions=cols, layout = c(5, 1), panel = function(...) {
              panel.fill(col = "light gray")
              panel.levelplot(...)
            }))  
  print(levelplot(SST.max~longitude*latitude|year_range, data=aggE, xlab=NULL, ylab=NULL, main=paste('SST.max_rcp', rcp, '  season:',season, '  Model:', modelrun[i]), at = seq(min(aggE$SST.max-.1),max(aggE$SST.max+.1), length.out=40), 
            col.regions=cols, layout = c(5, 1), panel = function(...) {
              panel.fill(col = "light gray")
              panel.levelplot(...)
            }))  
  print(levelplot(SBT.seasonal~longitude*latitude|year_range, data=aggW, xlab=NULL, ylab=NULL, main=paste('sbt.seasonal_rcp', rcp, '  season:',season, '  Model:', modelrun[i]), at = seq(min(aggW$SBT.seasonal-.1, na.rm=T),max(aggW$SBT.seasonal+.1, na.rm=T), length.out=40), 
            col.regions=cols, layout = c(5, 1), panel = function(...) {
              panel.fill(col = "light gray")
              panel.levelplot(...)
            }))  
  print(levelplot(SST.seasonal.mean~longitude*latitude|year_range, data=aggW, xlab=NULL, ylab=NULL, main=paste('SST.seasonal.mean_rcp', rcp, '  season:',season, '  Model:', modelrun[i]), at = seq(min(aggW$SST.seasonal.mean-.1),max(aggW$SST.seasonal.mean+.1), length.out=40), 
            col.regions=cols, layout = c(5, 1), panel = function(...) {
              panel.fill(col = "light gray")
              panel.levelplot(...)
            }))  
  print(levelplot(SBT.min~longitude*latitude|year_range, data=aggW, xlab=NULL, ylab=NULL, main=paste('SBT.min_rcp', rcp, '  season:',season, '  Model:', modelrun[i]), at = seq(min(aggW$SBT.min-.1, na.rm=T),max(aggW$SBT.min+.1, na.rm=T), length.out=40), 
            col.regions=cols, layout = c(5, 1), panel = function(...) {
              panel.fill(col = "light gray")
              panel.levelplot(...)
            }))  
  print(levelplot(SBT.max~longitude*latitude|year_range, data=aggW, xlab=NULL, ylab=NULL, main=paste('SBT.max_rcp', rcp, '  season:',season, '  Model:', modelrun[i]), at = seq(min(aggW$SBT.max-.1, na.rm=T),max(aggW$SBT.max+.1, na.rm=T), length.out=40), 
            col.regions=cols, layout = c(5, 1), panel = function(...) {
              panel.fill(col = "light gray")
              panel.levelplot(...)
            }))  
  print(levelplot(SST.max~longitude*latitude|year_range, data=aggW, xlab=NULL, ylab=NULL, main=paste('SST.max_rcp', rcp, '  season:',season, '  Model:', modelrun[i]), at = seq(min(aggW$SST.max-.1),max(aggW$SST.max+.1), length.out=40), 
            col.regions=cols, layout = c(5, 1), panel = function(...) {
              panel.fill(col = "light gray")
              panel.levelplot(...)
            }))  
  dev.off()

  filenameE <- paste('data/prediction_files_Nov2017/predictionEASTnoBias_rcp',rcp, '_', season, '_', modelrun[i], '.RData', sep='')
  filenameW <- paste('data/prediction_files_Nov2017/predictionWESTnoBias_rcp',rcp, '_', season, '_', modelrun[i], '.RData', sep='')
  save(pred.bathE, file=filenameE)
  save(pred.bathW, file=filenameW)
}
   
# To plot a specific area for troubleshooting_can delete this later ==================
#pred2 <- pred[pred$year > 2080,]
#print(levelplot(SBT.seasonal~lonClimgrid*latClimgrid|year, data=pred2, ylim=c(55,61), xlim=c(-172,-142), xlab=NULL, ylab=NULL, main=paste('sbt.seasonal_rcp', rcp, '  season:',season, '  Model:', modelrun[i]), at = seq(0,18, length.out=40), col.regions=cols))
 
#cols = colorRampPalette(colors = c('dark blue', 'blue', 'white', 'red', 'dark red'))
#depths <- proj.grid[proj.grid$longrid > -172 & proj.grid$longrid < -142 & proj.grid$latgrid > 55 & proj.grid$latgrid < 61,]
#levelplot(depth~longrid*latgrid, at = seq(min(depths$depth-.1),max(depths$depth+.1), length.out=40), col.regions=cols, data=depths)
 

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



# =================================================================================
# Below here is Malin's original code_was not used for final manuscript.
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
