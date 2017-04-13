# Read in temperature fields and models, then make range projections
  
## Set working directory
if(Sys.info()["nodename"] == "pinsky-macbookair"){
	setwd('~/Documents/Rutgers/Range projections/proj_ranges/')
	projfolder = '../CEmodels_proj'
	modfolder = '../CEModels'
	climgridfolder <- '../data/'
	numcorestouse <- 2
	}
if(Sys.info()["nodename"] == "amphiprion.deenr.rutgers.edu"){
	setwd('~/Documents/range_projections/')
	projfolder = 'CEmodels_proj'
	modfolder = 'CEmodels'
	climgridfolder <- 'data/'
	numcorestouse <- 12
	# .libPaths(new='~/R/x86_64-redhat-linux-gnu-library/3.1/') # so that it can find my old packages_Muted for Jim's use
}
if(Sys.info()["user"] == "jamesmorley"){
  setwd('/Users/jamesmorley/Documents/project_velocity')
  projfolder = 'output/CEmodels_proj/'
  modfolder <- 'output/CEmodels/'
  climgridfolder <- 'data/'
}

###################
## Load libraries
###################

require(mgcv)
require(Hmisc)
library(lattice)
library(RColorBrewer)
#library(ggplot2)
#library(reshape2)
# require(parallel) # for multi-core calculations

###############################################
# Choose the model fit and other flags to use
###############################################
rcp <- 85
#rcp <- 26

# different seasons to project with (letters represent the first letter of month)
# season <- 'jfm'
# season <- 'amj'
season <- 'jas'
# season <- 'ond'

#runtype <- 'test'
#runtype <- 'testseason'
#runtype <- 'testK6noSeas'
runtype <- 'fitallreg_2017'
stayinregion <- FALSE


#############################
# Choose species to project #
#############################
  
# Distribution models were done in a few parts, due to 4 errors that stopped the loop. So model diagnostics need to be merged together
# The species that got dropped due to problems with model fitting were: "ophidion holbrookii_Atl", "suberites ficus_Pac"; a cusk eel (non fishery) and a sponge, so no big deal about going back and fixing them
load(paste('output/modeldiag_', runtype, '.Rdata', sep='')) # model diagnostics
modeldiag1 <- modeldiag; rm(modeldiag)
load(paste('output/modeldiag_PART2_', runtype, '.Rdata', sep='')) # model diagnostics
modeldiag2 <- modeldiag; rm(modeldiag)
load(paste('output/modeldiag_PART3_', runtype, '.Rdata', sep='')) # model diagnostics
modeldiag3 <- modeldiag; rm(modeldiag)
load(paste('output/modeldiag_PART4_', runtype, '.Rdata', sep='')) # model diagnostics
modeldiag4 <- modeldiag; rm(modeldiag)
load(paste('output/modeldiag_PART5_', runtype, '.Rdata', sep='')) # model diagnostics
modeldiag5 <- modeldiag; rm(modeldiag)
modeldiag1 <- modeldiag1[!modeldiag1$sppocean=='calappa flammea_Atl',] # this was redone with modeldiag2
modeldiag2 <- modeldiag2[1:171,] # above row 171 was also on modeldiag3 (the file didn't close when the loop failed)
modeldiag <- rbind(modeldiag1, modeldiag2, modeldiag3, modeldiag4, modeldiag5)
modeldiag <- modeldiag[!is.na(modeldiag$sppocean),]
length(unique(modeldiag$sppocean)) # 702 total_and all rows are unique species

# 18 species did not fit the testing/training models successfully_I didn't recognize any of these as major fisheries
modeldiag <- modeldiag[!is.na(modeldiag$acc.tt),] # 684 species
#With additional criteria: ***from Elith et al.
modeldiag <- modeldiag[modeldiag$auc.tt >= 0.75,] #669 species
modeldiag <- modeldiag[((modeldiag$dev.pres - modeldiag$dev.pres.null > 0.05) | (modeldiag$dev.biomass - modeldiag$dev.biomass.null > 0.05)),] # down to 658 species

projspp <- modeldiag$sppocean 
rm(modeldiag1, modeldiag2, modeldiag3, modeldiag4, modeldiag5)

#boxplot <- data.frame(cbind(Presence_models = modeldiag$dev.pres, Biomass_models = modeldiag$dev.biomass)) 
#boxplot <- melt(boxplot)
#ggplot(boxplot, aes(x=factor(variable), y=value)) + geom_boxplot() + labs(x='', y='%deviance explained')

print(paste(length(projspp), 'models to project')) # number of species to project to

	# look at species not selected
#	hist(modeldiag$auc.tt)
#	hist(modeldiag$auc)
#	hist(modeldiag$dev.pres - modeldiag$dev.pres.null); abline(v=0.05, col='red')
#	hist(modeldiag$dev.biomass - modeldiag$dev.biomass.null); abline(v=0.05, col='red')
#	notselected <- modeldiag[!(modeldiag$sppocean %in% projspp), c('sppocean', 'auc.tt', 'dev.pres', 'dev.pres.null', 'dev.biomass', 'dev.biomass.null')]
#		sum(notselected$auc.tt < 0.75, na.rm=TRUE)
#		with(notselected, sum((dev.pres - dev.pres.null <= 0.05) & (dev.biomass - dev.biomass.null <= 0.05), na.rm=TRUE))
#		with(notselected, sum(auc.tt < 0.75 & (dev.pres - dev.pres.null <= 0.05) & (dev.biomass - dev.biomass.null <= 0.05), na.rm=TRUE))
#
#	# those that would be selected using auc instead of auc.tt
#	modeldiag$sppocean[modeldiag$auc.tt < 0.75 & modeldiag$auc >= 0.75 & !is.na(modeldiag$auc.tt) & !is.na(modeldiag$auc) & ((modeldiag$dev.pres - modeldiag$dev.pres.null > 0.05) | (modeldiag$dev.biomass - modeldiag$dev.biomass.null > 0.05))]

# find the files with these species for our chosen model fit
files <- list.files(modfolder)
files <- files[grepl(paste('_', runtype, '_', sep=''), files) & grepl(paste(gsub('/|\\(|\\)', '', projspp), collapse='|'), gsub('/|\\(|\\)', '', files))] # have to strip out parentheses and slashes from file and taxon names so that grep doesn't interpret them
length(files) # should match length of projspp

## Remove spp from planned projections IF the projection file already exists (OPTIONAL). 
#If this step is skipped, the existing files will be overwritten.
if(stayinregion) donefiles <- list.files(projfolder, pattern=paste(runtype, '_rcp', rcp, sep='')) # models made earlier
if(!stayinregion) donefiles <- list.files(projfolder, pattern=paste(runtype, '_xreg_rcp', rcp, sep='')) # models made earlier

# trim out prefix and suffix
donespp <- gsub(paste('summproj_', runtype, '_', sep=''), '', gsub('.Rdata', '', donefiles)) 
if(!stayinregion) donespp <- gsub('xreg_', '', donespp)
donespp <- gsub(paste('rcp', rcp, '_', sep=''), '', donespp)

# remove spp that we've already projected
if(length(donespp)>0){
	files <- files[!grepl(paste(gsub('/|\\(|\\)', '', donespp), collapse='|'), gsub('/|\\(|\\)', '', files))] # remove any models that we made earlier
	projspp <- projspp[!grepl(paste(gsub('/|\\(|\\)', '', donespp), collapse='|'), gsub('/|\\(|\\)', '', projspp))]
}
length(files)
length(projspp)

#################################
# Prep environmental data
#################################
   
load('data/projectionGrid_Feb24_2017.RData')# load projection grid to get lat/lon values
clim.grid <- proj.grid # rename as a different 'proj.grid' imported below with bathymetry
clim.grid$depth <- NULL
rm(proj.grid)

if(rcp == 85){
  pred.folder <- c('sst_rcp85/tos_Omon_','sbt_rcp85/temp_btm_1950_2100_')
} else{
  pred.folder <- c('sst_rcp26/tos_Omon_','sbt_rcp26/temp_btm_1950_2100_')
}

modelrun <- c('bcc-csm1-1-m','bcc-csm1-1','CanESM2','CCSM4','CESM1-CAM5','CNRM-CM5','GFDL-CM3','GFDL-ESM2M','GFDL-ESM2G','GISS-E2-R','GISS-E2-H','IPSL-CM5A-LR','IPSL-CM5A-MR','MIROC-ESM','MPI-ESM-LR','NorESM1-ME')
pred.metric <- c('max', 'min', 'mean')

# Malin's functions for reading in the data_need separate function for the min/max files as they have an extra column
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

# For reading in all models for a particular type of data
tempsmeansbt <- matrix(as.numeric(NA), nrow=1281878, ncol=length(modelrun))
tempsmeansst <- matrix(as.numeric(NA), nrow=1281878, ncol=length(modelrun))
tempsminsbt <- matrix(as.numeric(NA), nrow=1281878, ncol=length(modelrun))
tempsmaxsbt <- matrix(as.numeric(NA), nrow=1281878, ncol=length(modelrun))
# tempsminsst <- matrix(as.numeric(NA), nrow=1281878, ncol=length(modelrun)) # Not in habitat models
tempsmaxsst <- matrix(as.numeric(NA), nrow=1281878, ncol=length(modelrun))

# The loops below produce matrix where each column is a different climate model and the years are stacked in one vector
# seasonal mean sbt
for(i in 1:length(modelrun)){
  print(i)
  filename = paste('data/', pred.folder[2], modelrun[i], '_rcp', rcp, '_regrid.nc_fill.nc_2006_2100_', season,'_', pred.metric[3], '.txt', sep="")
  filein <- readLines(filename)
  filein <- t(sapply(filein, processlinesMean))
  rownames(filein) <- NULL
  tempsmeansbt[,i] <- as.vector(filein)
}
# seasonal mean sst
for(i in 1:length(modelrun)){
  print(i)
  filename = paste('data/', pred.folder[1], modelrun[i], '_rcp', rcp, '_r1i1p1_1950_2100.nc_regrid.nc_fill.nc_2006_2100_', season,'_', pred.metric[3], '.txt', sep="")
  filein <- readLines(filename)
  filein <- t(sapply(filein, processlinesMean))
  rownames(filein) <- NULL
  tempsmeansst[,i] <- as.vector(filein)
}
# min annual sbt
for(i in 1:length(modelrun)){
  print(i)
  filename = paste('data/', pred.folder[2], modelrun[i], '_rcp', rcp, '_regrid.nc_fill.nc_2006_2100_', season,'_', pred.metric[2], '.txt', sep="")
  filein <- readLines(filename)
  filein <- t(sapply(filein, processlinesMinMax))
  rownames(filein) <- NULL
  fileinVec <- as.vector(filein)
  fileinVec <- fileinVec[1:1281878]
  tempsminsbt[,i] <- fileinVec
}
# min annual sst
for(i in 1:length(modelrun)){
  print(i)
  filename = paste('data/', pred.folder[1], modelrun[i], '_rcp', rcp, '_r1i1p1_1950_2100.nc_regrid.nc_fill.nc_2006_2100_', season,'_', pred.metric[2], '.txt', sep="")
  filein <- readLines(filename)
  filein <- t(sapply(filein, processlinesMinMax))
  rownames(filein) <- NULL
  fileinVec <- as.vector(filein)
  fileinVec <- fileinVec[1:1281878]
  tempsminsst[,i] <- fileinVec
}
# max annual sbt
for(i in 1:length(modelrun)){
  print(i)
  filename = paste('data/', pred.folder[2], modelrun[i], '_rcp', rcp, '_regrid.nc_fill.nc_2006_2100_', season,'_', pred.metric[1], '.txt', sep="")
  filein <- readLines(filename)
  filein <- t(sapply(filein, processlinesMinMax))
  rownames(filein) <- NULL
  fileinVec <- as.vector(filein)
  fileinVec <- fileinVec[1:1281878]
  tempsmaxsbt[,i] <- fileinVec
}
# max annual sst
for(i in 1:length(modelrun)){
  print(i)
  filename = paste('data/', pred.folder[1], modelrun[i], '_rcp', rcp, '_r1i1p1_1950_2100.nc_regrid.nc_fill.nc_2006_2100_', season,'_', pred.metric[1], '.txt', sep="")
  filein <- readLines(filename)
  filein <- t(sapply(filein, processlinesMinMax))
  rownames(filein) <- NULL
  fileinVec <- as.vector(filein)
  fileinVec <- fileinVec[1:1281878]
  tempsmaxsst[,i] <- fileinVec
}
   
colnames(clim.grid) <- c('lonClimgrid', 'latClimgrid') 
clim.grid$index <- c(1:nrow(clim.grid))
load('data/ProjectionBathGrid_Feb27_2017.RData')# load bathymetry projection grid
proj.grid$depth <- NULL # To reduce file size
year.bin <- data.frame(year=2007:2100, bin=c(rep('2007-2020', 14), rep('2021-2040', 20), rep('2041-2060', 20), rep('2061-2080', 20), rep('2081-2100', 20))) # for making summary graphs in the loop
cols = colorRampPalette(colors = c('dark blue', 'blue', 'white', 'red', 'dark red'))

# Build and save all 16 prediction data sets_make some summary figures
# The following arrays are muted b/c my CPU couldn't handle the large files, so I save a prediction file for each projection model, split by east vs. west
# ....the array may work on amphiprion, but it is a huge file so I'm thinking this method will be faster
#pred_arrayEast <- array(numeric(), dim=c(7763648, 15, 16))
#pred_arrayWest <- array(numeric(), dim=c(6187644, 15, 16))
for(i in 1:length(modelrun)){ 
  print(paste('Beginning model number ', i, ': ', modelrun[i], sep=''))
  pred <- data.frame(cbind(clim.grid, year=rep(2007:2100, rep=94, each=13637), SBT.seasonal = tempsmeansbt[,i], SST.seasonal.mean = tempsmeansst[,i], SBT.min = tempsminsbt[,i], SBT.max = tempsmaxsbt[,i], SST.max = tempsmaxsst[,i]))
  pred <- merge(pred, year.bin, by='year', all.x=T)  # Bin into 20 year periods to make some plots
  pred.bathE <- merge(pred, proj.grid[proj.grid$lonClimgrid > -105,], by=c('lonClimgrid', 'latClimgrid'), all.y = T, sort=F)
  pred.bathW <- merge(pred, proj.grid[proj.grid$lonClimgrid < -105,], by=c('lonClimgrid', 'latClimgrid'), all.y = T, sort=F)
  pred.bathE <- pred.bathE[order(pred.bathE$year, pred.bathE$index),] # reorder by year and index
  pred.bathW <- pred.bathW[order(pred.bathW$year, pred.bathW$index),] # reorder by year and index
  
  aggE <- aggregate(list(sbt.seasonal = pred.bathE$SBT.seasonal, SST.seasonal.mean = pred.bathE$SST.seasonal.mean, SBT.min = pred.bathE$SBT.min,
            SBT.max = pred.bathE$SBT.max, SST.max = pred.bathE$SST.max), by=list(year_range=pred.bathE$bin, latitude=pred.bathE$latBathgrid, longitude=pred.bathE$lonBathgrid), FUN=mean)
  aggW <- aggregate(list(sbt.seasonal = pred.bathW$SBT.seasonal, SST.seasonal.mean = pred.bathW$SST.seasonal.mean, SBT.min = pred.bathW$SBT.min,
            SBT.max = pred.bathW$SBT.max, SST.max = pred.bathW$SST.max), by=list(year_range=pred.bathW$bin, latitude=pred.bathW$latBathgrid, longitude=pred.bathW$lonBathgrid), FUN=mean)
  
  # 20yr. summary figures for each model
  pdf(width=14, height=5, file=paste('figures/Temp_projections/', modelrun[i], '_20yr_tempProj.pdf', sep=''))
  print(levelplot(sbt.seasonal~longitude*latitude|year_range, data=aggE, xlab=NULL, ylab=NULL, main=paste('sbt.seasonal_rcp', rcp, '  season:',season, '  Model:', modelrun[i]), at = seq(min(aggE$sbt.seasonal-.1),max(aggE$sbt.seasonal+.1), length.out=40), 
            col.regions=cols, layout = c(5, 1), panel = function(...) {
              panel.fill(col = "light gray")
              panel.levelplot(...)
            }))  
  print(levelplot(SST.seasonal.mean~longitude*latitude|year_range, data=aggE, xlab=NULL, ylab=NULL, main=paste('SST.seasonal.mean_rcp', rcp, '  season:',season, '  Model:', modelrun[i]), at = seq(min(aggE$SST.seasonal.mean-.1),max(aggE$SST.seasonal.mean+.1), length.out=40), 
            col.regions=cols, layout = c(5, 1), panel = function(...) {
              panel.fill(col = "light gray")
              panel.levelplot(...)
            }))  
  print(levelplot(SBT.min~longitude*latitude|year_range, data=aggE, xlab=NULL, ylab=NULL, main=paste('SBT.min_rcp', rcp, '  season:',season, '  Model:', modelrun[i]), at = seq(min(aggE$SBT.min-.1),max(aggE$SBT.min+.1), length.out=40), 
            col.regions=cols, layout = c(5, 1), panel = function(...) {
              panel.fill(col = "light gray")
              panel.levelplot(...)
            }))  
  print(levelplot(SBT.max~longitude*latitude|year_range, data=aggE, xlab=NULL, ylab=NULL, main=paste('SBT.max_rcp', rcp, '  season:',season, '  Model:', modelrun[i]), at = seq(min(aggE$SBT.max-.1),max(aggE$SBT.max+.1), length.out=40), 
            col.regions=cols, layout = c(5, 1), panel = function(...) {
              panel.fill(col = "light gray")
              panel.levelplot(...)
            }))  
  print(levelplot(SST.max~longitude*latitude|year_range, data=aggE, xlab=NULL, ylab=NULL, main=paste('SST.max_rcp', rcp, '  season:',season, '  Model:', modelrun[i]), at = seq(min(aggE$SST.max-.1),max(aggE$SST.max+.1), length.out=40), 
            col.regions=cols, layout = c(5, 1), panel = function(...) {
              panel.fill(col = "light gray")
              panel.levelplot(...)
            }))  
  print(levelplot(sbt.seasonal~longitude*latitude|year_range, data=aggW, xlab=NULL, ylab=NULL, main=paste('sbt.seasonal_rcp', rcp, '  season:',season, '  Model:', modelrun[i]), at = seq(min(aggW$sbt.seasonal-.1),max(aggW$sbt.seasonal+.1), length.out=40), 
            col.regions=cols, layout = c(5, 1), panel = function(...) {
              panel.fill(col = "light gray")
              panel.levelplot(...)
            }))  
  print(levelplot(SST.seasonal.mean~longitude*latitude|year_range, data=aggW, xlab=NULL, ylab=NULL, main=paste('SST.seasonal.mean_rcp', rcp, '  season:',season, '  Model:', modelrun[i]), at = seq(min(aggW$SST.seasonal.mean-.1),max(aggW$SST.seasonal.mean+.1), length.out=40), 
            col.regions=cols, layout = c(5, 1), panel = function(...) {
              panel.fill(col = "light gray")
              panel.levelplot(...)
            }))  
  print(levelplot(SBT.min~longitude*latitude|year_range, data=aggW, xlab=NULL, ylab=NULL, main=paste('SBT.min_rcp', rcp, '  season:',season, '  Model:', modelrun[i]), at = seq(min(aggW$SBT.min-.1),max(aggW$SBT.min+.1), length.out=40), 
            col.regions=cols, layout = c(5, 1), panel = function(...) {
              panel.fill(col = "light gray")
              panel.levelplot(...)
            }))  
  print(levelplot(SBT.max~longitude*latitude|year_range, data=aggW, xlab=NULL, ylab=NULL, main=paste('SBT.max_rcp', rcp, '  season:',season, '  Model:', modelrun[i]), at = seq(min(aggW$SBT.max-.1),max(aggW$SBT.max+.1), length.out=40), 
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

  filenameE <- paste('data/prediction_files/predictionEAST_rcp',rcp, '_', season, '_', modelrun[i], '.RData', sep='')
  filenameW <- paste('data/prediction_files/predictionWEST_rcp',rcp, '_', season, '_', modelrun[i], '.RData', sep='')
  save(pred.bathE, file=filenameE)
  save(pred.bathW, file=filenameW)
  #pred_arrayEast[,,i] <- pred.bathE
  #pred_arrayWest[,,i] <- pred.bathW
}
  
# =======================================================================================
# CODE BELOW TO LOOK AT PROJECTION DATA A DIFFERENT WAY THAN THE ABOVE 20yr. SUMMARIES
# =======================================================================================
# Below are four blocks that can be recycled to bring in any of the data_adjust 'i' for model runs_this could be condensed if necessary
proj.clim <- unique(data.frame(lonClimgrid=proj.grid$lonClimgrid, latClimgrid=proj.grid$latClimgrid))
# seasonal mean sbt
  filename = paste('data/', pred.folder[2], modelrun[i], '_rcp', rcp, '_regrid.nc_fill.nc_2006_2100_', season,'_', pred.metric[3], '.txt', sep="")
  filein <- readLines(filename)
  filein <- t(sapply(filein, processlinesMean))
  rownames(filein) <- NULL
   
# seasonal mean sst
  filename = paste('data/', pred.folder[1], modelrun[i], '_rcp', rcp, '_r1i1p1_1950_2100.nc_regrid.nc_fill.nc_2006_2100_', season,'_', pred.metric[3], '.txt', sep="")
  filein <- readLines(filename)
  filein <- t(sapply(filein, processlinesMean))
  rownames(filein) <- NULL

# min or max sbt_for these two need to adjust 'i' but also 'pred.metric' can be 1 or 2 (max or min)
  filename = paste('data/', pred.folder[2], modelrun[i], '_rcp', rcp, '_regrid.nc_fill.nc_2006_2100_', season,'_', pred.metric[2], '.txt', sep="")
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

plot(latClimgrid~lonClimgrid, cex=.1, data=abc[!is.na(abc$X7),])
points(latClimgrid~lonClimgrid, col='blue', cex=.5, data=abc[abc$X7 < -5, ])
points(latClimgrid~lonClimgrid, col='red', data=abc[abc$X89 > 45, ])
levelplot(X89 ~ lonClimgrid * latClimgrid, data = abc)
# another plot option for showing just the outier areas well
forPlot <- abc$X89# choose a given year to work with the plots_this is only necessary for establishing the cutpoints
cutpts <- c(min(forPlot, na.rm=T)-1, -10, -5, 40, 80, max(forPlot, na.rm=T)+1) # add some cutpoints to show outliers_gray is relatively normal, with two levels of both negative and positive outliers
levelplot(X89 ~ lonClimgrid * latClimgrid, data = abc[!is.na(abc$X89),], at = cutpts, cuts = 6, pretty = T, col.regions = (rev(brewer.pal(9, "RdBu")))) 

# ===========================================
# Malin's code for looking at empty cells_this can probably be removed as all empty cells appear to have been filled
# ===========================================

# how much data do we have? (need array from previous code block)
modeldat <- apply(temps, MARGIN=c(1,3), FUN = function(x) all(!is.na(x))) # whether or not a model has data at each grid cell
colnames(modeldat) <- modelrun
t(t(colSums(modeldat))) # how many grid cells covered by each model

# how many models have data at each grid cell?
nmodeldat <- rowSums(modeldat)

hist(nmodeldat, breaks=seq(-0.5, 16.5, by=1))
sum(nmodeldat==16)
sum(nmodeldat>=13) # 9479
sum(nmodeldat>=10) # 11736

# how many models have data at each grid cell? W/OUT GFDL or Can
keep <- !grepl('GFDL', modelrun)
keep <- !grepl('GFDL|CanESM', modelrun)
sum(keep)
nmodeldat <- rowSums(modeldat[, keep])

hist(nmodeldat, breaks=seq(-0.5, 16.5, by=1))
sum(nmodeldat>=13)
sum(nmodeldat>=12) 
sum(nmodeldat>=10)

dim(temps)

# =======================================================================================
# The code below is still in progress
# =======================================================================================

load(filename)
load('output/CEmodels/CEmods_fitallreg_2017_archosargus probatocephalus_Atl.RData')
mygam1 <- mods[[1]]
mygam2 <- mods[[2]]
pred.bath$regionfact <- 'SCDNR_SEUS'
colnames(pred.bath) <- c('lonClimgrid', 'latClimgrid', 'index', 'year', 'SBT.seasonal', 'SST.seasonal.mean', 'SBT.min', 'SST.min', 'SBT.max', 'SST.max', 'latBathgrid', 'lonBathgrid', 'rugosity', 'GRAINSIZE', 'regionfact')
pred.bathE <- pred.bath[pred.bath$lonClimgrid > -105,]
pred.bathE <- pred.bathE[!is.na(pred.bathE$SBT.seasonal),]
pred.bathE <- pred.bathE[!is.na(pred.bathE$rugosity),]
#Change the na.action to just drop those rows? 
preds1 <- predict(mygam1, newdata = pred.bathE[pred.bathE$year == 2007,], type='response')
preds2 <- exp(predict(mygam2, newdata = pred.bathE[pred.bathE$year == 2007,], type='response')) 
preds <- preds1*preds2
pred.grid <- data.frame(cbind(pred.bathE[pred.bathE$year == 2007,], preds))

levelplot(SBT.seasonal ~ lonBathgrid*latBathgrid, data=pred.grid)

# =================================================================================
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
