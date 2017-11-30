setwd('/Users/jamesmorley/Documents/project_velocity')
setwd("/Users/abigail/Desktop/temp")
setwd("/Users/abigailporay/Desktop/temp")

require(Hmisc)
library(latticeExtra)
library(maps)
library(geosphere)
library(plotrix)
library(quantreg)
library(SDMTools)

modelrun <- c('bcc-csm1-1-m','bcc-csm1-1','CanESM2','CCSM4','CESM1-CAM5','CNRM-CM5','GFDL-CM3','GFDL-ESM2M','GFDL-ESM2G','GISS-E2-R','GISS-E2-H','IPSL-CM5A-LR','IPSL-CM5A-MR','MIROC-ESM','MPI-ESM-LR','NorESM1-ME')
rcp <- c(26,85)
projfolder = 'output/CEmodels_proj_May2017/'
#projfolder = 'temp/'
load('data/speciesProjectionList.RData')
load('speciesProjectionList.RData')
projspp <- projspp[-168] # For some reason this species didn't run rcp85
regionFreq <- regionFreq[-168] # For some reason this species didn't run rcp85
save(projspp, regionFreq, file='data/speciesProjectionList_N658.RData')
# ====================================================================================
# Make some figures of the 'raw' data for a random subset of species
# ====================================================================================
testSpp <- sort(sample(projspp, size=50, replace=F)) # random 50 species to look at
testSpp <- projspp[316]
for(i in 1:length(testSpp)){
  species = testSpp[i]
  print(paste("Beginning figures for ", species, sep=''))
  #filename <- paste(projfolder, species, '_rcp', rcp[1], '_jas_prediction_AGG.RData', sep='')
  filename <- paste('output/CEmodels_proj_May2017/', species, '_rcp', rcp[1], '_jas_prediction_AGG.RData', sep='')
  load(filename)
  pred.agg26 <- pred.agg
  
  #filename <- paste('output/CEmodels_proj_May2017/', species, '_rcp', rcp[2], '_jas_prediction_AGG.RData', sep='')
  filename <- paste('output/CEmodels_proj_May2017/', species, '_rcp', rcp[2], '_jas_prediction_AGG.RData', sep='')
  
  load(filename)
  pred.agg85 <- pred.agg
  rm(pred.agg)
  # Look at the model-level predictions
  #pdf(width=14, height=5, file=paste('figures/speciesProjections/MS/', species, '_projections.pdf', sep=''))
  pdf(width=14, height=5, file=paste('figures/', species, '_projections.pdf', sep=''))
  
  #tiff(width=14, height=5, units='in', res=500, file=paste('figures/speciesProjections/MS/', species, '_projections.tiff', sep=''))
  for(i in 4:19){
    print(i-3)
    
    scale26log = seq(0, max(log(pred.agg26[,i] + 1)) + .00001, length.out=20)
    scale85log = seq(0, max(log(pred.agg85[,i] + 1)) + .00001, length.out=20)
    scale26 = seq(0, max(pred.agg26[,i]) + .00001, length.out=20)
    scale85 = seq(0, max(pred.agg85[,i]) + .00001, length.out=20)
    cols = colorRampPalette(colors = c('gray90', 'blue', 'dark blue', 'black'))
    
    # STANDARDIZE PLOT DIMENSIONS FOR ALL FOUR PLOTS
    xlimit1 = c(min(pred.agg26$longitude[pred.agg26[,i] > scale26[2]/2] - 1),max(pred.agg26$longitude[pred.agg26[,i] > scale26[2]/2] + 1))
    ylimit1 = c(min(pred.agg26$latitude[pred.agg26[,i] > scale26[2]/2] - .25), max(pred.agg26$latitude[pred.agg26[,i] > scale26[2]/2] + .25))
    xlimit2 = c(min(pred.agg26$longitude[log(pred.agg26[,i] + 1) > scale26log[2]/2] - 1),max(pred.agg26$longitude[log(pred.agg26[,i] + 1) > scale26log[2]/2] + 1))
    ylimit2 = c(min(pred.agg26$latitude[log(pred.agg26[,i] + 1) > scale26log[2]/2] - .25), max(pred.agg26$latitude[log(pred.agg26[,i] + 1) > scale26log[2]/2] + .25))
    xlimit3 = c(min(pred.agg85$longitude[pred.agg85[,i] > scale85[2]/2] - 1),max(pred.agg85$longitude[pred.agg85[,i] > scale85[2]/2] + 1))
    ylimit3 = c(min(pred.agg85$latitude[pred.agg85[,i] > scale85[2]/2] - .25), max(pred.agg85$latitude[pred.agg85[,i] > scale85[2]/2] + .25))
    xlimit4 = c(min(pred.agg85$longitude[log(pred.agg85[,i] + 1) > scale85log[2]/2] - 1),max(pred.agg85$longitude[log(pred.agg85[,i] + 1) > scale85log[2]/2] + 1))
    ylimit4 = c(min(pred.agg85$latitude[log(pred.agg85[,i] + 1) > scale85log[2]/2] - .25), max(pred.agg85$latitude[log(pred.agg85[,i] + 1) > scale85log[2]/2] + .25))
    xlimit = c(min(xlimit1[1],xlimit2[1],xlimit3[1],xlimit4[1]),max(xlimit1[2],xlimit2[2],xlimit3[2],xlimit4[2]))
    ylimit = c(min(ylimit1[1],ylimit2[1],ylimit3[1],ylimit4[1]),max(ylimit1[2],ylimit2[2],ylimit3[2],ylimit4[2]))
    
    # RCP26 plots_one for normal, one for Logged+1
    Projmap <- map('world', xlim=xlimit, ylim=ylimit, plot=F) 
    Projmap <- data.frame(lon=Projmap$x, lat=Projmap$y)
    print(levelplot(pred.agg26[,i] ~ longitude*latitude|year_range, data=pred.agg26, xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, main=paste('rcp', rcp[1], 'summer', 'Model', modelrun[i-3], species, sep='_'), 
                    at = scale26, col.regions=cols, layout = c(5, 1)) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')) 
    
    Projmap <- map('world', xlim=xlimit, ylim=ylimit, plot=F) 
    Projmap <- data.frame(lon=Projmap$x, lat=Projmap$y)
    print(levelplot(log(pred.agg26[,i] + 1) ~ longitude*latitude|year_range, data=pred.agg26, xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, main=paste('rcp', rcp[1], 'summer', 'Model', modelrun[i-3], species, 'LOGGED', sep='_'), 
                    at = scale26log, col.regions=cols, layout = c(5, 1)) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')) 
     
    # RCP85 plots_one for normal, one for Logged+1
    Projmap <- map('world', xlim=xlimit, ylim=ylimit, plot=F) 
    Projmap <- data.frame(lon=Projmap$x, lat=Projmap$y)
    print(levelplot(pred.agg85[,i] ~ longitude*latitude|year_range, data=pred.agg85, xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, main=paste('rcp', rcp[2], 'summer', 'Model', modelrun[i-3], species, sep='_'), 
                    at = scale85, col.regions=cols, layout = c(5, 1)) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')) 
    
    Projmap <- map('world', xlim=xlimit, ylim=ylimit, plot=F) 
    Projmap <- data.frame(lon=Projmap$x, lat=Projmap$y)
    print(levelplot(log(pred.agg85[,i] + 1) ~ longitude*latitude|year_range, data=pred.agg85, xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, main=paste('rcp', rcp[2], 'summer', 'Model', modelrun[i-3], species, 'LOGGED', sep='_'), 
                    at = scale85log, col.regions=cols, layout = c(5, 1)) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')) 
  }
  dev.off()
}
  
# ============================================================================================================
# Compare centroid shifts with RCPs and logged vs. non-logged_b/c the visualizations looks better when logged
# Ideally, we don't want to log the centroid calculations
# ============================================================================================================
testSpp <- 'sebastes babcocki_Pac'
for(i in 1:length(testSpp)){
  species = testSpp[i]
  print(paste("Beginning figures for ", species, sep=''))
  filename <- paste(species, '_rcp', rcp[1], '_jas_prediction_AGG.RData', sep='')
  
  #filename <- paste(projfolder, species, '_rcp', rcp[1], '_jas_prediction_AGG.RData', sep='')
  load(filename)
  pred.agg26 <- pred.agg
  pred.agg26pre <- pred.agg26[pred.agg26$year_range == '2007-2020',]
  pred.agg26post <- pred.agg26[pred.agg26$year_range == '2081-2100',]
  
  #filename <- paste('output/CEmodels_proj_May2017/', species, '_rcp', rcp[2], '_jas_prediction_AGG.RData', sep='')
  filename <- paste(species, '_rcp', rcp[2], '_jas_prediction_AGG.RData', sep='')
  
  load(filename)
  pred.agg85 <- pred.agg
  pred.agg85pre <- pred.agg85[pred.agg85$year_range == '2007-2020',]
  pred.agg85post <- pred.agg85[pred.agg85$year_range == '2081-2100',]
  rm(pred.agg)
  
  #pdf(width=7, height=5, file=paste('figures/speciesProjections/MS/centroids/', species, 'centroid_projections.pdf', sep=''))
  pdf(width=7, height=5, file=paste('figures/', species, 'centroid_projections.pdf', sep=''))
  
  for(j in 4:19){
    print(j-3)
    cent1Lat26 = wtd.mean(pred.agg26pre$latitude, weights = pred.agg26pre[,j]) 
    cent2Lat26 = wtd.mean(pred.agg26post$latitude, weights = pred.agg26post[,j])
    cent1logLat26 = wtd.mean(pred.agg26pre$latitude, weights = log(pred.agg26pre[,j] + 1)) 
    cent2logLat26 = wtd.mean(pred.agg26post$latitude, weights = log(pred.agg26post[,j] + 1)) 
    
    cent1Lon26 = wtd.mean(pred.agg26pre$longitude, weights = pred.agg26pre[,j]) 
    cent2Lon26 = wtd.mean(pred.agg26post$longitude, weights = pred.agg26post[,j])
    cent1logLon26 = wtd.mean(pred.agg26pre$longitude, weights = log(pred.agg26pre[,j] + 1)) 
    cent2logLon26 = wtd.mean(pred.agg26post$longitude, weights = log(pred.agg26post[,j] + 1)) 
    
    cent1Lat85 = wtd.mean(pred.agg85pre$latitude, weights = pred.agg85pre[,j]) 
    cent2Lat85 = wtd.mean(pred.agg85post$latitude, weights = pred.agg85post[,j])
    cent1logLat85 = wtd.mean(pred.agg85pre$latitude, weights = log(pred.agg85pre[,j] + 1)) 
    cent2logLat85 = wtd.mean(pred.agg85post$latitude, weights = log(pred.agg85post[,j] + 1)) 
    
    cent1Lon85 = wtd.mean(pred.agg85pre$longitude, weights = pred.agg85pre[,j]) 
    cent2Lon85 = wtd.mean(pred.agg85post$longitude, weights = pred.agg85post[,j])
    cent1logLon85 = wtd.mean(pred.agg85pre$longitude, weights = log(pred.agg85pre[,j] + 1)) 
    cent2logLon85 = wtd.mean(pred.agg85post$longitude, weights = log(pred.agg85post[,j] + 1)) 
    
    # STANDARDIZE PLOT DIMENSIONS FOR ALL FOUR PLOTS
    scale26log = seq(0, max(log(pred.agg26[,j] + 1)) + .00001, length.out=20)
    scale85log = seq(0, max(log(pred.agg85[,j] + 1)) + .00001, length.out=20)
    scale26 = seq(0, max(pred.agg26[,j]) + .00001, length.out=20)
    scale85 = seq(0, max(pred.agg85[,j]) + .00001, length.out=20)
    
    xlimit1 = c(min(pred.agg26$longitude[pred.agg26[,j] > scale26[2]/2] - 1),max(pred.agg26$longitude[pred.agg26[,j] > scale26[2]/2] + 1))
    ylimit1 = c(min(pred.agg26$latitude[pred.agg26[,j] > scale26[2]/2] - .25), max(pred.agg26$latitude[pred.agg26[,j] > scale26[2]/2] + .25))
    xlimit2 = c(min(pred.agg26$longitude[log(pred.agg26[,j] + 1) > scale26log[2]/2] - 1),max(pred.agg26$longitude[log(pred.agg26[,j] + 1) > scale26log[2]/2] + 1))
    ylimit2 = c(min(pred.agg26$latitude[log(pred.agg26[,j] + 1) > scale26log[2]/2] - .25), max(pred.agg26$latitude[log(pred.agg26[,j] + 1) > scale26log[2]/2] + .25))
    xlimit3 = c(min(pred.agg85$longitude[pred.agg85[,j] > scale85[2]/2] - 1),max(pred.agg85$longitude[pred.agg85[,j] > scale85[2]/2] + 1))
    ylimit3 = c(min(pred.agg85$latitude[pred.agg85[,j] > scale85[2]/2] - .25), max(pred.agg85$latitude[pred.agg85[,j] > scale85[2]/2] + .25))
    xlimit4 = c(min(pred.agg85$longitude[log(pred.agg85[,j] + 1) > scale85log[2]/2] - 1),max(pred.agg85$longitude[log(pred.agg85[,j] + 1) > scale85log[2]/2] + 1))
    ylimit4 = c(min(pred.agg85$latitude[log(pred.agg85[,j] + 1) > scale85log[2]/2] - .25), max(pred.agg85$latitude[log(pred.agg85[,j] + 1) > scale85log[2]/2] + .25))
    xlimit = c(min(xlimit1[1],xlimit2[1],xlimit3[1],xlimit4[1]),max(xlimit1[2],xlimit2[2],xlimit3[2],xlimit4[2]))
    ylimit = c(min(ylimit1[1],ylimit2[1],ylimit3[1],ylimit4[1]),max(ylimit1[2],ylimit2[2],ylimit3[2],ylimit4[2]))
    
    Projmap <- map('world', xlim=xlimit, ylim=ylimit, plot=T) 
    points(x=cent1Lon26, y=cent1Lat26, col='light blue')
    arrows(x0=cent1Lon26, y0=cent1Lat26, x1=cent2Lon26, y1=cent2Lat26, length=.1, lwd=1, col='light blue')
    points(x=cent1logLon26, y=cent1logLat26, col='dark blue')
    arrows(x0=cent1logLon26, y0=cent1logLat26, x1=cent2logLon26, y1=cent2logLat26, length=.1, lwd=1, col='dark blue')
    points(x=cent1Lon85, y=cent1Lat85, col='pink')
    arrows(x0=cent1Lon85, y0=cent1Lat85, x1=cent2Lon85, y1=cent2Lat85, length=.1, lwd=1, col='pink')
    points(x=cent1logLon85, y=cent1logLat85, col='dark red')
    arrows(x0=cent1logLon85, y0=cent1logLat85, x1=cent2logLon85, y1=cent2logLat85, length=.1, lwd=1, col='dark red')
  }  
  dev.off()  
}     

# ============================================================================================================
# ============================================================================================================
# Loop that calculates the mean lat/lon of centroid start/end for each rcp
# Also will calculate the mean/SD of latitudinal CHANGE (- or +) and mean/sd of distance shifted (km)
# ============================================================================================================
# ============================================================================================================

# columns for init and final lat and lon, mean + SD in latitudinal change, mean + SD in absolute shift distance
# Ultimately it would be better to do some 'along coastline' shift distance, but nothing would work for every region or species
centroids26 <- data.frame(latPre=numeric(), lonPre=numeric(), latPost=numeric(), lonPost=numeric(), meanLat=numeric(), sdLat=numeric(), meanDist=numeric(), sdDist=numeric(), radius=numeric())
centroids85 <- data.frame(latPre=numeric(), lonPre=numeric(), latPost=numeric(), lonPost=numeric(), meanLat=numeric(), sdLat=numeric(), meanDist=numeric(), sdDist=numeric(), radius=numeric())
  
for(i in 1:length(projspp)){ # one loop for each species
  species = projspp[i]
  print(paste("Beginning figures for ", species, '_spp', i, sep=''))
  filename <- paste(projfolder, species, '_rcp', rcp[1], '_jas_prediction_AGG.RData', sep='')
  load(filename)
  pred.agg26 <- pred.agg
   
  filename <- paste(projfolder, species, '_rcp', rcp[2], '_jas_prediction_AGG.RData', sep='')
  load(filename)
  pred.agg85 <- pred.agg
  rm(pred.agg)
  # RCP26 calculations
  yearRange <- c('2007-2020', '2081-2100')
  cent_preds26lat <- matrix(data=NA, nrow=2, ncol=16)
  cent_preds26lon <- matrix(data=NA, nrow=2, ncol=16)
  for(j in 1:2){
    years = yearRange[j]
    preds <- pred.agg26[pred.agg26$year_range == years,]
    for(k in 4:19){
      valueLat = wtd.mean(preds$latitude, weights = preds[,k]) 
      valueLon = wtd.mean(preds$longitude, weights = preds[,k]) 
      cent_preds26lat[j,k-3] = valueLat
      cent_preds26lon[j,k-3] = valueLon
    }
  } 
  mean_Lat_pre26 <- apply(cent_preds26lat, 1, FUN=mean, na.rm=T)[1]
  mean_Lat_post26 <- apply(cent_preds26lat, 1, FUN=mean, na.rm=T)[2]
  mean_Lon_pre26 <- apply(cent_preds26lon, 1, FUN=mean, na.rm=T)[1]
  mean_Lon_post26 <- apply(cent_preds26lon, 1, FUN=mean, na.rm=T)[2]
  # Need to calculate change in latitude for each model_then get the SD of those values
  lat_diff26 <- apply(cent_preds26lat, 2, FUN=diff, na.rm=T)
  meanLat26 <- mean(lat_diff26)
  sdLat26 <- sd(lat_diff26)
  # Need to calculate distance and angle of shifts for each model_then means and SDs
  dist26 <- distHaversine(p1=matrix(cbind(cent_preds26lon[1,], cent_preds26lat[1,]), ncol=2), p2=matrix(cbind(cent_preds26lon[2,], cent_preds26lat[2,]), ncol=2)) / 1000# distance of vectors in km (I think in km anyway)
  meanDis26 <- mean(dist26)
  sdDis26 <- sd(dist26)
  bearing26 <- bearing(p1=matrix(cbind(cent_preds26lon[1,], cent_preds26lat[1,]), ncol=2), p2=matrix(cbind(cent_preds26lon[2,], cent_preds26lat[2,]), ncol=2))
  bearing26scaled <- ifelse(bearing26 < 0, bearing26 + 360, bearing26)
  bearing26rads <- bearing26scaled * pi/180
  xy26rads <- matrix(cbind(1*cos(bearing26rads), 1*sin(bearing26rads)), ncol=2, nrow=16)
  plot(xy26rads[,2]~xy26rads[,1], xlim=c(-2,2), ylim=c(-2,2))
  xyMean <- c(mean(xy26rads[,1]),mean(xy26rads[,2]))
  radius26 <- sqrt(xyMean[1]^2 + xyMean[2]^2)
  
  # RCP85 calculations
  cent_preds85lat <- matrix(data=NA, nrow=2, ncol=16)
  cent_preds85lon <- matrix(data=NA, nrow=2, ncol=16)
  for(j in 1:2){
    years = yearRange[j]
    preds <- pred.agg85[pred.agg85$year_range == years,]
    for(k in 4:19){
      valueLat = wtd.mean(preds$latitude, weights = preds[,k]) 
      valueLon = wtd.mean(preds$longitude, weights = preds[,k]) 
      cent_preds85lat[j,k-3] = valueLat
      cent_preds85lon[j,k-3] = valueLon
    }
  } 
  mean_Lat_pre85 <- apply(cent_preds85lat, 1, FUN=mean, na.rm=T)[1]
  mean_Lat_post85 <- apply(cent_preds85lat, 1, FUN=mean, na.rm=T)[2]
  mean_Lon_pre85 <- apply(cent_preds85lon, 1, FUN=mean, na.rm=T)[1]
  mean_Lon_post85 <- apply(cent_preds85lon, 1, FUN=mean, na.rm=T)[2]
  # Need to calculate change in latitude for each model_then get the SD of those values
  lat_diff85 <- apply(cent_preds85lat, 2, FUN=diff, na.rm=T)
  meanLat85 <- mean(lat_diff85)
  sdLat85 <- sd(lat_diff85)
  # Need to calculate distance shifted for each model
  dist85 <- distHaversine(p1=matrix(cbind(cent_preds85lon[1,], cent_preds85lat[1,]), ncol=2), p2=matrix(cbind(cent_preds85lon[2,], cent_preds85lat[2,]), ncol=2)) / 1000# distance of vectors in km (I think in km anyway)
  meanDis85 <- mean(dist85)
  sdDis85 <- sd(dist85)
  bearing85 <- bearing(p1=matrix(cbind(cent_preds85lon[1,], cent_preds85lat[1,]), ncol=2), p2=matrix(cbind(cent_preds85lon[2,], cent_preds85lat[2,]), ncol=2))
  bearing85scaled <- ifelse(bearing85 < 0, bearing85 + 360, bearing85)
  bearing85rads <- bearing85scaled * pi/180
  xy85rads <- matrix(cbind(1*cos(bearing85rads), 1*sin(bearing85rads)), ncol=2, nrow=16)
  plot(xy85rads[,2]~xy85rads[,1], xlim=c(-2,2), ylim=c(-2,2))
  points(y2~y4, col="pink")
  xyMean <- c(mean(xy85rads[,1]),mean(xy85rads[,2]))
  radius85 <- sqrt(xyMean[1]^2 + xyMean[2]^2)
    
  centroids26[i,] <- data.frame(latPre=mean_Lat_pre26, lonPre=mean_Lon_pre26, latPost=mean_Lat_post26, lonPost=mean_Lon_post26, meanLat=meanLat26, sdLat=sdLat26, meanDist=meanDis26, sdDist=sdDis26, radius=radius26)
  centroids85[i,] <- data.frame(latPre=mean_Lat_pre85, lonPre=mean_Lon_pre85, latPost=mean_Lat_post85, lonPost=mean_Lon_post85, meanLat=meanLat85, sdLat=sdLat85, meanDist=meanDis85, sdDist=sdDis85, radius=radius85)
}
save(centroids26, centroids85, file='output/centroid_summaries_June2017_update.RData')
  
centroids85$species <- projspp
centroids26$species <- projspp
# probSpp <- centroids85[centroids85$meanLat < 0 & centroids85$lonPre > -100 & centroids85$latPre > 40,]
regionFreqB <- ifelse(regionFreq=='AFSC_WCTri', 'NWFSC_WCAnn', regionFreq)
regionFreqB <- ifelse(regionFreqB=='VIMS_NEAMAP', 'NEFSC_NEUS', regionFreqB)
regionFreqB <- ifelse(regionFreqB=='DFO_SoGulf', 'DFO_ScotianShelf', regionFreqB)
centroids26$region <- regionFreqB
centroids85$region <- regionFreqB

# ===============================================================================================
# Uncertainty analysis for centroids=============================================================
# ===============================================================================================

# Uncertainty method #1_aggreement with projected distance and angle among models

uncer26 <- lm(log(sdDist)~log(meanDist), centroids26)
qqPlot(uncer26)
hist(resid(uncer26))
centroids26$resid <- resid(uncer26)

uncer85 <- lm(log(sdDist)~log(meanDist), centroids85)
qqPlot(uncer85)
hist(resid(uncer85))
centroids85$resid <- resid(uncer85)

rad26 <- quantile(centroids26$radius, probs=c(.25, .05))
sd26 <- quantile(centroids26$resid, probs=c(.75, .95))
rad85 <- quantile(centroids85$radius, probs=c(.25, .05))
sd85 <- quantile(centroids85$resid, probs=c(.75, .95))

lm26 <- lm(log(sdLat)~log(abs(meanLat)), data=centroids26)
qqPlot(lm26)
hist(resid(lm26))
quan95.26 <- rq(log(sdLat)~log(abs(meanLat)), tau=.95, data=centroids26)
quan75.26 <- rq(log(sdLat)~log(abs(meanLat)), tau=.75, data=centroids26)
centroids26$uncerCat2 <- ifelse(log(centroids26$sdLat) > (coef(quan95)[1] + coef(quan95)[2]*log(abs(centroids26$meanLat))), 'high uncer.',
              ifelse(log(centroids26$sdLat) < (coef(quan75)[1] + coef(quan75)[2]*log(abs(centroids26$meanLat))), 'low uncer.', 'med uncer.'))

lm85 <- lm(log(sdLat)~log(abs(meanLat)), data=centroids85)
qqPlot(lm85)
hist(resid(lm85))
quan95.85 <- rq(log(sdLat)~log(abs(meanLat)), tau=.95, data=centroids85)
quan75.85 <- rq(log(sdLat)~log(abs(meanLat)), tau=.75, data=centroids85)
centroids85$uncerCat2 <- ifelse(log(centroids85$sdLat) > (coef(quan95)[1] + coef(quan95)[2]*log(abs(centroids85$meanLat))), 'high uncer.',
              ifelse(log(centroids85$sdLat) < (coef(quan75)[1] + coef(quan75)[2]*log(abs(centroids85$meanLat))), 'low uncer.', 'med uncer.'))
#points(log(sdLat)~log(abs(meanLat)), col='red', data=centroids85[centroids85$uncerCat2 == 'high uncer.',])
#points(log(sdLat)~log(abs(meanLat)), col='green', data=centroids85[centroids85$uncerCat2 == 'med uncer.',])
#points(log(sdLat)~log(abs(meanLat)), col='blue', data=centroids85[centroids85$uncerCat2 == 'low uncer.',])

centroids26$uncerCat <- ifelse((centroids26$radius < rad26[2] | centroids26$resid > sd26[2]), "high unc.", 
              ifelse((centroids26$radius > rad26[1] & centroids26$resid < sd26[1]), "low unc.", "med unc."))  
centroids85$uncerCat <- ifelse((centroids85$radius < rad85[2] | centroids85$resid > sd85[2]), "high unc.", 
                               ifelse((centroids85$radius > rad85[1] & centroids85$resid < sd85[1]), "low unc.", "med unc."))  
save(centroids26, centroids85, file='output/centroid_summaries_June2017_update.RData')
  
# The plot of the above uncertainty stuff ===================================

pdf(width=7, height=7, file=paste('figures/', 'uncertainty_analysis_update.pdf', sep=''))
par(mfrow=c(3,2), xpd=F, tcl=-.1, pch=20, cex.axis=1.1, cex.lab=1.1, mgp=c(1.5,0.1,0), mar=c(2.7,3,.4,.5), oma=c(.1,.1,.1,.1))  

plot(log(sdDist)~log(meanDist), ylab='Model uncertainty (logged st. dev.)', xlab='Mean projected shift (logged km)', centroids26)
abline(coef(uncer26))
mtext('  A', side=3, line=-1.5, adj=0)
plot(log(sdDist)~log(meanDist), ylab='Model uncertainty (logged st. dev.)', xlab='Mean projected shift (logged km)', centroids85)
abline(coef(uncer85))
mtext('  B', side=3, line=-1.5, adj=0)

plot(resid~radius, xlab='Directional aggreement among models', ylab='Model uncertainty of shift distance', 
     col=ifelse(centroids26$uncerCat == 'low unc.', 'blue', ifelse(centroids26$uncerCat == 'med unc.', 'darkorange' ,'red2')), data=centroids26)
abline(h=sd26, v=rad26, col=c('blue', 'red2'))
mtext('  C', side=3, line=-1.5, adj=0)

plot(resid~radius, xlab='Directional aggreement among models', ylab='Model uncertainty of shift distance', 
     col=ifelse(centroids85$uncerCat == 'low unc.', 'blue', ifelse(centroids85$uncerCat == 'med unc.', 'darkorange' ,'red2')), data=centroids85)
abline(h=sd85, v=rad85, col=c('blue', 'red2'))
mtext('  D', side=3, line=-1.5, adj=0)
legend(x=.2, y=-0.22, xjust=.5, cex=1.4, legend=c('low', 'medium', 'high'), col=c('blue', 'darkorange', 'red2'), pt.cex=1, x.intersp=.75, text.col = "black", bg = "white", bty='n', inset=c(-0.2,0), pch=19)

plot(log(sdLat)~log(abs(meanLat)), xlab='Absolute value of mean shift (logged lat)', ylab='Model uncertainty (logged st. dev.)', col=centroids26$color2, data=centroids26)
abline(coef(lm26))
abline(coef(quan95.26), col='red2')
abline(coef(quan75.26), col='blue')
mtext('  E', side=3, line=-1.5, adj=0)

plot(log(sdLat)~log(abs(meanLat)), xlab='Absolute value of mean shift (logged lat)', ylab='Model uncertainty (logged st. dev.)', col=centroids85$color2, data=centroids85)
abline(coef(lm85))
abline(coef(quan95.85), col='red2')
abline(coef(quan75.85), col='blue')
mtext('  F', side=3, line=-1.5, adj=0)

dev.off()
 
# ============================================================================================
# SHOW THE ONE BIG MAP WITH ALL THE ACTUAL VECTORS IN PLACE
# Could do a chi-squared table to see if the two metrics are related to e/o
# ============================================================================================
  
#MAKE LONGITUDE WITH MORE TICKS_REDUCE SIZE OF PROJ_GRID SO IT DOESN'T OVERLAP LAND
#MAKE IT A TIFF WITH RIGHT DIMENSIONS_TEXT SIZE FOR AXES_ADD RCP VALUE
load('data/ProjectionBathGrid_Feb27_2017.RData')# load bathymetry projection grid

Projmap <- map('world', xlim=c(min(centroids26$lonPre),max(centroids26$lonPre)), ylim=c(min(centroids26$latPre),max(centroids26$latPost)),plot=F) 
Projmap <- data.frame(lon=Projmap$x, lat=Projmap$y)
pdf(width=14, height=8, file=paste('figures/', 'species_centroid_shifts26.pdf', sep=''))
#proj.grid2 <- proj.grid[!(proj.grid$latClimgrid > 60.5 & proj.grid$lonClimgrid > -64.62603),]
#proj.grid2 <- proj.grid2[!(proj.grid2$latClimgrid > 55 & proj.grid2$lonClimgrid < -64.6 & proj.grid2$lonClimgrid > -75),]
plot(latClimgrid~lonClimgrid, col='gray90', cex=.7, pch=16, xlab='Longitude', ylab='Latitude', main='Predicted centroid shifts for 658 species_RCP26', data=proj.grid) 
points(lat ~ lon, Projmap, type='l', lty=1, lwd=.7, col='gray50')
high26 <- centroids26[centroids26$uncerCat == 'high unc.',]
med26 <- centroids26[centroids26$uncerCat == 'med unc.',]
low26 <- centroids26[centroids26$uncerCat == 'low unc.',]
arrows(x0=low26$lonPre, y0=low26$latPre, x1=low26$lonPost, y1=low26$latPost, length=.05, lwd=.8, col='dark blue')
arrows(x0=med26$lonPre, y0=med26$latPre, x1=med26$lonPost, y1=med26$latPost, length=.05, lwd=.8, col='dark green')
arrows(x0=high26$lonPre, y0=high26$latPre, x1=high26$lonPost, y1=high26$latPost, length=.05, lwd=.8, col='dark red')
dev.off()

pdf(width=14, height=8, file=paste('figures/', 'species_centroid_shifts85.pdf', sep=''))
#proj.grid2 <- proj.grid[!(proj.grid$latClimgrid > 60.5 & proj.grid$lonClimgrid > -64.62603),]
#proj.grid2 <- proj.grid2[!(proj.grid2$latClimgrid > 55 & proj.grid2$lonClimgrid < -64.6 & proj.grid2$lonClimgrid > -75),]
plot(latClimgrid~lonClimgrid, col='gray90', cex=.7, pch=16, xlab='Longitude', ylab='Latitude', main='Predicted centroid shifts for 658 species_RCP85', data=proj.grid) 
points(lat ~ lon, Projmap, type='l', lty=1, lwd=.7, col='gray50')
high85 <- centroids85[centroids85$uncerCat == 'high unc.',]
med85 <- centroids85[centroids85$uncerCat == 'med unc.',]
low85 <- centroids85[centroids85$uncerCat == 'low unc.',]
arrows(x0=low85$lonPre, y0=low85$latPre, x1=low85$lonPost, y1=low85$latPost, length=.05, lwd=.8, col='dark blue')
arrows(x0=med85$lonPre, y0=med85$latPre, x1=med85$lonPost, y1=med85$latPost, length=.05, lwd=.8, col='dark green')
arrows(x0=high85$lonPre, y0=high85$latPre, x1=high85$lonPost, y1=high85$latPost, length=.05, lwd=.8, col='dark red')
dev.off()

# NOW WITH THE OTHER METHOD OF DETERMINING UNCERTAINTY
pdf(width=14, height=8, file=paste('figures/', 'species_centroid_shifts26_meth2.pdf', sep=''))
#proj.grid2 <- proj.grid[!(proj.grid$latClimgrid > 60.5 & proj.grid$lonClimgrid > -64.62603),]
#proj.grid2 <- proj.grid2[!(proj.grid2$latClimgrid > 55 & proj.grid2$lonClimgrid < -64.6 & proj.grid2$lonClimgrid > -75),]
plot(latClimgrid~lonClimgrid, col='gray90', cex=.7, pch=16, xlab='Longitude', ylab='Latitude', main='Predicted centroid shifts for 658 species_RCP26', data=proj.grid) 
points(lat ~ lon, Projmap, type='l', lty=1, lwd=.7, col='gray50')
high26 <- centroids26[centroids26$uncerCat2 == 'high uncer.',]
med26 <- centroids26[centroids26$uncerCat2 == 'med uncer.',]
low26 <- centroids26[centroids26$uncerCat2 == 'low uncer.',]
arrows(x0=low26$lonPre, y0=low26$latPre, x1=low26$lonPost, y1=low26$latPost, length=.05, lwd=.8, col='dark blue')
arrows(x0=med26$lonPre, y0=med26$latPre, x1=med26$lonPost, y1=med26$latPost, length=.05, lwd=.8, col='dark green')
arrows(x0=high26$lonPre, y0=high26$latPre, x1=high26$lonPost, y1=high26$latPost, length=.05, lwd=.8, col='dark red')
dev.off()

pdf(width=14, height=8, file=paste('figures/', 'species_centroid_shifts85_meth2.pdf', sep=''))
#proj.grid2 <- proj.grid[!(proj.grid$latClimgrid > 60.5 & proj.grid$lonClimgrid > -64.62603),]
#proj.grid2 <- proj.grid2[!(proj.grid2$latClimgrid > 55 & proj.grid2$lonClimgrid < -64.6 & proj.grid2$lonClimgrid > -75),]
plot(latClimgrid~lonClimgrid, col='gray90', cex=.7, pch=16, xlab='Longitude', ylab='Latitude', main='Predicted centroid shifts for 658 species_RCP85', data=proj.grid) 
points(lat ~ lon, Projmap, type='l', lty=1, lwd=.7, col='gray50')
high85 <- centroids85[centroids85$uncerCat2 == 'high uncer.',]
med85 <- centroids85[centroids85$uncerCat2 == 'med uncer.',]
low85 <- centroids85[centroids85$uncerCat2 == 'low uncer.',]
arrows(x0=low85$lonPre, y0=low85$latPre, x1=low85$lonPost, y1=low85$latPost, length=.05, lwd=.8, col='dark blue')
arrows(x0=med85$lonPre, y0=med85$latPre, x1=med85$lonPost, y1=med85$latPost, length=.05, lwd=.8, col='dark green')
arrows(x0=high85$lonPre, y0=high85$latPre, x1=high85$lonPost, y1=high85$latPost, length=.05, lwd=.8, col='dark red')
dev.off()

# ===================================================================================================
# COMPASS PLOTS=======================================================================================
# ===================================================================================================

centroids26$bearing <- bearing(p1=matrix(cbind(centroids26$lonPre, centroids26$latPre), ncol=2), p2=matrix(cbind(centroids26$lonPost, centroids26$latPost), ncol=2))
centroids26$bearing_scaled <- ifelse(centroids26$bearing < 0, centroids26$bearing + 360, centroids26$bearing)

centroids85$bearing <- bearing(p1=matrix(cbind(centroids85$lonPre, centroids85$latPre), ncol=2), p2=matrix(cbind(centroids85$lonPost, centroids85$latPost), ncol=2))
centroids85$bearing_scaled <- ifelse(centroids85$bearing < 0, centroids85$bearing + 360, centroids85$bearing)

#Make the a tiff ==============================================================================
regions <- names(table(regionFreqB))

centroids26$uncertainty <- ifelse(centroids26$uncerCat == 'high unc.' | centroids26$uncerCat2 == 'high uncer.', 'high', 
                          ifelse(centroids26$uncerCat == 'med unc.' | centroids26$uncerCat2 == 'med uncer.', 'medium', 'low'))
centroids85$uncertainty <- ifelse(centroids85$uncerCat == 'high unc.' | centroids85$uncerCat2 == 'high uncer.', 'high', 
                                  ifelse(centroids85$uncerCat == 'med unc.' | centroids85$uncerCat2 == 'med uncer.', 'medium', 'low'))
centroids26$npres <- modeldiag$npres
centroids85$npres <- modeldiag$npres

plot(centroids26$npres~as.factor(centroids26$uncertainty))
summary(centroids26$npres[centroids26$uncertainty=='low'])
anova(aov(centroids26$npres~centroids26$uncertainty))

plot(centroids85$npres~as.factor(centroids85$uncertainty))
summary(centroids85$npres[centroids85$uncertainty=='high'])
anova(aov(centroids85$npres~centroids85$uncertainty))
 
centroids26$color <- ifelse(centroids26$uncertainty == 'high', 'red2', ifelse(centroids26$uncertainty == 'medium', 'darkorange', 'blue'))
centroids85$color <- ifelse(centroids85$uncertainty == 'high', 'red2', ifelse(centroids85$uncertainty == 'medium', 'darkorange', 'blue'))
#centroids26$color2 <- ifelse(centroids26$uncerCat2 == 'high uncer.', 'red2', ifelse(centroids26$uncerCat2 == 'med uncer.', 'darkorange', 'blue'))
#centroids85$color2 <- ifelse(centroids85$uncerCat2 == 'high uncer.', 'red2', ifelse(centroids85$uncerCat2 == 'med uncer.', 'darkorange', 'blue'))
regions <- names(table(regionFreqB))

pdf(width=8, height=4, file=paste('figures/', 'compass_plots_update', '.pdf', sep=''))
par(mfrow=c(1,2), cex.main=.75)

for(i in 1:9){
  region = regions[i]
  data26 <- centroids26[centroids26$region == region,]
  data85 <- centroids85[centroids85$region == region,]

  polar.plot(lengths=data26$meanDist, polar.pos=data26$bearing_scaled, radial.lim=c(0,1500), grid.col='gray40' , radial.labels=c('','',''), lwd=2.5, line.col=data26$color, labels=c(), label.pos=c(0,90,180,270), boxed.radial=F, start=90, clockwise=T, main=paste(region, 'rcp26',sep='_'))
  polar.plot(lengths=data85$meanDist, polar.pos=data85$bearing_scaled, radial.lim=c(0,1500), grid.col='gray40' , radial.labels=c('','',''), lwd=2.5, line.col=data85$color, labels=c(), label.pos=c(0,90,180,270), boxed.radial=F, start=90, clockwise=T, main=paste(region, 'rcp85',sep='_'))
  #polar.plot(lengths=data26$meanDist, polar.pos=data26$bearing_scaled, radial.lim=c(0,1500), grid.col='gray80' , radial.labels=c('','',''), lwd=1.5, line.col=data26$color2, labels=c(), label.pos=c(0,90,180,270), boxed.radial=F, start=90, clockwise=T, main=paste(region, 'rcp26_B',sep='_'))
  #polar.plot(lengths=data85$meanDist, polar.pos=data85$bearing_scaled, radial.lim=c(0,1500), grid.col='gray80' , radial.labels=c('','',''), lwd=1.5, line.col=data85$color2, labels=c(), label.pos=c(0,90,180,270), boxed.radial=F, start=90, clockwise=T, main=paste(region, 'rcp85_B',sep='_'))
} 
dev.off()

load('ProjectionBathGrid_Feb27_2017.RData')# load bathymetry projection grid

# Make the background projection grid
pdf(width=6, height=6, file=paste('figures/', 'Map_template_East.pdf', sep=''))
par(xpd=F, tcl=-.2, mgp=c(1.8,.5,0), cex.axis=.8) 
Projmap <- map('world', xlim=c(-96,-45.5), ylim=c(26,60), plot=F) 
Projmap <- data.frame(lon=Projmap$x, lat=Projmap$y)
plot(latClimgrid~lonClimgrid, col='gray93', cex=.4, pch=16, xlim=c(-96,-45.5), ylim=c(26,60), xlab='Longitude', ylab='Latitude', data=proj.grid) 
points(lat ~ lon, Projmap, type='l', lty=1, lwd=.7, col='gray50')
dev.off()

pdf(width=6, height=6, file=paste('figures/', 'Map_template_West.pdf', sep=''))
par(xpd=F, tcl=-.2, mgp=c(1.8,.5,0), cex.axis=.8) 
Projmap <- map('world', xlim=c(-188,-117), ylim=c(32.5,62), plot=F) 
Projmap <- data.frame(lon=Projmap$x, lat=Projmap$y)
plot(latClimgrid~lonClimgrid, col='gray93', cex=.4, pch=16, xlim=c(-188,-117), ylim=c(32.5,62), xlab='Longitude', ylab='Latitude', data=proj.grid) 
points(lat ~ lon, Projmap, type='l', lty=1, lwd=.7, col='gray50')
dev.off()
 
#Map inset
pdf(width=9, height=6, file=paste('figures/', 'Map_legend.pdf', sep=''))
par(xpd=F, tcl=-.2, mgp=c(1.8,.5,0), cex.axis=.8) 
Projmap <- map('world', xlim=c(-185,-50), ylim=c(25,60), plot=F) 
Projmap <- data.frame(lon=Projmap$x, lat=Projmap$y)
#plot(latClimgrid~lonClimgrid, col='gray93', cex=.4, pch=16, xlim=c(-96,-45.5), ylim=c(26,60), xlab='Longitude', ylab='Latitude', data=proj.grid) 
plot(lat ~ lon, Projmap, type='l', lty=1, lwd=.7, col='gray40')
dev.off()


# make a legend
pdf(width=8, height=4, file=paste('figures/', 'compass_plots_LEGEND_update', '.pdf', sep=''))
par(mfrow=c(1,2), cex.lab=1.1)
polar.plot(lengths=0, polar.pos=0, radial.lim=c(0,1500), grid.col='gray40' , radial.labels=c('0','500','1000','1500 km'), labels=c(), show.grid.labels=1, label.pos=c(0,90,180,270), boxed.radial=T, start=90, clockwise=T, main='Legend')
dev.off()

# Omit species with high or medium uncertainty with either metric (n = 77) and do a mean mag and direction within regions
# Would be two compass plots, one for each RCP, but scale them so you can acually see rcp26, colorize by region
# Can I get a mean of vectors, incorporating direction and magnitude, probably easy to do and good way to compare regions and RCPs

trimmed26 <- centroids26[centroids26$uncertainty == 'low',] 
trimmed85 <- centroids85[centroids85$uncertainty == 'low',] 
#trimmed26 <- centroids26[!(centroids26$uncerCat == 'high unc.' | centroids26$uncerCat2 == 'high uncer.'),] 
#trimmed85 <- centroids85[!(centroids85$uncerCat == 'high unc.' | centroids85$uncerCat2 == 'high uncer.'),] 

compSumm <- data.frame(dis26=numeric(), direc26=numeric(), dis85=numeric(), direc85=numeric())
for(i in 1:9){
  region = regions[i]
  data26 <- trimmed26[trimmed26$region == region,]
  data85 <- trimmed85[trimmed85$region == region,]
  # Calculate mean bearing and magnitude
  dist26 <- vector.averaging(direction=data26$bearing_scaled, distance=data26$meanDist)[[1]]
  dir26 <- circular.averaging(direction=data26$bearing_scaled)
  #if(dir26 < 0){
   # dir26 <- dir26 + 360
  #}
  dist85 <- vector.averaging(direction=data85$bearing_scaled, distance=data85$meanDist)[[1]]
  dir85 <- circular.averaging(direction=data85$bearing_scaled)
  #if(dir85 < 0){
   # dir85 <- dir85 + 360
  #}
  compSumm[i,] <- data.frame(dis26=dist26, direc26=dir26, dis85=dist85, direc85=dir85)
} 
 
compSumm$region <- regions
compSumm$colors <- c('blue', 'darkorange', 'black', 'purple', 'cyan3', 'brown', 'green', 'gray60', 'red2')

# Need to rearrange the dataframe so that the colors lay right and show up better (figured this after much plotting trials)
abc <- c("DFO_Newfoundland","AFSC_Aleutians","DFO_ScotianShelf","NEFSC_NEUS","NWFSC_WCAnn","SCDNR_SEUS","AFSC_EBS","AFSC_GOA","SEFSC_GOMex")     
compSummB <- compSumm[match(abc, compSumm$region),]
legend.names <- c("Newfoundland","Aleutian Is.","Scotian Shelf","Northeast U.S.","West U.S.","Southeast U.S.","East. Bering S.","Gulf of Alaska","Gulf of Mexico")     

pdf(file='figures/compass_summary.pdf', width = 7, height = 4)
par(mfrow=c(1,2), cex.lab=1, cex.main=1.2, mar=c(.1,.1,.1,.1), oma=c(.1,.1,.1,.1))
polar.plot(lengths=compSummB$dis26, polar.pos=compSummB$direc26, radial.lim=c(0,1000), mar=c(.1,.5,.1,.5), grid.col='gray75', radial.labels=c('','200','400','600','800', '1000 km'), show.grid.labels=1, lwd=6, line.col=compSummB$colors, labels=c(), label.pos=c(0,90,180,270), boxed.radial=T, start=90, clockwise=T)
title('A) RCP 2.6', line=-1.1, adj=0, cex=.4)

polar.plot(lengths=compSummB$dis85, polar.pos=compSummB$direc85, radial.lim=c(0,1000), mar=c(.1,.5,.1,.5),grid.col='gray75', radial.labels=c('','','','','', ''), show.grid.labels=1, lwd=6, line.col=compSummB$colors, labels=c(), label.pos=c(0,90,180,270), boxed.radial=T, start=90, clockwise=T)
title('B) RCP 8.5', line=-1.1, adj=0)
legend(x=-650, y=-150, xjust=.5, cex=.8, legend=legend.names, col=compSummB$colors, lty=1, lwd=4, pt.cex=0, x.intersp=.75, text.col = "black", bg = "white", box.lwd=.5, box.col='gray75', inset=c(-0.2,0), pch=c(1,3))
dev.off()
    
save(centroids26, centroids85, file='centroid_summaries_June2017_update.RData')

# ============================================================================================================
# ============================================================================================================
# Thermal biomass analysis
# ============================================================================================================
# ============================================================================================================

# columns for init and final; mean + SD in biomass change
biomass26 <- data.frame(biomMeanChan=numeric(), sdChan=numeric())
biomass85 <- data.frame(biomMeanChan=numeric(), sdChan=numeric())
 
for(i in 1:length(projspp)){ # one loop for each species
  species = projspp[i]
  print(paste("Beginning figures for ", species, '_spp', i, sep=''))
  filename <- paste(projfolder, species, '_rcp', rcp[1], '_jas_prediction_AGG.RData', sep='')
  load(filename)
  pred.agg26 <- pred.agg
  
  filename <- paste(projfolder, species, '_rcp', rcp[2], '_jas_prediction_AGG.RData', sep='')
  load(filename)
  pred.agg85 <- pred.agg
  rm(pred.agg)
  # RCP26 calculations
  yearRange <- c('2007-2020', '2081-2100')
  biom_preds26 <- matrix(data=NA, nrow=2, ncol=16)
  for(j in 1:2){
    years = yearRange[j]
    preds <- pred.agg26[pred.agg26$year_range == years,]
    for(k in 4:19){
      valueBiom = sum(preds[,k]) 
      biom_preds26[j,k-3] = valueBiom
    }
  } 
  # Need to calculate change in biomass for each model_then get the mean/SD of those values
  biom_diff26 <- apply(biom_preds26, 2, FUN=diff, na.rm=T)
  biom_perc26 <- (biom_diff26/biom_preds26[1,]) *100
  meanBiom26 <- mean(biom_perc26)
  sdBiom26 <- sd(biom_perc26)

  # RCP85 calculations
  biom_preds85 <- matrix(data=NA, nrow=2, ncol=16)
  for(j in 1:2){
    years = yearRange[j]
    preds <- pred.agg85[pred.agg85$year_range == years,]
    for(k in 4:19){
      valueBiom = sum(preds[,k]) 
      biom_preds85[j,k-3] = valueBiom
    }
  } 
  # Need to calculate change in biomass for each model_then get the mean/SD of those values
  biom_diff85 <- apply(biom_preds85, 2, FUN=diff, na.rm=T)
  biom_perc85 <- (biom_diff85/biom_preds85[1,]) *100
  meanBiom85 <- mean(biom_perc85)
  sdBiom85 <- sd(biom_perc85)
  
  biomass26[i,] <- data.frame(biomMeanChan=meanBiom26, sdChan=sdBiom26)
  biomass85[i,] <- data.frame(biomMeanChan=meanBiom85, sdChan=sdBiom85)
}
 
biomass26$species <- projspp
biomass85$species <- projspp
# probSpp <- centroids85[centroids85$meanLat < 0 & centroids85$lonPre > -100 & centroids85$latPre > 40,]
biomass26$region <- regionFreqB
biomass85$region <- regionFreqB
 
#biomass26$uncerCat <- centroids26$uncerCat
#biomass26$uncerCat2 <- centroids26$uncerCat2
#biomass85$uncerCat <- centroids85$uncerCat
#biomass85$uncerCat2 <- centroids85$uncerCat2
#biomass26$dropspp <- ifelse((biomass26$uncerCat == 'high unc.' | biomass26$uncerCat2 == 'high uncer.'), 'drop', 'keep')
#biomass85$dropspp <- ifelse((biomass85$uncerCat == 'high unc.' | biomass85$uncerCat2 == 'high uncer.'), 'drop', 'keep')
biomass26$uncertainty <- centroids26$uncertainty
biomass85$uncertainty <- centroids85$uncertainty

save(biomass26, biomass85, file='Therm.Biomass_summaries_June2017.RData')

# ============================================================================================================
# Thermal biomass uncertainty analysis
# ============================================================================================================

load('Therm.Biomass_summaries_June2017.RData')
# a couple qaqc plots
plot(sdChan~biomMeanChan, xlim=c(-110,600), ylim=c(0,2000),data=biomass26[biomass26$dropspp == 'keep',])
plot(sdChan~biomMeanChan, xlim=c(-110,2200), ylim=c(0,2000),data=biomass85[biomass85$dropspp == 'keep',])
 
# Assign numbers to regions for the plot (based on prelim looks), so gives an ascending order
biomass26$xaxis <- ifelse(biomass26$region=='AFSC_Aleutians', 1, ifelse(biomass26$region=='AFSC_EBS', 2,
              ifelse(biomass26$region=='AFSC_GOA', 6, ifelse(biomass26$region=='DFO_Newfoundland', 3,
              ifelse(biomass26$region=='DFO_ScotianShelf', 4, ifelse(biomass26$region=='NEFSC_NEUS', 7, 
              ifelse(biomass26$region=='NWFSC_WCAnn', 9, ifelse(biomass26$region=='SCDNR_SEUS', 5,  8)))))))) 
biomass85$xaxis <- ifelse(biomass85$region=='AFSC_Aleutians', 1, ifelse(biomass85$region=='AFSC_EBS', 2,
              ifelse(biomass85$region=='AFSC_GOA', 6, ifelse(biomass85$region=='DFO_Newfoundland', 3,
              ifelse(biomass85$region=='DFO_ScotianShelf', 4, ifelse(biomass85$region=='NEFSC_NEUS', 7, 
              ifelse(biomass85$region=='NWFSC_WCAnn', 9, ifelse(biomass85$region=='SCDNR_SEUS', 5,  8)))))))) 

#legend.names <- c("Aleut. Is.","E. Bering S.","G. Alaska","Newfound.","Scot. Shelf","NE U.S.","West U.S.","SE U.S.","G. Mexico")     
legend.names <- c("Aleut. Is.","E. Bering S.","Newfound.","Scot. Shelf","SE U.S.","G. Alaska","NE U.S.","G. Mexico","West U.S.")     

#boxplot(biomMeanChan~region, ylim=c(-75, 200), staplewex=.25, plot=F, data=biomass26[biomass26$dropspp == 'keep',])
boxplot(biomMeanChan~xaxis, ylim=c(-125, 950), staplewex=.25, plot=F, data=biomass85[biomass85$dropspp == 'keep',])
     
pdf(file='thermBiomass_boxPlots_UPDATE.pdf', width = 7, height = 5) 
par(mfrow=c(1,2), cex.lab=1, cex.main=1.2, mar=c(5,2.5,1.5,.1), oma=c(.1,.1,.1,.1))
boxplot(biomMeanChan~xaxis, ylim=c(-75, 165), xaxt='n', yaxt='n', pch=20, staplewex=.25, names=legend.names, data=biomass26[biomass26$uncertainty == 'low',])
abline(a=0, b=0, lwd=.5, col='gray50')
boxplot(biomMeanChan~xaxis, ylim=c(-75, 165), xaxt='n', yaxt='n', pch=20, staplewex=.25, names=legend.names, col="white", data=biomass26[biomass26$uncertainty == 'low',], add=TRUE)
axis(side=1, at=seq(1:9), labels=F, tcl=-.2)
text(cex=1, x=c(1:9)-.4, y=-93, legend.names, pos=4, xpd=TRUE, srt=290)
text(cex=.7, x=c(1:9), y=-69, labels=c(17,56,37,6,16,32,25,56,79), pos=1, xpd=TRUE)
text(cex=.7, x=c(1:9), y=178, labels=c('','','','',2,1,'','',''), pos=1, xpd=TRUE)
axis(side=2,tcl=-.2, padj=1)
mtext('Change in thermal habitat (% of 2006-2020)', side=2, srt=290, padj=-2.7)
mtext('A) RCP 2.6', adj=0.05, line=0.2)

boxplot(biomMeanChan~xaxis, ylim=c(-115, 860), xaxt='n', yaxt='n', staplewex=.25, pch=20, names=legend.names, data=biomass85[biomass85$uncertainty == 'low',])
abline(a=0, b=0, lwd=.5, col='gray50')
boxplot(biomMeanChan~xaxis, ylim=c(-115, 860), xaxt='n', yaxt='n', staplewex=.25, col="white", pch=20, names=legend.names, data=biomass85[biomass85$uncertainty == 'low',], add=TRUE)
axis(side=1, at=seq(1:9), labels=F, tcl=-.2)
text(cex=1, x=seq(1:9)-.4, y=-185, legend.names, pos=4, xpd=TRUE, srt=290)
text(cex=.7, x=c(1:9), y=-90, labels=c(6,44,48,8,18,29,39,73,105), pos=1, xpd=TRUE)
text(cex=.7, x=c(1:9), y=915, labels=c('','','','','',1,1,7,1), pos=1, xpd=TRUE)
axis(side=2, at=c(-100, 0, 100, 300, 500, 700, 900), tcl=-.2, padj=1)
mtext('B) RCP 8.5', adj=0.05, line=0.2)
dev.off()
  
# ===========================================================================================================
# example graphs of biomass shift with specific species
# ===========================================================================================================
# EAST COAST 

  species1 = 'lutjanus griseus_Atl'
  species3 = 'gadus morhua_Atl'
  species2 = 'scomberomorus maculatus_Atl'
  
  filename <- paste(species1, '_rcp', rcp[2], '_jas_prediction_AGG.RData', sep='')
  load(filename)
  pred.agg85 <- as.matrix(pred.agg[,4:19], nrow=412960, ncol=16)
  # aggregate by lon lat and year range
  ensMean <- apply(pred.agg85, 1, FUN=mean, na.rm=T)
  pred.agg1 <- data.frame(cbind(pred.agg, ensMean=ensMean))
  
  filename <- paste(species2, '_rcp', rcp[2], '_jas_prediction_AGG.RData', sep='')
  load(filename)
  pred.agg85 <- as.matrix(pred.agg[,4:19], nrow=412960, ncol=16)
  # aggregate by lon lat and year range
  ensMean <- apply(pred.agg85, 1, FUN=mean, na.rm=T)
  pred.agg2 <- data.frame(cbind(pred.agg, ensMean=ensMean))

  filename <- paste(species3, '_rcp', rcp[2], '_jas_prediction_AGG.RData', sep='')
  load(filename)
  pred.agg85 <- as.matrix(pred.agg[,4:19], nrow=412960, ncol=16)
  # aggregate by lon lat and year range
  ensMean <- apply(pred.agg85, 1, FUN=mean, na.rm=T)
  pred.agg3 <- data.frame(cbind(pred.agg, ensMean=ensMean))
  
  pdf(width=7.5, height=7.75, file='figures/rcp85_ensemble_projections_ATLANTICspp.pdf')
    
    cols = colorRampPalette(colors = c('gray92', 'blue', 'dark blue', 'black'))
    
    scale85 = seq(0, max(pred.agg1$ensMean) + .00001, length.out=20)
    # STANDARDIZE PLOT DIMENSIONS FOR ALL PLOTS
    xlimit = c(min(pred.agg1$longitude[pred.agg1$ensMean > scale85[2]/2] - 1),max(pred.agg1$longitude[pred.agg1$ensMean > scale85[2]/2] + 1))
    ylimit = c(min(pred.agg1$latitude[pred.agg1$ensMean > scale85[2]/2] - .25), max(pred.agg1$latitude[pred.agg1$ensMean > scale85[2]/2] + .25))

    # RCP85 plots_one for normal, one for Logged+1
    Projmap <- map('world', xlim=xlimit, ylim=ylimit, plot=F) 
    Projmap <- data.frame(lon=Projmap$x, lat=Projmap$y)
    plot1 <- levelplot(ensMean ~ longitude*latitude|year_range, data=pred.agg1, index.cond=list(c(5,4,3,2,1)),  colorkey=F, par.settings=list(layout.heights=list(top.padding=-2, bottom.padding=-2),layout.widths=list(left.padding=-2, right.padding=-2.5)),
                    scales=list(tck=c(0,0), alternating=c(0,0)), strip=F, strip.left=TRUE,  aspect='fill', xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL,  
                    at = scale85, col.regions=cols, layout = c(1, 5)) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray',
                    par.settings=list(layout.heights=list(top.padding=0, bottom.padding=0)))
    
    scale85 = seq(0, max(pred.agg2$ensMean) + .00001, length.out=20)
    # STANDARDIZE PLOT DIMENSIONS FOR ALL PLOTS
    xlimit = c(min(pred.agg2$longitude[pred.agg2$ensMean > scale85[2]/2] - 1),max(pred.agg2$longitude[pred.agg2$ensMean > scale85[2]/2] + 1))
    ylimit = c(min(pred.agg2$latitude[pred.agg2$ensMean > scale85[2]/2] - .25), max(pred.agg2$latitude[pred.agg2$ensMean > scale85[2]/2] + .25))

    # RCP85 plots_one for normal, one for Logged+1
    Projmap <- map('world', xlim=xlimit, ylim=ylimit, plot=F) 
    Projmap <- data.frame(lon=Projmap$x, lat=Projmap$y)
    plot2 <-levelplot(ensMean ~ longitude*latitude|year_range, data=pred.agg2, index.cond=list(c(5,4,3,2,1)), colorkey=F, par.settings=list(layout.heights=list(top.padding=-2, bottom.padding=-2),layout.widths=list(left.padding=-2.25, right.padding=-2.25)),
                      scales=list(tck=c(0,0), alternating=c(0,0)), strip=F, aspect='fill', xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
                      at = scale85, col.regions=cols, layout = c(1, 5)) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
    
    scale85 = seq(0, max(pred.agg3$ensMean) + .00001, length.out=20)
    # STANDARDIZE PLOT DIMENSIONS FOR ALL PLOTS
    xlimit = c(min(pred.agg3$longitude[pred.agg3$ensMean > scale85[2]/2] - 1),max(pred.agg3$longitude[pred.agg3$ensMean > scale85[2]/2] + 1))
    ylimit = c(min(pred.agg3$latitude[pred.agg3$ensMean > scale85[2]/2] - .25), max(pred.agg3$latitude[pred.agg3$ensMean > scale85[2]/2] + .25))
    
    # RCP85 plots_one for normal, one for Logged+1
    Projmap <- map('world', xlim=xlimit, ylim=ylimit, plot=F) 
    Projmap <- data.frame(lon=Projmap$x, lat=Projmap$y)
    plot3 <-levelplot(ensMean ~ longitude*latitude|year_range, data=pred.agg3, index.cond=list(c(5,4,3,2,1)), colorkey=F, par.settings=list(layout.heights=list(top.padding=-2, bottom.padding=-2),layout.widths=list(left.padding=-2.5, right.padding=-2)),
                      scales=list(tck=c(0,0), alternating=c(0,0)), strip=F,  aspect='fill', xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
                      at = scale85, col.regions=cols, layout = c(1, 5)) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
    par(mar=c(.1,.1,.1,.1))
    print(plot1, split=c(1,1,3,1), more=TRUE) #position=c(0,0,.5,1), 
    print(plot2, split=c(2,1,3,1), more=TRUE)#position=c(.5,0,1,1))
    print(plot3, split=c(3,1,3,1))#position=c(.5,0,1,1))
    
  dev.off()

    
# WEST COAST
  species1 = 'theragra chalcogramma_Pac'
  species3 = 'sebastes pinniger_Pac'
  species2 = 'cancer productus_Pac'
    
  filename <- paste(species1, '_rcp', rcp[2], '_jas_prediction_AGG.RData', sep='')
  load(filename)
  pred.agg85 <- as.matrix(pred.agg[,4:19], nrow=329130, ncol=16)
  # aggregate by lon lat and year range
  ensMean <- apply(pred.agg85, 1, FUN=mean, na.rm=T)
  pred.agg1 <- data.frame(cbind(pred.agg, ensMean=ensMean))
  
  filename <- paste(species2, '_rcp', rcp[2], '_jas_prediction_AGG.RData', sep='')
  load(filename)
  pred.agg85 <- as.matrix(pred.agg[,4:19], nrow=329130, ncol=16)
  # aggregate by lon lat and year range
  ensMean <- apply(pred.agg85, 1, FUN=mean, na.rm=T)
  pred.agg2 <- data.frame(cbind(pred.agg, ensMean=ensMean))
  
  filename <- paste(species3, '_rcp', rcp[2], '_jas_prediction_AGG.RData', sep='')
  load(filename)
  pred.agg85 <- as.matrix(pred.agg[,4:19], nrow=329130, ncol=16)
  # aggregate by lon lat and year range
  ensMean <- apply(pred.agg85, 1, FUN=mean, na.rm=T)
  pred.agg3 <- data.frame(cbind(pred.agg, ensMean=ensMean))
  
  pdf(width=7.5, height=7.75, file='figures/rcp85_ensemble_projections_PACIFICspp.pdf')
  
  cols = colorRampPalette(colors = c('gray92', 'blue', 'dark blue', 'black'))
  
  scale85 = seq(0, max(pred.agg1$ensMean) + .00001, length.out=20)
  # STANDARDIZE PLOT DIMENSIONS FOR ALL PLOTS
  xlimit = c(min(pred.agg1$longitude[pred.agg1$ensMean > scale85[2]/2] - 1),max(pred.agg1$longitude[pred.agg1$ensMean > scale85[2]/2] + 1))
  ylimit = c(min(pred.agg1$latitude[pred.agg1$ensMean > scale85[2]/2] - .25), max(pred.agg1$latitude[pred.agg1$ensMean > scale85[2]/2] + .25))
  
  # RCP85 plots_one for normal, one for Logged+1
  Projmap <- map('world', xlim=xlimit, ylim=ylimit, plot=F) 
  Projmap <- data.frame(lon=Projmap$x, lat=Projmap$y)
  plot1 <- levelplot(ensMean ~ longitude*latitude|year_range, data=pred.agg1, index.cond=list(c(5,4,3,2,1)),  colorkey=F, par.settings=list(layout.heights=list(top.padding=-2, bottom.padding=-2),layout.widths=list(left.padding=-2, right.padding=-2.5)),
                     scales=list(tck=c(0,0), alternating=c(0,0)), strip=F, strip.left=TRUE,  aspect='fill', xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL,  
                     at = scale85, col.regions=cols, layout = c(1, 5)) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray',
                                                                                par.settings=list(layout.heights=list(top.padding=0, bottom.padding=0)))
  
  scale85 = seq(0, max(pred.agg2$ensMean) + .00001, length.out=20)
  # STANDARDIZE PLOT DIMENSIONS FOR ALL PLOTS
  xlimit = c(min(pred.agg2$longitude[pred.agg2$ensMean > scale85[2]/2] - 1),max(pred.agg2$longitude[pred.agg2$ensMean > scale85[2]/2] + 1))
  ylimit = c(min(pred.agg2$latitude[pred.agg2$ensMean > scale85[2]/2] - .25), max(pred.agg2$latitude[pred.agg2$ensMean > scale85[2]/2] + .25))
  
  # RCP85 plots_one for normal, one for Logged+1
  Projmap <- map('world', xlim=xlimit, ylim=ylimit, plot=F) 
  Projmap <- data.frame(lon=Projmap$x, lat=Projmap$y)
  plot2 <-levelplot(ensMean ~ longitude*latitude|year_range, data=pred.agg2, index.cond=list(c(5,4,3,2,1)), colorkey=F, par.settings=list(layout.heights=list(top.padding=-2, bottom.padding=-2),layout.widths=list(left.padding=-2.25, right.padding=-2.25)),
                    scales=list(tck=c(0,0), alternating=c(0,0)), strip=F, aspect='fill', xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
                    at = scale85, col.regions=cols, layout = c(1, 5)) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
  
  scale85 = seq(0, max(pred.agg3$ensMean) + .00001, length.out=20)
  # STANDARDIZE PLOT DIMENSIONS FOR ALL PLOTS
  xlimit = c(min(pred.agg3$longitude[pred.agg3$ensMean > scale85[2]/2] - 1),max(pred.agg3$longitude[pred.agg3$ensMean > scale85[2]/2] + 1))
  ylimit = c(min(pred.agg3$latitude[pred.agg3$ensMean > scale85[2]/2] - .25), max(pred.agg3$latitude[pred.agg3$ensMean > scale85[2]/2] + .25))
  
  # RCP85 plots_one for normal, one for Logged+1
  Projmap <- map('world', xlim=xlimit, ylim=ylimit, plot=F) 
  Projmap <- data.frame(lon=Projmap$x, lat=Projmap$y)
  plot3 <-levelplot(ensMean ~ longitude*latitude|year_range, data=pred.agg3, index.cond=list(c(5,4,3,2,1)), colorkey=F, par.settings=list(layout.heights=list(top.padding=-2, bottom.padding=-2),layout.widths=list(left.padding=-2.5, right.padding=-2)),
                    scales=list(tck=c(0,0), alternating=c(0,0)), strip=F,  aspect='fill', xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
                    at = scale85, col.regions=cols, layout = c(1, 5)) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
  par(mar=c(.1,.1,.1,.1))
  print(plot1, split=c(1,1,3,1), more=TRUE) #position=c(0,0,.5,1), 
  print(plot2, split=c(2,1,3,1), more=TRUE)#position=c(.5,0,1,1))
  print(plot3, split=c(3,1,3,1))#position=c(.5,0,1,1))
  
  dev.off()
  
# ====================================================================== 
# Redoing the loop from above for exploring species options
# ====================================================================== 

# I originally did this as a loop to look at all these potential species for example graphs
# now it's retooled to just make the graphs for species of interest.
testSpp <- c('caretta caretta_Atl','paralichthys albigutta_Atl','rhizoprionodon terraenovae_Atl',
             'rhinoptera bonasus_Atl','scomberomorus cavalla_Atl','scomberomorus maculatus_Atl',
             'stomolophus meleagris_Atl','trachinotus carolinus_Atl','myliobatis freminvillii_Atl',
             'gymnura micrura_Atl','menippe mercenaria_Atl','archosargus probatocephalus_Atl')
for(i in 1:length(testSpp)){
  species = testSpp[i]
  print(paste("Beginning figures for ", species, sep=''))
  filename <- paste(species, '_rcp', rcp[2], '_jas_prediction_AGG.RData', sep='')
  load(filename)
  pred.agg85 <- as.matrix(pred.agg[,4:19], nrow=329130, ncol=16)
  # aggregate by lon lat and year range
  ensMean <- apply(pred.agg85, 1, FUN=mean, na.rm=T)
  pred.agg85 <- data.frame(cbind(pred.agg, ensMean=ensMean))
  
  pdf(width=14, height=5, file=paste('figures/', species, '_rcp85_ensemble_projections.pdf', sep=''))

  #scale85log = seq(0, max(log(pred.agg1$ensMean + 1)) + .00001, length.out=20)
  scale85 = seq(0, max(pred.agg85$ensMean) + .00001, length.out=20)
  cols = colorRampPalette(colors = c('gray95', 'blue', 'dark blue', 'black'))
  
  # STANDARDIZE PLOT DIMENSIONS FOR ALL PLOTS
  xlimit = c(min(pred.agg85$longitude[pred.agg85$ensMean > scale85[2]/2] - 1),max(pred.agg85$longitude[pred.agg85$ensMean > scale85[2]/2] + 1))
  ylimit = c(min(pred.agg85$latitude[pred.agg85$ensMean > scale85[2]/2] - .25), max(pred.agg85$latitude[pred.agg85$ensMean > scale85[2]/2] + .25))
  #xlimitLog = c(min(pred.agg1$longitude[log(pred.agg1$ensMean + 1) > scale85log[2]/2] - 1),max(pred.agg1$longitude[log(pred.agg1$ensMean + 1) > scale85log[2]/2] + 1))
  #ylimitLog = c(min(pred.agg1$latitude[log(pred.agg1$ensMean + 1) > scale85log[2]/2] - .25), max(pred.agg1$latitude[log(pred.agg1$ensMean + 1) > scale85log[2]/2] + .25))
  #xlimit = c(min(xlimit1[1],xlimit2[1],xlimit3[1],xlimit4[1]),max(xlimit1[2],xlimit2[2],xlimit3[2],xlimit4[2]))
  #ylimit = c(min(ylimit1[1],ylimit2[1],ylimit3[1],ylimit4[1]),max(ylimit1[2],ylimit2[2],ylimit3[2],ylimit4[2]))
  
  # RCP85 plots_one for normal, one for Logged+1
  Projmap <- map('world', xlim=xlimit, ylim=ylimit, plot=F) 
  Projmap <- data.frame(lon=Projmap$x, lat=Projmap$y)
  print(levelplot(ensMean ~ longitude*latitude|year_range, data=pred.agg85, xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL,  
                     at = scale85, col.regions=cols, layout = c(5, 1)) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray'))
  
  dev.off()
}
 
# ==============================
# appendix
# ==============================

S1app <- data.frame(species = modeldiag$sppocean, region = biomass26$region, devPA = modeldiag$dev.pres, 
        devBiom = modeldiag$dev.biomass, uncertainty26 = centroids26$uncertainty,
        shift26 = centroids26$meanDist, sd26 = centroids26$sdDist, uncertainty85 = centroids85$uncertainty, 
        shift85 = centroids85$meanDist, sd85 = centroids85$sdDist, biomass26 = biomass26$biomMeanChan,
        biomsd26 = biomass26$sdChan, biomass85 = biomass85$biomMeanChan, biomsd85 = biomass85$sdChan)

save(S1app, file='Projections_S1app.RData')
write.csv(S1app, file='Projections_S1app.csv')

S1app <- data.frame(species = modeldiag$sppocean, region = biomass26$region, devPA = modeldiag$dev.pres, 
                    devBiom = modeldiag$dev.biomass, uncertainty = centroids26$uncertainty,
                    shift = centroids26$meanDist, sd = centroids26$sdDist, biomass = biomass26$biomMeanChan,
                    biomsd = biomass26$sdChan, stringsAsFactors=F)
S2app <- data.frame(species = modeldiag$sppocean, region = NA, devPA = NA, 
                    devBiom = NA, uncertainty = centroids85$uncertainty, 
                    shift = centroids85$meanDist, sd = centroids85$sdDist, biomass = biomass85$biomMeanChan, 
                    biomsd = biomass85$sdChan, stringsAsFactors=F)
# add a '2' to species names for S2app
S2app$species <- paste(S2app$species, '2', sep='')
Sapp <- rbind(S1app, S2app)
Sapp <- Sapp[order(Sapp$species),]
Sapp$species <- ifelse(grepl(2, Sapp$species), '', Sapp$species) 
Sapp$region <- ifelse(is.na(Sapp$region), '', Sapp$region) 
Sapp$devPA <- ifelse(is.na(Sapp$devPA), '', Sapp$devPA) 
Sapp$devBiom <- ifelse(is.na(Sapp$devBiom), '', Sapp$devBiom) 

save(Sapp, file='Projections_S1app.RData')
write.csv(Sapp, file='Projections_S1app.csv')




