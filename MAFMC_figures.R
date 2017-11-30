# Figures for MAFMC meeting_can be adapted to other regions
 
setwd('/Users/jamesmorley/Documents/project_velocity')
library(lattice)
library(latticeExtra)
library(RColorBrewer)
library(maps)
library(Hmisc)

load('data/speciesProjectionList.RData')
doneSpp <- sort(projspp[c(419,600,105,310,523,517,270,438,587,465,525,526,392,17,79,135,196,26,28,498,307,164,275,37,66,126,136,151,163,170,188,264)])
projspp <- projspp[!(projspp %in% doneSpp)]
 
rcp <- c(26,85)
modelrun <- c('bcc-csm1-1-m','bcc-csm1-1','CanESM2','CCSM4','CESM1-CAM5','CNRM-CM5','GFDL-CM3','GFDL-ESM2M','GFDL-ESM2G','GISS-E2-R','GISS-E2-H','IPSL-CM5A-LR','IPSL-CM5A-MR','MIROC-ESM','MPI-ESM-LR','NorESM1-ME')
yearRange <- c('2007-2020', '2021-2040', '2041-2060', '2061-2080', '2081-2100')
 
# Figures for predicted biomass
for(k in 207:215){
  species <- projspp[k]
  print(paste("Beginning figures for ", species, sep=''))
  filename <- paste('output/CEmodels_proj_May2017/', species, '_rcp', rcp[1], '_jas_prediction_AGG.RData', sep='')
  load(filename)
  # restrict region to east coast, mainly dropping GMex, but also that bit north of Labrador
  #pred.agg <- pred.agg[!pred.agg$longitude < -82,]
  #pred.agg <- pred.agg[!pred.agg$latitude < 26.6,]
  #pred.agg <- pred.agg[!(pred.agg$latitude > 55 & pred.agg$longitude < -64),]
  pred.agg26 <- pred.agg
   
  filename <- paste('output/CEmodels_proj_May2017/', species, '_rcp', rcp[2], '_jas_prediction_AGG.RData', sep='')
  load(filename)
  # restrict region to east coast, mainly dropping GMex, but also that bit north of Labrador
  #pred.agg <- pred.agg[!pred.agg$longitude < -82,]
  #pred.agg <- pred.agg[!pred.agg$latitude < 26.6,]
  #pred.agg <- pred.agg[!(pred.agg$latitude > 55 & pred.agg$longitude < -64),]
  pred.agg85 <- pred.agg
  rm(pred.agg)
      
  # Look at the model-level predictions
  pdf(width=14, height=5, file=paste('figures/speciesProjections/', species, '_projections_LOGGED.pdf', sep=''))
    for(i in 4:19){
      Projmap <- map('world', xlim=c(min(pred.agg26$longitude - 1),max(pred.agg26$longitude + 1)), ylim=c(min(pred.agg26$latitude),max(pred.agg26$latitude)), plot=FALSE) 
      Projmap <- data.frame(lon=Projmap$x, lat=Projmap$y)
      cols = colorRampPalette(colors = c('gray90', 'blue', 'dark blue', 'black'))
      #pred.agg$logPred <- log(pred.agg$pred.mean + 1)
      print(levelplot(log(pred.agg26[,i] + 1) ~ longitude*latitude|year_range, data=pred.agg26, xlab=NULL, ylab=NULL, main=paste('rcp', rcp[1], 'summer', 'Model', modelrun[i-3], species, sep='_'), 
            at = seq(0,max(log(pred.agg26[,i] + 1)) + .00001, length.out=20), col.regions=cols, layout = c(5, 1)) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')) 
      print(levelplot(log(pred.agg85[,i] + 1) ~ longitude*latitude|year_range, data=pred.agg85, xlab=NULL, ylab=NULL, main=paste('rcp', rcp[2], 'summer', 'Model', modelrun[i-3], species, sep='_'), 
            at = seq(0,max(log(pred.agg85[,i] + 1)) + .00001, length.out=20), col.regions=cols, layout = c(5, 1)) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')) 
    }
  dev.off()
} 

 
# Figure to plot centroid changes
meanSD <- data.frame(meanShift=numeric(), sd=numeric())
pdf(width=11, height=6, file='figures/MAFMC/centroid shifts_B.pdf')
par(mfrow=c(1,2), mar=c(3,3,1,1))
for(k in 1:length(mafmc)){
  species <- mafmc[k]
  print(paste("Beginning figures for ", species, sep=''))
  
  filename <- paste('output/CEmodels_proj_May2017/', species, '_rcp', rcp[1], '_jas_prediction_AGG.RData', sep='')
  load(filename)
  # restrict region to east coast, mainly dropping GMex, but also that bit north of Labrador
  pred.agg <- pred.agg[!pred.agg$longitude < -82,]
  pred.agg <- pred.agg[!pred.agg$latitude < 26.6,]
  pred.agg <- pred.agg[!(pred.agg$latitude > 55 & pred.agg$longitude < -64),]
  pred.agg26 <- pred.agg
  cent26preds <- matrix(data=NA, nrow=5, ncol=19)
  for(i in 1:5){
    years = yearRange[i]
    preds <- pred.agg26[pred.agg26$year_range == years,]
    for(j in 4:19){
      value = wtd.mean(preds$latitude, weights = preds[,j]) 
      cent26preds[i,j] = value
    }
  }

  filename <- paste('output/CEmodels_proj_May2017/', species, '_rcp', rcp[2], '_jas_prediction_AGG.RData', sep='')
  load(filename)
  # restrict region to east coast, mainly dropping GMex, but also that bit north of Labrador
  pred.agg <- pred.agg[!pred.agg$longitude < -82,]
  pred.agg <- pred.agg[!pred.agg$latitude < 26.6,]
  pred.agg <- pred.agg[!(pred.agg$latitude > 55 & pred.agg$longitude < -64),]
  pred.agg85 <- pred.agg
  rm(pred.agg)
  cent85preds <- matrix(data=NA, nrow=5, ncol=19)
  for(i in 1:5){
    years = yearRange[i]
    preds <- pred.agg85[pred.agg85$year_range == years,]
    for(j in 4:19){
      value = wtd.mean(preds$latitude, weights = preds[,j]) 
      cent85preds[i,j] = value
    }
  }
  mean26 <- apply(cent26preds, 1, FUN=mean, na.rm=T)
  mean85 <- apply(cent85preds, 1, FUN=mean, na.rm=T)
  sd85 <- apply(cent85preds, 1, FUN=sd, na.rm=T)
  #shift85 <- apply(cent85preds[,4:19], 2, FUN=max, na.rm=T) - apply(cent85preds[,4:19], 2, FUN=min, na.rm=T)
  shift85 <- cent85preds[5,4:19] - cent85preds[1,4:19]
  mean85shift <- mean(t(shift85))
  sd85shift <- sd(t(shift85))
  
  cent26preds <- data.frame(cent26preds)
  cent85preds <- data.frame(cent85preds)
  cent26preds$mean <- mean26
  cent85preds$mean <- mean85
  cent26preds$year <- c(2010, 2030, 2050, 2070, 2090); cent85preds$year <- c(2010, 2030, 2050, 2070, 2090)
  
  plot(X4~year, col='gray', lwd=1, type='l', ylab='Latitude of centroid', xlab='', ylim=c(min(cent26preds[,4:19], cent85preds[,4:19])-.1, max(cent26preds[,4:19], cent85preds[,4:19])+.1), cent26preds)
  mtext(paste('  ', species, sep=''), side=3, adj=0, line=-1)
  mtext('  rcp26', side=3, adj=0, line=-2)
  points(X5~year, col='gray', lwd=1, type='l',cent26preds)
  points(X6~year, col='gray', lwd=1, type='l',cent26preds)
  points(X7~year, col='gray', lwd=1, type='l',cent26preds)
  points(X8~year, col='gray', lwd=1, type='l',cent26preds)
  points(X9~year, col='gray', lwd=1, type='l',cent26preds)
  points(X10~year, col='gray', lwd=1, type='l',cent26preds)
  points(X11~year, col='gray', lwd=1, type='l',cent26preds)
  points(X12~year, col='gray', lwd=1, type='l',cent26preds)
  points(X13~year, col='gray', lwd=1, type='l',cent26preds)
  points(X14~year, col='gray', lwd=1, type='l',cent26preds)
  points(X15~year, col='gray', lwd=1, type='l',cent26preds)
  points(X16~year, col='gray', lwd=1, type='l',cent26preds)
  points(X17~year, col='gray', lwd=1, type='l',cent26preds)
  points(X18~year, col='gray', lwd=1, type='l',cent26preds)
  points(X19~year, col='gray', lwd=1, type='l',cent26preds)
  points(mean~year, col='dark red', lwd=4, type='l', cent26preds)
  
  plot(X4~year, col='gray', lwd=1, type='l', sub='rcp85', ylab='', xlab='', ylim=c(min(cent26preds[,4:19], cent85preds[,4:19])-.1, max(cent26preds[,4:19], cent85preds[,4:19])+.1), cent85preds)
  mtext('  rcp85', side=3, adj=0, line=-2)
  points(X5~year, col='gray', lwd=1, type='l',cent85preds)
  points(X6~year, col='gray', lwd=1, type='l',cent85preds)
  points(X7~year, col='gray', lwd=1, type='l',cent85preds)
  points(X8~year, col='gray', lwd=1, type='l',cent85preds)
  points(X9~year, col='gray', lwd=1, type='l',cent85preds)
  points(X10~year, col='gray', lwd=1, type='l',cent85preds)
  points(X11~year, col='gray', lwd=1, type='l',cent85preds)
  points(X12~year, col='gray', lwd=1, type='l',cent85preds)
  points(X13~year, col='gray', lwd=1, type='l',cent85preds)
  points(X14~year, col='gray', lwd=1, type='l',cent85preds)
  points(X15~year, col='gray', lwd=1, type='l',cent85preds)
  points(X16~year, col='gray', lwd=1, type='l',cent85preds)
  points(X17~year, col='gray', lwd=1, type='l',cent85preds)
  points(X18~year, col='gray', lwd=1, type='l',cent85preds)
  points(X19~year, col='gray', lwd=1, type='l',cent85preds)
  points(mean~year, col='dark red', lwd=4, type='l', cent85preds)
  meanSD[k,] = data.frame(meanShift=mean85shift, sd=sd85shift)
}
dev.off()

meanSD$species <- mafmc
meanSD <- meanSD[-c(1,2,12),] # drop the few species with funky projections

library(calibrate)
meanSD$meanShift[meanSD$species == 'decapterus macarellus_Atl'] <- .3174
mod1 <- lm(sd~meanShift, meanSD)
 
pdf(width=11, height=6, file='figures/MAFMC/centroid_uncertainty.pdf')
par(mfrow=c(1,2))
plot(sd~meanShift, ylab="Uncertainty of shift (SD)", xlab="Mean magnitude of shift (lat)", meanSD)# took out monkfish and st.anchovy as projections were a bit fishy
abline(a=0.78293, b=0.26914, lty=2)
plot(sd~meanShift, ylab="Uncertainty of shift (SD)", xlab="Mean magnitude of shift (lat)", meanSD)# took out monkfish and st.anchovy as projections were a bit fishy
abline(a=0.78293, b=0.26914, lty=2)
textxy(meanSD$meanShift, meanSD$sd, labs=meanSD$species, offset=0, m = c(0, 0))
dev.off()
# ============================================================================================

pdf(width=11, height=6, file='figures/MAFMC/thermal_abundance.pdf')
par(mfrow=c(1,2), mar=c(3,3,1,1))
for(k in 1:length(mafmc)){
  species <- mafmc[k]
  print(paste("Beginning figures for ", species, sep=''))
  filename <- paste('output/CEmodels_proj_May2017/', species, '_rcp', rcp[1], '_jas_prediction_AGG.RData', sep='')
  load(filename)
  
  # restrict region to east coast, mainly dropping GMex, but also that bit north of Labrador
  pred.agg <- pred.agg[!pred.agg$longitude < -82,]
  pred.agg <- pred.agg[!pred.agg$latitude < 26.6,]
  pred.agg <- pred.agg[!(pred.agg$latitude > 55 & pred.agg$longitude < -64),]
  pred.agg$region <- ifelse(pred.agg$latitude < 35.25, "Southeast", 
              ifelse(pred.agg$latitude > 35.25 & pred.agg$longitude < -71.4, "Mid-Atlantic", 
                     ifelse(pred.agg$longitude > -71.4 & pred.agg$longitude < -66.3 & pred.agg$latitude < 47, "New England", NA)))
  pred.agg <- pred.agg[!is.na(pred.agg$region),]
  #coords <- unique(pred.agg[,c(2:3,20)])
  #plot(latitude~longitude, cex=.1,col="grey",coords)
  #points(latitude~longitude, cex=.2,col='red',coords[coords$region=='Southeast',])
  #points(latitude~longitude, cex=.2,col='green',coords[coords$region=='Mid-Atlantic',])
  #points(latitude~longitude, cex=.2,col='blue',coords[coords$region=='New England',])
  #points(latitude~longitude, cex=.2,col='black',coords[is.na(coords$region),])
  pred.agg26 <- pred.agg
  
  #Need to make three of these matrices, with different loops, for each region.......!
  hab26se <- matrix(data=NA, nrow=5, ncol=19)
  hab26mab <- matrix(data=NA, nrow=5, ncol=19)
  hab26ne <- matrix(data=NA, nrow=5, ncol=19)
  for(i in 1:5){
    years = yearRange[i]
    preds <- pred.agg26[pred.agg26$year_range == years,]
    for(j in 4:19){
      valueSE = sum(preds[preds$region == "Southeast",][,j]) 
      valueMAB = sum(preds[preds$region == "Mid-Atlantic",][,j]) 
      valueNE = sum(preds[preds$region == "New England",][,j]) 
      hab26se[i,j] = valueSE
      hab26mab[i,j] = valueMAB
      hab26ne[i,j] = valueNE
    }
  }
  mean26se <- apply(hab26se, 1, FUN=mean, na.rm=T)
  sd26se <- apply(hab26se, 1, FUN=sd, na.rm=T)
  mean26mab <- apply(hab26mab, 1, FUN=mean, na.rm=T)
  sd26mab <- apply(hab26mab, 1, FUN=sd, na.rm=T)
  mean26ne <- apply(hab26ne, 1, FUN=mean, na.rm=T)
  sd26ne <- apply(hab26ne, 1, FUN=sd, na.rm=T)
  
  preds26 <- data.frame(year = c(2010, 2030, 2050, 2070, 2090), meanSE=mean26se, sdSE=sd26se, meanMAB=mean26mab, sdMAB=sd26mab, meanNE=mean26ne, sdNE=sd26ne)

  
  filename <- paste('output/CEmodels_proj_May2017/', species, '_rcp', rcp[2], '_jas_prediction_AGG.RData', sep='')
  load(filename)
  # restrict region to east coast, mainly dropping GMex, but also that bit north of Labrador
  pred.agg <- pred.agg[!pred.agg$longitude < -82,]
  pred.agg <- pred.agg[!pred.agg$latitude < 26.6,]
  pred.agg <- pred.agg[!(pred.agg$latitude > 55 & pred.agg$longitude < -64),]
  pred.agg$region <- ifelse(pred.agg$latitude < 35.25, "Southeast", 
                            ifelse(pred.agg$latitude > 35.25 & pred.agg$longitude < -71.4, "Mid-Atlantic", 
                                   ifelse(pred.agg$longitude > -71.4 & pred.agg$longitude < -66.3 & pred.agg$latitude < 47, "New England", NA)))
  pred.agg <- pred.agg[!is.na(pred.agg$region),]
  pred.agg85 <- pred.agg
  rm(pred.agg)
   
  #Need to make three of these matrices, with different loops, for each region.......!
  hab85se <- matrix(data=NA, nrow=5, ncol=19)
  hab85mab <- matrix(data=NA, nrow=5, ncol=19)
  hab85ne <- matrix(data=NA, nrow=5, ncol=19)
  for(i in 1:5){
    years = yearRange[i]
    preds <- pred.agg85[pred.agg85$year_range == years,]
    for(j in 4:19){
      valueSE = sum(preds[preds$region == "Southeast",][,j]) 
      valueMAB = sum(preds[preds$region == "Mid-Atlantic",][,j]) 
      valueNE = sum(preds[preds$region == "New England",][,j]) 
      hab85se[i,j] = valueSE
      hab85mab[i,j] = valueMAB
      hab85ne[i,j] = valueNE
    }
  }
  mean85se <- apply(hab85se, 1, FUN=mean, na.rm=T)
  sd85se <- apply(hab85se, 1, FUN=sd, na.rm=T)
  mean85mab <- apply(hab85mab, 1, FUN=mean, na.rm=T)
  sd85mab <- apply(hab85mab, 1, FUN=sd, na.rm=T)
  mean85ne <- apply(hab85ne, 1, FUN=mean, na.rm=T)
  sd85ne <- apply(hab85ne, 1, FUN=sd, na.rm=T)
  
  preds85 <- data.frame(year = c(2010, 2030, 2050, 2070, 2090), meanSE=mean85se, sdSE=sd85se, meanMAB=mean85mab, sdMAB=sd85mab, meanNE=mean85ne, sdNE=sd85ne)
  
  # SHOULD I LOG THIS FIRST? TO DAMP DOWN THE UNREALISTIC HIGHS, MAYBE DO ONE LOGGED AND ONE W/O
  # ADD THE RCP26 PLOT_THEN DO THE LOOP
  plot(meanSE~year, type='l', lwd=4, col='red', ylab="", xlab="", ylim=c(min(preds26[,c(2,4,6)], preds85[,c(2,4,6)]), max(preds26[,c(2,4,6)], preds85[,c(2,4,6)])), preds26)
  errbar(x=preds26$year, y=preds26$meanSE, yplus=preds26$meanSE + preds26$sdSE, yminus=preds26$meanSE - preds26$sdSE, errbar.col='red', add=T)
  points(meanSE~year, col='red', preds26) 
  points(meanMAB~year, type='l', lwd=4, col='dark green', preds26)
  errbar(x=preds26$year, y=preds26$meanMAB, yplus=preds26$meanMAB + preds26$sdMAB, yminus=preds26$meanMAB - preds26$sdMAB, errbar.col='dark green', add=T)
  points(meanMAB~year, col='green', preds26)
  points(meanNE~year, type='l', lwd=4, col='blue', preds26)
  errbar(x=preds26$year, y=preds26$meanNE, yplus=preds26$meanNE + preds26$sdNE, yminus=preds26$meanNE - preds26$sdNE, errbar.col='blue', add=T)
  points(meanNE~year, col='blue', preds26) 
  mtext(paste('  ', species, sep=''), side=3, adj=0, line=-1)
  mtext('  rcp26', side=3, adj=0, line=-2)
  
  
  plot(meanSE~year, type='l', lwd=4, col='red', ylab="", xlab="", ylim=c(min(preds85[,c(2,4,6)], preds85[,c(2,4,6)]), max(preds85[,c(2,4,6)], preds85[,c(2,4,6)])), preds85)
  errbar(x=preds85$year, y=preds85$meanSE, yplus=preds85$meanSE + preds85$sdSE, yminus=preds85$meanSE - preds85$sdSE, errbar.col='red', add=T)
  points(meanSE~year, col='red', preds85) 
  points(meanMAB~year, type='l', lwd=4, col='dark green', preds85)
  errbar(x=preds85$year, y=preds85$meanMAB, yplus=preds85$meanMAB + preds85$sdMAB, yminus=preds85$meanMAB - preds85$sdMAB, errbar.col='dark green', add=T)
  points(meanMAB~year, col='green', preds85)
  points(meanNE~year, type='l', lwd=4, col='blue', preds85)
  errbar(x=preds85$year, y=preds85$meanNE, yplus=preds85$meanNE + preds85$sdNE, yminus=preds85$meanNE - preds85$sdNE, errbar.col='blue', add=T)
  points(meanNE~year, col='blue', preds85) 
  #mtext(paste('  ', species, sep=''), side=3, adj=0, line=-1)
  mtext('  rcp85', side=3, adj=0, line=-2)
  
  #meanSD[k,] = data.frame(meanShift=cent85preds$mean[5] - cent85preds$mean[1], sd=sd85[5])
}
dev.off()

# ===============================================================
# NOW DO THE SAME FIGURE AS JUST ABOVE, BUT LOG THE PREDICTIONS
# ===============================================================

meanSDhab <- data.frame(meanShift=numeric(), sd=numeric())

pdf(width=11, height=6, file='figures/MAFMC/thermal_abundance_LOGGED.pdf')
par(mfrow=c(1,2), mar=c(3,3,2,1))
for(k in 1:length(mafmc)){
  species <- mafmc[k]
  print(paste("Beginning figures for ", species, sep=''))
  filename <- paste('output/CEmodels_proj_May2017/', species, '_rcp', rcp[1], '_jas_prediction_AGG.RData', sep='')
  load(filename)
  
  # restrict region to east coast, mainly dropping GMex, but also that bit north of Labrador
  pred.agg <- pred.agg[!pred.agg$longitude < -82,]
  pred.agg <- pred.agg[!pred.agg$latitude < 26.6,]
  pred.agg <- pred.agg[!(pred.agg$latitude > 55 & pred.agg$longitude < -64),]
  pred.agg$region <- ifelse(pred.agg$latitude < 35.25, "Southeast", 
                            ifelse(pred.agg$latitude > 35.25 & pred.agg$longitude < -71.4, "Mid-Atlantic", 
                                   ifelse(pred.agg$longitude > -71.4 & pred.agg$longitude < -66.3 & pred.agg$latitude < 47, "New England", NA)))
  pred.agg <- pred.agg[!is.na(pred.agg$region),]
  #coords <- unique(pred.agg[,c(2:3,20)])
  #plot(latitude~longitude, cex=.1,col="grey",coords)
  #points(latitude~longitude, cex=.2,col='red',coords[coords$region=='Southeast',])
  #points(latitude~longitude, cex=.2,col='green',coords[coords$region=='Mid-Atlantic',])
  #points(latitude~longitude, cex=.2,col='blue',coords[coords$region=='New England',])
  #points(latitude~longitude, cex=.2,col='black',coords[is.na(coords$region),])
  pred.agg26 <- pred.agg
  
  #Need to make three of these matrices, with different loops, for each region.......!
  hab26se <- matrix(data=NA, nrow=5, ncol=19)
  hab26mab <- matrix(data=NA, nrow=5, ncol=19)
  hab26ne <- matrix(data=NA, nrow=5, ncol=19)
  for(i in 1:5){
    years = yearRange[i]
    preds <- pred.agg26[pred.agg26$year_range == years,]
    for(j in 4:19){
      valueSE = sum(log(preds[preds$region == "Southeast",][,j] + 1)) 
      valueMAB = sum(log(preds[preds$region == "Mid-Atlantic",][,j] + 1)) 
      valueNE = sum(log(preds[preds$region == "New England",][,j] + 1)) 
      hab26se[i,j] = valueSE
      hab26mab[i,j] = valueMAB
      hab26ne[i,j] = valueNE
    }
  }
  mean26se <- apply(hab26se, 1, FUN=mean, na.rm=T)
  sd26se <- apply(hab26se, 1, FUN=sd, na.rm=T)
  mean26mab <- apply(hab26mab, 1, FUN=mean, na.rm=T)
  sd26mab <- apply(hab26mab, 1, FUN=sd, na.rm=T)
  mean26ne <- apply(hab26ne, 1, FUN=mean, na.rm=T)
  sd26ne <- apply(hab26ne, 1, FUN=sd, na.rm=T)
  
  preds26 <- data.frame(year = c(2010, 2030, 2050, 2070, 2090), meanSE=mean26se, sdSE=sd26se, meanMAB=mean26mab, sdMAB=sd26mab, meanNE=mean26ne, sdNE=sd26ne)
  
  
  filename <- paste('output/CEmodels_proj_May2017/', species, '_rcp', rcp[2], '_jas_prediction_AGG.RData', sep='')
  load(filename)
  # restrict region to east coast, mainly dropping GMex, but also that bit north of Labrador
  pred.agg <- pred.agg[!pred.agg$longitude < -82,]
  pred.agg <- pred.agg[!pred.agg$latitude < 26.6,]
  pred.agg <- pred.agg[!(pred.agg$latitude > 55 & pred.agg$longitude < -64),]
  pred.agg$region <- ifelse(pred.agg$latitude < 35.25, "Southeast", 
                            ifelse(pred.agg$latitude > 35.25 & pred.agg$longitude < -71.4, "Mid-Atlantic", 
                                   ifelse(pred.agg$longitude > -71.4 & pred.agg$longitude < -66.3 & pred.agg$latitude < 47, "New England", NA)))
  pred.agg <- pred.agg[!is.na(pred.agg$region),]
  pred.agg85 <- pred.agg
  rm(pred.agg)
  
  #Need to make three of these matrices, with different loops, for each region.......!
  hab85se <- matrix(data=NA, nrow=5, ncol=19)
  hab85mab <- matrix(data=NA, nrow=5, ncol=19)
  hab85ne <- matrix(data=NA, nrow=5, ncol=19)
  for(i in 1:5){
    years = yearRange[i]
    preds <- pred.agg85[pred.agg85$year_range == years,]
    for(j in 4:19){
      valueSE = sum(log(preds[preds$region == "Southeast",][,j] + 1)) 
      valueMAB = sum(log(preds[preds$region == "Mid-Atlantic",][,j] + 1)) 
      valueNE = sum(log(preds[preds$region == "New England",][,j] + 1)) 
      hab85se[i,j] = valueSE
      hab85mab[i,j] = valueMAB
      hab85ne[i,j] = valueNE
    }
  }
  mean85se <- apply(hab85se, 1, FUN=mean, na.rm=T)
  sd85se <- apply(hab85se, 1, FUN=sd, na.rm=T)
  mean85mab <- apply(hab85mab, 1, FUN=mean, na.rm=T)
  sd85mab <- apply(hab85mab, 1, FUN=sd, na.rm=T)
  mean85ne <- apply(hab85ne, 1, FUN=mean, na.rm=T)
  sd85ne <- apply(hab85ne, 1, FUN=sd, na.rm=T)
  
  preds85 <- data.frame(year = c(2010, 2030, 2050, 2070, 2090), meanSE=mean85se, sdSE=sd85se, meanMAB=mean85mab, sdMAB=sd85mab, meanNE=mean85ne, sdNE=sd85ne)
  
  shift85 <- ((hab85mab[5,4:19] - hab85mab[1,4:19])/hab85mab[1,4:19]) * 100
  mean85shift <- mean(shift85)
  sd85shift <- sd(shift85) 
  
  plot(meanSE~year, type='l', lwd=4, col='red', ylab="", xlab="", ylim=c(min(preds26[,c(2,4,6)], preds85[,c(2,4,6)]), max(preds26[,c(2,4,6)], preds85[,c(2,4,6)])), preds26)
  errbar(x=preds26$year, y=preds26$meanSE, yplus=preds26$meanSE + preds26$sdSE, yminus=preds26$meanSE - preds26$sdSE, errbar.col='red', add=T)
  points(meanSE~year, col='red', preds26) 
  points(meanMAB~year, type='l', lwd=4, col='dark green', preds26)
  errbar(x=preds26$year, y=preds26$meanMAB, yplus=preds26$meanMAB + preds26$sdMAB, yminus=preds26$meanMAB - preds26$sdMAB, errbar.col='dark green', add=T)
  points(meanMAB~year, col='green', preds26)
  points(meanNE~year, type='l', lwd=4, col='blue', preds26)
  errbar(x=preds26$year, y=preds26$meanNE, yplus=preds26$meanNE + preds26$sdNE, yminus=preds26$meanNE - preds26$sdNE, errbar.col='blue', add=T)
  points(meanNE~year, col='blue', preds26) 
  mtext(species, side=3, adj=0, line=0.25)

  plot(meanSE~year, type='l', lwd=4, col='red', ylab="", xlab="", ylim=c(min(preds85[,c(2,4,6)], preds85[,c(2,4,6)]), max(preds85[,c(2,4,6)], preds85[,c(2,4,6)])), preds85)
  errbar(x=preds85$year, y=preds85$meanSE, yplus=preds85$meanSE + preds85$sdSE, yminus=preds85$meanSE - preds85$sdSE, errbar.col='red', add=T)
  points(meanSE~year, col='red', preds85) 
  points(meanMAB~year, type='l', lwd=4, col='dark green', preds85)
  errbar(x=preds85$year, y=preds85$meanMAB, yplus=preds85$meanMAB + preds85$sdMAB, yminus=preds85$meanMAB - preds85$sdMAB, errbar.col='dark green', add=T)
  points(meanMAB~year, col='green', preds85)
  points(meanNE~year, type='l', lwd=4, col='blue', preds85)
  errbar(x=preds85$year, y=preds85$meanNE, yplus=preds85$meanNE + preds85$sdNE, yminus=preds85$meanNE - preds85$sdNE, errbar.col='blue', add=T)
  points(meanNE~year, col='blue', preds85) 

  meanSDhab[k,] = data.frame(meanShift=mean85shift, sd=sd85shift)
}
dev.off()
 
meanSDhab$species <- mafmc
par(mfrow=c(1,1))
plot(sd~meanShift, ylab="Uncertainty of shift (coef of Var)", xlab="Mean percent magnitude of change (biomass)", meanSDhab)# took out monkfish and st.anchovy as projections were a bit fishy
textxy(meanSDhab$meanShift, meanSDhab$sd, labs=meanSDhab$species, cx = 0.5, dcol = "black", offset=0, m = c(0, 0))
habPlot <- meanSDhab[c(5,9,11,12,14,15,16,17,19,20,21,22,23),]
habPlot$names <- c('Black sea bass', 'Illex squid', 'Longfin squid', 'Monkfish', 'Summer flounder', 'Butterfish', 
          'Bluefish', 'Cobia', 'At. mackerel', 'King mackerel', 'Spanish mackerel', 'Spiny dogfish', 'Scup')
habPlot <- habPlot[order(habPlot$meanShift),]
habPlot$order <- c(1:nrow(habPlot))
par(mfrow=c(1,1))
library(ggplot2)

pdf(width=8, height=7, file='figures/MAFMC/shift_summary.pdf')
par(mfrow=c(1,1))

ggplot(habPlot, aes(x=order, y=meanShift)) + geom_point() + geom_hline(yintercept=0, col="dark gray", linetype="dashed") + scale_x_continuous(breaks=habPlot$order, labels=habPlot$names) +  scale_y_continuous(breaks=c(-100,0,100,200,300,400,500)) +
  theme(axis.text.x = element_text(angle=310, hjust=0, vjust=1, size=17), axis.text.y = element_text(size=15)) + geom_errorbar(aes(ymin=meanShift-sd, ymax=meanShift+sd), width=0) +
  geom_point(col="dark blue", size=6) 
dev.off()




  
  
  

