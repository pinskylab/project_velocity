library(latticeExtra)
library(maps)
library(reshape2)
library(dplyr)

setwd('/Users/jamesmorley/Documents/project_velocity')
modfolder <- 'output/CEmodels_proj_PresAbs_Jun2018/'
rcp <- c(26, 45, 85)
modelrun <- c('bcc-csm1-1-m','bcc-csm1-1','CanESM2','CCSM4','CESM1-CAM5','CNRM-CM5','GFDL-CM3','GFDL-ESM2M','GFDL-ESM2G','GISS-E2-R','GISS-E2-H','HadGEM2-ES','IPSL-CM5A-LR','IPSL-CM5A-MR','MIROC-ESM','MIROC5','MPI-ESM-LR','NorESM1-ME')

 
# =======================================================================================
# First block produces a bunch of figures to take a look at familiar species results ====
# =======================================================================================

load('data/hauls_catch_Dec2017.RData') # hauls and catch data
load('data/speciesProjectionList.RData')
test <- sort(projspp[c(14,27,29,32,42,52,62,64,81,82,93,95,96,97,104,108,113,120,122,138,162,177,215,216,220,263,265,277,282,319,
                  325,326,327,315,316,351,358,361,367,352,370,372,395,396,432,433,434,454,475,481,515,517,
                  528,545,546,565,573,590,607,623,610,652,658,671)])
  
for(i in 1:length(test)){ 
  #for loop begin here
  species <- 'homarus americanus_Atl' # test[i]
  print(paste('Beginning  ', species, i))
  pdf(height=5, width=12, file=paste('figures/speciesProjections_PRES-ABS_Jun2018/', species, '_Pres-Abs_projections', '.pdf', sep=''))
   
  mydat <- dat[dat$sppocean==species,] 
  spdata <- merge(hauls, mydat, by='haulid', all.x = T, sort=F) # Add haul data and empty hauls
  myocean <- head(spdata$ocean[spdata$presfit == TRUE])[1] # identify if this is a Pacific or Atlantic species
  spdata <- spdata[spdata$ocean == myocean,] # trim master hauls file to the ocean of interest 
  spdata$presfit[is.na(spdata$presfit)] <- FALSE
  
  # make map of species occurences first
  plot(lat~lon, pch=20, cex=.2, spdata[spdata$presfit==T,])
  
      # rcp 26
      filename <- paste(modfolder, species, '_rcp', rcp[1], '_jas_PresAbs_projection_AGG.RData', sep='')
      load(filename)
      if(grepl('_Atl', species)){
        pred.agg26 <- as.matrix(pred.agg[,4:21], nrow=412960, ncol=18)
      }
      if(grepl('_Pac', species)){
        pred.agg26 <- as.matrix(pred.agg[,4:21], nrow=329130, ncol=18)
      }
      # aggregate by lon lat and year range
      ensMean26 <- apply(pred.agg26, 1, FUN=mean, na.rm=T)
      ensSD26 <- apply(pred.agg26, 1, FUN=sd, na.rm=T)
        
      # rcp 45
      filename <- paste(modfolder, species, '_rcp', rcp[2], '_jas_PresAbs_projection_AGG.RData', sep='')
      load(filename)
      if(grepl('_Atl', species)){
        pred.agg45 <- as.matrix(pred.agg[,4:21], nrow=412960, ncol=18)
      }
      if(grepl('_Pac', species)){
        pred.agg45 <- as.matrix(pred.agg[,4:21], nrow=329130, ncol=18)
      }
      # aggregate by lon lat and year range
      ensMean45 <- apply(pred.agg45, 1, FUN=mean, na.rm=T)
      ensSD45 <- apply(pred.agg45, 1, FUN=sd, na.rm=T)
      
      # rcp 85
      filename <- paste(modfolder, species, '_rcp', rcp[3], '_jas_PresAbs_projection_AGG.RData', sep='')
      load(filename)
      if(grepl('_Atl', species)){
        pred.agg85 <- as.matrix(pred.agg[,4:21], nrow=412960, ncol=18)
      }
      if(grepl('_Pac', species)){
        pred.agg85 <- as.matrix(pred.agg[,4:21], nrow=329130, ncol=18)
      }
      # aggregate by lon lat and year range
      ensMean85 <- apply(pred.agg85, 1, FUN=mean, na.rm=T)
      ensSD85 <- apply(pred.agg85, 1, FUN=sd, na.rm=T)
      
  pred_26_45_85 <- data.frame(cbind(pred.agg[,1:3], ensMean26=ensMean26, ensSD26=ensSD26, ensMean45=ensMean45, ensSD45=ensSD45, ensMean85=ensMean85, ensSD85=ensSD85))
  

xlimit = c(min(pred_26_45_85$longitude),max(pred_26_45_85$longitude))
ylimit = c(min(pred_26_45_85$latitude), max(pred_26_45_85$latitude))
scale85 = seq(0, max(pred_26_45_85[,c(4,6,8)]), length.out=20)
cols = colorRampPalette(colors = c('gray90', 'blue', 'dark blue', 'black'))

Projmap <- map('world', xlim=xlimit, ylim=ylimit, plot=F) 
Projmap <- data.frame(lon=Projmap$x, lat=Projmap$y)
print(levelplot(ensMean26 ~ longitude*latitude|year_range, data=pred_26_45_85, xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
                at = scale85, col.regions=cols, layout = c(5, 1)) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')) 
print(levelplot(ensMean45 ~ longitude*latitude|year_range, data=pred_26_45_85, xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
                at = scale85, col.regions=cols, layout = c(5, 1)) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')) 
print(levelplot(ensMean85 ~ longitude*latitude|year_range, data=pred_26_45_85, xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
                at = scale85, col.regions=cols, layout = c(5, 1)) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')) 
dev.off()
}


# Figure out what species need to be represented where
regions <- data.frame(West_coast=FALSE, Alaska=FALSE, Gulf_of_Mexico=FALSE, East_coast=FALSE)

for(i in 1:length(projspp)){
  species <- projspp[i]
  print(paste('Beginning  ', species, i))
  
  mydat <- dat[dat$sppocean==species,] 
  spdata <- merge(hauls, mydat, by='haulid', all.x = T, sort=F) # Add haul data and empty hauls
  myocean <- spdata$ocean[spdata$presfit == TRUE][1] # identify if this is a Pacific or Atlantic species
  spdata <- spdata[spdata$ocean == myocean,] # trim master hauls file to the ocean of interest 
  spdata$presfit[is.na(spdata$presfit)] <- FALSE
  
  spdata <- droplevels(spdata) # drop the west or east coast 'regionfact' levels
  myregions <- table(spdata$region[spdata$presfit==TRUE])  
  
  if(grepl('_Atl', species)){
    reg_numb <- length(intersect(names(myregions), c("SCDNR_SEUS","NEFSC_NEUS","VIMS_NEAMAP","DFO_ScotianShelf")))
    regions[i,] <- data.frame(West_coast=FALSE, Alaska=FALSE, Gulf_of_Mexico=any(names(myregions)=='SEFSC_GOMex'), East_coast=reg_numb > 0)
  }
  if(grepl('_Pac', species)){
    reg_numb <- length(intersect(names(myregions), c("AFSC_WCTri","NWFSC_WCAnn")))
    regions[i,] <- data.frame(West_coast=reg_numb > 0, Alaska=TRUE, Gulf_of_Mexico=FALSE, East_coast=FALSE)
  }
}
regions$species <- projspp
regions_mod <- melt(regions, id.vars='species', variable.name='region', value.name='include')
regions_mod <- regions_mod[!regions_mod$include==FALSE,]
regions_mod <- regions_mod[order(regions_mod$species),]
regions_mod$include <- NULL
regions_mod <- regions_mod[!regions_mod$species == 'squalus suckleyi_Pac',] # the 'squalus acanthias_Pac' is better model_these two species need to be merged and redone


# ==================================================================================================
# Calculate % change in habitat ====================================================================
# ==================================================================================================
# This is presently modified to only give results for one GCM =========================================
# ==================================================================================================

# load the U.S. EEZ reference vector
load('data/EEZ_grid_east.RData')
load('data/EEZ_grid_west.RData')

#HabChang_fullEns <- data.frame(Mean26_t1=numeric(), sd26_t1=numeric(), Mean26_t2=numeric(), sd26_t2=numeric(), Mean26_t3=numeric(), sd26_t3=numeric(), Mean26_t4=numeric(), sd26_t4=numeric(), 
 #                              Mean45_t1=numeric(), sd45_t1=numeric(), Mean45_t2=numeric(), sd45_t2=numeric(), Mean45_t3=numeric(), sd45_t3=numeric(), Mean45_t4=numeric(), sd45_t4=numeric(),
  #                             Mean85_t1=numeric(), sd85_t1=numeric(), Mean85_t2=numeric(), sd85_t2=numeric(), Mean85_t3=numeric(), sd85_t3=numeric(), Mean85_t4=numeric(), sd85_t4=numeric())

# Don't need the sd b/c only extracting on GCM
HabChang_fullEns <- data.frame(Mean26_t1=numeric(), Mean26_t2=numeric(), Mean26_t3=numeric(), Mean26_t4=numeric(), 
                               Mean45_t1=numeric(), Mean45_t2=numeric(), Mean45_t3=numeric(), Mean45_t4=numeric(), 
                               Mean85_t1=numeric(), Mean85_t2=numeric(), Mean85_t3=numeric(), Mean85_t4=numeric())

for(i in 1:nrow(regions_mod)){ 
  species <- regions_mod$species[i]
  
  print(paste('Beginning  ', species, ' in region:', regions_mod$region[i], i))
    
  # rcp 26
  filename <- paste(modfolder, species, '_rcp', rcp[1], '_jas_PresAbs_projection_AGG.RData', sep='')
  load(filename)
  if(grepl('_Atl', species)){
    pred.agg$EEZ <- EEZ_east
    pred.agg <- pred.agg[pred.agg$EEZ == TRUE,]
    if(regions_mod$region[i] == 'East_coast'){
      pred.agg <- pred.agg[!pred.agg$longitude < -82,]
      pred.agg26 <- pred.agg[!(pred.agg$longitude < -80.75 & pred.agg$latitude < 27),]
    }
    if(regions_mod$region[i] == 'Gulf_of_Mexico'){
      pred.agg <- pred.agg[!pred.agg$longitude > -80.75,]
      pred.agg26 <- pred.agg[!(pred.agg$longitude > -81.5 & pred.agg$latitude > 28),]
    }
  }
  if(grepl('_Pac', species)){
    pred.agg$EEZ <- EEZ_west
    pred.agg <- pred.agg[pred.agg$EEZ == TRUE,]
    if(regions_mod$region[i] == 'West_coast'){
      pred.agg26 <- pred.agg[!pred.agg$latitude > 50,]
    }
    if(regions_mod$region[i] == 'Alaska'){
      pred.agg26 <- pred.agg[!pred.agg$latitude < 50,]
    }
  }
  #Aggregate w/n year bins to get sum habitat suitability per time period
  #agg26 <- aggregate(pred.agg26[,4:21], by=list(year_range=pred.agg26$year_range), FUN=sum)
  #time_1 <- as.vector((agg26[2,2:19] - agg26[1,2:19])/agg26[1,2:19])
  #time_2 <- as.vector((agg26[3,2:19] - agg26[1,2:19])/agg26[1,2:19])
  #time_3 <- as.vector((agg26[4,2:19] - agg26[1,2:19])/agg26[1,2:19])
  #time_4 <- as.vector((agg26[5,2:19] - agg26[1,2:19])/agg26[1,2:19])
  
  # For the 5 GCMs the EPA want run this block instead of the preceeding 
  #agg26 <- aggregate(pred.agg26[,c(6,7,13,15,19)], by=list(year_range=pred.agg26$year_range), FUN=sum)
  #time_1 <- as.vector((agg26[2,2:6] - agg26[1,2:6])/agg26[1,2:6])
  #time_2 <- as.vector((agg26[3,2:6] - agg26[1,2:6])/agg26[1,2:6])
  #time_3 <- as.vector((agg26[4,2:6] - agg26[1,2:6])/agg26[1,2:6])
  #time_4 <- as.vector((agg26[5,2:6] - agg26[1,2:6])/agg26[1,2:6])
  
  # For the 1 GCM (GFDL-CM3) run this code instead of the preceeding
  agg26 <- aggregate(pred.agg26[,c(10)], by=list(year_range=pred.agg26$year_range), FUN=sum)
  time_1 <- as.vector((agg26[2,2] - agg26[1,2])/agg26[1,2])
  time_2 <- as.vector((agg26[3,2] - agg26[1,2])/agg26[1,2])
  time_3 <- as.vector((agg26[4,2] - agg26[1,2])/agg26[1,2])
  time_4 <- as.vector((agg26[5,2] - agg26[1,2])/agg26[1,2])
  
  agg26_chan <- rbind(time_1, time_2, time_3, time_4)
  
  # Don't need these running when only looking at one GCM
  #Mean26_t1 <- mean(as.numeric(agg26_chan[1,]))
  #sd26_t1 <- sd(as.numeric(agg26_chan[1,]))
  #Mean26_t2 <- mean(as.numeric(agg26_chan[2,]))
  #sd26_t2 <- sd(as.numeric(agg26_chan[2,]))
  #Mean26_t3 <- mean(as.numeric(agg26_chan[3,]))
  #sd26_t3 <- sd(as.numeric(agg26_chan[3,]))
  #Mean26_t4 <- mean(as.numeric(agg26_chan[4,]))
  #sd26_t4 <- sd(as.numeric(agg26_chan[4,]))
  
  
  # rcp 45
  filename <- paste(modfolder, species, '_rcp', rcp[2], '_jas_PresAbs_projection_AGG.RData', sep='')
  load(filename)
  if(grepl('_Atl', species)){
    pred.agg$EEZ <- EEZ_east
    pred.agg <- pred.agg[pred.agg$EEZ == TRUE,]
    if(regions_mod$region[i] == 'East_coast'){
      pred.agg <- pred.agg[!pred.agg$longitude < -82,]
      pred.agg45 <- pred.agg[!(pred.agg$longitude < -80.75 & pred.agg$latitude < 27),]
    }
    if(regions_mod$region[i] == 'Gulf_of_Mexico'){
      pred.agg <- pred.agg[!pred.agg$longitude > -80.75,]
      pred.agg45 <- pred.agg[!(pred.agg$longitude > -81.5 & pred.agg$latitude > 28),]
    }
  }
  if(grepl('_Pac', species)){
    pred.agg$EEZ <- EEZ_west
    pred.agg <- pred.agg[pred.agg$EEZ == TRUE,]
    if(regions_mod$region[i] == 'West_coast'){
      pred.agg45 <- pred.agg[!pred.agg$latitude > 50,]
    }
    if(regions_mod$region[i] == 'Alaska'){
      pred.agg45 <- pred.agg[!pred.agg$latitude < 50,]
    }
  }
  #Aggregate w/n year bins to get sum habitat suitability per time period
  #agg45 <- aggregate(pred.agg45[,4:21], by=list(year_range=pred.agg45$year_range), FUN=sum)
  #time_1 <- as.vector((agg45[2,2:19] - agg45[1,2:19])/agg45[1,2:19])
  #time_2 <- as.vector((agg45[3,2:19] - agg45[1,2:19])/agg45[1,2:19])
  #time_3 <- as.vector((agg45[4,2:19] - agg45[1,2:19])/agg45[1,2:19])
  #time_4 <- as.vector((agg45[5,2:19] - agg45[1,2:19])/agg45[1,2:19])
  
  # For the 5 GCMs the EPA want run this block instead of the preceeding 
  #agg45 <- aggregate(pred.agg45[,c(6,7,13,15,19)], by=list(year_range=pred.agg45$year_range), FUN=sum)
  #time_1 <- as.vector((agg45[2,2:6] - agg45[1,2:6])/agg45[1,2:6])
  #time_2 <- as.vector((agg45[3,2:6] - agg45[1,2:6])/agg45[1,2:6])
  #time_3 <- as.vector((agg45[4,2:6] - agg45[1,2:6])/agg45[1,2:6])
  #time_4 <- as.vector((agg45[5,2:6] - agg45[1,2:6])/agg45[1,2:6])
  
  # For the 1 GCM (GFDL-CM3) run this code instead of the preceeding
  agg45 <- aggregate(pred.agg45[,c(10)], by=list(year_range=pred.agg45$year_range), FUN=sum)
  time_1 <- as.vector((agg45[2,2] - agg45[1,2])/agg45[1,2])
  time_2 <- as.vector((agg45[3,2] - agg45[1,2])/agg45[1,2])
  time_3 <- as.vector((agg45[4,2] - agg45[1,2])/agg45[1,2])
  time_4 <- as.vector((agg45[5,2] - agg45[1,2])/agg45[1,2])
  
  agg45_chan <- rbind(time_1, time_2, time_3, time_4)
  
  #Mean45_t1 <- mean(as.numeric(agg45_chan[1,]))
  #sd45_t1 <- sd(as.numeric(agg45_chan[1,]))
  #Mean45_t2 <- mean(as.numeric(agg45_chan[2,]))
  #sd45_t2 <- sd(as.numeric(agg45_chan[2,]))
  #Mean45_t3 <- mean(as.numeric(agg45_chan[3,]))
  #sd45_t3 <- sd(as.numeric(agg45_chan[3,]))
  #Mean45_t4 <- mean(as.numeric(agg45_chan[4,]))
  #sd45_t4 <- sd(as.numeric(agg45_chan[4,]))
  

  # rcp 85
  filename <- paste(modfolder, species, '_rcp', rcp[3], '_jas_PresAbs_projection_AGG.RData', sep='')
  load(filename)
  if(grepl('_Atl', species)){
    pred.agg$EEZ <- EEZ_east
    pred.agg <- pred.agg[pred.agg$EEZ == TRUE,]
    if(regions_mod$region[i] == 'East_coast'){
      pred.agg <- pred.agg[!pred.agg$longitude < -82,]
      pred.agg85 <- pred.agg[!(pred.agg$longitude < -80.75 & pred.agg$latitude < 27),]
    }
    if(regions_mod$region[i] == 'Gulf_of_Mexico'){
      pred.agg <- pred.agg[!pred.agg$longitude > -80.75,]
      pred.agg85 <- pred.agg[!(pred.agg$longitude > -81.5 & pred.agg$latitude > 28),]
    }
  }
  if(grepl('_Pac', species)){
    pred.agg$EEZ <- EEZ_west
    pred.agg <- pred.agg[pred.agg$EEZ == TRUE,]
    if(regions_mod$region[i] == 'West_coast'){
      pred.agg85 <- pred.agg[!pred.agg$latitude > 50,]
      #pred.aggMod <- pred.agg[!pred.agg$latitude > 50,]
    }
    if(regions_mod$region[i] == 'Alaska'){
      pred.agg85 <- pred.agg[!pred.agg$latitude < 50,]
      pred.aggMod <- pred.agg[!pred.agg$latitude < 50,]
    }
  }
  
    #pred.agg85 <- as.matrix(pred.aggMod[,4:21], nrow=288755, ncol=18)
    # aggregate by lon lat and year range
    #ensMean85 <- apply(pred.agg85, 1, FUN=mean, na.rm=T)
    #pred.agg85 <- cbind(pred.aggMod[,1:3], ensMean85 = ensMean85)
  
    #xlimit = c(min(pred.agg85$longitude),max(pred.agg85$longitude))
    #ylimit = c(min(pred.agg85$latitude), max(pred.agg85$latitude))
    #scale85 = seq(0, max(pred.agg85$ensMean85), length.out=20)
    #cols = colorRampPalette(colors = c('gray90', 'blue', 'dark blue', 'black'))
    
    #Projmap <- map('world', xlim=xlimit, ylim=ylimit, plot=F) 
    #Projmap <- data.frame(lon=Projmap$x, lat=Projmap$y)
    #pdf(height=5, width=12, file=paste('figures/speciesProjections_PRES-ABS_Jun2018/', species, '_Pres-Abs_projections_DEMO2', '.pdf', sep=''))
    #print(levelplot(ensMean85 ~ longitude*latitude|year_range, data=pred.agg85, xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
     #       at = scale85, col.regions=cols, layout = c(5, 1)) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray'))
    #dev.off()
    
    
  #Aggregate w/n year bins to get sum habitat suitability per time period
  #agg85 <- aggregate(pred.agg85[,4:21], by=list(year_range=pred.agg85$year_range), FUN=sum)
  #time_1 <- as.vector((agg85[2,2:19] - agg85[1,2:19])/agg85[1,2:19])
  #time_2 <- as.vector((agg85[3,2:19] - agg85[1,2:19])/agg85[1,2:19])
  #time_3 <- as.vector((agg85[4,2:19] - agg85[1,2:19])/agg85[1,2:19])
  #time_4 <- as.vector((agg85[5,2:19] - agg85[1,2:19])/agg85[1,2:19])
  
  # For the 5 GCMs the EPA want run this block instead of the preceeding 
  #agg85 <- aggregate(pred.agg85[,c(6,7,13,15,19)], by=list(year_range=pred.agg85$year_range), FUN=sum)
  #time_1 <- as.vector((agg85[2,2:6] - agg85[1,2:6])/agg85[1,2:6])
  #time_2 <- as.vector((agg85[3,2:6] - agg85[1,2:6])/agg85[1,2:6])
  #time_3 <- as.vector((agg85[4,2:6] - agg85[1,2:6])/agg85[1,2:6])
  #time_4 <- as.vector((agg85[5,2:6] - agg85[1,2:6])/agg85[1,2:6])
  
  # For the 1 GCM (GFDL-CM3) run this code instead of the preceeding
  agg85 <- aggregate(pred.agg85[,c(10)], by=list(year_range=pred.agg85$year_range), FUN=sum)
  time_1 <- as.vector((agg85[2,2] - agg85[1,2])/agg85[1,2])
  time_2 <- as.vector((agg85[3,2] - agg85[1,2])/agg85[1,2])
  time_3 <- as.vector((agg85[4,2] - agg85[1,2])/agg85[1,2])
  time_4 <- as.vector((agg85[5,2] - agg85[1,2])/agg85[1,2])
  
  agg85_chan <- rbind(time_1, time_2, time_3, time_4)
   
  #Mean85_t1 <- mean(as.numeric(agg85_chan[1,]))
  #sd85_t1 <- sd(as.numeric(agg85_chan[1,]))
  #Mean85_t2 <- mean(as.numeric(agg85_chan[2,]))
  #sd85_t2 <- sd(as.numeric(agg85_chan[2,]))
  #Mean85_t3 <- mean(as.numeric(agg85_chan[3,]))
  #sd85_t3 <- sd(as.numeric(agg85_chan[3,]))
  #Mean85_t4 <- mean(as.numeric(agg85_chan[4,]))
  #sd85_t4 <- sd(as.numeric(agg85_chan[4,]))
  
  
  #HabChang_fullEns[i,] <- data.frame(Mean26_t1=Mean26_t1, sd26_t1=sd26_t1, Mean26_t2=Mean26_t2, sd26_t2=sd26_t2, Mean26_t3=Mean26_t3, sd26_t3=sd26_t3, Mean26_t4=Mean26_t4, sd26_t4=sd26_t4, 
   #                              Mean45_t1=Mean45_t1, sd45_t1=sd45_t1, Mean45_t2=Mean45_t2, sd45_t2=sd45_t2, Mean45_t3=Mean45_t3, sd45_t3=sd45_t3, Mean45_t4=Mean45_t4, sd45_t4=sd45_t4,
    #                             Mean85_t1=Mean85_t1, sd85_t1=sd85_t1, Mean85_t2=Mean85_t2, sd85_t2=sd85_t2, Mean85_t3=Mean85_t3, sd85_t3=sd85_t3, Mean85_t4=Mean85_t4, sd85_t4=sd85_t4)
  HabChang_fullEns[i,] <- data.frame(Mean26_t1=agg26_chan[1], Mean26_t2=agg26_chan[2], Mean26_t3=agg26_chan[3], Mean26_t4=agg26_chan[4], 
                                     Mean45_t1=agg45_chan[1], Mean45_t2=agg45_chan[2], Mean45_t3=agg45_chan[3], Mean45_t4=agg45_chan[4],
                                     Mean85_t1=agg85_chan[1], Mean85_t2=agg85_chan[2], Mean85_t3=agg85_chan[3], Mean85_t4=agg85_chan[4])
  
}
   
HabChang_fullEns$species <- regions_mod$species
HabChang_fullEns$region <- regions_mod$region  
HabChang_fullEns$species_Mod <- substr(HabChang_fullEns$species, 1, nchar(HabChang_fullEns$species)-4)

#save(HabChang_fullEns, regions_mod, file='output/Prop_habitat_change_Jun2018.RData')
#write.csv(HabChang_fullEns, file='output/Prop_habitat_change_Jun2018.csv')
#save(HabChang_fullEns, regions_mod, file='output/Prop_habitat_change_5GCM_Jun2018.RData')
#write.csv(HabChang_fullEns, file='output/Prop_habitat_change_5GCM_Jun2018.csv')
save(HabChang_fullEns, regions_mod, file='output/Prop_habitat_change_GFDLCM3_Mar2019.RData')
write.csv(HabChang_fullEns, file='output/Prop_habitat_change_GFDLCM3_Mar2019.csv')

#plot(log(sd26_t4)~log(abs(Mean26_t4)), HabChang_fullEns)

#plot(latitude~longitude, cex=.01, data=pred.agg[pred.agg$year_range=='2021-2040',])
#points(latitude~longitude, cex=.01, col='green', data=pred.agg[pred.agg$year_range=='2021-2040' & pred.agg$EEZ==TRUE,])
#points(latitude~longitude, cex=.01, col='green', data=pred.agg[pred.agg$year_range=='2021-2040' & pred.agg$latitude>50,])
