# This script gets the mean projections in order for all of the different niche models.
  # It produces 'centroid_final' and 'biomass_final' which are used in one of the main figures for the uncertainty paper.
  # The main figure is produced at the bottom of the script.
 
library(latticeExtra)
library(maps)
library(Hmisc)
 
setwd('/Users/jamesmorley/Documents/project_velocity')
projfolder <- 'output/CEmodels_Uncertainty_Proj_MEAN_2018/'
  
NICHE <- c('GLM', 'GAM', 'BRT')
RCP <- c(26,45,85)
SPECIES <- c('homarus americanus_Atl','paralichthys dentatus_Atl','centropristis striata_Atl','anoplopoma fimbria_Pac','doryteuthis opalescens_Pac','hippoglossus stenolepis_Pac','sebastes alutus_Pac')
modelrun <- c('bcc-csm1-1-m','bcc-csm1-1','CanESM2','CCSM4','CESM1-CAM5','CNRM-CM5','GFDL-CM3','GFDL-ESM2M','GFDL-ESM2G','GISS-E2-R','GISS-E2-H','HadGEM2-ES','IPSL-CM5A-LR','IPSL-CM5A-MR','MIROC-ESM','MIROC5','MPI-ESM-LR','NorESM1-ME')
YEARS <- c('2007-2020','2021-2040','2041-2060','2061-2080','2081-2100')
load('grid_area_correction_lm.RData')

centroid <- data.frame(species=character(), nicheMod=character(), rcp=character(), years=character(), centroid_biom=numeric(), centroid_PA=numeric(), sd_biom=numeric(), sd_PA=numeric(), stringsAsFactors = F)
centroid_final <- NULL
   
for(i in 1:length(SPECIES)){
  species = SPECIES[i]
  print(paste('Beginning', species))
  for(j in 1:3){ 
    rcp = RCP[j]
    pdf(width=9, height=11, file=paste('figures/Uncertainty_MS_2018/Projections/', species, '_rcp', rcp, '.pdf', sep=''))
    for(k in 1:3){
      niche = NICHE[k]
      load(paste(projfolder, species, '_rcp', rcp, '_jas_', niche, '_proj_AGG_2018.RData', sep=''))
      # ADD A MEAN PROJECTIONS FIGURE AS WELL_ONE FOR EACH SPECIES-RCP COMBINATION
      if(grepl('_Atl', species)){
        pred_means_PA <- as.matrix(pred.agg_PA[,4:21], nrow=412960, ncol=18)
        pred_means_BIOM <- as.matrix(pred.agg_biom[,4:21], nrow=412960, ncol=18)
        ensMean_PA <- apply(pred_means_PA, 1, FUN=mean, na.rm=T)
        ensMean_BIOM <- apply(pred_means_BIOM, 1, FUN=mean, na.rm=T)
      }
      if(grepl('_Pac', species)){
        pred_means_PA <- as.matrix(pred.agg_PA[,4:21], nrow=329130, ncol=18)
        pred_means_BIOM <- as.matrix(pred.agg_biom[,4:21], nrow=329130, ncol=18)
        ensMean_PA <- apply(pred_means_PA, 1, FUN=mean, na.rm=T)
        ensMean_BIOM <- apply(pred_means_BIOM, 1, FUN=mean, na.rm=T)
      }
       
      pred.agg_PA <- data.frame(cbind(pred.agg_PA, ensMean_PA=ensMean_PA))
      pred.agg_biom <- data.frame(cbind(pred.agg_biom, ensMean_BIOM=ensMean_BIOM))
      
      # ADD THE RELATIVE (OR ABSOLUTE) GRID AREA CORRECTIONS TO THE DATA HERE
      pred.agg_biom$area <- coef(area_mod)[1] + coef(area_mod)[2]*pred.agg_biom$latitude + coef(area_mod)[3]*(pred.agg_biom$latitude^2)
      pred.agg_PA$area <- coef(area_mod)[1] + coef(area_mod)[2]*pred.agg_PA$latitude + coef(area_mod)[3]*(pred.agg_PA$latitude^2)
      pred.agg_PA$area_rel <- pred.agg_PA$area/max(pred.agg_PA$area)
      
      # STANDARDIZE PLOT DIMENSIONS
      xlimit = c(min(pred.agg_PA$longitude) - 1, max(pred.agg_PA$longitude) + 1)
      ylimit = c(min(pred.agg_PA$latitude) - .25, max(pred.agg_PA$latitude) + .25)
      cols = colorRampPalette(colors = c('gray92', 'blue', 'dark blue', 'black'))
      Projmap <- map('world', xlim=xlimit, ylim=ylimit, plot=F) 
      Projmap <- data.frame(lon=Projmap$x, lat=Projmap$y)
       
      scale85 <- seq(0, max(pred.agg_PA$ensMean_PA) + .0001, length.out=20)
      plot1 <- levelplot(ensMean_PA ~ longitude*latitude|year_range, data=pred.agg_PA, index.cond=list(c(5,4,3,2,1)),  colorkey=F, par.settings=list(layout.heights=list(top.padding=-2, bottom.padding=-2),layout.widths=list(left.padding=-2, right.padding=-2.5)),
                  scales=list(tck=c(0,0), alternating=c(0,0)), strip=F, strip.left=TRUE,  aspect='fill', xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL,  
                  at = scale85, col.regions=cols, layout = c(1, 5)) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray',
                  par.settings=list(layout.heights=list(top.padding=0, bottom.padding=0)))
       
      scale85 <- seq(0, log(max(pred.agg_biom$ensMean_BIOM) + 1) + .01, length.out=20)
      plot2 <- levelplot(log(pred.agg_biom$ensMean_BIOM + 1) ~ longitude*latitude|year_range, data=pred.agg_biom, index.cond=list(c(5,4,3,2,1)),  colorkey=F, par.settings=list(layout.heights=list(top.padding=-2, bottom.padding=-2),layout.widths=list(left.padding=-2, right.padding=-2.5)),
                  scales=list(tck=c(0,0), alternating=c(0,0)), strip=F, strip.left=TRUE,  aspect='fill', xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL,  
                  at = scale85, col.regions=cols, layout = c(1, 5)) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray',
                  par.settings=list(layout.heights=list(top.padding=0, bottom.padding=0)))
      
      par(mar=c(.1,.1,.1,.1))
      print(plot1, split=c(1,1,2,1), more=TRUE) #position=c(0,0,.5,1), 
      print(plot2, split=c(2,1,2,1))#position=c(.5,0,1,1))
      
      # Need to carve out Gmex for Atlantic species
      if(grepl('_Atl', species)){
        pred.agg_PA <- pred.agg_PA[!pred.agg_PA$longitude < -82,]
        pred.agg_PA <- pred.agg_PA[!(pred.agg_PA$longitude < -80.75 & pred.agg_PA$latitude < 27),]
        pred.agg_biom <- pred.agg_biom[!pred.agg_biom$longitude < -82,]
        pred.agg_biom <- pred.agg_biom[!(pred.agg_biom$longitude < -80.75 & pred.agg_biom$latitude < 27),]
      }
      
      for(l in 1:5){
        years = YEARS[l]
        # CALCULATE SEPARATE CENTROID VALUES FOR EACH GCM AND YEAR RANGE_CALCULATE MEANS AND SD FROM THOSE
        predPA_sub <- pred.agg_PA[pred.agg_PA$year_range==years,]
        predBIOM_sub <- pred.agg_biom[pred.agg_biom$year_range==years,]
        
        PA_cents <- c(wtd.mean(predPA_sub$latitude, weights=(predPA_sub$mean1*predPA_sub$area_rel)), wtd.mean(predPA_sub$latitude, weights=(predPA_sub$mean2*predPA_sub$area_rel)),
                  wtd.mean(predPA_sub$latitude, weights=(predPA_sub$mean3*predPA_sub$area_rel)), wtd.mean(predPA_sub$latitude, weights=(predPA_sub$mean4*predPA_sub$area_rel)),
                  wtd.mean(predPA_sub$latitude, weights=(predPA_sub$mean5*predPA_sub$area_rel)), wtd.mean(predPA_sub$latitude, weights=(predPA_sub$mean6*predPA_sub$area_rel)),
                  wtd.mean(predPA_sub$latitude, weights=(predPA_sub$mean7*predPA_sub$area_rel)), wtd.mean(predPA_sub$latitude, weights=(predPA_sub$mean8*predPA_sub$area_rel)),
                  wtd.mean(predPA_sub$latitude, weights=(predPA_sub$mean9*predPA_sub$area_rel)), wtd.mean(predPA_sub$latitude, weights=(predPA_sub$mean10*predPA_sub$area_rel)),
                  wtd.mean(predPA_sub$latitude, weights=(predPA_sub$mean11*predPA_sub$area_rel)), wtd.mean(predPA_sub$latitude, weights=(predPA_sub$mean12*predPA_sub$area_rel)),
                  wtd.mean(predPA_sub$latitude, weights=(predPA_sub$mean13*predPA_sub$area_rel)), wtd.mean(predPA_sub$latitude, weights=(predPA_sub$mean14*predPA_sub$area_rel)),
                  wtd.mean(predPA_sub$latitude, weights=(predPA_sub$mean15*predPA_sub$area_rel)), wtd.mean(predPA_sub$latitude, weights=(predPA_sub$mean16*predPA_sub$area_rel)),
                  wtd.mean(predPA_sub$latitude, weights=(predPA_sub$mean17*predPA_sub$area_rel)), wtd.mean(predPA_sub$latitude, weights=(predPA_sub$mean18*predPA_sub$area_rel)))
        BIOM_cents <- c(wtd.mean(predBIOM_sub$latitude, weights=(predBIOM_sub$mean1*predBIOM_sub$area)), wtd.mean(predBIOM_sub$latitude, weights=(predBIOM_sub$mean2*predBIOM_sub$area)),
                  wtd.mean(predBIOM_sub$latitude, weights=(predBIOM_sub$mean3*predBIOM_sub$area)), wtd.mean(predBIOM_sub$latitude, weights=(predBIOM_sub$mean4*predBIOM_sub$area)),
                  wtd.mean(predBIOM_sub$latitude, weights=(predBIOM_sub$mean5*predBIOM_sub$area)), wtd.mean(predBIOM_sub$latitude, weights=(predBIOM_sub$mean6*predBIOM_sub$area)),
                  wtd.mean(predBIOM_sub$latitude, weights=(predBIOM_sub$mean7*predBIOM_sub$area)), wtd.mean(predBIOM_sub$latitude, weights=(predBIOM_sub$mean8*predBIOM_sub$area)),
                  wtd.mean(predBIOM_sub$latitude, weights=(predBIOM_sub$mean9*predBIOM_sub$area)), wtd.mean(predBIOM_sub$latitude, weights=(predBIOM_sub$mean10*predBIOM_sub$area)),
                  wtd.mean(predBIOM_sub$latitude, weights=(predBIOM_sub$mean11*predBIOM_sub$area)), wtd.mean(predBIOM_sub$latitude, weights=(predBIOM_sub$mean12*predBIOM_sub$area)),
                  wtd.mean(predBIOM_sub$latitude, weights=(predBIOM_sub$mean13*predBIOM_sub$area)), wtd.mean(predBIOM_sub$latitude, weights=(predBIOM_sub$mean14*predBIOM_sub$area)),
                  wtd.mean(predBIOM_sub$latitude, weights=(predBIOM_sub$mean15*predBIOM_sub$area)), wtd.mean(predBIOM_sub$latitude, weights=(predBIOM_sub$mean16*predBIOM_sub$area)),
                  wtd.mean(predBIOM_sub$latitude, weights=(predBIOM_sub$mean17*predBIOM_sub$area)), wtd.mean(predBIOM_sub$latitude, weights=(predBIOM_sub$mean18*predBIOM_sub$area)))
         
        centroid[l,] <- data.frame(species=species, nicheMod=niche, rcp=rcp, years=years, centroid_biom=mean(BIOM_cents), centroid_PA=mean(PA_cents), sd_biom=sd(BIOM_cents), sd_PA=sd(PA_cents), stringsAsFactors = F)
      }
      centroid_final <- rbind(centroid_final, centroid)
    }
    dev.off()  
  }
}

centroid_final$xAxis <- 1 
centroid_final$xAxis[centroid_final$years=='2021-2040'] <- 2
centroid_final$xAxis[centroid_final$years=='2041-2060'] <- 3
centroid_final$xAxis[centroid_final$years=='2061-2080'] <- 4
centroid_final$xAxis[centroid_final$years=='2081-2100'] <- 5
centroid_final$xAxis_plus <- centroid_final$xAxis + .1 
centroid_final$xAxis_min <- centroid_final$xAxis - .1 

save(centroid_final, file='output/Uncertainty_centroidMEAN_summaries_Dec2018.RData')
#load('output/Uncertainty_centroidMEAN_summaries_Dec2018.RData')

# =========================================================================================================
# Biomass by region ==========================================================================================
# =========================================================================================================

load('data/EEZ_grid_east.RData') # load the EEZ for U.S.
load('data/EEZ_grid_west.RData') # load the EEZ for U.S.
YEARS <- c('2021-2040','2041-2060','2061-2080','2081-2100') # drop first time period, as everything is compared to that (i.e. change not absolute)

biomass <- data.frame(species=character(), nicheMod=character(), rcp=character(), years=character(), north_biom=numeric(), north_PA=numeric(), south_biom=numeric(), south_PA=numeric(), sdNorth_biom=numeric(), sdNorth_PA=numeric(), sdSouth_biom=numeric(), sdSouth_PA=numeric(), stringsAsFactors = F)
biomass_final <- NULL

for(i in 1:length(SPECIES)){
  species = SPECIES[i]
  print(paste('Beginning', species))
  for(j in 1:3){ 
    rcp = RCP[j]
    for(k in 1:3){
      niche = NICHE[k]
      load(paste(projfolder, species, '_rcp', rcp, '_jas_', niche, '_proj_AGG_2018.RData', sep=''))
      
      # ADD THE RELATIVE (OR ABSOLUTE) GRID AREA CORRECTIONS TO THE DATA HERE
      pred.agg_biom$area <- coef(area_mod)[1] + coef(area_mod)[2]*pred.agg_biom$latitude + coef(area_mod)[3]*(pred.agg_biom$latitude^2)
      pred.agg_PA$area <- coef(area_mod)[1] + coef(area_mod)[2]*pred.agg_PA$latitude + coef(area_mod)[3]*(pred.agg_PA$latitude^2)
      pred.agg_PA$area_rel <- pred.agg_PA$area/max(pred.agg_PA$area)
      
      #Add EEZ identifier
      if(grepl('_Atl', species)){
        pred.agg_PA$EEZ <- EEZ_east
        pred.agg_biom$EEZ <- EEZ_east
      }
      if(grepl('_Pac', species)){
        pred.agg_PA$EEZ <- EEZ_west
        pred.agg_biom$EEZ <- EEZ_west
      }
      
      # Carve out Gmexico
      if(grepl('_Atl', species)){
        pred.agg_PA <- pred.agg_PA[!pred.agg_PA$longitude < -82,]
        pred.agg_PA <- pred.agg_PA[!(pred.agg_PA$longitude < -80.75 & pred.agg_PA$latitude < 27),]
        pred.agg_biom <- pred.agg_biom[!pred.agg_biom$longitude < -82,]
        pred.agg_biom <- pred.agg_biom[!(pred.agg_biom$longitude < -80.75 & pred.agg_biom$latitude < 27),]
      }
       
      #ADD THE NORTH SOUTH DESIGNATIONS
      if(species == 'paralichthys dentatus_Atl' | species == 'centropristis striata_Atl'){ # For a Cape Hatteras split (BSB)
        pred.agg_PA$region <- ifelse(pred.agg_PA$latitude > 35.25, 'North', 'South')
        pred.agg_PA$region[pred.agg_PA$EEZ == FALSE] <- NA # drop non U.S. waters
        pred.agg_biom$region <- ifelse(pred.agg_biom$latitude > 35.25, 'North', 'South')
        pred.agg_biom$region[pred.agg_biom$EEZ == FALSE] <- NA # drop non U.S. waters
      }
      if(species == 'homarus americanus_Atl'){ # For a Gulf of Maine/Georges Bank versus Southern NE/MAB split (Lobster)
        pred.agg_PA$region <- ifelse(pred.agg_PA$longitude > -70, 'North', 'South')
        pred.agg_PA$region[pred.agg_PA$region == 'South' & pred.agg_PA$latitude > 41.75] <- 'North'
        pred.agg_PA$region[pred.agg_PA$latitude < 35.25] <- NA # drop south of Cape Hatteras
        pred.agg_PA$region[pred.agg_PA$EEZ == FALSE] <- NA # drop non U.S. waters
        
        pred.agg_biom$region <- ifelse(pred.agg_biom$longitude > -70, 'North', 'South')
        pred.agg_biom$region[pred.agg_biom$region == 'South' & pred.agg_biom$latitude > 41.75] <- 'North'
        pred.agg_biom$region[pred.agg_biom$latitude < 35.25] <- NA # drop south of Cape Hatteras
        pred.agg_biom$region[pred.agg_biom$EEZ == FALSE] <- NA # drop non U.S. waters
      }
      if(species == 'anoplopoma fimbria_Pac' | species == 'hippoglossus stenolepis_Pac' | species == 'sebastes alutus_Pac'){ # For an Alaska-U.S. split, won't include Canada for habitat analysis
        pred.agg_PA$region <- ifelse(pred.agg_PA$latitude > 50, 'North', 'South')
        pred.agg_PA$region[pred.agg_PA$EEZ == FALSE] <- NA # drop non U.S. waters
        pred.agg_biom$region <- ifelse(pred.agg_biom$latitude > 50, 'North', 'South')
        pred.agg_biom$region[pred.agg_biom$EEZ == FALSE] <- NA # drop non U.S. waters
      }
      if(species == 'doryteuthis opalescens_Pac'){ # For U.S. west coast vs. West Canada split
        pred.agg_PA$region <- ifelse(pred.agg_PA$longitude > -130 & pred.agg_PA$EEZ == TRUE, 'South', NA)
        pred.agg_PA$region[pred.agg_PA$latitude > 45 & pred.agg_PA$latitude < 58 & pred.agg_PA$EEZ == FALSE] <- 'North'
        pred.agg_biom$region <- ifelse(pred.agg_biom$longitude > -130 & pred.agg_biom$EEZ == TRUE, 'South', NA)
        pred.agg_biom$region[pred.agg_biom$latitude > 45 & pred.agg_biom$latitude < 58 & pred.agg_biom$EEZ == FALSE] <- 'North'
      }
      
      # CALCULATE total habitat for the baseline time period
      predPA_sub <- pred.agg_PA[pred.agg_PA$year_range=='2007-2020',]
      predBIOM_sub <- pred.agg_biom[pred.agg_biom$year_range=='2007-2020',]
      
      # Calculate habitat area for each region, using the grid area scalar correction
      south_biom <- predBIOM_sub[!is.na(predBIOM_sub$region) & predBIOM_sub$region == 'South',]
      south_biomMat <- as.matrix(south_biom[,c(4:21)], nrow=nrow(south_biom), ncol=18)  
      south_biomMat <- south_biomMat*south_biom$area
      south_biom <- apply(south_biomMat, 2, FUN=sum, na.rm=T)
        
      north_biom <- predBIOM_sub[!is.na(predBIOM_sub$region) & predBIOM_sub$region == 'North',]
      north_biomMat <- as.matrix(north_biom[,c(4:21)], nrow=nrow(north_biom), ncol=18)  
      north_biomMat <- north_biomMat*north_biom$area
      north_biom <- apply(north_biomMat, 2, FUN=sum, na.rm=T)
        
      south_biomPA <- predPA_sub[!is.na(predPA_sub$region) & predPA_sub$region == 'South',]
      south_biomPAMat <- as.matrix(south_biomPA[,c(4:21)], nrow=nrow(south_biomPA), ncol=18)  
      south_biomPAMat <- south_biomPAMat*south_biomPA$area
      south_biomPA <- apply(south_biomPAMat, 2, FUN=sum, na.rm=T)
        
      north_biomPA <- predPA_sub[!is.na(predPA_sub$region) & predPA_sub$region == 'North',]
      north_biomPAMat <- as.matrix(north_biomPA[,c(4:21)], nrow=nrow(north_biomPA), ncol=18)  
      north_biomPAMat <- north_biomPAMat*north_biomPA$area
      north_biomPA <- apply(north_biomPAMat, 2, FUN=sum, na.rm=T)
      
      for(l in 1:4){
        years = YEARS[l]
        predPA_sub <- pred.agg_PA[pred.agg_PA$year_range==years,]
        predBIOM_sub <- pred.agg_biom[pred.agg_biom$year_range==years,]
        
        south_biomB <- predBIOM_sub[!is.na(predBIOM_sub$region) & predBIOM_sub$region == 'South',]
        south_biomMatB <- as.matrix(south_biomB[,c(4:21)], nrow=nrow(south_biomB), ncol=18)  
        south_biomMatB <- south_biomMatB*south_biomB$area
        south_biomB <- apply(south_biomMatB, 2, FUN=sum, na.rm=T)
        south_biomB <- south_biomB - south_biom
        south_biomB <- (south_biomB/south_biom)*100
        
        north_biomB <- predBIOM_sub[!is.na(predBIOM_sub$region) & predBIOM_sub$region == 'North',]
        north_biomMatB <- as.matrix(north_biomB[,c(4:21)], nrow=nrow(north_biomB), ncol=18)  
        north_biomMatB <- north_biomMatB*north_biomB$area
        north_biomB <- apply(north_biomMatB, 2, FUN=sum, na.rm=T)
        north_biomB <- north_biomB - north_biom
        north_biomB <- (north_biomB/north_biom)*100
        
        south_biomPAB <- predPA_sub[!is.na(predPA_sub$region) & predPA_sub$region == 'South',]
        south_biomPAMatB <- as.matrix(south_biomPAB[,c(4:21)], nrow=nrow(south_biomPAB), ncol=18)  
        south_biomPAMatB <- south_biomPAMatB*south_biomPAB$area
        south_biomPAB <- apply(south_biomPAMatB, 2, FUN=sum, na.rm=T)
        south_biomPAB <- south_biomPAB - south_biomPA
        south_biomPAB <- (south_biomPAB/south_biomPA)*100
        
        north_biomPAB <- predPA_sub[!is.na(predPA_sub$region) & predPA_sub$region == 'North',]
        north_biomPAMatB <- as.matrix(north_biomPAB[,c(4:21)], nrow=nrow(north_biomPAB), ncol=18)  
        north_biomPAMatB <- north_biomPAMatB*north_biomPAB$area
        north_biomPAB <- apply(north_biomPAMatB, 2, FUN=sum, na.rm=T)
        north_biomPAB <- north_biomPAB - north_biomPA
        north_biomPAB <- (north_biomPAB/north_biomPA)*100
        
        biomass[l,] <- data.frame(species=species, nicheMod=niche, rcp=rcp, years=years, north_biom=mean(north_biomB), north_PA=mean(north_biomPAB), south_biom=mean(south_biomB), south_PA=mean(south_biomPAB), sdNorth_biom=sd(north_biomB), sdNorth_PA=sd(north_biomPAB), sdSouth_biom=sd(south_biomB), sdSouth_PA=sd(south_biomPAB), stringsAsFactors = F)
      }
      biomass_final <- rbind(biomass_final, biomass)
    }
  }
}

biomass_final$xAxis <- 1 
biomass_final$xAxis[biomass_final$years=='2041-2060'] <- 2
biomass_final$xAxis[biomass_final$years=='2061-2080'] <- 3
biomass_final$xAxis[biomass_final$years=='2081-2100'] <- 4
biomass_final$xAxis_plus <- biomass_final$xAxis + .1 
biomass_final$xAxis_min <- biomass_final$xAxis - .1 

save(biomass_final, file='output/Uncertainty_habitatMEAN_summaries_Dec2018.RData')
#load('output/Uncertainty_habitatMEAN_summaries_Dec2018.RData')






# =========================================================================================================
# Plot up summaries ============================================================================================
# =========================================================================================================

# First load centroid_final and biomass_final from code above_the rest should run as is.

centroid_final85 <- centroid_final[centroid_final$rcp == 85,] 
biomass_final85 <- biomass_final[biomass_final$rcp == 85,] 
# Read in fish graphics
library(png)
#rock <- readPNG('Fish_graphics/widow.png')

pdf(width=6.5, height=8, file='figures/Uncertainty_MS_2018/NicheMod_comparison_Dec2018.pdf')
#plot.new()
plot_matrix <- rbind(c(16:22),c(23,1,24,5,25,9,26),c(13,2,14,6,15,10,27),c(13,3,14,7,15,11,28),c(29,4,30,8,31,12,32),c(33:39))
layout(mat=matrix(plot_matrix, nrow=6, ncol=7), widths=c(1.2,4,1.2,4,.6,4,.5), heights=c(0.1,4,4,4,4,1.6))# 1.2, Need to adjust this to fit the three figures
#layout.show(n=15)
 
# HALIBUT centroid
i=6
par(mar=c(0,0,0,0), oma=c(0,0,0,0), tcl=-.25, mgp=c(2.25,0.5,0), cex.lab=1.2, xpd = NA)
  plot(centroid_biom~xAxis_plus, xlim=c(0.9,5.1), ylim=c(56.5,58.25), xaxt='n', yaxt='n', xlab='', ylab='Centroid (latitude)', type='l', lwd=2, col='dodgerblue2', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM',])
  errbar(x=c(1.05,2.05,3.05,4.05,5.05)+0.05, y=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'], yplus=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'] + centroid_final85$sd_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'], yminus=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'] - centroid_final85$sd_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'], errbar.col='dodgerblue2', lwd=2, cap=0, add=T)
  points(centroid_biom~xAxis_plus, type='p', pch=19, col='dodgerblue2', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM',])
  points(centroid_PA~xAxis_plus, type='l', lwd=2, lty=2, col='dodgerblue2', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM',])
  errbar(x=c(1.05,2.05,3.05,4.05,5.05)+0.05, y=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'], yplus=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'] + centroid_final85$sd_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'], yminus=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'] - centroid_final85$sd_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'], errbar.col='dodgerblue2', lwd=2, cap=0, add=T)
  points(centroid_PA~xAxis_plus, type='p', pch=19, col='dodgerblue2', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM',])
  axis(2,at=c(56.5,57,57.5,58), lwd=1.3, cex.axis=0.95, las=1)
  mtext('(a)',line=-1.25,adj=0.01,font=2)
  
  points(centroid_biom~xAxis, type='l', lwd=2, col='gray20', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM',])
  errbar(x=c(1:5), y=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'], yplus=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'] + centroid_final85$sd_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'], yminus=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'] - centroid_final85$sd_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'], errbar.col='gray20', lwd=2, cap=0, add=T)
  points(centroid_biom~xAxis, type='p', pch=19, col='gray20', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM',])
  points(centroid_PA~xAxis, type='l', lwd=2, lty=2, col='gray20', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM',])
  errbar(x=c(1:5), y=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'], yplus=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'] + centroid_final85$sd_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'], yminus=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'] - centroid_final85$sd_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'], errbar.col='gray20', lwd=2, cap=0, add=T)
  points(centroid_PA~xAxis, type='p', pch=19, col='gray20', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM',])

  points(centroid_biom~xAxis_min, type='l', lwd=2, col='red', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT',])
  errbar(x=c(0.95,1.95,2.95,3.95,4.95)-0.05, y=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'], yplus=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'] + centroid_final85$sd_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'], yminus=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'] - centroid_final85$sd_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'], errbar.col='red', lwd=2, cap=0, add=T)
  points(centroid_biom~xAxis_min, type='p', pch=19, col='red', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT',])
  points(centroid_PA~xAxis_min, type='l', lwd=2, lty=2, col='red', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT',])
  errbar(x=c(0.95,1.95,2.95,3.95,4.95)-0.05, y=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'], yplus=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'] + centroid_final85$sd_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'], yminus=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'] - centroid_final85$sd_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'], errbar.col='red', lwd=2, cap=0, add=T)
  points(centroid_PA~xAxis_min, type='p', pch=19, col='red', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT',])
  #rasterImage(rock, .5,57.5,2,58) # placeholder for adding the fish graphics
  
# PERCH centroid
i=7
  plot(centroid_biom~xAxis_plus, xlim=c(0.9,5.1), ylim=c(52.25,56.75), xaxt='n', yaxt='n', xlab='', ylab='Centroid (latitude)', type='l', lwd=2, col='dodgerblue2', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM',])
  errbar(x=c(1.05,2.05,3.05,4.05,5.05)+0.05, y=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'], yplus=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'] + centroid_final85$sd_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'], yminus=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'] - centroid_final85$sd_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'], errbar.col='dodgerblue2', lwd=2, cap=0, add=T)
  points(centroid_biom~xAxis_plus, type='p', pch=19, col='dodgerblue2', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM',])
  points(centroid_PA~xAxis_plus, type='l', lwd=2, lty=2, col='dodgerblue2', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM',])
  errbar(x=c(1.05,2.05,3.05,4.05,5.05)+0.05, y=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'], yplus=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'] + centroid_final85$sd_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'], yminus=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'] - centroid_final85$sd_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'], errbar.col='dodgerblue2', lwd=2, cap=0, add=T)
  points(centroid_PA~xAxis_plus, type='p', pch=19, col='dodgerblue2', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM',])
  axis(2,at=c(53,54,55,56), lwd=1.3, cex.axis=0.95, las=1)
  mtext('(b)',line=-1.25,adj=0.01,font=2)
  
  points(centroid_biom~xAxis, type='l', lwd=2, col='gray20', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM',])
  errbar(x=c(1:5), y=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'], yplus=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'] + centroid_final85$sd_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'], yminus=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'] - centroid_final85$sd_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'], errbar.col='gray20', lwd=2, cap=0, add=T)
  points(centroid_biom~xAxis, type='p', pch=19, col='gray20', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM',])
  points(centroid_PA~xAxis, type='l', lwd=2, lty=2, col='gray20', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM',])
  errbar(x=c(1:5), y=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'], yplus=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'] + centroid_final85$sd_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'], yminus=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'] - centroid_final85$sd_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'], errbar.col='gray20', lwd=2, cap=0, add=T)
  points(centroid_PA~xAxis, type='p', pch=19, col='gray20', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM',])
  
  points(centroid_biom~xAxis_min, type='l', lwd=2, col='red', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT',])
  errbar(x=c(0.95,1.95,2.95,3.95,4.95)-0.05, y=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'], yplus=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'] + centroid_final85$sd_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'], yminus=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'] - centroid_final85$sd_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'], errbar.col='red', lwd=2, cap=0, add=T)
  points(centroid_biom~xAxis_min, type='p', pch=19, col='red', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT',])
  points(centroid_PA~xAxis_min, type='l', lwd=2, lty=2, col='red', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT',])
  errbar(x=c(0.95,1.95,2.95,3.95,4.95)-0.05, y=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'], yplus=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'] + centroid_final85$sd_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'], yminus=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'] - centroid_final85$sd_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'], errbar.col='red', lwd=2, cap=0, add=T)
  points(centroid_PA~xAxis_min, type='p', pch=19, col='red', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT',])

# SUMMER FLOUNDER centroid
i=2 # 
  plot(centroid_biom~xAxis_plus, xlim=c(0.9,5.1), ylim=c(36,43.75), xaxt='n', yaxt='n', xlab='', ylab='Centroid (latitude)', type='l', lwd=2, col='dodgerblue2', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM',])
  errbar(x=c(1.05,2.05,3.05,4.05,5.05)+0.05, y=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'], yplus=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'] + centroid_final85$sd_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'], yminus=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'] - centroid_final85$sd_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'], errbar.col='dodgerblue2', lwd=2, cap=0, add=T)
  points(centroid_biom~xAxis_plus, type='p', pch=19, col='dodgerblue2', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM',])
  points(centroid_PA~xAxis_plus, type='l', lwd=2, lty=2, col='dodgerblue2', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM',])
  errbar(x=c(1.05,2.05,3.05,4.05,5.05)+0.05, y=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'], yplus=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'] + centroid_final85$sd_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'], yminus=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'] - centroid_final85$sd_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'], errbar.col='dodgerblue2', lwd=2, cap=0, add=T)
  points(centroid_PA~xAxis_plus, type='p', pch=19, col='dodgerblue2', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM',])
  axis(2,at=c(37,39,41,43), lwd=1.3, cex.axis=0.95, las=1)
  mtext('(c)',line=-1.25,adj=0.01,font=2)
  
  points(centroid_biom~xAxis, type='l', lwd=2, col='gray20', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM',])
  errbar(x=c(1:5), y=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'], yplus=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'] + centroid_final85$sd_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'], yminus=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'] - centroid_final85$sd_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'], errbar.col='gray20', lwd=2, cap=0, add=T)
  points(centroid_biom~xAxis, type='p', pch=19, col='gray20', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM',])
  points(centroid_PA~xAxis, type='l', lwd=2, lty=2, col='gray20', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM',])
  errbar(x=c(1:5), y=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'], yplus=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'] + centroid_final85$sd_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'], yminus=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'] - centroid_final85$sd_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'], errbar.col='gray20', lwd=2, cap=0, add=T)
  points(centroid_PA~xAxis, type='p', pch=19, col='gray20', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM',])
  
  points(centroid_biom~xAxis_min, type='l', lwd=2, col='red', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT',])
  errbar(x=c(0.95,1.95,2.95,3.95,4.95)-0.05, y=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'], yplus=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'] + centroid_final85$sd_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'], yminus=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'] - centroid_final85$sd_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'], errbar.col='red', lwd=2, cap=0, add=T)
  points(centroid_biom~xAxis_min, type='p', pch=19, col='red', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT',])
  points(centroid_PA~xAxis_min, type='l', lwd=2, lty=2, col='red', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT',])
  errbar(x=c(0.95,1.95,2.95,3.95,4.95)-0.05, y=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'], yplus=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'] + centroid_final85$sd_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'], yminus=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'] - centroid_final85$sd_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'], errbar.col='red', lwd=2, cap=0, add=T)
  points(centroid_PA~xAxis_min, type='p', pch=19, col='red', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT',])

# LOBSTER centroid
i=1 # 
  plot(centroid_biom~xAxis_plus, xlim=c(0.9,5.1), ylim=c(41.5,54.5), xaxt='n', yaxt='n', xlab='', ylab='Centroid (latitude)', type='l', lwd=2, col='dodgerblue2', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM',])
  errbar(x=c(1.05,2.05,3.05,4.05,5.05)+0.05, y=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'], yplus=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'] + centroid_final85$sd_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'], yminus=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'] - centroid_final85$sd_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'], errbar.col='dodgerblue2', lwd=2, cap=0, add=T)
  points(centroid_biom~xAxis_plus, type='p', pch=19, col='dodgerblue2', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM',])
  points(centroid_PA~xAxis_plus, type='l', lwd=2, lty=2, col='dodgerblue2', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM',])
  errbar(x=c(1.05,2.05,3.05,4.05,5.05)+0.05, y=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'], yplus=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'] + centroid_final85$sd_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'], yminus=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'] - centroid_final85$sd_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'], errbar.col='dodgerblue2', lwd=2, cap=0, add=T)
  points(centroid_PA~xAxis_plus, type='p', pch=19, col='dodgerblue2', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM',])
  axis(2,at=c(42,46,50,54), lwd=1.3, cex.axis=0.95, las=1)
  axis(1,at=c(1,2,3,4,5), labels=F, lwd=1.3, cex.axis=.9, las=1)
  text(x=c(1:5+0.44), y=39.5, c('2007-2020','2021-2040','2041-2060','2061-2080','2081-2100'), srt = 305, pos = 1, cex=1)
  text(x=3, y=37.1, 'Years', pos = 1, cex=1.2)
  mtext('(d)',line=-1.25,adj=0.01,font=2)
  
  points(centroid_biom~xAxis, type='l', lwd=2, col='gray20', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM',])
  errbar(x=c(1:5), y=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'], yplus=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'] + centroid_final85$sd_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'], yminus=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'] - centroid_final85$sd_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'], errbar.col='gray20', lwd=2, cap=0, add=T)
  points(centroid_biom~xAxis, type='p', pch=19, col='gray20', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM',])
  points(centroid_PA~xAxis, type='l', lwd=2, lty=2, col='gray20', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM',])
  errbar(x=c(1:5), y=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'], yplus=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'] + centroid_final85$sd_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'], yminus=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'] - centroid_final85$sd_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'], errbar.col='gray20', lwd=2, cap=0, add=T)
  points(centroid_PA~xAxis, type='p', pch=19, col='gray20', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM',])
  
  points(centroid_biom~xAxis_min, type='l', lwd=2, col='red', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT',])
  par(xpd=F)
  errbar(x=c(0.95,1.95,2.95,3.95,4.95)-0.05, y=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'], yplus=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'] + centroid_final85$sd_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'], yminus=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'] - centroid_final85$sd_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'], errbar.col='red', lwd=2, cap=0, add=T)
  par(xpd=NA)
  points(centroid_biom~xAxis_min, type='p', pch=19, col='red', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT',])
  points(centroid_PA~xAxis_min, type='l', lwd=2, lty=2, col='red', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT',])
  errbar(x=c(0.95,1.95,2.95,3.95,4.95)-0.05, y=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'], yplus=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'] + centroid_final85$sd_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'], yminus=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'] - centroid_final85$sd_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'], errbar.col='red', lwd=2, cap=0, add=T)
  points(centroid_PA~xAxis_min, type='p', pch=19, col='red', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT',])
  
# HALIBUT north biomass =================================================================================================================
i=6
  plot(north_biom~xAxis_plus, xlim=c(0.9,4.1), ylim=c(-10,57), xaxt='n', yaxt='n', xlab='', ylab='Change in habitat (%)', type='l', lwd=2, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
  errbar(x=c(1.05,2.05,3.05,4.05)+0.05, y=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yplus=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] + biomass_final85$sdNorth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yminus=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] - biomass_final85$sdNorth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], errbar.col='dodgerblue2', lwd=2, cap=0, add=T)
  points(north_biom~xAxis_plus, type='p', pch=19, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
  points(north_PA~xAxis_plus, type='l', lwd=2, lty=2, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
  errbar(x=c(1.05,2.05,3.05,4.05)+0.05, y=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yplus=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] + biomass_final85$sdNorth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yminus=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] - biomass_final85$sdNorth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], errbar.col='dodgerblue2', lwd=2, cap=0, add=T)
  points(north_PA~xAxis_plus, type='p', pch=19, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
  abline(a=0, b=0, lwd=.75, col='gray50', lty=2, xpd=F)
  axis(2,at=c(-10,0,10,20,30,40,50), lwd=1.3, cex.axis=0.95, las=1)
  mtext('(e)',line=-1.25,adj=0.01,font=2)
  # 
  points(north_biom~xAxis, type='l', lwd=2, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
  errbar(x=c(1:4), y=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yplus=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] + biomass_final85$sdNorth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yminus=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] - biomass_final85$sdNorth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], errbar.col='gray20', lwd=2, cap=0, add=T)
  points(north_biom~xAxis, type='p', pch=19, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
  points(north_PA~xAxis, type='l', lwd=2, lty=2, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
  errbar(x=c(1:4), y=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yplus=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] + biomass_final85$sdNorth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yminus=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] - biomass_final85$sdNorth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], errbar.col='gray20', lwd=2, cap=0, add=T)
  points(north_PA~xAxis, type='p', pch=19, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
  
  points(north_biom~xAxis_min, type='l', lwd=2, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
  errbar(x=c(0.95,1.95,2.95,3.95)-0.05, y=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yplus=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] + biomass_final85$sdNorth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yminus=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] - biomass_final85$sdNorth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], errbar.col='red', lwd=2, cap=0, add=T)
  points(north_biom~xAxis_min, type='p', pch=19, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
  points(north_PA~xAxis_min, type='l', lwd=2, lty=2, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
  errbar(x=c(0.95,1.95,2.95,3.95)-0.05, y=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yplus=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] + biomass_final85$sdNorth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yminus=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] - biomass_final85$sdNorth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], errbar.col='red', lwd=2, cap=0, add=T)
  points(north_PA~xAxis_min, type='p', pch=19, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
  legend(x=1.06, y=61, legend=c('GLM pres-abs','GLM biom','GAM pres-abs','GAM biom','BRT pres-abs','BRT biom'),col=c('dodgerblue2','dodgerblue2','gray20','gray20','red','red'),lty=c(2,1,2,1,2,1), cex=0.9, lwd=2, bty='n')

# Perch north biomass
i=7
  plot(north_biom~xAxis_plus, xlim=c(0.9,4.1), ylim=c(-5,75), xaxt='n', yaxt='n', xlab='', ylab='Change in habitat (%)', type='l', lwd=2, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
  errbar(x=c(1.05,2.05,3.05,4.05)+0.05, y=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yplus=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] + biomass_final85$sdNorth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yminus=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] - biomass_final85$sdNorth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], errbar.col='dodgerblue2', lwd=2, cap=0, add=T)
  points(north_biom~xAxis_plus, type='p', pch=19, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
  points(north_PA~xAxis_plus, type='l', lwd=2, lty=2, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
  errbar(x=c(1.05,2.05,3.05,4.05)+0.05, y=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yplus=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] + biomass_final85$sdNorth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yminus=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] - biomass_final85$sdNorth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], errbar.col='dodgerblue2', lwd=2, cap=0, add=T)
  points(north_PA~xAxis_plus, type='p', pch=19, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
  abline(a=0, b=0, lwd=.75, col='gray50', lty=2, xpd=F)
  axis(2,at=c(0,20,40,60), lwd=1.3, cex.axis=0.95, las=1)
  mtext('(f)',line=-1.25,adj=0.01,font=2)
  # 
  points(north_biom~xAxis, type='l', lwd=2, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
  errbar(x=c(1:4), y=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yplus=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] + biomass_final85$sdNorth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yminus=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] - biomass_final85$sdNorth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], errbar.col='gray20', lwd=2, cap=0, add=T)
  points(north_biom~xAxis, type='p', pch=19, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
  points(north_PA~xAxis, type='l', lwd=2, lty=2, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
  errbar(x=c(1:4), y=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yplus=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] + biomass_final85$sdNorth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yminus=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] - biomass_final85$sdNorth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], errbar.col='gray20', lwd=2, cap=0, add=T)
  points(north_PA~xAxis, type='p', pch=19, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
  
  points(north_biom~xAxis_min, type='l', lwd=2, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
  errbar(x=c(0.95,1.95,2.95,3.95)-0.05, y=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yplus=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] + biomass_final85$sdNorth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yminus=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] - biomass_final85$sdNorth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], errbar.col='red', lwd=2, cap=0, add=T)
  points(north_biom~xAxis_min, type='p', pch=19, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
  points(north_PA~xAxis_min, type='l', lwd=2, lty=2, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
  errbar(x=c(0.95,1.95,2.95,3.95)-0.05, y=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yplus=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] + biomass_final85$sdNorth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yminus=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] - biomass_final85$sdNorth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], errbar.col='red', lwd=2, cap=0, add=T)
  points(north_PA~xAxis_min, type='p', pch=19, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
  
# SUMMER FLOUNDER north biomass
i=2
  plot(north_biom~xAxis_plus, xlim=c(0.9,4.1), ylim=c(-40,270), xaxt='n', yaxt='n', xlab='', ylab='Change in habitat (%)', type='l', lwd=2, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
  par(xpd=T)
  errbar(x=c(1.05,2.05,3.05,4.05)+0.05, y=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yplus=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] + biomass_final85$sdNorth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yminus=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] - biomass_final85$sdNorth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], errbar.col='dodgerblue2', lwd=2, cap=0, add=T)
  par(xpd=NA)
  points(north_biom~xAxis_plus, type='p', pch=19, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
  points(north_PA~xAxis_plus, type='l', lwd=2, lty=2, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
  errbar(x=c(1.05,2.05,3.05,4.05)+0.05, y=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yplus=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] + biomass_final85$sdNorth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yminus=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] - biomass_final85$sdNorth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], errbar.col='dodgerblue2', lwd=2, cap=0, add=T)
  points(north_PA~xAxis_plus, type='p', pch=19, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
  abline(a=0, b=0, lwd=.75, col='gray50', lty=2, xpd=F)
  axis(2,at=c(-30,30,90,150,210,270), lwd=1.3, cex.axis=0.95, las=1)
  mtext('(g)',line=-1.25,adj=0.01,font=2)
  # 
  points(north_biom~xAxis, type='l', lwd=2, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
  par(xpd=T)
  errbar(x=c(1:4), y=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yplus=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] + biomass_final85$sdNorth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yminus=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] - biomass_final85$sdNorth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], errbar.col='gray20', lwd=2, cap=0, add=T)
  par(xpd=NA)
  points(north_biom~xAxis, type='p', pch=19, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
  points(north_PA~xAxis, type='l', lwd=2, lty=2, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
  errbar(x=c(1:4), y=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yplus=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] + biomass_final85$sdNorth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yminus=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] - biomass_final85$sdNorth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], errbar.col='gray20', lwd=2, cap=0, add=T)
  points(north_PA~xAxis, type='p', pch=19, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
  
  points(north_biom~xAxis_min, type='l', lwd=2, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
  errbar(x=c(0.95,1.95,2.95,3.95)-0.05, y=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yplus=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] + biomass_final85$sdNorth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yminus=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] - biomass_final85$sdNorth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], errbar.col='red', lwd=2, cap=0, add=T)
  points(north_biom~xAxis_min, type='p', pch=19, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
  points(north_PA~xAxis_min, type='l', lwd=2, lty=2, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
  errbar(x=c(0.95,1.95,2.95,3.95)-0.05, y=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yplus=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] + biomass_final85$sdNorth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yminus=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] - biomass_final85$sdNorth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], errbar.col='red', lwd=2, cap=0, add=T)
  points(north_PA~xAxis_min, type='p', pch=19, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
  
# LOBSTER north biomass
i=1
  plot(north_biom~xAxis_plus, xlim=c(0.9,4.1), ylim=c(-100,145), xaxt='n', yaxt='n', xlab='', ylab='Change in habitat (%)', type='l', lwd=2, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
  errbar(x=c(1.05,2.05,3.05,4.05)+0.05, y=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yplus=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] + biomass_final85$sdNorth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yminus=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] - biomass_final85$sdNorth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], errbar.col='dodgerblue2', lwd=2, cap=0, add=T)
  points(north_biom~xAxis_plus, type='p', pch=19, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
  points(north_PA~xAxis_plus, type='l', lwd=2, lty=2, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
  errbar(x=c(1.05,2.05,3.05,4.05)+0.05, y=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yplus=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] + biomass_final85$sdNorth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yminus=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] - biomass_final85$sdNorth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], errbar.col='dodgerblue2', lwd=2, cap=0, add=T)
  points(north_PA~xAxis_plus, type='p', pch=19, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
  abline(a=0, b=0, lwd=.75, col='gray50', lty=2, xpd=F)
  axis(2,at=c(-60,-20,20,60,100,140), lwd=1.3, cex.axis=0.95, las=1)
  axis(1,at=c(1,2,3,4), labels=F, lwd=1.3, cex.axis=.9, las=1)
  text(x=c(1:4+0.34), y=-139.4, c('2021-2040','2041-2060','2061-2080','2081-2100'), srt = 305, pos = 1, cex=1)
  text(x=2.55, y=-190, 'Years', pos = 1, cex=1.2)
  mtext('(h)',line=-1.25,adj=0.01,font=2)
  # 
  points(north_biom~xAxis, type='l', lwd=2, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
  errbar(x=c(1:4), y=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yplus=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] + biomass_final85$sdNorth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yminus=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] - biomass_final85$sdNorth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], errbar.col='gray20', lwd=2, cap=0, add=T)
  points(north_biom~xAxis, type='p', pch=19, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
  points(north_PA~xAxis, type='l', lwd=2, lty=2, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
  errbar(x=c(1:4), y=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yplus=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] + biomass_final85$sdNorth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yminus=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] - biomass_final85$sdNorth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], errbar.col='gray20', lwd=2, cap=0, add=T)
  points(north_PA~xAxis, type='p', pch=19, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
  
  points(north_biom~xAxis_min, type='l', lwd=2, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
  errbar(x=c(0.95,1.95,2.95,3.95)-0.05, y=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yplus=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] + biomass_final85$sdNorth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yminus=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] - biomass_final85$sdNorth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], errbar.col='red', lwd=2, cap=0, add=T)
  points(north_biom~xAxis_min, type='p', pch=19, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
  points(north_PA~xAxis_min, type='l', lwd=2, lty=2, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
  errbar(x=c(0.95,1.95,2.95,3.95)-0.05, y=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yplus=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] + biomass_final85$sdNorth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yminus=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] - biomass_final85$sdNorth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], errbar.col='red', lwd=2, cap=0, add=T)
  points(north_PA~xAxis_min, type='p', pch=19, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
  
# HALIBUT south biomass
i=6
  plot(south_biom~xAxis_plus, xlim=c(0.9,4.1), ylim=c(-77,5), xaxt='n', yaxt='n', xlab='', ylab='', type='l', lwd=2, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
  errbar(x=c(1.05,2.05,3.05,4.05)+0.05, y=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yplus=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] + biomass_final85$sdSouth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yminus=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] - biomass_final85$sdSouth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], errbar.col='dodgerblue2', lwd=2, cap=0, add=T)
  points(south_biom~xAxis_plus, type='p', pch=19, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
  points(south_PA~xAxis_plus, type='l', lwd=2, lty=2, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
  errbar(x=c(1.05,2.05,3.05,4.05)+0.05, y=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yplus=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] + biomass_final85$sdSouth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yminus=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] - biomass_final85$sdSouth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], errbar.col='dodgerblue2', lwd=2, cap=0, add=T)
  points(south_PA~xAxis_plus, type='p', pch=19, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
  abline(a=0, b=0, lwd=.75, col='gray50', lty=2, xpd=F)
  axis(2,at=c(-60,-40,-20,0), lwd=1.3, cex.axis=0.95, las=1)
  mtext('(i)',line=-1.25,adj=0.01,font=2)
  # 
  points(south_biom~xAxis, type='l', lwd=2, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
  errbar(x=c(1:4), y=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yplus=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] + biomass_final85$sdSouth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yminus=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] - biomass_final85$sdSouth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], errbar.col='gray20', lwd=2, cap=0, add=T)
  points(south_biom~xAxis, type='p', pch=19, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
  points(south_PA~xAxis, type='l', lwd=2, lty=2, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
  errbar(x=c(1:4), y=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yplus=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] + biomass_final85$sdSouth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yminus=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] - biomass_final85$sdSouth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], errbar.col='gray20', lwd=2, cap=0, add=T)
  points(south_PA~xAxis, type='p', pch=19, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
  
  points(south_biom~xAxis_min, type='l', lwd=2, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
  errbar(x=c(0.95,1.95,2.95,3.95)-0.05, y=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yplus=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] + biomass_final85$sdSouth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yminus=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] - biomass_final85$sdSouth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], errbar.col='red', lwd=2, cap=0, add=T)
  points(south_biom~xAxis_min, type='p', pch=19, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
  points(south_PA~xAxis_min, type='l', lwd=2, lty=2, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
  errbar(x=c(0.95,1.95,2.95,3.95)-0.05, y=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yplus=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] + biomass_final85$sdSouth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yminus=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] - biomass_final85$sdSouth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], errbar.col='red', lwd=2, cap=0, add=T)
  points(south_PA~xAxis_min, type='p', pch=19, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
  
# PERCH south biomass
i=7
  plot(south_biom~xAxis_plus, xlim=c(0.9,4.1), ylim=c(-85,8), xaxt='n', yaxt='n', xlab='', ylab='', type='l', lwd=2, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
  errbar(x=c(1.05,2.05,3.05,4.05)+0.05, y=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yplus=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] + biomass_final85$sdSouth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yminus=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] - biomass_final85$sdSouth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], errbar.col='dodgerblue2', lwd=2, cap=0, add=T)
  points(south_biom~xAxis_plus, type='p', pch=19, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
  points(south_PA~xAxis_plus, type='l', lwd=2, lty=2, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
  errbar(x=c(1.05,2.05,3.05,4.05)+0.05, y=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yplus=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] + biomass_final85$sdSouth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yminus=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] - biomass_final85$sdSouth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], errbar.col='dodgerblue2', lwd=2, cap=0, add=T)
  points(south_PA~xAxis_plus, type='p', pch=19, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
  abline(a=0, b=0, lwd=.75, col='gray50', lty=2, xpd=F)
  axis(2,at=c(-80,-60,-40,-20,0), lwd=1.3, cex.axis=0.95, las=1)
  mtext('(j)',line=-1.25,adj=0.01,font=2)
  # 
  points(south_biom~xAxis, type='l', lwd=2, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
  errbar(x=c(1:4), y=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yplus=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] + biomass_final85$sdSouth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yminus=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] - biomass_final85$sdSouth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], errbar.col='gray20', lwd=2, cap=0, add=T)
  points(south_biom~xAxis, type='p', pch=19, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
  points(south_PA~xAxis, type='l', lwd=2, lty=2, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
  errbar(x=c(1:4), y=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yplus=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] + biomass_final85$sdSouth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yminus=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] - biomass_final85$sdSouth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], errbar.col='gray20', lwd=2, cap=0, add=T)
  points(south_PA~xAxis, type='p', pch=19, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
  
  points(south_biom~xAxis_min, type='l', lwd=2, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
  errbar(x=c(0.95,1.95,2.95,3.95)-0.05, y=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yplus=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] + biomass_final85$sdSouth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yminus=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] - biomass_final85$sdSouth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], errbar.col='red', lwd=2, cap=0, add=T)
  points(south_biom~xAxis_min, type='p', pch=19, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
  points(south_PA~xAxis_min, type='l', lwd=2, lty=2, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
  errbar(x=c(0.95,1.95,2.95,3.95)-0.05, y=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yplus=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] + biomass_final85$sdSouth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yminus=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] - biomass_final85$sdSouth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], errbar.col='red', lwd=2, cap=0, add=T)
  points(south_PA~xAxis_min, type='p', pch=19, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
  
# SUMMER FLOUNDER south biomass
i=2
  plot(south_biom~xAxis_plus, xlim=c(0.9,4.1), ylim=c(-100,5), xaxt='n', yaxt='n', xlab='', ylab='', type='l', lwd=2, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
  errbar(x=c(1.05,2.05,3.05,4.05)+0.05, y=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yplus=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] + biomass_final85$sdSouth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yminus=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] - biomass_final85$sdSouth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], errbar.col='dodgerblue2', lwd=2, cap=0, add=T)
  points(south_biom~xAxis_plus, type='p', pch=19, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
  points(south_PA~xAxis_plus, type='l', lwd=2, lty=2, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
  errbar(x=c(1.05,2.05,3.05,4.05)+0.05, y=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yplus=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] + biomass_final85$sdSouth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yminus=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] - biomass_final85$sdSouth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], errbar.col='dodgerblue2', lwd=2, cap=0, add=T)
  points(south_PA~xAxis_plus, type='p', pch=19, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
  abline(a=0, b=0, lwd=.75, col='gray50', lty=2, xpd=F)
  axis(2,at=c(-90,-70,-50,-30,-10), lwd=1.3, cex.axis=0.95, las=1)
  mtext('(k)',line=-1.25,adj=0.01,font=2)
  # 
  points(south_biom~xAxis, type='l', lwd=2, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
  errbar(x=c(1:4), y=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yplus=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] + biomass_final85$sdSouth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yminus=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] - biomass_final85$sdSouth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], errbar.col='gray20', lwd=2, cap=0, add=T)
  points(south_biom~xAxis, type='p', pch=19, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
  points(south_PA~xAxis, type='l', lwd=2, lty=2, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
  errbar(x=c(1:4), y=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yplus=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] + biomass_final85$sdSouth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yminus=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] - biomass_final85$sdSouth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], errbar.col='gray20', lwd=2, cap=0, add=T)
  points(south_PA~xAxis, type='p', pch=19, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
  
  points(south_biom~xAxis_min, type='l', lwd=2, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
  errbar(x=c(0.95,1.95,2.95,3.95)-0.05, y=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yplus=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] + biomass_final85$sdSouth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yminus=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] - biomass_final85$sdSouth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], errbar.col='red', lwd=2, cap=0, add=T)
  points(south_biom~xAxis_min, type='p', pch=19, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
  points(south_PA~xAxis_min, type='l', lwd=2, lty=2, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
  errbar(x=c(0.95,1.95,2.95,3.95)-0.05, y=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yplus=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] + biomass_final85$sdSouth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yminus=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] - biomass_final85$sdSouth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], errbar.col='red', lwd=2, cap=0, add=T)
  points(south_PA~xAxis_min, type='p', pch=19, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
  
# LOBSTER south biomass
i=1
  plot(south_biom~xAxis_plus, xlim=c(0.9,4.1), ylim=c(-85,16), xaxt='n', yaxt='n', xlab='', ylab='', type='l', lwd=2, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
  errbar(x=c(1.05,2.05,3.05,4.05)+0.05, y=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yplus=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] + biomass_final85$sdSouth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yminus=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] - biomass_final85$sdSouth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], errbar.col='dodgerblue2', lwd=2, cap=0, add=T)
  points(south_biom~xAxis_plus, type='p', pch=19, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
  points(south_PA~xAxis_plus, type='l', lwd=2, lty=2, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
  errbar(x=c(1.05,2.05,3.05,4.05)+0.05, y=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yplus=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] + biomass_final85$sdSouth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yminus=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] - biomass_final85$sdSouth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], errbar.col='dodgerblue2', lwd=2, cap=0, add=T)
  points(south_PA~xAxis_plus, type='p', pch=19, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
  abline(a=0, b=0, lwd=.75, col='gray50', lty=2, xpd=F)
  axis(2,at=c(-80,-60,-40,-20,0), lwd=1.3, cex.axis=.95, las=1)
  axis(1,at=c(1,2,3,4), labels=F, lwd=1.3, cex.axis=.9, las=1)
  text(x=c(1:4+0.33), y=-100.6, c('2021-2040','2041-2060','2061-2080','2081-2100'), srt = 305, pos = 1, cex=1)
  text(x=2.5, y=-120, 'Years', pos = 1, cex=1.2)
  mtext('(l)',line=-1.25,adj=0.01,font=2)
  
  points(south_biom~xAxis, type='l', lwd=2, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
  par(xpd=T)
  errbar(x=c(1:4), y=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yplus=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] + biomass_final85$sdSouth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yminus=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] - biomass_final85$sdSouth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], errbar.col='gray20', lwd=2, cap=0, add=T)
  par(xpd=NA)
  points(south_biom~xAxis, type='p', pch=19, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
  points(south_PA~xAxis, type='l', lwd=2, lty=2, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
  errbar(x=c(1:4), y=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yplus=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] + biomass_final85$sdSouth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yminus=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] - biomass_final85$sdSouth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], errbar.col='gray20', lwd=2, cap=0, add=T)
  points(south_PA~xAxis, type='p', pch=19, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
  
  points(south_biom~xAxis_min, type='l', lwd=2, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
  errbar(x=c(0.95,1.95,2.95,3.95)-0.05, y=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yplus=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] + biomass_final85$sdSouth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yminus=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] - biomass_final85$sdSouth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], errbar.col='red', lwd=2, cap=0, add=T)
  points(south_biom~xAxis_min, type='p', pch=19, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
  points(south_PA~xAxis_min, type='l', lwd=2, lty=2, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
  errbar(x=c(0.95,1.95,2.95,3.95)-0.05, y=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yplus=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] + biomass_final85$sdSouth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yminus=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] - biomass_final85$sdSouth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], errbar.col='red', lwd=2, cap=0, add=T)
  points(south_PA~xAxis_min, type='p', pch=19, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
  
dev.off()
  



# =====================================================================================================
# APPENDIX figure ======================================================================
# =====================================================================================================
  
pdf(width=6.5, height=6.3, file='figures/Uncertainty_MS_2018/NicheMod_comparison_Appendix_Dec2018.pdf')
#plot.new()
plot_matrix <- rbind(c(10:16),c(17,1,20,4,21,7,22),c(17,2,18,5,19,8,23),c(17,3,18,6,19,9,24),c(25:31))
layout(mat=matrix(plot_matrix, nrow=5, ncol=7), widths=c(1.2,4,1.2,4,.6,4,.5), heights=c(0.1,4,4,4,1.6))# 1.2, Need to adjust this to fit the three figures
#layout.show(n=9)

# SABLEFISH centroid
i=4
par(mar=c(0,0,0,0), oma=c(0,0,0,0), tcl=-.25, mgp=c(2.25,0.5,0), cex.lab=1.2, xpd = NA)
plot(centroid_biom~xAxis_plus, xlim=c(0.9,5.1), ylim=c(54.4,57.1), yaxt='n', xaxt='n', xlab='', ylab='Centroid (latitude)', type='l', lwd=2, col='dodgerblue2', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM',])
errbar(x=c(1.05,2.05,3.05,4.05,5.05)+0.05, y=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'], yplus=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'] + centroid_final85$sd_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'], yminus=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'] - centroid_final85$sd_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'], errbar.col='dodgerblue2', lwd=2, cap=0, add=T)
points(centroid_biom~xAxis_plus, type='p', pch=19, col='dodgerblue2', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM',])
points(centroid_PA~xAxis_plus, type='l', lwd=2, lty=2, col='dodgerblue2', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM',])
errbar(x=c(1.05,2.05,3.05,4.05,5.05)+0.05, y=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'], yplus=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'] + centroid_final85$sd_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'], yminus=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'] - centroid_final85$sd_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'], errbar.col='dodgerblue2', lwd=2, cap=0, add=T)
points(centroid_PA~xAxis_plus, type='p', pch=19, col='dodgerblue2', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM',])
axis(2,at=c(55,56,57), lwd=1.3, cex.axis=0.95, las=1)
mtext('(a)',line=-1.25,adj=0.01,font=2)
#  
points(centroid_biom~xAxis, type='l', lwd=2, col='gray20', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM',])
errbar(x=c(1:5), y=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'], yplus=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'] + centroid_final85$sd_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'], yminus=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'] - centroid_final85$sd_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'], errbar.col='gray20', lwd=2, cap=0, add=T)
points(centroid_biom~xAxis, type='p', pch=19, col='gray20', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM',])
points(centroid_PA~xAxis, type='l', lwd=2, lty=2, col='gray20', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM',])
errbar(x=c(1:5), y=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'], yplus=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'] + centroid_final85$sd_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'], yminus=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'] - centroid_final85$sd_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'], errbar.col='gray20', lwd=2, cap=0, add=T)
points(centroid_PA~xAxis, type='p', pch=19, col='gray20', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM',])

points(centroid_biom~xAxis_min, type='l', lwd=2, col='red', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT',])
errbar(x=c(0.95,1.95,2.95,3.95,4.95)-0.05, y=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'], yplus=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'] + centroid_final85$sd_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'], yminus=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'] - centroid_final85$sd_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'], errbar.col='red', lwd=2, cap=0, add=T)
points(centroid_biom~xAxis_min, type='p', pch=19, col='red', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT',])
points(centroid_PA~xAxis_min, type='l', lwd=2, lty=2, col='red', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT',])
errbar(x=c(0.95,1.95,2.95,3.95,4.95)-0.05, y=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'], yplus=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'] + centroid_final85$sd_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'], yminus=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'] - centroid_final85$sd_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'], errbar.col='red', lwd=2, cap=0, add=T)
points(centroid_PA~xAxis_min, type='p', pch=19, col='red', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT',])
#rasterImage(rock, .5,57.5,2,58) # placeholder for adding the fish graphics

# SQUID centroid
i=5
plot(centroid_biom~xAxis_plus, xlim=c(0.9,5.1), ylim=c(36,54), yaxt='n', xaxt='n', xlab='', ylab='Centroid (latitude)', type='l', lwd=2, col='dodgerblue2', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM',])
errbar(x=c(1.05,2.05,3.05,4.05,5.05)+0.05, y=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'], yplus=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'] + centroid_final85$sd_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'], yminus=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'] - centroid_final85$sd_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'], errbar.col='dodgerblue2', lwd=2, cap=0, add=T)
points(centroid_biom~xAxis_plus, type='p', pch=19, col='dodgerblue2', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM',])
points(centroid_PA~xAxis_plus, type='l', lwd=2, lty=2, col='dodgerblue2', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM',])
errbar(x=c(1.05,2.05,3.05,4.05,5.05)+0.05, y=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'], yplus=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'] + centroid_final85$sd_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'], yminus=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'] - centroid_final85$sd_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'], errbar.col='dodgerblue2', lwd=2, cap=0, add=T)
points(centroid_PA~xAxis_plus, type='p', pch=19, col='dodgerblue2', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM',])
axis(2,at=c(40,45,50), lwd=1.3, cex.axis=0.95, las=1)
mtext('(b)',line=-1.25,adj=0.01,font=2)
#
points(centroid_biom~xAxis, type='l', lwd=2, col='gray20', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM',])
errbar(x=c(1:5), y=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'], yplus=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'] + centroid_final85$sd_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'], yminus=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'] - centroid_final85$sd_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'], errbar.col='gray20', lwd=2, cap=0, add=T)
points(centroid_biom~xAxis, type='p', pch=19, col='gray20', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM',])
points(centroid_PA~xAxis, type='l', lwd=2, lty=2, col='gray20', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM',])
errbar(x=c(1:5), y=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'], yplus=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'] + centroid_final85$sd_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'], yminus=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'] - centroid_final85$sd_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'], errbar.col='gray20', lwd=2, cap=0, add=T)
points(centroid_PA~xAxis, type='p', pch=19, col='gray20', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM',])

points(centroid_biom~xAxis_min, type='l', lwd=2, col='red', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT',])
errbar(x=c(0.95,1.95,2.95,3.95,4.95)-0.05, y=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'], yplus=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'] + centroid_final85$sd_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'], yminus=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'] - centroid_final85$sd_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'], errbar.col='red', lwd=2, cap=0, add=T)
points(centroid_biom~xAxis_min, type='p', pch=19, col='red', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT',])
points(centroid_PA~xAxis_min, type='l', lwd=2, lty=2, col='red', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT',])
errbar(x=c(0.95,1.95,2.95,3.95,4.95)-0.05, y=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'], yplus=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'] + centroid_final85$sd_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'], yminus=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'] - centroid_final85$sd_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'], errbar.col='red', lwd=2, cap=0, add=T)
points(centroid_PA~xAxis_min, type='p', pch=19, col='red', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT',])

# BSB centroid
i=3 # 
plot(centroid_biom~xAxis_plus, xlim=c(0.9,5.1), ylim=c(29,42), yaxt='n', xaxt='n', xlab='', ylab='Centroid (latitude)', type='l', lwd=2, col='dodgerblue2', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM',])
errbar(x=c(1.05,2.05,3.05,4.05,5.05)+0.05, y=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'], yplus=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'] + centroid_final85$sd_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'], yminus=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'] - centroid_final85$sd_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'], errbar.col='dodgerblue2', lwd=2, cap=0, add=T)
points(centroid_biom~xAxis_plus, type='p', pch=19, col='dodgerblue2', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM',])
points(centroid_PA~xAxis_plus, type='l', lwd=2, lty=2, col='dodgerblue2', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM',])
errbar(x=c(1.05,2.05,3.05,4.05,5.05)+0.05, y=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'], yplus=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'] + centroid_final85$sd_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'], yminus=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'] - centroid_final85$sd_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM'], errbar.col='dodgerblue2', lwd=2, cap=0, add=T)
points(centroid_PA~xAxis_plus, type='p', pch=19, col='dodgerblue2', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GLM',])
axis(2,at=c(30,34,38,42), lwd=1.3, cex.axis=0.95, las=1)
axis(1,at=c(1,2,3,4,5), labels=F, lwd=1.3, cex.axis=.9, las=1)
text(x=c(1:5+0.44), y=26.85, c('2007-2020','2021-2040','2041-2060','2061-2080','2081-2100'), srt = 305, pos = 1, cex=1)
text(x=3, y=24.3, 'Years', pos = 1, cex=1.2)
mtext('(c)',line=-1.25,adj=0.01,font=2)
# 
points(centroid_biom~xAxis, type='l', lwd=2, col='gray20', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM',])
errbar(x=c(1:5), y=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'], yplus=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'] + centroid_final85$sd_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'], yminus=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'] - centroid_final85$sd_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'], errbar.col='gray20', lwd=2, cap=0, add=T)
points(centroid_biom~xAxis, type='p', pch=19, col='gray20', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM',])
points(centroid_PA~xAxis, type='l', lwd=2, lty=2, col='gray20', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM',])
errbar(x=c(1:5), y=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'], yplus=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'] + centroid_final85$sd_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'], yminus=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'] - centroid_final85$sd_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM'], errbar.col='gray20', lwd=2, cap=0, add=T)
points(centroid_PA~xAxis, type='p', pch=19, col='gray20', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='GAM',])

points(centroid_biom~xAxis_min, type='l', lwd=2, col='red', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT',])
par(xpd=F)
errbar(x=c(0.95,1.95,2.95,3.95,4.95)-0.05, y=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'], yplus=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'] + centroid_final85$sd_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'], yminus=centroid_final85$centroid_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'] - centroid_final85$sd_biom[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'], errbar.col='red', lwd=2, cap=0, add=T)
par(xpd=NA)
points(centroid_biom~xAxis_min, type='p', pch=19, col='red', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT',])
points(centroid_PA~xAxis_min, type='l', lwd=2, lty=2, col='red', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT',])
errbar(x=c(0.95,1.95,2.95,3.95,4.95)-0.05, y=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'], yplus=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'] + centroid_final85$sd_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'], yminus=centroid_final85$centroid_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'] - centroid_final85$sd_PA[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT'], errbar.col='red', lwd=2, cap=0, add=T)
points(centroid_PA~xAxis_min, type='p', pch=19, col='red', centroid_final85[centroid_final85$species==SPECIES[i] & centroid_final85$nicheMod=='BRT',])

# SABLEFISH north biomass =================================================================================================================
i=4
par(xpd=NA)
plot(north_biom~xAxis_plus, xlim=c(0.9,4.1), ylim=c(-26,105), yaxt='n', xaxt='n', xlab='', ylab='Change in habitat (%)', type='l', lwd=2, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
errbar(x=c(1.05,2.05,3.05,4.05)+0.05, y=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yplus=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] + biomass_final85$sdNorth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yminus=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] - biomass_final85$sdNorth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], errbar.col='dodgerblue2', lwd=2, cap=0, add=T)
points(north_biom~xAxis_plus, type='p', pch=19, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
points(north_PA~xAxis_plus, type='l', lwd=2, lty=2, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
errbar(x=c(1.05,2.05,3.05,4.05)+0.05, y=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yplus=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] + biomass_final85$sdNorth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yminus=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] - biomass_final85$sdNorth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], errbar.col='dodgerblue2', lwd=2, cap=0, add=T)
points(north_PA~xAxis_plus, type='p', pch=19, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
abline(a=0, b=0, lwd=.75, col='gray50', lty=2, xpd=F)
axis(2,at=c(-20,0,20,40,60,80), lwd=1.3, cex.axis=0.95, las=1)
mtext('(d)',line=-1.25,adj=0.01,font=2)
# 
points(north_biom~xAxis, type='l', lwd=2, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
errbar(x=c(1:4), y=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yplus=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] + biomass_final85$sdNorth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yminus=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] - biomass_final85$sdNorth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], errbar.col='gray20', lwd=2, cap=0, add=T)
points(north_biom~xAxis, type='p', pch=19, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
points(north_PA~xAxis, type='l', lwd=2, lty=2, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
errbar(x=c(1:4), y=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yplus=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] + biomass_final85$sdNorth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yminus=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] - biomass_final85$sdNorth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], errbar.col='gray20', lwd=2, cap=0, add=T)
points(north_PA~xAxis, type='p', pch=19, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])

points(north_biom~xAxis_min, type='l', lwd=2, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
errbar(x=c(0.95,1.95,2.95,3.95)-0.05, y=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yplus=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] + biomass_final85$sdNorth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yminus=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] - biomass_final85$sdNorth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], errbar.col='red', lwd=2, cap=0, add=T)
points(north_biom~xAxis_min, type='p', pch=19, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
points(north_PA~xAxis_min, type='l', lwd=2, lty=2, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
errbar(x=c(0.95,1.95,2.95,3.95)-0.05, y=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yplus=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] + biomass_final85$sdNorth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yminus=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] - biomass_final85$sdNorth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], errbar.col='red', lwd=2, cap=0, add=T)
points(north_PA~xAxis_min, type='p', pch=19, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
legend(x=1.06, y=115, legend=c('GLM pres-abs','GLM biom','GAM pres-abs','GAM biom','BRT pres-abs','BRT biom'),col=c('dodgerblue2','dodgerblue2','gray20','gray20','red','red'),lty=c(2,1,2,1,2,1), cex=0.9, lwd=2, bty='n')

# SQUID north biomass
i=5
par(xpd=NA)
plot(north_PA~xAxis_plus, xlim=c(0.9,4.1), ylim=c(-1,375), yaxt='n', xaxt='n', xlab='', ylab='Change in habitat (%)', type='l', lwd=2, lty=2, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
errbar(x=c(1.05,2.05,3.05,4.05)+0.05, y=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yplus=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] + biomass_final85$sdNorth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yminus=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] - biomass_final85$sdNorth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], errbar.col='dodgerblue2', lwd=2, cap=0, add=T)
points(north_PA~xAxis_plus, type='p', pch=19, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
par(xpd=T)
points(north_biom~xAxis_plus, type='l', lwd=2, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
errbar(x=c(1.05,2.05,3.05,4.05)+0.05, y=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yplus=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] + biomass_final85$sdNorth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yminus=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] - biomass_final85$sdNorth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], errbar.col='dodgerblue2', lwd=2, cap=0, add=T)
points(north_biom~xAxis_plus, type='p', pch=19, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
abline(a=0, b=0, lwd=.75, col='gray50', lty=2, xpd=F)
axis(2,at=c(0,100,200,300), lwd=1.3, cex.axis=0.95, las=1)
mtext('(e)',line=-1.25,adj=0.01,font=2)
#  
par(xpd=T)
points(north_biom~xAxis, type='l', lwd=2, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
errbar(x=c(1:4), y=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yplus=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] + biomass_final85$sdNorth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yminus=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] - biomass_final85$sdNorth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], errbar.col='gray20', lwd=2, cap=0, add=T)
points(north_biom~xAxis, type='p', pch=19, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
points(north_PA~xAxis, type='l', lwd=2, lty=2, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
errbar(x=c(1:4), y=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yplus=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] + biomass_final85$sdNorth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yminus=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] - biomass_final85$sdNorth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], errbar.col='gray20', lwd=2, cap=0, add=T)
points(north_PA~xAxis, type='p', pch=19, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])

points(north_biom~xAxis_min, type='l', lwd=2, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
errbar(x=c(0.95,1.95,2.95,3.95)-0.05, y=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yplus=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] + biomass_final85$sdNorth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yminus=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] - biomass_final85$sdNorth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], errbar.col='red', lwd=2, cap=0, add=T)
points(north_biom~xAxis_min, type='p', pch=19, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
points(north_PA~xAxis_min, type='l', lwd=2, lty=2, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
errbar(x=c(0.95,1.95,2.95,3.95)-0.05, y=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yplus=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] + biomass_final85$sdNorth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yminus=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] - biomass_final85$sdNorth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], errbar.col='red', lwd=2, cap=0, add=T)
points(north_PA~xAxis_min, type='p', pch=19, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])

# BSB north biomass
i=3
par(xpd=NA)
plot(north_PA~xAxis_plus, xlim=c(0.9,4.1), ylim=c(-45,900), yaxt='n', xaxt='n', xlab='', ylab='Change in habitat (%)', type='l', lwd=2, lty=2, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
errbar(x=c(1.05,2.05,3.05,4.05)+0.05, y=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yplus=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] + biomass_final85$sdNorth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yminus=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] - biomass_final85$sdNorth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], errbar.col='dodgerblue2', lwd=2, cap=0, add=T)
points(north_PA~xAxis_plus, type='p', pch=19, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
par(xpd=T)
points(north_biom~xAxis_plus, type='l', lwd=2, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
errbar(x=c(1.05,2.05,3.05,4.05)+0.05, y=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yplus=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] + biomass_final85$sdNorth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yminus=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] - biomass_final85$sdNorth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], errbar.col='dodgerblue2', lwd=2, cap=0, add=T)
points(north_biom~xAxis_plus, type='p', pch=19, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
abline(a=0, b=0, lwd=.75, col='gray50', lty=2, xpd=F)
par(xpd=NA)
axis(2,at=c(0,200,400,600,800), lwd=1.3, cex.axis=0.95, las=1)
axis(1,at=c(1,2,3,4), labels=F, lwd=1.3, cex.axis=.9, las=1)
text(x=c(1:4+0.34), y=-200, c('2021-2040','2041-2060','2061-2080','2081-2100'), srt = 305, pos = 1, cex=1)
text(x=2.55, y=-386, 'Years', pos = 1, cex=1.2)
mtext('(f)',line=-1.25,adj=0.01,font=2)
#  
par(xpd=T)
points(north_biom~xAxis, type='l', lwd=2, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
errbar(x=c(1:4), y=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yplus=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] + biomass_final85$sdNorth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yminus=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] - biomass_final85$sdNorth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], errbar.col='gray20', lwd=2, cap=0, add=T)
points(north_biom~xAxis, type='p', pch=19, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
points(north_PA~xAxis, type='l', lwd=2, lty=2, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
errbar(x=c(1:4), y=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yplus=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] + biomass_final85$sdNorth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yminus=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] - biomass_final85$sdNorth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], errbar.col='gray20', lwd=2, cap=0, add=T)
points(north_PA~xAxis, type='p', pch=19, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])

points(north_biom~xAxis_min, type='l', lwd=2, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
errbar(x=c(0.95,1.95,2.95,3.95)-0.05, y=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yplus=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] + biomass_final85$sdNorth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yminus=biomass_final85$north_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] - biomass_final85$sdNorth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], errbar.col='red', lwd=2, cap=0, add=T)
points(north_biom~xAxis_min, type='p', pch=19, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
points(north_PA~xAxis_min, type='l', lwd=2, lty=2, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
errbar(x=c(0.95,1.95,2.95,3.95)-0.05, y=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yplus=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] + biomass_final85$sdNorth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yminus=biomass_final85$north_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] - biomass_final85$sdNorth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], errbar.col='red', lwd=2, cap=0, add=T)
points(north_PA~xAxis_min, type='p', pch=19, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])

# SABLEFISH south biomass
i=4
par(xpd=NA)
plot(south_biom~xAxis_plus, xlim=c(0.9,4.1), ylim=c(-88,5), yaxt='n', xaxt='n',xlab='', ylab='', type='l', lwd=2, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
errbar(x=c(1.05,2.05,3.05,4.05)+0.05, y=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yplus=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] + biomass_final85$sdSouth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yminus=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] - biomass_final85$sdSouth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], errbar.col='dodgerblue2', lwd=2, cap=0, add=T)
points(south_biom~xAxis_plus, type='p', pch=19, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
points(south_PA~xAxis_plus, type='l', lwd=2, lty=2, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
errbar(x=c(1.05,2.05,3.05,4.05)+0.05, y=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yplus=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] + biomass_final85$sdSouth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yminus=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] - biomass_final85$sdSouth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], errbar.col='dodgerblue2', lwd=2, cap=0, add=T)
points(south_PA~xAxis_plus, type='p', pch=19, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
abline(a=0, b=0, lwd=.75, col='gray50', lty=2, xpd=F)
axis(2,at=c(-80,-60,-40,-20,0), lwd=1.3, cex.axis=0.95, las=1)
mtext('(g)',line=-1.25,adj=0.01,font=2)
# 
points(south_biom~xAxis, type='l', lwd=2, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
errbar(x=c(1:4), y=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yplus=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] + biomass_final85$sdSouth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yminus=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] - biomass_final85$sdSouth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], errbar.col='gray20', lwd=2, cap=0, add=T)
points(south_biom~xAxis, type='p', pch=19, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
points(south_PA~xAxis, type='l', lwd=2, lty=2, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
errbar(x=c(1:4), y=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yplus=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] + biomass_final85$sdSouth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yminus=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] - biomass_final85$sdSouth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], errbar.col='gray20', lwd=2, cap=0, add=T)
points(south_PA~xAxis, type='p', pch=19, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])

points(south_biom~xAxis_min, type='l', lwd=2, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
errbar(x=c(0.95,1.95,2.95,3.95)-0.05, y=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yplus=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] + biomass_final85$sdSouth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yminus=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] - biomass_final85$sdSouth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], errbar.col='red', lwd=2, cap=0, add=T)
points(south_biom~xAxis_min, type='p', pch=19, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
points(south_PA~xAxis_min, type='l', lwd=2, lty=2, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
errbar(x=c(0.95,1.95,2.95,3.95)-0.05, y=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yplus=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] + biomass_final85$sdSouth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yminus=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] - biomass_final85$sdSouth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], errbar.col='red', lwd=2, cap=0, add=T)
points(south_PA~xAxis_min, type='p', pch=19, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])

# SQUID south biomass
i=5
par(xpd=T)
plot(south_biom~xAxis_plus, xlim=c(0.9,4.1), ylim=c(-75,85), yaxt='n', xaxt='n',xlab='', ylab='', type='l', lwd=2, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
errbar(x=c(1.05,2.05,3.05,4.05)+0.05, y=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yplus=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] + biomass_final85$sdSouth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yminus=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] - biomass_final85$sdSouth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], errbar.col='dodgerblue2', lwd=2, cap=0, add=T)
points(south_biom~xAxis_plus, type='p', pch=19, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
points(south_PA~xAxis_plus, type='l', lwd=2, lty=2, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
errbar(x=c(1.05,2.05,3.05,4.05)+0.05, y=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yplus=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] + biomass_final85$sdSouth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yminus=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] - biomass_final85$sdSouth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], errbar.col='dodgerblue2', lwd=2, cap=0, add=T)
points(south_PA~xAxis_plus, type='p', pch=19, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
abline(a=0, b=0, lwd=.75, col='gray50', lty=2, xpd=F)
par(xpd=NA)
axis(2,at=c(-50,0,50), lwd=1.3, cex.axis=0.95, las=1)
mtext('(h)',line=-1.25,adj=0.01,font=2)
# 
par(xpd=T)
points(south_biom~xAxis, type='l', lwd=2, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
errbar(x=c(1:4), y=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yplus=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] + biomass_final85$sdSouth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yminus=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] - biomass_final85$sdSouth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], errbar.col='gray20', lwd=2, cap=0, add=T)
points(south_biom~xAxis, type='p', pch=19, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
points(south_PA~xAxis, type='l', lwd=2, lty=2, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
errbar(x=c(1:4), y=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yplus=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] + biomass_final85$sdSouth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yminus=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] - biomass_final85$sdSouth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], errbar.col='gray20', lwd=2, cap=0, add=T)
points(south_PA~xAxis, type='p', pch=19, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])

points(south_biom~xAxis_min, type='l', lwd=2, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
errbar(x=c(0.95,1.95,2.95,3.95)-0.05, y=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yplus=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] + biomass_final85$sdSouth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yminus=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] - biomass_final85$sdSouth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], errbar.col='red', lwd=2, cap=0, add=T)
points(south_biom~xAxis_min, type='p', pch=19, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
points(south_PA~xAxis_min, type='l', lwd=2, lty=2, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
errbar(x=c(0.95,1.95,2.95,3.95)-0.05, y=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yplus=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] + biomass_final85$sdSouth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yminus=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] - biomass_final85$sdSouth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], errbar.col='red', lwd=2, cap=0, add=T)
points(south_PA~xAxis_min, type='p', pch=19, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])

# BSB south biomass
i=3
par(xpd=T)
plot(south_biom~xAxis_plus, xlim=c(0.9,4.1), ylim=c(-78,125), yaxt='n', xaxt='n', xlab='', ylab='', type='l', lwd=2, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
errbar(x=c(1.05,2.05,3.05,4.05)+0.05, y=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yplus=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] + biomass_final85$sdSouth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yminus=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] - biomass_final85$sdSouth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], errbar.col='dodgerblue2', lwd=2, cap=0, add=T)
points(south_biom~xAxis_plus, type='p', pch=19, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
points(south_PA~xAxis_plus, type='l', lwd=2, lty=2, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
errbar(x=c(1.05,2.05,3.05,4.05)+0.05, y=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yplus=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] + biomass_final85$sdSouth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], yminus=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'] - biomass_final85$sdSouth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM'], errbar.col='dodgerblue2', lwd=2, cap=0, add=T)
points(south_PA~xAxis_plus, type='p', pch=19, col='dodgerblue2', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GLM',])
abline(a=0, b=0, lwd=.75, col='gray50', lty=2, xpd=F)
par(xpd=NA)
axis(2,at=c(-50,0,50,100), lwd=1.3, cex.axis=.95, las=1)
axis(1,at=c(1,2,3,4), labels=F, lwd=1.3, cex.axis=.9, las=1)
text(x=c(1:4+0.33), y=-111.5, c('2021-2040','2041-2060','2061-2080','2081-2100'), srt = 305, pos = 1, cex=1)
text(x=2.5, y=-151, 'Years', pos = 1, cex=1.2)
mtext('(i)',line=-1.25,adj=0.01,font=2)
#  
par(xpd=T)
points(south_biom~xAxis, type='l', lwd=2, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
errbar(x=c(1:4), y=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yplus=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] + biomass_final85$sdSouth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yminus=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] - biomass_final85$sdSouth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], errbar.col='gray20', lwd=2, cap=0, add=T)
points(south_biom~xAxis, type='p', pch=19, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
points(south_PA~xAxis, type='l', lwd=2, lty=2, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])
errbar(x=c(1:4), y=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yplus=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] + biomass_final85$sdSouth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], yminus=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'] - biomass_final85$sdSouth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM'], errbar.col='gray20', lwd=2, cap=0, add=T)
points(south_PA~xAxis, type='p', pch=19, col='gray20', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='GAM',])

points(south_biom~xAxis_min, type='l', lwd=2, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
errbar(x=c(0.95,1.95,2.95,3.95)-0.05, y=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yplus=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] + biomass_final85$sdSouth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yminus=biomass_final85$south_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] - biomass_final85$sdSouth_biom[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], errbar.col='red', lwd=2, cap=0, add=T)
points(south_biom~xAxis_min, type='p', pch=19, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
points(south_PA~xAxis_min, type='l', lwd=2, lty=2, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])
errbar(x=c(0.95,1.95,2.95,3.95)-0.05, y=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yplus=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] + biomass_final85$sdSouth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], yminus=biomass_final85$south_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'] - biomass_final85$sdSouth_PA[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT'], errbar.col='red', lwd=2, cap=0, add=T)
points(south_PA~xAxis_min, type='p', pch=19, col='red', biomass_final85[biomass_final85$species==SPECIES[i] & biomass_final85$nicheMod=='BRT',])

dev.off()
