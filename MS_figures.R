setwd('/Users/jamesmorley/Documents/project_velocity')

require(Hmisc)
library(latticeExtra)
library(maps)
library(geosphere)
library(plotrix)
library(quantreg)
library(SDMTools)
library(ggplot2)
 
modelrun <- c('bcc-csm1-1-m','bcc-csm1-1','CanESM2','CCSM4','CESM1-CAM5','CNRM-CM5','GFDL-CM3','GFDL-ESM2M','GFDL-ESM2G','GISS-E2-R','GISS-E2-H','IPSL-CM5A-LR','IPSL-CM5A-MR','MIROC-ESM','MPI-ESM-LR','NorESM1-ME')
rcp <- c(26,85)
projfolder = 'output/CEmodels_proj_Nov2017/'
#projfolder = 'temp/'
load('data/speciesProjectionList.RData')

# ====================================================================================
# Make vector with maximum observed depth for each spp.
# ====================================================================================

load('data/hauls_catch_Dec2017.RData') 
minDepth <- rep(NA, length=length(projspp))
maxDepth <- rep(NA, length=length(projspp))

for(i in 1:length(projspp)){
  print(i)
  spp <- projspp[i]  
  spp.data <- dat[dat$sppocean == spp,]
  spp.data <- merge(spp.data, hauls, by='haulid', all.x = T, sort=F) # Add empty hauls
  minDepth[i] <- min(spp.data$depth, na.rm = TRUE)
  maxDepth[i] <- max(spp.data$depth, na.rm = TRUE)
}
save(projspp, regionFreq, minDepth, maxDepth, file='data/MinMaxDepths_spp_Dec2017.RData')
 
abc <- data.frame(spp=projspp, minDepth=minDepth, maxDepth=maxDepth, regionFreq=regionFreq)
# Do most species occur up to 400m on the west coast? 
PacAtl <- grepl('_Pac', abc$spp)
summary(abc[PacAtl,])
# Plot up max depth by region (even though that leaves out some info, could be informative)
ggplot(abc, aes(x=maxDepth)) + geom_histogram(binwidth=20) + facet_wrap(~regionFreq)

# ====================================================================================
# Make some figures of the 'raw' data for a random subset of species
# ====================================================================================

testSpp <- sort(sample(projspp_adj$spp_adj, size=30, replace=F)) # random 50 species to look at
# testSpp <- projspp[138]
for(i in 1:length(testSpp)){
  species = testSpp[i]
  species = 'centropristis striata_Atl'
  speciesMod = gsub('_Gmex', '_Atl', species)
  print(paste("Beginning figures for ", species, sep=''))
  #filename <- paste(projfolder, species, '_rcp', rcp[1], '_jas_prediction_AGG.RData', sep='')
  filename <- paste('output/CEmodels_proj_Nov2017/', speciesMod, '_rcp', rcp[1], '_jas_prediction_AGG.RData', sep='')
  load(filename)
  pred.agg26 <- pred.agg
  
  #filename <- paste('output/CEmodels_proj_May2017/', species, '_rcp', rcp[2], '_jas_prediction_AGG.RData', sep='')
  filename <- paste('output/CEmodels_proj_Nov2017/', speciesMod, '_rcp', rcp[2], '_jas_prediction_AGG.RData', sep='')
  
  load(filename)
  pred.agg85 <- pred.agg
  rm(pred.agg)
  
  if(grepl('_Atl', species)){
    pred.agg26 <- pred.agg26[!pred.agg26$longitude < -82,]
    pred.agg26 <- pred.agg26[!(pred.agg26$longitude < -80.75 & pred.agg26$latitude < 27),]
    pred.agg85 <- pred.agg85[!pred.agg85$longitude < -82,]
    pred.agg85 <- pred.agg85[!(pred.agg85$longitude < -80.75 & pred.agg85$latitude < 27),]
  }
  if(grepl('_Gmex', species)){
    pred.agg26 <- pred.agg26[!pred.agg26$longitude > -80.75,]
    pred.agg26 <- pred.agg26[!(pred.agg26$longitude > -81.5 & pred.agg26$latitude > 28),]
    pred.agg85 <- pred.agg85[!pred.agg85$longitude > -80.75,]
    pred.agg85 <- pred.agg85[!(pred.agg85$longitude > -81.5 & pred.agg85$latitude > 28),]
  }
  
  # Look at the model-level predictions
  pdf(width=14, height=5, file=paste('figures/speciesProjections_Dec2017/', species, '_projections.pdf', sep=''))

  for(j in 4:19){
    print(j-3)
    
    #scale26log = seq(0, max(log(pred.agg26[,j] + 1)), length.out=20)
    scale85log = seq(0, max(log(pred.agg85[,j] + 1)), length.out=20)
    #scale26 = seq(0, max(pred.agg26[,j]), length.out=20)
    scale85 = seq(0, max(pred.agg85[,j]), length.out=20)
    cols = colorRampPalette(colors = c('gray90', 'blue', 'dark blue', 'black'))
    
    # STANDARDIZE PLOT DIMENSIONS FOR ALL FOUR PLOTS
    #xlimit1 = c(min(pred.agg26$longitude[pred.agg26[,j] > scale26[2]/2] - 1),max(pred.agg26$longitude[pred.agg26[,j] > scale26[2]/2] + 1))
    #ylimit1 = c(min(pred.agg26$latitude[pred.agg26[,j] > scale26[2]/2] - .25), max(pred.agg26$latitude[pred.agg26[,j] > scale26[2]/2] + .25))
    #xlimit2 = c(min(pred.agg26$longitude[log(pred.agg26[,j] + 1) > scale26log[2]/2] - 1),max(pred.agg26$longitude[log(pred.agg26[,j] + 1) > scale26log[2]/2] + 1))
    #ylimit2 = c(min(pred.agg26$latitude[log(pred.agg26[,j] + 1) > scale26log[2]/2] - .25), max(pred.agg26$latitude[log(pred.agg26[,j] + 1) > scale26log[2]/2] + .25))
    xlimit3 = c(min(pred.agg85$longitude[pred.agg85[,j] > scale85[2]/2] - 1),max(pred.agg85$longitude[pred.agg85[,j] > scale85[2]/2] + 1))
    ylimit3 = c(min(pred.agg85$latitude[pred.agg85[,j] > scale85[2]/2] - .25), max(pred.agg85$latitude[pred.agg85[,j] > scale85[2]/2] + .25))
    xlimit4 = c(min(pred.agg85$longitude[log(pred.agg85[,j] + 1) > scale85log[2]/2] - 1),max(pred.agg85$longitude[log(pred.agg85[,j] + 1) > scale85log[2]/2] + 1))
    ylimit4 = c(min(pred.agg85$latitude[log(pred.agg85[,j] + 1) > scale85log[2]/2] - .25), max(pred.agg85$latitude[log(pred.agg85[,j] + 1) > scale85log[2]/2] + .25))
    #xlimit = c(min(xlimit1[1],xlimit2[1],xlimit3[1],xlimit4[1]),max(xlimit1[2],xlimit2[2],xlimit3[2],xlimit4[2]))
    #ylimit = c(min(ylimit1[1],ylimit2[1],ylimit3[1],ylimit4[1]),max(ylimit1[2],ylimit2[2],ylimit3[2],ylimit4[2]))
    xlimit = c(min(xlimit3[1],xlimit4[1]),max(xlimit3[2],xlimit4[2]))
    ylimit = c(min(ylimit3[1],ylimit4[1]),max(ylimit3[2],ylimit4[2]))
    
    # RCP26 plots_one for normal, one for Logged+1
    #Projmap <- map('world', xlim=xlimit, ylim=ylimit, plot=F) 
    #Projmap <- data.frame(lon=Projmap$x, lat=Projmap$y)
    #print(levelplot(pred.agg26[,j] ~ longitude*latitude|year_range, data=pred.agg26, xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, main=paste('rcp', rcp[1], 'summer', 'Model', modelrun[j-3], species, sep='_'), 
     #               at = scale26, col.regions=cols, layout = c(5, 1)) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')) 
    
    #Projmap <- map('world', xlim=xlimit, ylim=ylimit, plot=F) 
    #Projmap <- data.frame(lon=Projmap$x, lat=Projmap$y)
    #print(levelplot(log(pred.agg26[,j] + 1) ~ longitude*latitude|year_range, data=pred.agg26, xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, main=paste('rcp', rcp[1], 'summer', 'Model', modelrun[j-3], species, 'LOGGED', sep='_'), 
     #               at = scale26log, col.regions=cols, layout = c(5, 1)) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')) 
     
    # RCP85 plots_one for normal, one for Logged+1
    Projmap <- map('world', xlim=xlimit, ylim=ylimit, plot=F) 
    Projmap <- data.frame(lon=Projmap$x, lat=Projmap$y)
    print(levelplot(pred.agg85[,j] ~ longitude*latitude|year_range, data=pred.agg85, xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, main=paste('rcp', rcp[2], 'summer', 'Model', modelrun[j-3], species, sep='_'), 
                    at = scale85, col.regions=cols, layout = c(5, 1)) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')) 
    
    Projmap <- map('world', xlim=xlimit, ylim=ylimit, plot=F) 
    Projmap <- data.frame(lon=Projmap$x, lat=Projmap$y)
    print(levelplot(log(pred.agg85[,j] + 1) ~ longitude*latitude|year_range, data=pred.agg85, xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, main=paste('rcp', rcp[2], 'summer', 'Model', modelrun[j-3], species, 'LOGGED', sep='_'), 
                    at = scale85log, col.regions=cols, layout = c(5, 1)) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')) 
  }
  dev.off()
}
 
# ============================================================================================================
# Partition species based on GMex as separate body
# ============================================================================================================

# Make two column logic vectors one saying should it be in GMEX and the other for Atlantic side (Pac spp should be both FALSE obviously)
load('data/hauls_catch_Dec2017.RData') 
gmex <- rep(NA, length=length(projspp))
atl <- rep(NA, length=length(projspp))

for(i in 1:length(projspp)){
  print(i)
  spp <- projspp[i]  
  spp.data <- dat[dat$sppocean == spp,]
  spp.data <- merge(spp.data, hauls, by='haulid', all.x = T, sort=F) # Add empty hauls
  regs <- names(table(spp.data$region))
  
  gmex[i] <- ifelse(length(regs[regs=='SEFSC_GOMex'])>0, TRUE, FALSE)
  if(grepl('_Pac', spp)){
    atl[i] <- FALSE
  }else{
    atl[i] <- ifelse(length(regs)>1, TRUE, ifelse(length(regs[regs=='SEFSC_GOMex'])>0, FALSE, TRUE))
  }
}
save(atl, gmex, projspp, regionFreq, file='data/Gmex_divide_spp_Dec2017.RData')
abc <- data.frame(spp=projspp, region=regionFreq, gmex=gmex, atl=atl, stringsAsFactors = FALSE)

# If Gmex is T and Atl is F, the last bit needs to be changed to Gmex_i.e. only project into the Gmex
Gmex_only <- abc[abc$gmex==TRUE & abc$atl==FALSE,]
abc <- abc[!(abc$gmex==TRUE & abc$atl==FALSE),]
abc$spp_adj <- abc$spp
Gmex_only$spp_adj <- gsub('_Atl', '_Gmex', Gmex_only$spp)

# If they are both T, a new row needs to be added with the same spp and Gmex
Gmex_Atl_both <- abc[abc$gmex==TRUE & abc$atl==TRUE,]
Gmex_Atl_both$spp_adj <- gsub('_Atl', '_Gmex', Gmex_Atl_both$spp)

projspp_adj <- rbind(abc, Gmex_only, Gmex_Atl_both)
projspp_adj <- projspp_adj[order(projspp_adj$spp_adj),] 
save(projspp_adj, file='data/projspp_wGmex.RData')

# POTENTIAL ERROR OBSERVATIONS_SPECIES THAT NEED TO BE DROPPED FROM GMEX PROJECTIONS
# get number of presences for each GMEX species in the GMEX_I noticed some species that really shouldn't be there (maybe misIDs)
Gmex_spp <- centroids26$species[grepl('_Gmex', centroids26$species)]
Gmex_spp <- gsub('_Gmex', '_Atl', Gmex_spp)
# Now loop to figure out how many Gmex occurences for each of these species
Gmex_freq <- rep(NA, length=length(Gmex_spp))
for(i in 1:length(Gmex_spp)){
  print(i)
  spp <- Gmex_spp[i]  
  spp.data <- dat[dat$sppocean == spp,]
  spp.data <- merge(spp.data, hauls, by='haulid', all.x = T, sort=F) # Add haul info
  spp.data <- spp.data[spp.data$region=='SEFSC_GOMex',]
  spp.data <- spp.data[!is.na(spp.data$sppocean),]
  Gmex_freq[i] <- nrow(spp.data)
}
Gmex_freq <- data.frame(cbind(spp=Gmex_spp, freq=Gmex_freq), stringsAsFactors = F)

# These species are most likely mis-IDs in the Gmex data set and shouldn't have projections in that region_they were all rare in observations
drops <- c('alosa pseudoharengus_Gmex', 'brevoortia tyrannus_Gmex', 'clupea harengus_Gmex', 'dipturus laevis_Gmex', 'paralichthys dentatus_Gmex',
           'hippoglossoides platessoides_Gmex', 'pseudopleuronectes americanus_Gmex', 'scomber scombrus_Gmex', 'cryptacanthodes maculatus_Gmex',
           'echinarachnius parma_Gmex', 'illex illecebrosus_Gmex', 'melanostigma atlanticum_Gmex', 'menidia menidia_Gmex', 'ovalipes ocellatus_Gmex','placopecten magellanicus_Gmex')

projspp_adj <- projspp_adj[!(projspp_adj$spp_adj %in% drops), ]
rownames(projspp_adj) <- NULL
save(projspp_adj, file='data/projspp_wGmex.RData')

# ============================================================================================================
# Compare centroid shifts with RCPs and logged vs. non-logged_b/c the visualizations looks better when logged
# Ideally, we don't want to log the centroid calculations
# ============================================================================================================
 
# testSpp <- 'sebastes babcocki_Pac'
for(i in 1:nrow(projspp_adj)){
  species = testSpp[i]
  print(paste("Beginning figures for ", species, sep=''))
  filename <- paste('output/CEmodels_proj_Nov2017/', species, '_rcp', rcp[1], '_jas_prediction_AGG.RData', sep='')
  load(filename)
  pred.agg26 <- pred.agg
  pred.agg26pre <- pred.agg26[pred.agg26$year_range == '2007-2020',]
  pred.agg26post <- pred.agg26[pred.agg26$year_range == '2081-2100',]
  
  filename <- paste('output/CEmodels_proj_Nov2017/', species, '_rcp', rcp[2], '_jas_prediction_AGG.RData', sep='')
  load(filename)
  pred.agg85 <- pred.agg
  pred.agg85pre <- pred.agg85[pred.agg85$year_range == '2007-2020',]
  pred.agg85post <- pred.agg85[pred.agg85$year_range == '2081-2100',]
  rm(pred.agg)
  
  pdf(width=7, height=5, file=paste('figures/speciesProjections_Dec2017/', species, '_centroid_projections.pdf', sep=''))
  
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
  
# FIRST create an area column to correct for the biomass weighting due to converging longitude (i.e. convert predicted cpue to biomass in a cell)
library(raster)
load('data/ProjectionBathGrid_Feb27_2017.RData')# load bathymetry projection grid

grid_rast <- raster(ymn=min(proj.grid$latBathgrid), ymx=max(proj.grid$latBathgrid), xmn=min(proj.grid$lonBathgrid), xmx=max(proj.grid$lonBathgrid), resolution=1/20)
grid_rast[] <- area(grid_rast)[] # in km2_PROB SHOULD CONVERT THIS TO HECTARES? Whatever is most common
grid_area <- data.frame(rasterToPoints(grid_rast), stringsAsFactors = FALSE) # all the lat/lon expanded and all three values in the rasterStack included
grid_area <- unique(data.frame(lat = grid_area$y, area = grid_area$layer))

# now get the latitude values to match up as sig-figs altered in raster
lats <- unique(proj.grid$latBathgrid)
grid_area$latitude<- NA

for(i in 1:nrow(grid_area)){
  print(i)
  lat <- grid_area$lat[i]
  lat_adj <- data.frame(lats=lats, diff=lats-lat)
  lat_adj$min <- abs(lat_adj$diff)
  grid_area$latitude[i] <- lat_adj$lats[lat_adj$min==min(lat_adj$min)][1]
}
grid_area$lat <- NULL
rownames(grid_area) <- NULL
grid_area$latitude <- round(grid_area$latitude, digits=3)# to prevent rounding error in the loop
grid_area <- rbind(grid_area, c(14.413, 62.125))# this latitude row got knocked off in raster
save(grid_area, file='data/lat_area_corrections.RData')

# columns for init and final lat and lon, mean + SD in latitudinal change, mean + SD in absolute shift distance
centroids26 <- data.frame(latPre=numeric(), lonPre=numeric(), latPost=numeric(), lonPost=numeric(), meanLat=numeric(), sdLat=numeric(), meanDist=numeric(), sdDist=numeric(), radius=numeric())
centroids85 <- data.frame(latPre=numeric(), lonPre=numeric(), latPost=numeric(), lonPost=numeric(), meanLat=numeric(), sdLat=numeric(), meanDist=numeric(), sdDist=numeric(), radius=numeric())
options(warn=0) 
 
for(i in 1:nrow(projspp_adj)){ # one loop for each species
  species <- projspp_adj$spp_adj[i]
  print(paste("Beginning figures for ", species, '_spp', i, sep=''))
  filename <- paste(projfolder, projspp_adj$spp[i], '_rcp', rcp[1], '_jas_prediction_AGG.RData', sep='')
  load(filename)
  pred.agg26 <- pred.agg
   
  filename <- paste(projfolder, projspp_adj$spp[i], '_rcp', rcp[2], '_jas_prediction_AGG.RData', sep='')
  load(filename)
  pred.agg85 <- pred.agg
  rm(pred.agg)
  
  # Need to trim the grid based on if '_Atl' or '_Gmex' is in the adjusted name  
  if(grepl('_Atl', species)){
    pred.agg26 <- pred.agg26[!pred.agg26$longitude < -82,]
    pred.agg26 <- pred.agg26[!(pred.agg26$longitude < -80.75 & pred.agg26$latitude < 27),]
    pred.agg85 <- pred.agg85[!pred.agg85$longitude < -82,]
    pred.agg85 <- pred.agg85[!(pred.agg85$longitude < -80.75 & pred.agg85$latitude < 27),]
  }
  if(grepl('_Gmex', species)){
    pred.agg26 <- pred.agg26[!pred.agg26$longitude > -80.75,]
    pred.agg26 <- pred.agg26[!(pred.agg26$longitude > -81.5 & pred.agg26$latitude > 28),]
    pred.agg85 <- pred.agg85[!pred.agg85$longitude > -80.75,]
    pred.agg85 <- pred.agg85[!(pred.agg85$longitude > -81.5 & pred.agg85$latitude > 28),]
  }
  
  # Need to correct some rounding issues before merging grid_area
  pred.agg26$latitude <- round(pred.agg26$latitude, digits=3)
  pred.agg85$latitude <- round(pred.agg85$latitude, digits=3)
  # Add the grid area corrections
  pred.agg26 <- merge(pred.agg26, grid_area, by='latitude', all.x=T)
  pred.agg85 <- merge(pred.agg85, grid_area, by='latitude', all.x=T)
  #levelplot(area ~ longitude*latitude, data=pred.agg85[pred.agg85$year_range=='2081-2100',])
  
  # RCP26 calculations
  yearRange <- c('2007-2020', '2081-2100')
  cent_preds26lat <- matrix(data=NA, nrow=2, ncol=16)
  cent_preds26lon <- matrix(data=NA, nrow=2, ncol=16)
  for(j in 1:2){
    years = yearRange[j]
    preds <- pred.agg26[pred.agg26$year_range == years,]
    for(k in 4:19){
      weights_adj <- preds[,k] * (preds$area) #*100) # the '*100' is in case I want to do hectares_it shouldn't matter
      valueLat = wtd.mean(preds$latitude, weights = weights_adj, na.rm=FALSE) 
      valueLon = wtd.mean(preds$longitude, weights = weights_adj, na.rm=FALSE) 
      if(is.na(valueLat)){
        print('Problem')
      }
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
  #plot(xy26rads[,2]~xy26rads[,1], xlim=c(-2,2), ylim=c(-2,2))
  xyMean <- c(mean(xy26rads[,1]),mean(xy26rads[,2]))
  radius26 <- sqrt(xyMean[1]^2 + xyMean[2]^2)
  
  # RCP85 calculations
  cent_preds85lat <- matrix(data=NA, nrow=2, ncol=16)
  cent_preds85lon <- matrix(data=NA, nrow=2, ncol=16)
  for(j in 1:2){
    years = yearRange[j]
    preds <- pred.agg85[pred.agg85$year_range == years,]
    for(k in 4:19){
      weights_adj <- preds[,k] * (preds$area) #*100) # the '*100' is in case I want to do hectares_it shouldn't matter
      valueLat = wtd.mean(preds$latitude, weights = weights_adj, na.rm=FALSE) 
      valueLon = wtd.mean(preds$longitude, weights = weights_adj, na.rm=FALSE) 
      if(is.na(valueLat)){
        print('Problem')
      }
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
  #plot(xy85rads[,2]~xy85rads[,1], xlim=c(-2,2), ylim=c(-2,2))
  xyMean <- c(mean(xy85rads[,1]),mean(xy85rads[,2]))
  radius85 <- sqrt(xyMean[1]^2 + xyMean[2]^2)
    
  centroids26[i,] <- data.frame(latPre=mean_Lat_pre26, lonPre=mean_Lon_pre26, latPost=mean_Lat_post26, lonPost=mean_Lon_post26, meanLat=meanLat26, sdLat=sdLat26, meanDist=meanDis26, sdDist=sdDis26, radius=radius26)
  centroids85[i,] <- data.frame(latPre=mean_Lat_pre85, lonPre=mean_Lon_pre85, latPost=mean_Lat_post85, lonPost=mean_Lon_post85, meanLat=meanLat85, sdLat=sdLat85, meanDist=meanDis85, sdDist=sdDis85, radius=radius85)
}

centroids85$species <- projspp_adj$spp_adj
centroids26$species <- projspp_adj$spp_adj
save(centroids26, centroids85, file='output/centroid_summaries_Dec2017.RData')

# GROUPS SPECIES INTO REGIONS FOR PLOTTING_I don't use this method anymore, but the following one instead
load('data/hauls_catch_Dec2017.RData') 
haulsMod <- hauls
haulsMod$regionGroup <- ifelse(haulsMod$region=='AFSC_WCTri', 'NWFSC_WCAnn', haulsMod$region)
haulsMod$regionGroup <- ifelse(haulsMod$regionGroup=='VIMS_NEAMAP', 'NEFSC_NEUS', haulsMod$regionGroup)
haulsMod$regionGroup <- ifelse(haulsMod$regionGroup=='DFO_SoGulf', 'DFO_ScotianShelf', haulsMod$regionGroup)

regions <- character(length=length(centroids85$species))
for(i in 1:length(regions)){
  if(grepl('_Gmex', centroids85$species[i])){
    regions[i] <- 'SEFSC_GOMex'
  }else{ 
    spdata <- dat[dat$sppocean==centroids85$species[i],]
    spdata <- merge(spdata, haulsMod, by='haulid', all.x=T)
    regs <- data.frame(table(spdata$regionGroup))
    abc <- as.character(regs$Var1[regs$Freq == max(regs$Freq)[1]])
    if(abc == 'SEFSC_GOMex'){
      regs <- regs[!regs$Var1=='SEFSC_GOMex',]
      abc <- as.character(regs$Var1[regs$Freq == max(regs$Freq)[1]])
    }
    regions[i] <- abc
  }
}
 
centroids26$region <- regions
centroids85$region <- regions


# SEEMS BETTER TO GROUP INTO REGIONS BASED ON INITIAL CENTROID LOCATION ===============
centroids85$regionsB <- ifelse(grepl('_Atl', centroids85$species) & centroids85$latPre < 35.5, 'SCDNR_SEUS', centroids85$region) 
centroids85$regionsB <- ifelse((grepl('_Atl', centroids85$species) & centroids85$latPre > 35.5 & centroids85$latPre < 48.2 & centroids85$lonPre < -66.2), 'NEFSC_NEUS', centroids85$regionsB) 
centroids85$regionsB <- ifelse((grepl('_Atl', centroids85$species) & (centroids85$latPre > 48.2 | centroids85$lonPre > -66.2)), 'DFO_Newfoundland', centroids85$regionsB)
centroids85$regionsB <- ifelse(grepl('_Pac', centroids85$species) & centroids85$latPre < 49, 'NWFSC_WCAnn', centroids85$regionsB)
centroids85$regionsB <- ifelse((grepl('_Pac', centroids85$species) & centroids85$latPre > 49 & centroids85$lonPre > -159), 'AFSC_GOA', centroids85$regionsB)
centroids85$regionsB <- ifelse(grepl('_Pac', centroids85$species) & centroids85$lonPre < -159, 'AFSC_EBS', centroids85$regionsB)
centroids26$regionsB <- centroids85$regionsB

# Plot this up and compare with previous run ===========================================
# add colors for each region_useful for some plotting diagnostics 
regs <- names(table(centroids85$regionsB))
regs <- data.frame(regionsB=regs, colors=c('blue', 'green', 'blue', 'green', 'red', 'black', 'red'), stringsAsFactors = F) 
centroids85 <- merge(centroids85, regs, by='regionsB', all.x=T)
centroids26 <- merge(centroids26, regs, by='regionsB', all.x=T)

load('data/ProjectionBathGrid_Feb27_2017.RData')# load bathymetry projection grid
Projmap <- map('world', xlim=c(min(centroids26$lonPre),max(centroids26$lonPre)), ylim=c(min(centroids26$latPre),max(centroids26$latPost)),plot=F) 
Projmap <- data.frame(lon=Projmap$x, lat=Projmap$y)

pdf(width=14, height=8, file=paste('figures/', 'species_centroid_shifts85_Dec2017.pdf', sep=''))
plot(latClimgrid~lonClimgrid, col='gray90', cex=.7, pch=16, xlab='Longitude', ylab='Latitude', main='Predicted centroid shifts for 686 species_RCP26', data=proj.grid) 
points(lat ~ lon, Projmap, type='l', lty=1, lwd=.7, col='gray50')
#points(latPre ~ lonPre, cex=.3, pch=20, col=centroids85$colors, centroids85)
arrows(x0=centroids85$lonPre, y0=centroids85$latPre, x1=centroids85$lonPost, y1=centroids85$latPost, length=.05, lwd=.8, col='dark blue')
dev.off()

save(centroids26, centroids85, file='output/centroid_summaries_Dec2017.RData')

# ===============================================================================================
# Uncertainty analysis for centroids=============================================================
# ===============================================================================================

# Uncertainty method #1_aggreement with projected distance and angle among models
library('car')
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
centroids26$uncerCat2 <- ifelse(log(centroids26$sdLat) > (coef(quan95.26)[1] + coef(quan95.26)[2]*log(abs(centroids26$meanLat))), 'high uncer.',
              ifelse(log(centroids26$sdLat) < (coef(quan75.26)[1] + coef(quan75.26)[2]*log(abs(centroids26$meanLat))), 'low uncer.', 'med uncer.'))

lm85 <- lm(log(sdLat)~log(abs(meanLat)), data=centroids85)
qqPlot(lm85)
hist(resid(lm85))
quan95.85 <- rq(log(sdLat)~log(abs(meanLat)), tau=.95, data=centroids85)
quan75.85 <- rq(log(sdLat)~log(abs(meanLat)), tau=.75, data=centroids85)
centroids85$uncerCat2 <- ifelse(log(centroids85$sdLat) > (coef(quan95.85)[1] + coef(quan95.85)[2]*log(abs(centroids85$meanLat))), 'high uncer.',
              ifelse(log(centroids85$sdLat) < (coef(quan75.85)[1] + coef(quan75.85)[2]*log(abs(centroids85$meanLat))), 'low uncer.', 'med uncer.'))

centroids26$uncerCat <- ifelse((centroids26$radius < rad26[2] | centroids26$resid > sd26[2]), "high unc.", 
              ifelse((centroids26$radius > rad26[1] & centroids26$resid < sd26[1]), "low unc.", "med unc."))  
centroids85$uncerCat <- ifelse((centroids85$radius < rad85[2] | centroids85$resid > sd85[2]), "high unc.", 
              ifelse((centroids85$radius > rad85[1] & centroids85$resid < sd85[1]), "low unc.", "med unc."))  

centroids26$uncertainty <- ifelse(centroids26$uncerCat == 'high unc.' | centroids26$uncerCat2 == 'high uncer.', 'high', 
                                  ifelse(centroids26$uncerCat == 'med unc.' | centroids26$uncerCat2 == 'med uncer.', 'medium', 'low'))
centroids85$uncertainty <- ifelse(centroids85$uncerCat == 'high unc.' | centroids85$uncerCat2 == 'high uncer.', 'high', 
                                  ifelse(centroids85$uncerCat == 'med unc.' | centroids85$uncerCat2 == 'med uncer.', 'medium', 'low'))

save(centroids26, centroids85, file='output/centroid_summaries_Dec2017.RData')
  
# The plot of the above uncertainty stuff ===================================

pdf(width=7, height=7, file=paste('figures/', 'uncertainty_analysis_Dec2017.pdf', sep=''))
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
legend(x=.23, y=2.1, xjust=.5, cex=1.4, legend=c('low', 'medium', 'high'), col=c('blue', 'darkorange', 'red2'), pt.cex=1, x.intersp=.75, text.col = "black", bg = "white", bty='n', inset=c(-0.2,0), pch=19)

plot(log(sdLat)~log(abs(meanLat)), xlab='Absolute value of mean shift (logged lat)', ylab='Model uncertainty (logged st. dev.)', 
     col=ifelse(centroids26$uncerCat2 == 'low uncer.', 'blue', ifelse(centroids26$uncerCat2 == 'med uncer.', 'darkorange' ,'red2')), data=centroids26)
abline(coef(lm26))
abline(coef(quan95.26), col='red2')
abline(coef(quan75.26), col='blue')
mtext('  E', side=3, line=-1.5, adj=0)

plot(log(sdLat)~log(abs(meanLat)), xlab='Absolute value of mean shift (logged lat)', ylab='Model uncertainty (logged st. dev.)', 
     col=ifelse(centroids85$uncerCat2 == 'low uncer.', 'blue', ifelse(centroids85$uncerCat2 == 'med uncer.', 'darkorange' ,'red2')), data=centroids85)
abline(coef(lm85))
abline(coef(quan95.85), col='red2')
abline(coef(quan75.85), col='blue')
mtext('  F', side=3, line=-1.5, adj=0)

dev.off()
  
# ============================================================================================
# SHOW THE ONE BIG MAP WITH ALL THE ACTUAL VECTORS IN PLACE
# Could do a chi-squared table to see if the two metrics are related to e/o
# ============================================================================================
   
load('data/ProjectionBathGrid_Feb27_2017.RData')# load bathymetry projection grid

Projmap <- map('world', xlim=c(min(centroids26$lonPre),max(centroids26$lonPre)), ylim=c(min(centroids26$latPre),max(centroids26$latPost)),plot=F) 
Projmap <- data.frame(lon=Projmap$x, lat=Projmap$y)
#Projmap <- Projmap[Projmap$lon > -110,]
pdf(width=14, height=8, file=paste('figures/', 'species_centroid_shifts26_uncerCol_Dec2017.pdf', sep=''))
plot(latClimgrid~lonClimgrid, col='gray90', cex=.7, pch=16, xlab='Longitude', ylab='Latitude', main='Predicted centroid shifts for 686 species_RCP26', data=proj.grid) 
points(lat ~ lon, Projmap, type='l', lty=1, lwd=.7, col='gray50')
high26 <- centroids26[centroids26$uncertainty == 'high',]
med26 <- centroids26[centroids26$uncertainty == 'medium',]
low26 <- centroids26[centroids26$uncertainty == 'low',]
arrows(x0=high26$lonPre, y0=high26$latPre, x1=high26$lonPost, y1=high26$latPost, length=.05, lwd=.8, col='red2')
arrows(x0=med26$lonPre, y0=med26$latPre, x1=med26$lonPost, y1=med26$latPost, length=.05, lwd=.8, col='darkorange')
arrows(x0=low26$lonPre, y0=low26$latPre, x1=low26$lonPost, y1=low26$latPost, length=.05, lwd=.8, col='blue')
dev.off()

pdf(width=14, height=8, file=paste('figures/', 'species_centroid_shifts85_LOWuncer_Dec2017.pdf', sep=''))
plot(latClimgrid~lonClimgrid, col='gray90', cex=.7, pch=16, xlab='Longitude', ylab='Latitude', main='Predicted centroid shifts for 686 species_RCP85', data=proj.grid) 
points(lat ~ lon, Projmap, type='l', lty=1, lwd=.7, col='gray50')
#low85 <- centroids85[centroids85$uncertainty == 'low' & centroids85$lonPre > -110,]
low85 <- centroids85[centroids85$uncertainty == 'low',]
med85 <- centroids85[centroids85$uncertainty == 'medium',]
high85 <- centroids85[centroids85$uncertainty == 'high',]
arrows(x0=high85$lonPre, y0=high85$latPre, x1=high85$lonPost, y1=high85$latPost, length=.05, lwd=.8, col='red2')
arrows(x0=med85$lonPre, y0=med85$latPre, x1=med85$lonPost, y1=med85$latPost, length=.05, lwd=.8, col='darkorange')
arrows(x0=low85$lonPre, y0=low85$latPre, x1=low85$lonPost, y1=low85$latPost, length=.05, lwd=.8, col='blue')
dev.off()

 
# ===================================================================================================
# COMPASS PLOTS=======================================================================================
# ===================================================================================================

centroids26$bearing <- bearing(p1=matrix(cbind(centroids26$lonPre, centroids26$latPre), ncol=2), p2=matrix(cbind(centroids26$lonPost, centroids26$latPost), ncol=2))
centroids26$bearing_scaled <- ifelse(centroids26$bearing < 0, centroids26$bearing + 360, centroids26$bearing)

centroids85$bearing <- bearing(p1=matrix(cbind(centroids85$lonPre, centroids85$latPre), ncol=2), p2=matrix(cbind(centroids85$lonPost, centroids85$latPost), ncol=2))
centroids85$bearing_scaled <- ifelse(centroids85$bearing < 0, centroids85$bearing + 360, centroids85$bearing)

#Make the a tiff ==============================================================================
regions <- names(table(centroids85$regionsB))
load('output/modeldiag_Nov2017_fitallreg_2017.Rdata')
npres <- data.frame(species_orig=modeldiag$sppocean, npres=modeldiag$npres, stringsAsFactors = F)
centroids26$species_orig <- gsub('_Gmex', '_Atl', centroids26$species)
centroids85$species_orig <- gsub('_Gmex', '_Atl', centroids85$species)

centroids26 <- merge(centroids26, npres, by='species_orig', all.x=T)
centroids85 <- merge(centroids85, npres, by='species_orig', all.x=T)

library(car)
plot(centroids26$npres~as.factor(centroids26$uncertainty))
summary(centroids26$npres[centroids26$uncertainty=='high'])
abc <- aov(npres~uncertainty, data=centroids26)
Anova(lm(npres~uncertainty, data=centroids26), type=3)
TukeyHSD(abc, which=c('uncertainty'))

plot(centroids85$npres~as.factor(centroids85$uncertainty))
summary(centroids85$npres[centroids85$uncertainty=='high'])
abc <- aov(npres~uncertainty, data=centroids85)
Anova(lm(npres~uncertainty, data=centroids85), type=3)

centroids26$colors <- centroids85$colors <- NULL  
centroids26$color <- ifelse(centroids26$uncertainty == 'high', 'red2', ifelse(centroids26$uncertainty == 'medium', 'darkorange', 'blue'))
centroids85$color <- ifelse(centroids85$uncertainty == 'high', 'red2', ifelse(centroids85$uncertainty == 'medium', 'darkorange', 'blue'))
regions <- names(table(centroids85$regionsB))

pdf(width=8, height=4, file=paste('figures/', 'compass_plots_Dec2017', '.pdf', sep=''))
par(mfrow=c(1,2), cex.main=.75)
 
for(i in 1:7){
  region = regions[i]
  data26 <- centroids26[centroids26$regionsB == region,]
  data85 <- centroids85[centroids85$regionsB == region,]
  if(region %in% c('AFSC_EBS', 'AFSC_GOA', 'NWFSC_WCAnn')){
  polar.plot(lengths=data26$meanDist, polar.pos=data26$bearing_scaled, radial.lim=c(0,1500), grid.col='gray40' , radial.labels=c('','',''), lwd=2.5, line.col=data26$color, labels=c(), label.pos=c(0,90,180,270), boxed.radial=F, start=90, clockwise=T, main=paste(region, 'rcp26',sep='_'))
  polar.plot(lengths=data85$meanDist, polar.pos=data85$bearing_scaled, radial.lim=c(0,1500), grid.col='gray40' , radial.labels=c('','',''), lwd=2.5, line.col=data85$color, labels=c(), label.pos=c(0,90,180,270), boxed.radial=F, start=90, clockwise=T, main=paste(region, 'rcp85',sep='_'))
  #polar.plot(lengths=data26$meanDist, polar.pos=data26$bearing_scaled, radial.lim=c(0,1500), grid.col='gray80' , radial.labels=c('','',''), lwd=1.5, line.col=data26$color2, labels=c(), label.pos=c(0,90,180,270), boxed.radial=F, start=90, clockwise=T, main=paste(region, 'rcp26_B',sep='_'))
  #polar.plot(lengths=data85$meanDist, polar.pos=data85$bearing_scaled, radial.lim=c(0,1500), grid.col='gray80' , radial.labels=c('','',''), lwd=1.5, line.col=data85$color2, labels=c(), label.pos=c(0,90,180,270), boxed.radial=F, start=90, clockwise=T, main=paste(region, 'rcp85_B',sep='_'))
  } else{
    polar.plot(lengths=data26$meanDist, polar.pos=data26$bearing_scaled, radial.lim=c(0,800), grid.col='gray40' , radial.labels=c('','','',''), lwd=2.5, line.col=data26$color, labels=c(), label.pos=c(0,90,180,270), boxed.radial=F, start=90, clockwise=T, main=paste(region, 'rcp26',sep='_'))
    polar.plot(lengths=data85$meanDist, polar.pos=data85$bearing_scaled, radial.lim=c(0,800), grid.col='gray40' , radial.labels=c('','','',''), lwd=2.5, line.col=data85$color, labels=c(), label.pos=c(0,90,180,270), boxed.radial=F, start=90, clockwise=T, main=paste(region, 'rcp85',sep='_'))
  }
}
dev.off()

load('ProjectionBathGrid_Feb27_2017.RData')# load bathymetry projection grid

# Make the background projection grid
pdf(width=6, height=6, file=paste('figures/', 'Map_template_East.pdf', sep=''))
par(xpd=F, tcl=-.2, mgp=c(1.8,.5,0), cex.axis=.8, las=1) 
Projmap <- map('world', xlim=c(-96,-45.5), ylim=c(26,60), plot=F) 
Projmap <- data.frame(lon=Projmap$x, lat=Projmap$y)
plot(latClimgrid~lonClimgrid, col='gray93', cex=.4, pch=16, xlim=c(-96,-45.5), ylim=c(26,60), xlab='Longitude', ylab='Latitude', data=proj.grid) 
points(lat ~ lon, Projmap, type='l', lty=1, lwd=.7, col='gray50')
dev.off()

pdf(width=6, height=6, file=paste('figures/', 'Map_template_West.pdf', sep=''))
par(xpd=F, tcl=-.2, mgp=c(1.8,.5,0), cex.axis=.8, las=1) 
Projmap <- map('world', xlim=c(-188,-117), ylim=c(32.5,62), plot=F) 
Projmap <- data.frame(lon=Projmap$x, lat=Projmap$y)
plot(latClimgrid~lonClimgrid, col='gray93', cex=.4, pch=16, xlim=c(-188,-117), ylim=c(32.5,62), xlab='Longitude', ylab='Latitude', data=proj.grid) 
points(lat ~ lon, Projmap, type='l', lty=1, lwd=.7, col='gray50')
dev.off()

# make a legend
pdf(width=8, height=4, file=paste('figures/', 'compass_plots_LEGEND', '.pdf', sep=''))
par(mfrow=c(1,2), cex.lab=1.1)
polar.plot(lengths=0, polar.pos=0, radial.lim=c(0,1500), grid.col='gray40' , radial.labels=c('0','500','1000','1500 km'), labels=c(), show.grid.labels=1, label.pos=c(0,90,180,270), boxed.radial=T, start=90, clockwise=T, main='Legend')
polar.plot(lengths=0, polar.pos=0, radial.lim=c(0,800), grid.col='gray40' , radial.labels=c('0','200','400','600', '800 km'), labels=c(), show.grid.labels=1, label.pos=c(0,90,180,270), boxed.radial=T, start=90, clockwise=T, main='Legend')
dev.off()

# Omit species with high or medium uncertainty with either metric and do a mean mag and direction within regions
 
trimmed26 <- centroids26[centroids26$uncertainty == 'low',] 
trimmed85 <- centroids85[centroids85$uncertainty == 'low',] 
  
compSumm <- data.frame(dis26=numeric(), direc26=numeric(), dis85=numeric(), direc85=numeric())
for(i in 1:7){
  region = regions[i]
  data26 <- trimmed26[trimmed26$regionsB == region,]
  data85 <- trimmed85[trimmed85$regionsB == region,]
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
write.csv(compSumm, file='/Users/jamesmorley/Desktop/Assemblage_shifts.csv') # For Pew

compSumm$colors <- c('darkorange', 'gray60', 'black', 'cyan3', 'brown', 'blue', 'red2')

# Need to rearrange the dataframe so that the colors lay right and show up better (figured this after much plotting trials)
abc <- c("DFO_Newfoundland","NEFSC_NEUS","NWFSC_WCAnn","SCDNR_SEUS","AFSC_EBS","AFSC_GOA","SEFSC_GOMex")     
compSummB <- compSumm[match(abc, compSumm$region),]
legend.names <- c("E. Canada","Northeast U.S.","West U.S.","Southeast U.S.","E. Bering Sea","G. Alaska-W. Canada","Gulf of Mexico")     

pdf(file='figures/compass_summary_Dec2017.pdf', width = 7, height = 4)
par(mfrow=c(1,2), cex.lab=1, cex.main=1.2, mar=c(.1,.1,.1,.1), oma=c(.1,.1,.1,.1))
polar.plot(lengths=compSummB$dis26, polar.pos=compSummB$direc26, radial.lim=c(0,1000), mar=c(.1,.5,.1,.5), grid.col='gray75', radial.labels=c('','200','400','600','800', '1000 km'), show.grid.labels=1, lwd=6, line.col=compSummB$colors, labels=c(), label.pos=c(0,90,180,270), boxed.radial=T, start=90, clockwise=T)
title('A) RCP 2.6', line=-1.1, adj=0, cex=.4)

polar.plot(lengths=compSummB$dis85, polar.pos=compSummB$direc85, radial.lim=c(0,1000), mar=c(.1,.5,.1,.5),grid.col='gray75', radial.labels=c('','','','','', ''), show.grid.labels=1, lwd=6, line.col=compSummB$colors, labels=c(), label.pos=c(0,90,180,270), boxed.radial=T, start=90, clockwise=T)
title('B) RCP 8.5', line=-1.1, adj=0)
legend(x=-545, y=-90, xjust=.5, cex=.8, legend=legend.names, col=compSummB$colors, lty=1, lwd=4, pt.cex=0, x.intersp=.75, text.col = "black", bg = "white", box.lwd=.5, box.col='gray75', inset=c(-0.2,0), pch=c(1,3))
dev.off()
    
save(projspp_adj, centroids26, centroids85, file='output/centroid_summaries_Dec2017.RData')
  
# ============================================================================================================
# ============================================================================================================
# Thermal biomass analysis
# ============================================================================================================
# ============================================================================================================
  
# columns for init and final; mean + SD in biomass change
biomass26 <- data.frame(biomMeanChan=numeric(), sdChan=numeric())
biomass85 <- data.frame(biomMeanChan=numeric(), sdChan=numeric())
  
for(i in 1:nrow(projspp_adj)){ # one loop for each species
  species = projspp_adj$spp_adj[i]
  print(paste("Beginning figures for ", species, '_spp', i, sep=''))
  filename <- paste(projfolder, projspp_adj$spp[i], '_rcp', rcp[1], '_jas_prediction_AGG.RData', sep='')
  load(filename)
  pred.agg26 <- pred.agg
  
  filename <- paste(projfolder, projspp_adj$spp[i], '_rcp', rcp[2], '_jas_prediction_AGG.RData', sep='')
  load(filename)
  pred.agg85 <- pred.agg
  rm(pred.agg)
  
  # Need to trim the grid based on if '_Atl' or '_Gmex' is in the adjusted name  
  if(grepl('_Atl', species)){
    pred.agg26 <- pred.agg26[!pred.agg26$longitude < -82,]
    pred.agg26 <- pred.agg26[!(pred.agg26$longitude < -80.75 & pred.agg26$latitude < 27),]
    pred.agg85 <- pred.agg85[!pred.agg85$longitude < -82,]
    pred.agg85 <- pred.agg85[!(pred.agg85$longitude < -80.75 & pred.agg85$latitude < 27),]
  }
  if(grepl('_Gmex', species)){
    pred.agg26 <- pred.agg26[!pred.agg26$longitude > -80.75,]
    pred.agg26 <- pred.agg26[!(pred.agg26$longitude > -81.5 & pred.agg26$latitude > 28),]
    pred.agg85 <- pred.agg85[!pred.agg85$longitude > -80.75,]
    pred.agg85 <- pred.agg85[!(pred.agg85$longitude > -81.5 & pred.agg85$latitude > 28),]
  }
  
  # Need to correct some rounding issues before merging grid_area
  pred.agg26$latitude <- round(pred.agg26$latitude, digits=3)
  pred.agg85$latitude <- round(pred.agg85$latitude, digits=3)
  # Add the grid area corrections
  pred.agg26 <- merge(pred.agg26, grid_area, by='latitude', all.x=T)
  pred.agg85 <- merge(pred.agg85, grid_area, by='latitude', all.x=T)
  
  # RCP26 calculations
  yearRange <- c('2007-2020', '2081-2100')
  biom_preds26 <- matrix(data=NA, nrow=2, ncol=16)
  for(j in 1:2){
    years = yearRange[j]
    preds <- pred.agg26[pred.agg26$year_range == years,]
    for(k in 4:19){
      valueBiom = sum(preds[,k]*preds$area, na.rm=FALSE) 
      if(is.na(valueBiom)){
        print('Problem')
      }
      biom_preds26[j,k-3] = valueBiom
    }
  } 
  # Need to calculate change in biomass for each model_then get the mean/SD of those values
  biom_diff26 <- apply(biom_preds26, 2, FUN=diff, na.rm=FALSE)
  biom_perc26 <- (biom_diff26/biom_preds26[1,]) *100
  meanBiom26 <- mean(biom_perc26)
  sdBiom26 <- sd(biom_perc26)

  # RCP85 calculations
 # pred.agg85 <- pred.agg85[pred.agg85$latitude < 45 & pred.agg85$longitude < -66,]
  biom_preds85 <- matrix(data=NA, nrow=2, ncol=16)
  for(j in 1:2){
    years = yearRange[j]
    preds <- pred.agg85[pred.agg85$year_range == years,]
    for(k in 4:19){
      valueBiom = sum(preds[,k]*preds$area, na.rm=FALSE) 
      if(is.na(valueBiom)){
        print('Problem')
      }
      biom_preds85[j,k-3] = valueBiom
    }
  }  
  # Need to calculate change in biomass for each model_then get the mean/SD of those values
  biom_diff85 <- apply(biom_preds85, 2, FUN=diff, na.rm=FALSE)
  biom_perc85 <- (biom_diff85/biom_preds85[1,])*100
  meanBiom85 <- mean(biom_perc85)
  sdBiom85 <- sd(biom_perc85)
  
  biomass26[i,] <- data.frame(biomMeanChan=meanBiom26, sdChan=sdBiom26)
  biomass85[i,] <- data.frame(biomMeanChan=meanBiom85, sdChan=sdBiom85)
}
 
biomass26$species <- projspp_adj$spp_adj
biomass85$species <- projspp_adj$spp_adj

regs_spp <- data.frame(species=centroids85$species, regionsB=centroids85$regionsB, stringsAsFactors = FALSE)
biomass26 <- merge(biomass26, regs_spp, by='species', all.x=T)
biomass85 <- merge(biomass85, regs_spp, by='species', all.x=T)
 
uncer_spp26 <- data.frame(species=centroids26$species, uncertainty=centroids26$uncertainty, stringsAsFactors = FALSE)
uncer_spp85 <- data.frame(species=centroids85$species, uncertainty=centroids85$uncertainty, stringsAsFactors = FALSE)
biomass26 <- merge(biomass26, uncer_spp26, by='species', all.x=T)
biomass85 <- merge(biomass85, uncer_spp85, by='species', all.x=T)

save(biomass26, biomass85, file='output/Therm.Biomass_summaries_Dec2017.RData')

# ============================================================================================================
# Thermal biomass figure
# ============================================================================================================
  
# a couple qaqc plots
plot(sdChan~biomMeanChan, data=biomass26[biomass26$uncertainty == 'low',])
plot(sdChan~biomMeanChan, xlim=c(-110,2200), ylim=c(0,2000),data=biomass85[biomass85$uncertainty == 'low',])
 
# Assign numbers to regions for the plot (based on prelim looks), so gives an ascending order
regs <- data.frame(regionsB=c('AFSC_EBS','AFSC_GOA','DFO_Newfoundland','NEFSC_NEUS','NWFSC_WCAnn','SCDNR_SEUS','SEFSC_GOMex'), 
                   xaxis=c(1,5,2,4,6,7,3), stringsAsFactors = FALSE)
biomass26 <- merge(biomass26, regs, by='regionsB', all.x=TRUE)
biomass85 <- merge(biomass85, regs, by='regionsB', all.x=TRUE)

legend.names <- c("E. Bering S.","E. Canada","G. Mexico","NE U.S.","G. Alaska","West U.S.","SE U.S.")     


pdf(file='figures/thermBiomass_boxPlots_Dec2016.pdf', width = 6.5, height = 5) 
layout(mat=matrix(c(4,1,5,2,3,6), 1, 6), widths=c(1,4.9,1,4,.7,.5))# Need to adjust this to fit the three figures
layout.show(n=3)
par(mar=c(6,.1,1.5,.1), oma=c(.1,.1,.1,.1), cex.axis=1.2, bty='u')
boxplot(biomMeanChan~xaxis, xlim=c(.5,7.5), ylim=c(-75, 260), axes=F, xaxt='n', yaxt='n', width=c(1,1,1,1,1,1,1), pch=20, staplewex=.25, names=legend.names, data=biomass26[biomass26$uncertainty == 'low',])
abline(a=0, b=0, lwd=.5, col='gray50')
boxplot(biomMeanChan~xaxis, xlim=c(.5,7.5), ylim=c(-75, 260), axes=F, xaxt='n', yaxt='n', width=c(1,1,1,1,1,1,1), pch=20, staplewex=.25, col='white',names=legend.names, data=biomass26[biomass26$uncertainty == 'low',], add=TRUE)
axis(side=1, at=seq(1:7), labels=F, tcl=-.2, bty='l')
axis(side=2,tcl=-.2, padj=1)
text(cex=1.2, x=c(1:7)-.3, y=-93, legend.names, pos=4, xpd=TRUE, srt=280)
#text(cex=.7, x=c(1:7), y=-69, labels=c('84','67','76','17','71','34','80'), pos=1, xpd=TRUE)
text(cex=.9, x=c(1:7), y=273, labels=c('','1','1','','','','1'), pos=1, xpd=TRUE)
mtext('Change in thermal habitat (% of 2006-2020)', side=2, srt=290, padj=-2.7)
mtext('A) RCP 2.6', adj=0.05, line=0.2)
box(bty = "L") 

abc <- biomass85[!biomass85$regionsB=='SCDNR_SEUS',]
boxplot(biomMeanChan~xaxis, axes=F, xlim=c(.5,6.25), ylim=c(-115, 770), width=c(1,1,1,1,1,1), xaxt='n', staplewex=.25, pch=20, names=legend.names[-7], data=abc[abc$uncertainty == 'low',])
abline(a=0, b=0, lwd=.5, col='gray50')
boxplot(biomMeanChan~xaxis, axes=F, xlim=c(.5,6.25), ylim=c(-115, 770), xaxt='n', staplewex=.25, pch=20, col='white', names=legend.names[-7], data=abc[abc$uncertainty == 'low',], add=TRUE)
axis(1, at=c(1:6), labels=F, tcl=-.2)
axis(2, at=c(-100, 0, 100, 200, 300, 400, 500, 600, 700, 800))
text(cex=1.2, x=c(1:6)-.3, y=-165, legend.names[-7], pos=4, xpd=TRUE, srt=280)
text(cex=.9, x=c(1:6), y=805, labels=c('','','5','','',''), pos=1, xpd=TRUE)
mtext('B) RCP 8.5', adj=0.05, line=0.2)
box(bty = "L") 

boxplot(biomMeanChan~xaxis, axes=F, ylim=c(-115, 4500), width=c(1), at=7, xaxt='n', staplewex=.25, pch=20, names=legend.names[7], data=biomass85[biomass85$uncertainty == 'low' & biomass85$regionsB=='SCDNR_SEUS',])
abline(a=0, b=0, lwd=.5, col='gray50')
boxplot(biomMeanChan~xaxis, axes=F, ylim=c(-115, 4500), width=c(1), at=7, xaxt='n', staplewex=.25, pch=20, names=legend.names[7], data=biomass85[biomass85$uncertainty == 'low' & biomass85$regionsB=='SCDNR_SEUS',], add=T)
axis(1, at=c(6.5,7,7.5),labels=F, tcl=-.2)
text(cex=1.2, x=c(7)-.3, y=-360, legend.names[7], pos=4, xpd=TRUE, srt=280)
text(cex=.9, x=c(7), y=4700, labels=c('13'), pos=1, xpd=TRUE)
axis(4)
dev.off()

# ===========================================================================================================
# example graphs of biomass shift with specific species
# ===========================================================================================================
load('data/MinMaxDepths_spp_Dec2017.RData')
depths <- data.frame(spp = projspp, maxDepth = maxDepth, minDepth = minDepth, stringsAsFactors = F)
# load the bathygrid, fine scale
load('data/ProjectionBathGrid_Feb27_2017.RData')
depth.grid <- data.frame(latitude = proj.grid$latBathgrid, longitude = proj.grid$lonBathgrid, depth = -1*proj.grid$depth, stringsAsFactors = F)
# EAST COAST 
  
  species1 = 'sebastes fasciatus_Atl'
  species2 = 'lutjanus campechanus_Atl'
  species3 = 'seriola dumerili_Atl'

  filename <- paste(projfolder, species1, '_rcp', rcp[2], '_jas_prediction_AGG.RData', sep='')
  load(filename)
  #pred.agg85 <- as.matrix(pred.agg[,4:19], nrow=412960, ncol=16)
  pred.agg85 <- as.matrix(log(pred.agg[,4:19] + 1), nrow=412960, ncol=16)
  # aggregate by lon lat and year range
  ensMean <- apply(pred.agg85, 1, FUN=mean, na.rm=T)
  pred.agg1 <- data.frame(cbind(pred.agg, ensMean=ensMean))
  #pred.agg1 <- pred.agg1[pred.agg1$latitude > 31,]
  #pred.agg1 <- merge(pred.agg1, depth.grid, by=c('latitude', 'longitude'), all.x=T)
  #summary(pred.agg1$depth)
  #depthLim <- depths$maxDepth[depths$spp == species1]
  #pred.agg1 <- pred.agg1[pred.agg1$depth < depthLim + 1,]
  
  filename <- paste(projfolder, species2, '_rcp', rcp[2], '_jas_prediction_AGG.RData', sep='')
  load(filename)
  #pred.agg85 <- as.matrix(pred.agg[,4:19], nrow=412960, ncol=16)
  pred.agg85 <- as.matrix(log(pred.agg[,4:19] + 1), nrow=412960, ncol=16)
  # aggregate by lon lat and year range
  ensMean <- apply(pred.agg85, 1, FUN=mean, na.rm=T)
  pred.agg2 <- data.frame(cbind(pred.agg, ensMean=ensMean))
  #pred.agg2 <- pred.agg2[pred.agg2$latitude < 45.5,]
  #pred.agg2 <- merge(pred.agg2, depth.grid, by=c('latitude', 'longitude'), all.x=T)
  #summary(pred.agg2$depth)
  #depthLim <- depths$maxDepth[depths$spp == species2]
  #pred.agg2 <- pred.agg2[pred.agg2$depth < depthLim + 1,]
  
  filename <- paste(projfolder, species3, '_rcp', rcp[2], '_jas_prediction_AGG.RData', sep='')
  load(filename)
  #pred.agg85 <- as.matrix(pred.agg[,4:19], nrow=412960, ncol=16)
  pred.agg85 <- as.matrix(log(pred.agg[,4:19] + 1), nrow=412960, ncol=16)
  # aggregate by lon lat and year range
  ensMean <- apply(pred.agg85, 1, FUN=mean, na.rm=T)
  pred.agg3 <- data.frame(cbind(pred.agg, ensMean=ensMean))
  #pred.agg3 <- pred.agg3[!pred.agg3$latitude < 31,]
  #pred.agg3 <- pred.agg3[!(pred.agg3$longitude < -80.75 & pred.agg3$latitude < 27),]
  
  #pred.agg3 <- pred.agg3[pred.agg3$latitude > 34,]
  #pred.agg3 <- pred.agg3[pred.agg3$longitude < -58,]
  #pred.agg3 <- pred.agg3[pred.agg3$latitude < 49,]
  #pred.agg3 <- merge(pred.agg3, depth.grid, by=c('latitude', 'longitude'), all.x=T)
  #pred.agg3 <- pred.agg3[pred.agg3$depth < 333,]
  #pred.agg3 <- merge(pred.agg3, depth.grid, by=c('latitude', 'longitude'), all.x=T)
  #summary(pred.agg3$depth)
  #depthLim <- depths$maxDepth[depths$spp == species3]
  #pred.agg3 <- pred.agg3[pred.agg3$depth < depthLim + 1,]
  
  #pdf(width=7.5, height=7.75, file='figures/rcp85_ensemble_projections_LOGGED_K.Mack-Cobia-At.Menh_Jan2018.pdf')
  #pdf(width=7.5, height=7.75, file='/Users/jamesmorley/Desktop/rcp85_ensemble_projections_LOG_pinfi-sandlan-stAnch-MODS_Jan2018.pdf')
  pdf(width=7.5, height=7.75, file='/Users/jamesmorley/Desktop/rcp85_ensemble_projections_LOG_RedFish-RedSnap-GAmber_Jan2018.pdf')
  #pdf(width=8, height=11, file='/Users/jamesmorley/Desktop/rcp85_ensemble_projections_lobster-NShrimp-NPinkShrim.pdf')
  
    cols = colorRampPalette(colors = c('gray92', 'blue', 'dark blue', 'black'))
    
    scale85 = seq(0, max(pred.agg1$ensMean) + .000001, length.out=20)
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
    
    scale85 = seq(0, max(pred.agg2$ensMean) + .000001, length.out=20)
    # STANDARDIZE PLOT DIMENSIONS FOR ALL PLOTS
    xlimit = c(min(pred.agg2$longitude[pred.agg2$ensMean > scale85[2]/2] - 1),max(pred.agg2$longitude[pred.agg2$ensMean > scale85[2]/2] + 1))
    ylimit = c(min(pred.agg2$latitude[pred.agg2$ensMean > scale85[2]/2] - .25), max(pred.agg2$latitude[pred.agg2$ensMean > scale85[2]/2] + .25))

    # RCP85 plots_one for normal, one for Logged+1
    Projmap <- map('world', xlim=xlimit, ylim=ylimit, plot=F) 
    Projmap <- data.frame(lon=Projmap$x, lat=Projmap$y)
    plot2 <-levelplot(ensMean ~ longitude*latitude|year_range, data=pred.agg2, index.cond=list(c(5,4,3,2,1)), colorkey=F, par.settings=list(layout.heights=list(top.padding=-2, bottom.padding=-2),layout.widths=list(left.padding=-2.25, right.padding=-2.25)),
                      scales=list(tck=c(0,0), alternating=c(0,0)), strip=F, aspect='fill', xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
                      at = scale85, col.regions=cols, layout = c(1, 5)) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
    
    scale85 = seq(0, max(pred.agg3$ensMean) + .000001, length.out=20)
    # STANDARDIZE PLOT DIMENSIONS FOR ALL PLOTS
    xlimit = c(min(pred.agg3$longitude[pred.agg3$ensMean > scale85[2]/2] - 1),max(pred.agg3$longitude[pred.agg3$ensMean > scale85[2]/2] + 1))
    ylimit = c(min(pred.agg3$latitude[pred.agg3$ensMean > scale85[2]/2] - .25), max(pred.agg3$latitude[pred.agg3$ensMean > scale85[2]/2] + .25))
    
    # RCP85 plots_one for normal, one for Logged+1
    Projmap <- map('world', xlim=xlimit, ylim=ylimit, plot=F) 
    Projmap <- data.frame(lon=Projmap$x, lat=Projmap$y)
    plot3 <-levelplot(ensMean ~ longitude*latitude|year_range, data=pred.agg3, index.cond=list(c(5,4,3,2,1)), colorkey=F, par.settings=list(layout.heights=list(top.padding=-2, bottom.padding=-2),layout.widths=list(left.padding=-2.5, right.padding=-2)),
                      scales=list(tck=c(0,0), alternating=c(0,0)), strip=F,  aspect='fill', xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
                      at = scale85, col.regions=cols, layout = c(1, 5)) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
   #plot3 <-levelplot(ensMean ~ longitude*latitude|year_range, data=pred.agg3, index.cond=list(c(1,2,3,4,5)), colorkey=F, par.settings=list(layout.heights=list(top.padding=-2, bottom.padding=-2),layout.widths=list(left.padding=-2.5, right.padding=-2)),
    #                  scales=list(tck=c(0,0), alternating=c(0,0)), strip=F,  aspect='fill', xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
     #                 at = scale85, col.regions=cols, layout = c(5, 1)) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
    #plot3 <-levelplot(ensMean ~ longitude*latitude, data=pred.agg3[pred.agg3$year_range=='2081-2100',], colorkey=F, par.settings=list(layout.heights=list(top.padding=-2, bottom.padding=-2),layout.widths=list(left.padding=-2.5, right.padding=-2)),
     #                 scales=list(tck=c(0,0), alternating=c(0,0)), strip=F,  aspect='fill', xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
      #                at = scale85, col.regions=cols) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray') +
      #xyplot(c(38.63375, 40.78177) ~ c(-73.81012, -70.97968), type='l', lty=1, lwd=4, col='firebrick1') 

    par(mar=c(.1,.1,.1,.1))
    print(plot1, split=c(1,1,3,1), more=TRUE) #position=c(0,0,.5,1), 
    print(plot2, split=c(2,1,3,1), more=TRUE)#position=c(.5,0,1,1))
    print(plot3, split=c(3,1,3,1))#position=c(.5,0,1,1))
    #print(plot3)      
  dev.off()

    
# WEST COAST
  species1 = 'hexagrammos stelleri_Pac'
  species2 = 'hippoglossus stenolepis_Pac'
  species3 = 'sebastes alutus_Pac'
  
  filename <- paste(projfolder, species1, '_rcp', rcp[2], '_jas_prediction_AGG.RData', sep='')
  load(filename)
  pred.agg85 <- as.matrix(pred.agg[,4:19], nrow=329130, ncol=16)
  #pred.agg85 <- as.matrix(log(pred.agg[,4:19] + 1), nrow=329130, ncol=16)
  # aggregate by lon lat and year range
  ensMean <- apply(pred.agg85, 1, FUN=mean, na.rm=T)
  pred.agg1 <- data.frame(cbind(pred.agg, ensMean=ensMean))
  
  filename <- paste(projfolder, species2, '_rcp', rcp[2], '_jas_prediction_AGG.RData', sep='')
  load(filename)
  pred.agg85 <- as.matrix(pred.agg[,4:19], nrow=329130, ncol=16)
  #pred.agg85 <- as.matrix(log(pred.agg[,4:19] + 1), nrow=329130, ncol=16)
  # aggregate by lon lat and year range
  ensMean <- apply(pred.agg85, 1, FUN=mean, na.rm=T)
  pred.agg2 <- data.frame(cbind(pred.agg, ensMean=ensMean))
  
  filename <- paste(projfolder, species3, '_rcp', rcp[2], '_jas_prediction_AGG.RData', sep='')
  load(filename)
  pred.agg85 <- as.matrix(pred.agg[,4:19], nrow=329130, ncol=16)
  #pred.agg85 <- as.matrix(log(pred.agg[,4:19] + 1), nrow=329130, ncol=16)
  # aggregate by lon lat and year range
  ensMean <- apply(pred.agg85, 1, FUN=mean, na.rm=T)
  pred.agg3 <- data.frame(cbind(pred.agg, ensMean=ensMean))
  
  #pdf(width=7.5, height=7.75, file='figures/rcp85_ensemble_projections_PACIFICspp_Jan2018.pdf')
  pdf(width=7.5, height=7.75, file='/Users/jamesmorley/Desktop/rcp85_ensemble_projections_WSGreen-Halibut-P.OceanPerch_Jan2018.pdf')
  
  cols = colorRampPalette(colors = c('gray92', 'blue', 'dark blue', 'black'))
  
  scale85 = seq(0, max(pred.agg1$ensMean) + .000001, length.out=20)
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
  
  scale85 = seq(0, max(pred.agg2$ensMean) + .000001, length.out=20)
  # STANDARDIZE PLOT DIMENSIONS FOR ALL PLOTS
  xlimit = c(min(pred.agg2$longitude[pred.agg2$ensMean > scale85[2]/2] - 1),max(pred.agg2$longitude[pred.agg2$ensMean > scale85[2]/2] + 1))
  ylimit = c(min(pred.agg2$latitude[pred.agg2$ensMean > scale85[2]/2] - .25), max(pred.agg2$latitude[pred.agg2$ensMean > scale85[2]/2] + .25))
  
  # RCP85 plots_one for normal, one for Logged+1
  Projmap <- map('world', xlim=xlimit, ylim=ylimit, plot=F) 
  Projmap <- data.frame(lon=Projmap$x, lat=Projmap$y)
  plot2 <-levelplot(ensMean ~ longitude*latitude|year_range, data=pred.agg2, index.cond=list(c(5,4,3,2,1)), colorkey=F, par.settings=list(layout.heights=list(top.padding=-2, bottom.padding=-2),layout.widths=list(left.padding=-2.25, right.padding=-2.25)),
                    scales=list(tck=c(0,0), alternating=c(0,0)), strip=F, aspect='fill', xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
                    at = scale85, col.regions=cols, layout = c(1, 5)) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
  
  scale85 = seq(0, max(pred.agg3$ensMean) + .000001, length.out=20)
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
 
#mafmc <- read.csv("/Users/jamesmorley/Desktop/MAFMC_species.csv", stringsAsFactors = FALSE)  
# Add '_Atl' then merge with centroids85$species to get a species list
#mafmc <- paste(mafmc$species, '_Atl', sep='')
#testSpp <- intersect(mafmc, centroids85$species)  
  
# I originally did this as a loop to look at all these potential species for example graphs
# now it's retooled to just make the graphs for species of interest.
testSpp <- c('caretta caretta_Atl','paralichthys albigutta_Atl','rhizoprionodon terraenovae_Atl',
             'rhinoptera bonasus_Atl','scomberomorus cavalla_Atl','scomberomorus maculatus_Atl',
             'stomolophus meleagris_Atl','trachinotus carolinus_Atl','myliobatis freminvillii_Atl',
             'gymnura micrura_Atl','menippe mercenaria_Atl','archosargus probatocephalus_Atl')
for(i in 1:length(testSpp)){ 
  species = testSpp[i] 
  print(paste("Beginning figures for ", species, sep=''))
  filename <- paste(projfolder, species, '_rcp', rcp[2], '_jas_prediction_AGG.RData', sep='')
  load(filename)
  
  pred.agg85 <- as.matrix(log(pred.agg[,4:19]+1), nrow=412960, ncol=16)
  # aggregate by lon lat and year range
  ensMean <- apply(pred.agg85, 1, FUN=mean, na.rm=T)
  pred.agg85 <- data.frame(cbind(pred.agg, ensMean=ensMean))
  if(grepl('_Atl', species)){
    pred.agg85 <- pred.agg85[!pred.agg85$longitude < -82,]
    pred.agg85 <- pred.agg85[!pred.agg85$latitude < 27,]
  }
  
  pdf(width=14, height=5, file=paste('figures/speciesProjections_Dec2017/ensemble_proj/', species, '_rcp85_ensemble_projections_Dec2017.pdf', sep=''))

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
library('stringr')
load('output/modeldiag_Nov2017_fitallreg_2017.Rdata')
diags <- data.frame(Species_orig = modeldiag$sppocean, devPA = modeldiag$dev.pres, devBiom = modeldiag$dev.biomass, stringsAsFactors = F)

load('output/centroid_summaries_Dec2017.RData')
load('output/Therm.Biomass_summaries_Dec2017.RData')
# sort biomass and centroid files so they have the same order and are all four identical_then make the following table and do a merger based on 'species_orig'
biomass26 <- biomass26[order(biomass26$species),]
biomass85 <- biomass85[order(biomass85$species),]
centroids26 <- centroids26[order(centroids26$species),]
centroids85 <- centroids85[order(centroids85$species),]

S1app <- data.frame(Species = centroids26$species, Species_orig = centroids26$species_orig, Region = centroids26$regionsB, Uncertainty = centroids26$uncertainty,
                    Shift = centroids26$meanDist, sd_shift = centroids26$sdDist, Biomass = biomass26$biomMeanChan,
                    sd_biom = biomass26$sdChan, stringsAsFactors=F)
S1app <- merge(S1app, diags, by='Species_orig', all.x=T)

S2app <- data.frame(Species_orig = centroids85$species_orig, Species = centroids85$species, Region = '', Uncertainty = centroids85$uncertainty,
                    Shift = centroids85$meanDist, sd_shift = centroids85$sdDist, Biomass = biomass85$biomMeanChan,
                    sd_biom = biomass85$sdChan, devPA = NA, devBiom = NA, stringsAsFactors=F)
# add a '2' to species names for S2app
S2app$Species <- paste(S2app$Species, '2', sep='')
Sapp <- rbind(S1app, S2app)
Sapp <- Sapp[order(Sapp$Species),]
Sapp$drops <- duplicated(Sapp$Species_orig)
Sapp$devPA <- ifelse(Sapp$drops==TRUE, '', Sapp$devPA)
Sapp$devBiom <- ifelse(Sapp$drops==TRUE, '', Sapp$devBiom)
Sapp$RCP <- ifelse(grepl(2, Sapp$Species), 85, 26)

sppNames <- str_split_fixed(Sapp$Species, '_', n=2)
Sapp$Species <- as.character(sppNames[,1])
Sapp$drops <- duplicated(Sapp$Species)
Sapp$Species_orig <- NULL
Sapp$Species <- ifelse(Sapp$drops==TRUE, '', Sapp$Species)
Sapp$drops <- NULL
regions <- data.frame(Region=c('','AFSC_EBS','AFSC_GOA','DFO_Newfoundland','NEFSC_NEUS',
            'NWFSC_WCAnn','SCDNR_SEUS','SEFSC_GOMex'),region=c('','E. Bering S.','G. Alaska',
            'E. Canada','NE U.S.','West U.S.','SE U.S.','G. Mexico'), stringsAsFactors=F)
 
Sapp$index <- c(1:nrow(Sapp))
Sapp <- merge(Sapp, regions, by='Region', all.x=T)
Sapp <- Sapp[order(Sapp$index),]

# The rest was modified in Excel (e.g., moving columns around)
save(Sapp, file='/Users/jamesmorley/Desktop/Projections_S1app.RData')
write.csv(Sapp, file='/Users/jamesmorley/Desktop/Projections_S1app.csv')
write.csv(Sapp, file='/Users/jamesmorley/Desktop/Projections_S1app_FILLED.csv')



