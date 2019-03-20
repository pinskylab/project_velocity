setwd('/Users/jamesmorley/Documents/project_velocity')

require(Hmisc)
library(latticeExtra)
library(maps)
library(geosphere)
library(plotrix)
library(quantreg)
library(SDMTools)
library(ggplot2)
library(reshape2)
library(stringr)
library(taxize)
library(data.table)

modelrun <- c('bcc-csm1-1-m','bcc-csm1-1','CanESM2','CCSM4','CESM1-CAM5','CNRM-CM5','GFDL-CM3','GFDL-ESM2M','GFDL-ESM2G','GISS-E2-R','GISS-E2-H','IPSL-CM5A-LR','IPSL-CM5A-MR','MIROC-ESM','MPI-ESM-LR','NorESM1-ME')
rcp <- c(26,85)
projfolder = 'output/CEmodels_proj_Nov2017/'
load('data/speciesProjectionList.RData')
load('data/hauls_catch_Dec2017.RData') 
load('data/projspp_wGmex.RData')

# Get common names of species for the website
sci_names <- gsub("_.*", "", projspp)
sci_namesMod <- data.frame(spp=paste(toupper(substr(sci_names, 1,1)), substr(sci_names,2, nchar(sci_names)), sep=''), stringsAsFactors = F)

comm_names <- sci2comm(sci_namesMod$spp, db='itis', simplify=F)
save(comm_names, file='data/SciNames_OceanAd_May2018.RData')

OceanAd_taxa <- data.frame(sci_name=as.character(), common=as.character(), stringsAsFactors = FALSE)

for(i in 1:length(comm_names)){
  species <- sci_namesMod[i,1]
  print(paste('Running species:', species))
  taxa <- data.frame(comm_names[[i]])
  taxaName <- taxa$commonName[taxa$language=='English'][1]
  if(is.null(taxaName)){
    taxaName <- 'No common name'
  }
  if(is.na(taxaName)){
    taxaName <- taxa$commonName[taxa$language=='Spanish'][1]
  }
  OceanAd_taxa[i,] <- data.frame(sci_name=as.character(species), common=as.character(taxaName), stringsAsFactors = F)
}

 
# May not need this, below.....
taxa <- data.table(read.csv("Jim/spp.key.csv", stringsAsFactors=F)) # Brings in taxonomy file from Ryan's trawlData package
taxa <- taxa[,c(2,16)]
taxa <- unique(taxa)
taxaMod <- merge(sci_namesMod, taxa, by='spp', all.x=T)
OceanAd_taxa$common_Alt <- taxaMod$common
write.csv(OceanAd_taxa, file='data/OceanAd_taxa.csv')

OceanAd_taxa$common_AB <- ifelse(OceanAd_taxa$common == 'No common name', OceanAd_taxa$common_Alt, OceanAd_taxa$common)
write.csv(OceanAd_taxa, file='data/OceanAd_taxa.csv')


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

 
# columns for init and final lat and lon, mean + SD in latitudinal change, mean + SD in absolute shift distance
centroids26 <- data.frame(lat_2007_2020=numeric(), lat_2021_2040=numeric(), lat_2041_2060=numeric(), lat_2061_2080=numeric(), lat_2081_2100=numeric(), sd_2007_2020=numeric(), sd_2021_2040=numeric(), sd_2041_2060=numeric(), sd_2061_2080=numeric(), sd_2081_2100=numeric())
centroids85 <- data.frame(lat_2007_2020=numeric(), lat_2021_2040=numeric(), lat_2041_2060=numeric(), lat_2061_2080=numeric(), lat_2081_2100=numeric(), sd_2007_2020=numeric(), sd_2021_2040=numeric(), sd_2041_2060=numeric(), sd_2061_2080=numeric(), sd_2081_2100=numeric())
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
  yearRange <- c('2007-2020', '2021-2040', '2041-2060', '2061-2080','2081-2100')
  cent_preds26lat <- matrix(data=NA, nrow=5, ncol=16)
  for(j in 1:5){
    years = yearRange[j]
    preds <- pred.agg26[pred.agg26$year_range == years,]
    for(k in 4:19){
      weights_adj <- log(preds[,k] + 1) * (preds$area) #*100) # the '*100' is in case I want to do hectares_it shouldn't matter
      valueLat = wtd.mean(preds$latitude, weights = weights_adj, na.rm=FALSE) 
      if(is.na(valueLat)){
        print('Problem')
      }
      cent_preds26lat[j,k-3] = valueLat
    }
  }  
  mean_Lat_26 <- apply(cent_preds26lat, 1, FUN=mean, na.rm=T)
  sd_Lat_26 <- apply(cent_preds26lat, 1, FUN=sd, na.rm=T)

  # RCP85 calculations
  cent_preds85lat <- matrix(data=NA, nrow=5, ncol=16)
  for(j in 1:5){
    years = yearRange[j]
    preds <- pred.agg85[pred.agg85$year_range == years,]
    for(k in 4:19){
      weights_adj <- log(preds[,k] + 1) * (preds$area) #*100) # the '*100' is in case I want to do hectares_it shouldn't matter
      valueLat = wtd.mean(preds$latitude, weights = weights_adj, na.rm=FALSE) 
      if(is.na(valueLat)){
        print('Problem')
      }
      cent_preds85lat[j,k-3] = valueLat
    }
  } 
  mean_Lat_85 <- apply(cent_preds85lat, 1, FUN=mean, na.rm=T)
  sd_Lat_85 <- apply(cent_preds85lat, 1, FUN=sd, na.rm=T)
  
  centroids26[i,] <- data.frame(lat_2007_2020=mean_Lat_26[1], lat_2021_2040=mean_Lat_26[2], lat_2041_2060=mean_Lat_26[3], lat_2061_2080=mean_Lat_26[4], lat_2081_2100=mean_Lat_26[5], sd_2007_2020=sd_Lat_26[1], sd_2021_2040=sd_Lat_26[2], sd_2041_2060=sd_Lat_26[3], sd_2061_2080=sd_Lat_26[4], sd_2081_2100=sd_Lat_26[5])
  centroids85[i,] <- data.frame(lat_2007_2020=mean_Lat_85[1], lat_2021_2040=mean_Lat_85[2], lat_2041_2060=mean_Lat_85[3], lat_2061_2080=mean_Lat_85[4], lat_2081_2100=mean_Lat_85[5], sd_2007_2020=sd_Lat_85[1], sd_2021_2040=sd_Lat_85[2], sd_2041_2060=sd_Lat_85[3], sd_2061_2080=sd_Lat_85[4], sd_2081_2100=sd_Lat_85[5])
}

centroids26$species <- projspp_adj$spp_adj
centroids85$species <- projspp_adj$spp_adj
centroids26$rcp <- 26
centroids85$rcp <- 85
centroids <- rbind(centroids26, centroids85)

# NEED TO MELT DATA_PERHAPS DO MEAN AND SD SEPARATE THEN MERGE_KEEP LOGICS (etc.) WITH ONE OR THE OTHER
centroids_mean <- centroids[,c(11, 12, 1:5)]
centroids_sd <- centroids[,c(11, 12, 6:10)] 
centroids_mean <- melt(centroids_mean, id.vars=c('species', 'rcp'), variable.name='lat_year_range', value.name='mean_latitude')
centroids_sd <- melt(centroids_sd, id.vars=c('species', 'rcp'), variable.name='lat_year_range', value.name='st.dev_latitude')
centroids_mean$year_range <- substr(centroids_mean$lat_year_range, 5,13)
centroids_sd$year_range <- substr(centroids_sd$lat_year_range, 4,12)
centroids_mean$lat_year_range <- NULL
centroids_sd$lat_year_range <- NULL

centroids <- merge(centroids_mean, centroids_sd, by=c('species', 'rcp', 'year_range'), all=TRUE)
centroids$region <- ifelse(grepl('_Atl', centroids$species), 'Atlantic', ifelse(grepl('_Gmex', centroids$species), 'Gulf of Mexico', 'Pacific'))
centroids$spp <- gsub("_.*","",centroids$species)
centroids$name <- paste(toupper(substr(centroids$spp, 1,1)), substr(centroids$spp,2, nchar(centroids$spp)), sep='')
centroids$spp <- NULL

sppList <- projspp_adj[,c(3:5)]
colnames(sppList)[3] <- 'species'
centroids <- merge(centroids, sppList, by='species', all.x=TRUE)

save(centroids, file='output/centroid_summaries_OceanAdapt2018.RData')
 
