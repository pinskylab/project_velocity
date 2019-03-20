library(latticeExtra)
library(maps)

setwd('Documents/project_velocity')
load('data/speciesProjectionList.RData')
  
rcp <- c(26,85)
 
for(i in 1:length(projspp)){
  species <- projspp[i]
  load(paste('output/CEmodels_proj_Nov2017/',species,'_rcp26_jas_prediction_AGG.RData', sep=''))
  proj26 <- pred.agg
  load(paste('output/CEmodels_proj_Nov2017/',species,'_rcp85_jas_prediction_AGG.RData', sep=''))
  proj85 <- pred.agg  
  rm(pred.agg)
  
  if(grepl('_Atl', species)){
    pred.agg26_log <- as.matrix(log(proj26[,4:19] + 1), nrow=412960, ncol=16)
    pred.agg85_log <- as.matrix(log(proj85[,4:19] + 1), nrow=412960, ncol=16)
  
    # aggregate by lon lat and year range
    ensMean26_log <- apply(pred.agg26_log, 1, FUN=mean, na.rm=T)
    ensMean85_log <- apply(pred.agg85_log, 1, FUN=mean, na.rm=T)
    
    pred_26_85 <- data.frame(cbind(proj85[,1:3], ensMean26_logged=ensMean26_log, ensMean85_logged=ensMean85_log))
  }
  if(grepl('_Pac', species)){
    pred.agg26_log <- as.matrix(log(proj26[,4:19] + 1), nrow=329130, ncol=16)
    pred.agg85_log <- as.matrix(log(proj85[,4:19] + 1), nrow=329130, ncol=16)
    
    # aggregate by lon lat and year range
    ensMean26_log <- apply(pred.agg26_log, 1, FUN=mean, na.rm=T)
    ensMean85_log <- apply(pred.agg85_log, 1, FUN=mean, na.rm=T)
    
    pred_26_85 <- data.frame(cbind(proj85[,1:3], ensMean26_logged=ensMean26_log, ensMean85_logged=ensMean85_log))
  }
  
  if(grepl('_Atl', species)){
    star_name <- toupper(substr(species, 1,1))
    end_name <- substr(species, 2, nchar(species))
    speciesB <- paste(star_name, end_name, 'antic', sep='')
  }
  if(grepl('_Pac', species)){
    star_name <- toupper(substr(species, 1,1))
    end_name <- substr(species, 2, nchar(species))
    speciesB <- paste(star_name, end_name, 'ific', sep='')
  }  
  print(paste('Finishing files for', speciesB))
  save(pred_26_85, file=paste('output/CEmodels_proj_OceanAd/', speciesB,'_OceanAd_projection.RData', sep=''))
}

# =====================================================================
# Plot the ocean adapt figures along with trawl observations figure 
# =====================================================================
species <- projspp[396]
if(grepl('_Atl', species)){
  star_name <- toupper(substr(species, 1,1))
  end_name <- substr(species, 2, nchar(species))
  speciesB <- paste(star_name, end_name, 'antic', sep='')
}
if(grepl('_Pac', species)){
  star_name <- toupper(substr(species, 1,1))
  end_name <- substr(species, 2, nchar(species))
  speciesB <- paste(star_name, end_name, 'ific', sep='')
}  

load(paste('output/CEmodels_proj_OceanAd/', speciesB,'_OceanAd_projection.RData', sep=''))

load('data/hauls_catch_Dec2017.RData') # hauls and catch data
mydat <- dat[dat$sppocean==species,] 
spdata <- merge(hauls, mydat, by='haulid', all.x = T, sort=F) # Add haul data and empty hauls
myocean <- head(spdata$ocean[spdata$presfit == TRUE])[1] # identify if this is a Pacific or Atlantic species
spdata <- spdata[spdata$ocean == myocean,] # trim master hauls file to the ocean of interest 
spdata$presfit[is.na(spdata$presfit)] <- FALSE

plot(lat~lon, pch=20, cex=.2, spdata[spdata$presfit==T,])

xlimit = c(min(pred_26_85$longitude),max(pred_26_85$longitude))
ylimit = c(min(pred_26_85$latitude), max(pred_26_85$latitude))
scale85 = seq(0, max(pred_26_85$ensMean85_logged), length.out=20)
cols = colorRampPalette(colors = c('gray90', 'blue', 'dark blue', 'black'))

Projmap <- map('world', xlim=xlimit, ylim=ylimit, plot=F) 
Projmap <- data.frame(lon=Projmap$x, lat=Projmap$y)
print(levelplot(ensMean85_logged ~ longitude*latitude|year_range, data=pred_26_85, xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, layout = c(5, 1)) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')) 




