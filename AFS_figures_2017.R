# Example of projected biomass_canary rockfish

load('output/centroid_summaries_June2017_update.RData')

species3 = 'sebastes pinniger_Pac'
filename <- paste(projfolder, species3, '_rcp', rcp[2], '_jas_prediction_AGG.RData', sep='')
load(filename)
pred.agg85 <- as.matrix(pred.agg[,4:19], nrow=329130, ncol=16)
# aggregate by lon lat and year range
ensMean <- apply(pred.agg85, 1, FUN=mean, na.rm=T)
pred.agg2 <- data.frame(cbind(pred.agg, ensMean=ensMean))

pdf(width=7.5, height=7.75, file='figures/AFS_figure_rockfish_centroid.pdf')

cols = colorRampPalette(colors = c('gray92', 'blue', 'dark blue', 'black'))
scale85 = seq(0, max(pred.agg2$ensMean) + .00001, length.out=20)
# STANDARDIZE PLOT DIMENSIONS FOR ALL PLOTS
xlimit = c(min(pred.agg2$longitude[pred.agg2$ensMean > scale85[2]/2] - 1),max(pred.agg2$longitude[pred.agg2$ensMean > scale85[2]/2] + 1))
ylimit = c(min(pred.agg2$latitude[pred.agg2$ensMean > scale85[2]/2] - .25), max(pred.agg2$latitude[pred.agg2$ensMean > scale85[2]/2] + .25))

# RCP85 plots_one for normal, one for Logged+1
Projmap <- map('world', xlim=xlimit, ylim=ylimit, plot=F) 
Projmap <- data.frame(lon=Projmap$x, lat=Projmap$y)
coords <- centroids85[centroids85$species==species3,]
plot1 <- levelplot(ensMean ~ longitude*latitude|year_range, data=pred.agg2[pred.agg2$year_range == '2081-2100',],  colorkey=F, par.settings=list(layout.heights=list(top.padding=-2, bottom.padding=-2),layout.widths=list(left.padding=-2, right.padding=-2.5)),
          scales=list(tck=c(0,0), alternating=c(0,0)), strip=F, strip.left=TRUE,  aspect='fill', xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL,  
          at = scale85, col.regions=cols) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray',
          par.settings=list(layout.heights=list(top.padding=0, bottom.padding=0))) +
          xyplot(c(52.78846, 56.4445) ~ c(-134.5157, -158.1879), type='l', lty=1, lwd=2, col='black')#,
par(mar=c(.1,.1,.1,.1))
print(plot1) #position=c(0,0,.5,1), 

dev.off()

# FULL MAP CENTROID SHIFT VECTORS
load('data/ProjectionBathGrid_Feb27_2017.RData')# load bathymetry projection grid
Projmap <- map('world', xlim=c(min(centroids26$lonPre),max(centroids26$lonPre)), ylim=c(min(centroids26$latPre),max(centroids26$latPost)),plot=F) 
Projmap <- data.frame(lon=Projmap$x, lat=Projmap$y)

pdf(width=14, height=8, file=paste('figures/', 'species_centroid_shifts85_nouncer.pdf', sep=''))
#proj.grid2 <- proj.grid[!(proj.grid$latClimgrid > 60.5 & proj.grid$lonClimgrid > -64.62603),]
#proj.grid2 <- proj.grid2[!(proj.grid2$latClimgrid > 55 & proj.grid2$lonClimgrid < -64.6 & proj.grid2$lonClimgrid > -75),]
plot(latClimgrid~lonClimgrid, col='gray90', cex=.7, pch=16, xlab='Longitude', ylab='Latitude', main='Predicted centroid shifts for 658 species_RCP85', data=proj.grid) 
points(lat ~ lon, Projmap, type='l', lty=1, lwd=.7, col='gray50')
arrows(x0=centroids85$lonPre, y0=centroids85$latPre, x1=centroids85$lonPost, y1=centroids85$latPost, length=.05, lwd=.8, col='dark blue')
dev.off()

# uncertainty compass plots =============================

#for(i in 1:length(projspp)){ # one loop for each species
  #species = 'sebastes pinniger_Pac'
  species = 'farfantepenaeus duorarum_Atl'
  filename <- paste(projfolder, species, '_rcp', rcp[2], '_jas_prediction_AGG.RData', sep='')
  load(filename)
  pred.agg85 <- pred.agg
  rm(pred.agg)
  yearRange <- c('2007-2020', '2081-2100')
  
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
  # Need to calculate distance shifted for each model
  dist85 <- distHaversine(p1=matrix(cbind(cent_preds85lon[1,], cent_preds85lat[1,]), ncol=2), p2=matrix(cbind(cent_preds85lon[2,], cent_preds85lat[2,]), ncol=2)) / 1000# distance of vectors in km (I think in km anyway)
  bearing85 <- bearing(p1=matrix(cbind(cent_preds85lon[1,], cent_preds85lat[1,]), ncol=2), p2=matrix(cbind(cent_preds85lon[2,], cent_preds85lat[2,]), ncol=2))
  bearing85scaled <- ifelse(bearing85 < 0, bearing85 + 360, bearing85)
  #bearing85rads <- bearing85scaled * pi/180
  #xy85rads <- matrix(cbind(1*cos(bearing85rads), 1*sin(bearing85rads)), ncol=2, nrow=16)
  #plot(xy85rads[,2]~xy85rads[,1], xlim=c(-2,2), ylim=c(-2,2))
  #xyMean <- c(mean(xy85rads[,1]),mean(xy85rads[,2]))
  #radius85 <- sqrt(xyMean[1]^2 + xyMean[2]^2)
  
# Need to make compass plot with dist85 and bearing85(scaled) vectors
  pdf(width=8, height=8, file=paste('figures/', 'pinkShrimp compass.pdf', sep=''))
  par(cex.lab=2.2)
  polar.plot(lengths=dist85, polar.pos=bearing85scaled, radial.lim=c(0,1200),   radial.labels=c('','200','400', '600', '800', '1000', '1200 km'), grid.col='gray40',  boxed.radial=T, lwd=2.5, line.col='darkblue', labels=c(), label.pos=c(0,90,180,270), show.grid.labels=1, start=90, clockwise=T)
  dev.off()
   
 
  
  
  
  
  
  
  
  
  
  
  
  
  
  
