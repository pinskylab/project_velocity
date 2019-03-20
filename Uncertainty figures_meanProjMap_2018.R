library(latticeExtra)
library(maps)
library(reshape2)
  
setwd('/Users/jamesmorley/Documents/project_velocity')
projfolder <- 'output/CEmodels_Uncertainty_Proj_MEAN_2018/'

NICHE <- c('GLM', 'GAM', 'BRT')
RCP <- c(26,45,85)
SPECIES <- c('homarus americanus_Atl','paralichthys dentatus_Atl','centropristis striata_Atl','anoplopoma fimbria_Pac','doryteuthis opalescens_Pac','hippoglossus stenolepis_Pac','sebastes alutus_Pac')
modelrun <- c('bcc-csm1-1-m','bcc-csm1-1','CanESM2','CCSM4','CESM1-CAM5','CNRM-CM5','GFDL-CM3','GFDL-ESM2M','GFDL-ESM2G','GISS-E2-R','GISS-E2-H','HadGEM2-ES','IPSL-CM5A-LR','IPSL-CM5A-MR','MIROC-ESM','MIROC5','MPI-ESM-LR','NorESM1-ME')
YEARS <- c('2007-2020','2021-2040','2041-2060','2061-2080','2081-2100')

species = SPECIES[3]
rcp = RCP[3]
# Load the model diagnostics to help choose the niche model to use for the figure
load(paste('output/modeldiag_Uncer_', species, 'Feb2018.Rdata', sep=''))

niche = NICHE[2]
load(paste(projfolder, species, '_rcp', rcp, '_jas_', niche, '_proj_AGG_2018.RData', sep=''))

# Add depth data_for some tests on limiting depths and eliminating no-analogue issues
#load('data/ProjectionBathGrid_Feb27_2017.RData')
#proj.grid <- proj.grid[,c(1,2,5)]


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
    
# STANDARDIZE PLOT DIMENSIONS
xlimit = c(min(pred.agg_PA$longitude), max(pred.agg_PA$longitude) + 0.5)
ylimit = c(min(pred.agg_PA$latitude), max(pred.agg_PA$latitude))
cols = colorRampPalette(colors = c('gray92', 'blue', 'dark blue', 'black'))
Projmap <- map('world', xlim=xlimit, ylim=ylimit, plot=F) 
Projmap <- data.frame(lon=Projmap$x, lat=Projmap$y)



# First the plot to check things over and decide on the niche model to use in the MS ===========================================================

pdf(width=9, height=11, file=paste('figures/Uncertainty_MS_2018/Projections/', species, '_rcp', rcp, '_', niche, 'noLog.pdf', sep=''))

scale85 <- seq(0, max(pred.agg_PA$ensMean_PA) + .0001, length.out=20)
plot1 <- levelplot(ensMean_PA ~ longitude*latitude|year_range, data=pred.agg_PA, index.cond=list(c(5,4,3,2,1)),  colorkey=F, par.settings=list(layout.heights=list(top.padding=-2, bottom.padding=-2),layout.widths=list(left.padding=-2, right.padding=-2.5)),
        scales=list(tck=c(0,0), alternating=c(0,0)), strip=F, strip.left=TRUE,  aspect='fill', xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL,  
        at = scale85, col.regions=cols, layout = c(1, 5)) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray',
        par.settings=list(layout.heights=list(top.padding=0, bottom.padding=0)))

scale85 <- seq(0, max(pred.agg_biom$ensMean_BIOM) + .01, length.out=20)
plot2 <- levelplot(ensMean_BIOM ~ longitude*latitude|year_range, data=pred.agg_biom, index.cond=list(c(5,4,3,2,1)),  colorkey=F, par.settings=list(layout.heights=list(top.padding=-2, bottom.padding=-2),layout.widths=list(left.padding=-2, right.padding=-2.5)),
        scales=list(tck=c(0,0), alternating=c(0,0)), strip=F, strip.left=TRUE,  aspect='fill', xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL,  
        at = scale85, col.regions=cols, layout = c(1, 5)) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray',
        par.settings=list(layout.heights=list(top.padding=0, bottom.padding=0)))

par(mar=c(.1,.1,.1,.1))
print(plot1, split=c(1,1,2,1), more=TRUE) #position=c(0,0,.5,1), 
print(plot2, split=c(2,1,2,1))#position=c(.5,0,1,1))
dev.off()


 
# ===========================================================
# NOW TO MESS AROUND WITH DEPTH ON THE SCALING ==========================================================
# ===========================================================

#pred.agg_biom <- pred.agg_biom[,c(1,2,3,22)]
#colnames(proj.grid) <- c('latitude','longitude','depth')
#pred.agg_depth <- merge(pred.agg_biom, proj.grid, by=c('latitude', 'longitude'), all.x=TRUE)

#pred.agg_depthMod <- pred.agg_depth[pred.agg_depth$depth > -25,]
#scale85 <- seq(0, max(pred.agg_depthMod$ensMean_BIOM) + .01, length.out=20)

#pdf(height=6, width=12, file="/Users/jamesmorley/Desktop/BSB_depth_tweaking_25.pdf")
#levelplot(ensMean_BIOM ~ longitude*latitude|year_range, data=pred.agg_depthMod, index.cond=list(c(1,2,3,4,5)),  colorkey=F, par.settings=list(layout.heights=list(top.padding=-2, bottom.padding=-2),layout.widths=list(left.padding=-2, right.padding=-2.5)),
 #         scales=list(tck=c(0,0), alternating=c(0,0)), strip=T, strip.left=FALSE,  aspect='fill', xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL,  
  #        at = scale85, col.regions=cols, layout = c(5, 1)) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray',
   #       par.settings=list(layout.heights=list(top.padding=0, bottom.padding=0)))
#dev.off()




# ===========================================================================================================
# The main plot showing present and proportional change in distribution ===========================================
# Use this if using a biomass model_if PA, go to next block down
# STANDARDIZE PLOT DIMENSIONS

xlimit = c(min(pred.agg_biom$longitude), max(pred.agg_biom$longitude) + 0.5)
ylimit = c(min(pred.agg_biom$latitude), max(pred.agg_biom$latitude))
cols = colorRampPalette(colors = c('gray92', 'blue', 'dark blue', 'black'))
Projmap <- map('world', xlim=xlimit, ylim=ylimit, plot=F) 
Projmap <- data.frame(lon=Projmap$x, lat=Projmap$y)

scale85 <- seq(0, max(pred.agg_biom$ensMean_BIOM[pred.agg_biom$year_range=='2007-2020']) + .001, length.out=20)
plot1 <-  levelplot(ensMean_BIOM ~ longitude*latitude, data=pred.agg_biom[pred.agg_biom$year_range=='2007-2020',],  colorkey=F, par.settings=list(layout.heights=list(top.padding=-3.5, bottom.padding=-3.5),layout.widths=list(left.padding=-2.5, right.padding=-2.5)),
          scales=list(tck=c(0,0), alternating=c(0,0)), strip=F, strip.left=FALSE,  aspect='fill', xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL,  
          at = scale85, col.regions=cols) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=0.2, col='gray50')

# panel = function(...) { # this is a function that can make the background black for the levelplot
  # panel.fill(col = "black")
  # panel.levelplot(...)}

# Calculate the proportion change from beginning to end of century
prop_change <- pred.agg_biom[,c(1:3, 22)]
prop_change <- prop_change[prop_change$year_range=='2007-2020' | prop_change$year_range=='2081-2100',]
prop_change <- dcast(prop_change, latitude+longitude~year_range, value.var='ensMean_BIOM')
prop_change$diff <- prop_change$`2081-2100`-prop_change$`2007-2020`
prop_change$prop_change <- prop_change$diff/max(abs(prop_change$diff))

scale85_prop <- seq(-1.001, 1.001, length.out=20)
cols_prop = colorRampPalette(colors = c('black', 'dark red', 'red', 'gray92', 'blue', 'dark blue', 'black'))

plot2 <- levelplot(prop_change ~ longitude*latitude, data=prop_change,  colorkey=F, par.settings=list(layout.heights=list(top.padding=-3.5, bottom.padding=-3.5),layout.widths=list(left.padding=-2.5, right.padding=-2.5)),
          scales=list(tck=c(0,0), alternating=c(0,0)), strip=F, strip.left=FALSE,  aspect='fill', xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL,  
          at = scale85_prop, col.regions=cols_prop)  + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=0.2, col='gray50')

pdf(width=6.5, height=1.95, file=paste('/Users/jamesmorley/Desktop/', species, '_rcp', rcp, '_', niche, 'noLog_TEST.pdf', sep=''))
  print(plot1, split=c(1,1,2,1), more=TRUE) #position=c(0,0,.5,1), 
  print(plot2, split=c(2,1,2,1))#position=c(.5,0,1,1))
dev.off()

save(plot1, plot2, prop_change, pred.agg_biom, pred.agg_PA, niche, file=paste('figures/Uncertainty_MS_2018/Final_proj_figure_', species, '.RData', sep=''))



# The main plot showing present and proportional change in distribution ===========================================
# Use this if using a PA model, biomass models use the block above

if(species == 'homarus americanus_Atl'){ # drop south of NC and north of Newfoundland
  pred.agg_PA <- pred.agg_PA[pred.agg_PA$latitude < 52,]
  pred.agg_PA <- pred.agg_PA[pred.agg_PA$latitude > 34,]
  xlimit = c(min(pred.agg_PA$longitude), max(pred.agg_PA$longitude))
  ylimit = c(min(pred.agg_PA$latitude), max(pred.agg_PA$latitude))
  Projmap <- map('world', xlim=xlimit, ylim=ylimit, plot=F) 
  Projmap <- data.frame(lon=Projmap$x, lat=Projmap$y)
}
if(species == 'paralichthys dentatus_Atl'){
  pred.agg_PA <- pred.agg_PA[pred.agg_PA$latitude < 51,]
  pred.agg_PA <- pred.agg_PA[!pred.agg_PA$longitude < -82,]
  pred.agg_PA <- pred.agg_PA[!(pred.agg_PA$longitude < -80.75 & pred.agg_PA$latitude < 27),]
  xlimit = c(min(pred.agg_PA$longitude)-.5, max(pred.agg_PA$longitude))
  ylimit = c(min(pred.agg_PA$latitude), max(pred.agg_PA$latitude))
  Projmap <- map('world', xlim=xlimit, ylim=ylimit, plot=F) 
  Projmap <- data.frame(lon=Projmap$x, lat=Projmap$y)
}
if(species == 'centropristis striata_Atl'){
  pred.agg_PA <- pred.agg_PA[pred.agg_PA$latitude < 51,]
  xlimit = c(min(pred.agg_PA$longitude), max(pred.agg_PA$longitude))
  ylimit = c(min(pred.agg_PA$latitude), max(pred.agg_PA$latitude))
  Projmap <- map('world', xlim=xlimit, ylim=ylimit, plot=F) 
  Projmap <- data.frame(lon=Projmap$x, lat=Projmap$y)
}
if(species == 'anoplopoma fimbria_Pac'){
  xlimit = c(min(pred.agg_PA$longitude), max(pred.agg_PA$longitude) + 0.5)
  ylimit = c(min(pred.agg_PA$latitude), max(pred.agg_PA$latitude))
  Projmap <- map('world', xlim=xlimit, ylim=ylimit, plot=F) 
  Projmap <- data.frame(lon=Projmap$x, lat=Projmap$y)
}
if(species == 'doryteuthis opalescens_Pac'){
  xlimit = c(min(pred.agg_PA$longitude), max(pred.agg_PA$longitude) + 0.5)
  ylimit = c(min(pred.agg_PA$latitude), max(pred.agg_PA$latitude))
  Projmap <- map('world', xlim=xlimit, ylim=ylimit, plot=F) 
  Projmap <- data.frame(lon=Projmap$x, lat=Projmap$y)
}

scale85 <- seq(0, max(pred.agg_PA$ensMean_PA[pred.agg_PA$year_range=='2007-2020']) + .001, length.out=20)
cols = colorRampPalette(colors = c('gray92', 'blue', 'dark blue', 'black'))

plot1 <-  levelplot(ensMean_PA ~ longitude*latitude, data=pred.agg_PA[pred.agg_PA$year_range=='2007-2020',],  colorkey=F, par.settings=list(layout.heights=list(top.padding=-3.5, bottom.padding=-3.5),layout.widths=list(left.padding=-2.5, right.padding=-2.5)),
        scales=list(tck=c(0,0), alternating=c(0,0)), strip=F, strip.left=FALSE,  aspect='fill', xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL,  
        at = scale85, col.regions=cols) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=0.2, col='gray50')


# Calculate the proportion change from beginning to end of century
prop_change <- pred.agg_PA[,c(1:3, 22)]
prop_change <- prop_change[prop_change$year_range=='2007-2020' | prop_change$year_range=='2081-2100',]
prop_change <- dcast(prop_change, latitude+longitude~year_range, value.var='ensMean_PA')
prop_change$diff <- prop_change$`2081-2100`-prop_change$`2007-2020`
prop_change$prop_change <- prop_change$diff/max(abs(prop_change$diff))
  
scale85_prop <- seq(-1.001, 1.001, length.out=20)
cols_prop = colorRampPalette(colors = c('black', 'dark red', 'red', 'gray92', 'blue', 'dark blue', 'black'))

plot2 <- levelplot(prop_change ~ longitude*latitude, data=prop_change,  colorkey=F, par.settings=list(layout.heights=list(top.padding=-3.5, bottom.padding=-3.5),layout.widths=list(left.padding=-2.5, right.padding=-2.5)),
        scales=list(tck=c(0,0), alternating=c(0,0)), strip=F, strip.left=FALSE,  aspect='fill', xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL,  
        at = scale85_prop, col.regions=cols_prop)  + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=0.2, col='gray50')

pdf(width=6.5, height=1.95, file=paste('/Users/jamesmorley/Desktop/', species, '_rcp', rcp, '_', niche, 'noLog_TEST.pdf', sep=''))
  print(plot1, split=c(1,1,2,1), more=TRUE) #position=c(0,0,.5,1), 
  print(plot2, split=c(2,1,2,1))#position=c(.5,0,1,1))
dev.off()
 
save(plot1, plot2, prop_change, pred.agg_PA, pred.agg_biom, niche, file=paste('figures/Uncertainty_MS_2018/Final_proj_figure_', species, '.RData', sep=''))
 



# PUT THE FIGURES TOGETHER INTO THE MAIN FIGURE =================================================================
    # first the main MS figure then the appendix
 
pdf(width=6.5, height=8, file=paste('figures/Uncertainty_MS_2018/Projections/Final projection maps_Dec2018.pdf', sep=''))
load(paste('figures/Uncertainty_MS_2018/Final_proj_figure_', SPECIES[6], '.RData', sep=''))
print(plot1, split=c(1,1,2,4), more=TRUE) #position=c(0,0,.5,1), 
print(plot2, split=c(2,1,2,4), more=TRUE)#position=c(.5,0,1,1))
load(paste('figures/Uncertainty_MS_2018/Final_proj_figure_', SPECIES[7], '.RData', sep=''))
print(plot1, split=c(1,2,2,4), more=TRUE) #position=c(0,0,.5,1), 
print(plot2, split=c(2,2,2,4), more=TRUE)#position=c(.5,0,1,1))
load(paste('figures/Uncertainty_MS_2018/Final_proj_figure_', SPECIES[2], '.RData', sep=''))
print(plot1, split=c(1,3,2,4), more=TRUE) #position=c(0,0,.5,1), 
print(plot2, split=c(2,3,2,4), more=TRUE)#position=c(.5,0,1,1))
load(paste('figures/Uncertainty_MS_2018/Final_proj_figure_', SPECIES[1], '.RData', sep=''))
print(plot1, split=c(1,4,2,4), more=TRUE) #position=c(0,0,.5,1), 
print(plot2, split=c(2,4,2,4))#position=c(.5,0,1,1))
dev.off()

# Appendix
pdf(width=6.5, height=6, file=paste('figures/Uncertainty_MS_2018/Projections/Final projection maps_Appendix_Dec2018.pdf', sep=''))
load(paste('figures/Uncertainty_MS_2018/Final_proj_figure_', SPECIES[4], '.RData', sep=''))
print(plot1, split=c(1,1,2,3), more=TRUE) #position=c(0,0,.5,1), 
print(plot2, split=c(2,1,2,3), more=TRUE)#position=c(.5,0,1,1))
load(paste('figures/Uncertainty_MS_2018/Final_proj_figure_', SPECIES[5], '.RData', sep=''))
print(plot1, split=c(1,2,2,3), more=TRUE) #position=c(0,0,.5,1), 
print(plot2, split=c(2,2,2,3), more=TRUE)#position=c(.5,0,1,1))
load(paste('figures/Uncertainty_MS_2018/Final_proj_figure_', SPECIES[3], '.RData', sep=''))
print(plot1, split=c(1,3,2,3), more=TRUE) #position=c(0,0,.5,1), 
print(plot2, split=c(2,3,2,3), more=TRUE)#position=c(.5,0,1,1))
dev.off()

 
# ===============================================================================================
# Trying to add the graphics =================================================================================================
# Below is messing around with the 'magick' package. But couldn't get it to work well...loss of resolution
    # So can basically disregard all that is below

library(magick)

fig <- image_graph(width = 1000, height = 600) # 
  print(plot1, split=c(1,1,2,1), more=TRUE) #position=c(0,0,.5,1), 
  print(plot2, split=c(2,1,2,1))#position=c(.5,0,1,1))
dev.off()

rock <- image_read('Fish_graphics/widow.png')
rock <- image_scale(rock, geometry='200x')
out <- image_composite(fig, rock, offset = "+130+250")

image_write(out, path='/Users/jamesmorley/Desktop/test.pdf', format='pdf')

