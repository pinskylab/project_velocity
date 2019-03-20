library(ggplot2)
library(latticeExtra)
# bring in catch and hauls to look at the observed distribution of temperatures occupied
load('data/hauls_catch_Dec2017.RData')
species <- 'morone saxatilis_Atl'
catch <- dat[dat$sppocean==species,]
catch <- merge(hauls, catch, by='haulid', all=T)
catch$presfit[is.na(catch$presfit)] <- FALSE
regions <- names(table(catch$region[catch$presfit==T]))
catch <- catch[catch$region %in% regions,]

# compare temps at observations vs. what is present in the region
pdf(height=6, width=7, file="/Users/jamesmorley/Desktop/bias_investigate.pdf")
ggplot(catch, aes(x=bottemp, fill=presfit)) + geom_histogram(binwidth = 1) + facet_wrap(~region, scales='free_y')
ggplot(catch, aes(x=SBT.seasonal, fill=presfit)) + geom_histogram(binwidth = 1) + facet_wrap(~region, scales='free_y')
ggplot(catch, aes(x=SST.seasonal.mean, fill=presfit)) + geom_histogram(binwidth = 1) + facet_wrap(~region, scales='free_y')
ggplot(catch, aes(x=SBT.min, fill=presfit)) + geom_histogram(binwidth = 1) + facet_wrap(~region, scales='free_y')
ggplot(catch, aes(x=SBT.max, fill=presfit)) + geom_histogram(binwidth = 1) + facet_wrap(~region, scales='free_y')
ggplot(catch, aes(x=SST.max, fill=presfit)) + geom_histogram(binwidth = 1) + facet_wrap(~region, scales='free_y')
dev.off()

# this can investigate one model at a time =================================================================
  filename <- paste('data/prediction_files_Nov2017/predictionEASTnoBias_rcp85_jas_', modelrun[i], '.RData', sep='')
  load(filename)
  cde <- pred.bathE
  seafloor <- pred.bathE[,c('latBathgrid', 'lonBathgrid', 'GRAINSIZE', 'rugosity')]
  seafloor <- unique(seafloor)
  
  pred.bathE <- pred.bathE[pred.bathE$bin=='2081-2100',]
  
  pred.bathE <- pred.bathE[pred.bathE$latBathgrid > 36,]
  pred.bathE <- pred.bathE[!pred.bathE$lonBathgrid < -82,]
  pred.bathE <- pred.bathE[!pred.bathE$latBathgrid < 27,]
  pred.bathE <- pred.bathE[!pred.bathE$lonBathgrid > -68,]
  
  #pred.bathE <- pred.bathE[pred.bathE$lonBathgrid < -68,]
  abc <- aggregate(list(SBT.seasonal=pred.bathE$SBT.seasonal, SBT.min=pred.bathE$SBT.min, SBT.max=pred.bathE$SBT.max,
                SST.seasonal.mean=pred.bathE$SST.seasonal.mean, SST.max=pred.bathE$SST.max), 
                by=list(latBathgrid=pred.bathE$latBathgrid, lonBathgrid=pred.bathE$lonBathgrid), FUN=mean)
  ggplot(abc, aes(x=SBT.seasonal)) + geom_histogram(binwidth = 1)
  ggplot(abc, aes(x=SBT.min)) + geom_histogram(binwidth = 1)
  ggplot(abc, aes(x=SBT.max)) + geom_histogram(binwidth = 1)
  ggplot(abc, aes(x=SST.seasonal.mean)) + geom_histogram(binwidth = 1)
  ggplot(abc, aes(x=SST.max)) + geom_histogram(binwidth = 1)
  ggplot(abc, aes(x=GRAINSIZE)) + geom_histogram(binwidth = 1)
  
  cde <- cde[cde$bin=='2007-2020',]
  cde <- cde[!cde$lonBathgrid < -82,]
  cde <- cde[!cde$latBathgrid < 27,]
  cde <- aggregate(list(SBT.seasonal=cde$SBT.seasonal, SBT.min=cde$SBT.min, SBT.max=cde$SBT.max,
        SST.seasonal.mean=cde$SST.seasonal.mean, SST.max=cde$SST.max), 
        by=list(latBathgrid=cde$latBathgrid, lonBathgrid=cde$lonBathgrid), FUN=mean)
# =======================================================================================================  
  
model_gam <- paste('output/CEmodels_Nov2017/CEmods_Nov2017_fitallreg_2017_', species, '.RData', sep='')
load(model_gam)
mygam1 <- mods[[1]]
mygam2 <- mods[[2]]

abc$regionfact <- 'VIMS_NEAMAP'
abc <- merge(abc, seafloor, by=c('latBathgrid', 'lonBathgrid'), all.x=T)
load('data/ProjectionBathGrid_Feb27_2017.RData')# load bathymetry projection grid
proj.grid <- proj.grid[,c(1,2,5)]
abc <- merge(abc, proj.grid, by=c('latBathgrid', 'lonBathgrid'), all.x=T)
abc$depth <- abc$depth * -1
#abc$rugosity <- 1
#abc$GRAINSIZE <- 2
#abc$SST.seasonal.mean <- 15
#abc$SST.max <- 25


predPA <- predict(mygam1FULL, newdata = abc, type='response', na.action='na.pass')
abc$predPA <- predPA
predBIOM <- predict(mygam2FULL, newdata = abc, type='response', na.action='na.pass')
abc$predBIOM <- exp(predBIOM)

#levelplot(rugosity~lonBathgrid*latBathgrid, abc)
#levelplot(GRAINSIZE~lonBathgrid*latBathgrid, abc)
#levelplot(SST.max~lonBathgrid*latBathgrid, abc)

levelplot(predPA~lonBathgrid*latBathgrid, abc)
#levelplot(predPA~lonBathgrid*latBathgrid, abc[abc$predPA > .18,])
levelplot(log(predBIOM+1)~lonBathgrid*latBathgrid, abc)
levelplot(log((predBIOM*predPA)+1)~lonBathgrid*latBathgrid, abc)

cde <- merge(cde, seafloor, by=c('latBathgrid', 'lonBathgrid'), all.x=T)
cde$regionfact <- 'VIMS_NEAMAP'
predPA <- predict(mygam1, newdata = cde, type='response', na.action='na.pass')
cde$predPA <- predPA
predBIOM <- predict(mygam2, newdata = cde, type='response', na.action='na.pass')
cde$predBIOM <- exp(predBIOM)

levelplot(rugosity~lonBathgrid*latBathgrid, cde)
levelplot(GRAINSIZE~lonBathgrid*latBathgrid, cde)
levelplot(predPA~lonBathgrid*latBathgrid, cde)
levelplot(predBIOM~lonBathgrid*latBathgrid, cde)
levelplot(log((predBIOM*predPA)+1)~lonBathgrid*latBathgrid, cde)
