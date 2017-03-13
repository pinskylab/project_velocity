# This script is used in conjunction with the main data files and some output of script '5.5_sdm_uncertainty.R

# This figure plots color coded survey effort
pdf(width=15, height=7, file='/Users/jamesmorley/Desktop/haul_coords_bycolor.pdf')
par(mai=c(0.7,0.7,0.3, 0.1), las=1, mgp=c(2,1,0))
map(database='world', fill=TRUE, col='dark gray', border=FALSE, xlim=range(hauls$lon), ylim=range(hauls$lat), xlab='Longitude', ylab='Latitude')
points(hauls$lon[hauls$region=='AFSC_Aleutians'], hauls$lat[hauls$region=='AFSC_Aleutians'], col='red', pch=20, cex=0.1)
points(hauls$lon[hauls$region=='AFSC_EBS'], hauls$lat[hauls$region=='AFSC_EBS'], col='blue', pch=20, cex=0.1)
points(hauls$lon[hauls$region=='AFSC_GOA'], hauls$lat[hauls$region=='AFSC_GOA'], col='green', pch=20, cex=0.1)
points(hauls$lon[hauls$region=='AFSC_WCTri'], hauls$lat[hauls$region=='AFSC_WCTri'], col='orange', pch=20, cex=0.1)
points(hauls$lon[hauls$region=='DFO_Newfoundland'], hauls$lat[hauls$region=='DFO_Newfoundland'], col='red', pch=20, cex=0.1)
points(hauls$lon[hauls$region=='DFO_ScotianShelf'], hauls$lat[hauls$region=='DFO_ScotianShelf'], col='blue', pch=20, cex=0.1)
points(hauls$lon[hauls$region=='DFO_SoGulf'], hauls$lat[hauls$region=='DFO_SoGulf'], col='green', pch=20, cex=0.1)
points(hauls$lon[hauls$region=='NEFSC_NEUS'], hauls$lat[hauls$region=='NEFSC_NEUS'], col='orange', pch=20, cex=0.1)
points(hauls$lon[hauls$region=='NWFSC_WCAnn'], hauls$lat[hauls$region=='NWFSC_WCAnn'], col='orange', pch=20, cex=0.1)
points(hauls$lon[hauls$region=='SCDNR_SEUS'], hauls$lat[hauls$region=='SCDNR_SEUS'], col='black', pch=20, cex=0.1)
points(hauls$lon[hauls$region=='SEFSC_GOMex'], hauls$lat[hauls$region=='SEFSC_GOMex'], col='red', pch=20, cex=0.1)
points(hauls$lon[hauls$region=='VIMS_NEAMAP'], hauls$lat[hauls$region=='VIMS_NEAMAP'], col='blue', pch=20, cex=0.1)
dev.off()
 
# ==============================================================================================================================

# This plot produces habitat preference curves using fitted gam models_currently this is tailored for "menippe mercenaria_Atl"
newdat <- data.frame(SBT.seasonal=seq(-1, 29, by=0.1), SST.seasonal.mean=rep(mean(haulsMod$SST.seasonal.mean[haulsMod$presfit==T]), 301), SBT.min=rep(mean(haulsMod$SBT.min[haulsMod$presfit==T]), 301), SBT.max=rep(mean(haulsMod$SBT.max[haulsMod$presfit==T]),301), SST.max=rep(mean(haulsMod$SST.max[haulsMod$presfit==T]), 301), rugosity=rep(mean(haulsMod$rugosity[haulsMod$presfit==T]),301), GRAINSIZE=rep(mean(haulsMod$GRAINSIZE[haulsMod$presfit==T]),301), regionfact='SCDNR_SEUS', habitatFact='soft')
SBT.pred <- data.frame(cbind(newdat$SBT.seasonal, exp(predict(mod2, newdata = newdat, type='response')), predict(mod1, newdata=newdat, type='response')))
SBT.pred$curve <- SBT.pred$X2 * SBT.pred$X3

newdat <- data.frame(SBT.seasonal=rep(mean(haulsMod$SBT.seasonal[haulsMod$presfit==T]), 251), SST.seasonal.mean=rep(mean(haulsMod$SST.seasonal.mean[haulsMod$presfit==T]), 251), SBT.min=seq(-2, 23, by=0.1), SBT.max=rep(mean(haulsMod$SBT.max[haulsMod$presfit==T]),251), SST.max=rep(mean(haulsMod$SST.max[haulsMod$presfit==T]), 251), rugosity=rep(mean(haulsMod$rugosity[haulsMod$presfit==T]),251), GRAINSIZE=rep(mean(haulsMod$GRAINSIZE[haulsMod$presfit==T]),251), regionfact='SCDNR_SEUS', habitatFact='soft')
SBTmin.pred <- data.frame(cbind(newdat$SBT.min, exp(predict(mod2, newdata = newdat, type='response')), predict(mod1, newdata=newdat, type='response')))
SBTmin.pred$curve <- SBTmin.pred$X2 * SBTmin.pred$X3

newdat <- data.frame(SBT.seasonal=rep(mean(haulsMod$SBT.seasonal[haulsMod$presfit==T]), 331), SST.seasonal.mean=rep(mean(haulsMod$SST.seasonal.mean[haulsMod$presfit==T]), 331), SBT.min=rep(mean(haulsMod$SBT.min[haulsMod$presfit==T]), 331), SBT.max=seq(-1, 32, by=0.1), SST.max=rep(mean(haulsMod$SST.max[haulsMod$presfit==T]), 331), rugosity=rep(mean(haulsMod$rugosity[haulsMod$presfit==T]),331), GRAINSIZE=rep(mean(haulsMod$GRAINSIZE[haulsMod$presfit==T]),331), regionfact='SCDNR_SEUS', habitatFact='soft')
SBTmax.pred <- data.frame(cbind(newdat$SBT.max, exp(predict(mod2, newdata = newdat, type='response')), predict(mod1, newdata=newdat, type='response')))
SBTmax.pred$curve <- SBTmax.pred$X2 * SBTmax.pred$X3
 
newdat <- data.frame(SBT.seasonal=rep(mean(haulsMod$SBT.seasonal[haulsMod$presfit==T]), 111), SST.seasonal.mean=rep(mean(haulsMod$SST.seasonal.mean[haulsMod$presfit==T]), 111), SBT.min=rep(mean(haulsMod$SBT.min[haulsMod$presfit==T]), 111), SBT.max=rep(mean(haulsMod$SBT.max[haulsMod$presfit==T]),111), SST.max=rep(mean(haulsMod$SST.max[haulsMod$presfit==T]), 111), rugosity=seq(0, 5.5, by=0.05), GRAINSIZE=rep(mean(haulsMod$GRAINSIZE[haulsMod$presfit==T]),111), regionfact='SCDNR_SEUS', habitatFact='soft')
SBTrug.pred <- data.frame(cbind(newdat$rugosity, exp(predict(mod2, newdata = newdat, type='response')), predict(mod1, newdata=newdat, type='response')))
SBTrug.pred$curve <- SBTrug.pred$X2 * SBTrug.pred$X3

newdat <- data.frame(SBT.seasonal=rep(mean(haulsMod$SBT.seasonal[haulsMod$presfit==T]), 142), SST.seasonal.mean=rep(mean(haulsMod$SST.seasonal.mean[haulsMod$presfit==T]), 142), SBT.min=rep(mean(haulsMod$SBT.min[haulsMod$presfit==T]), 142), SBT.max=rep(mean(haulsMod$SBT.max[haulsMod$presfit==T]),142), SST.max=rep(mean(haulsMod$SST.max[haulsMod$presfit==T]), 142), rugosity=rep(mean(haulsMod$rugosity[haulsMod$presfit==T]),142), GRAINSIZE=seq(-5.1, 9, by=0.1), regionfact='SCDNR_SEUS', habitatFact='soft')
SBTgrain.pred <- data.frame(cbind(newdat$GRAINSIZE, exp(predict(mod2, newdata = newdat, type='response')), predict(mod1, newdata=newdat, type='response')))
SBTgrain.pred$curve <- SBTgrain.pred$X2 * SBTgrain.pred$X3

pdf(width=10, height=7, file='/Users/jamesmorley/Desktop/spiny dogs_curves.pdf')
par(mfrow=c(2,3), mai=c(0.7,0.7,0.3, 0.1), las=0, mgp=c(2,1,0), cex=1.1)
plot(curve~X1, ylim=c(0,7), col='dark blue', xlab='Temperature at capture', ylab='Predicted CPUE', data=SBT.pred)
plot(curve~X1, ylim=c(0,7), col='dark blue', xlab='Minimum annual temperature', ylab='Predicted CPUE', data=SBTmin.pred)
plot(curve~X1, ylim=c(0,7), col='dark blue', xlab='Maximum annual temperature', ylab='Predicted CPUE', data=SBTmax.pred)
plot(curve~X1, ylim=c(0,7), col='dark blue', xlab='Seafloor rugosity', ylab='Predicted CPUE', data=SBTrug.pred)
plot(curve~X1, ylim=c(0,7), col='dark blue', xlab='Sediment grain size (inverse)', ylab='Predicted CPUE', data=SBTgrain.pred)
dev.off()

