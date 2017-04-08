## A script to test whether thermal envelopes change through time (non-stationarity)
## Only for pres/abs model
## Works with Jim Morley's Feb 2017 data and model format
## Run after nonstationarity_projection1.r, which fits models to first decade of dataset and projects them later to see if the fit is still good
## This script examines summary stats on models to first decade of dataset and projected later

## Set working directories depending on computer
if(Sys.info()["nodename"] == "pinsky-macbookair"){
	setwd('~/Documents/Rutgers/Range projections/project_velocity/')
	modfolder <- '../CEmodels_nonstationaritytest/'
	nthreads=2
	}
if(Sys.info()["nodename"] == "amphiprion.deenr.rutgers.edu"){
	setwd('~/Documents/range_projections/')
	.libPaths(new='~/R/library/') # so that it can find my old packages
	modfolder <- 'CEmodels_nonstationaritytest/'
	nthreads=10
	}

# Load model evaluations
modeldiag <- read.csv('output/modeldiag_nonstationarity_projection2_1960+1970.csv')

# Trim out species with no models, and poor models
modeldiag <- modeldiag[!is.na(modeldiag$ntot),] # no model was fit
keepspp <- modeldiag$sppocean[modeldiag$decade==1960 & modeldiag$auc >= 0.7 & modeldiag$tss>=0.3]
length(sort(unique(modeldiag$sppocean))) # 79 spp with models fit
length(keepspp) # 23 spp to keep (1960), 26 (1960+1970)
modeldiag <- modeldiag[modeldiag$sppocean %in% keepspp,]

##########################################
# Plot stats over decades for each spp
##########################################
spp <- sort(unique(modeldiag$sppocean[!is.na(modeldiag$auc)]))

# On individual points
col=rgb(0, 0, 0, 0.2)
quartz(width=12, height=9)
# pdf(width=12, height=9, file='figures/nonstationarity_projection2_1960+1970.pdf')
par(mfrow=c(3,4))
	# auc
with(modeldiag[modeldiag$sppocean==spp[1],], plot(decade, auc, type='l', col=col, ylim=c(0,1)))
for(i in 2:length(spp)) with(modeldiag[modeldiag$sppocean==spp[i],], lines(decade, auc, col=col))
lines(aggregate(modeldiag$auc, by=list(modeldiag$decade), FUN=mean, na.rm=TRUE), lwd=2)
	# tss
with(modeldiag[modeldiag$sppocean==spp[1],], plot(decade, tss, type='l', col=col, ylim=c(0,1)))
for(i in 2:length(spp)) with(modeldiag[modeldiag$sppocean==spp[i],], lines(decade, tss, col=col))
lines(aggregate(modeldiag$tss, by=list(modeldiag$decade), FUN=mean, na.rm=TRUE), lwd=2)
	# acc
with(modeldiag[modeldiag$sppocean==spp[1],], plot(decade, acc, type='l', col=col, ylim=c(0,1)))
for(i in 2:length(spp)) with(modeldiag[modeldiag$sppocean==spp[i],], lines(decade, acc, col=col))
lines(aggregate(modeldiag$acc, by=list(modeldiag$decade), FUN=mean, na.rm=TRUE), lwd=2)
	# sens
with(modeldiag[modeldiag$sppocean==spp[1],], plot(decade, sens, type='l', col=col, ylim=c(0,1)))
for(i in 2:length(spp)) with(modeldiag[modeldiag$sppocean==spp[i],], lines(decade, sens, col=col))
lines(aggregate(modeldiag$sens, by=list(modeldiag$decade), FUN=mean, na.rm=TRUE), lwd=2)
	# spec
with(modeldiag[modeldiag$sppocean==spp[1],], plot(decade, spec, type='l', col=col, ylim=c(0,1)))
for(i in 2:length(spp)) with(modeldiag[modeldiag$sppocean==spp[i],], lines(decade, spec, col=col))
lines(aggregate(modeldiag$spec, by=list(modeldiag$decade), FUN=mean, na.rm=TRUE), lwd=2)
	# kappa
with(modeldiag[modeldiag$sppocean==spp[1],], plot(decade, kappa, type='l', col=col, ylim=c(0,1)))
for(i in 2:length(spp)) with(modeldiag[modeldiag$sppocean==spp[i],], lines(decade, kappa, col=col))
lines(aggregate(modeldiag$kappa, by=list(modeldiag$decade), FUN=mean, na.rm=TRUE), lwd=2)
	# tssmax
with(modeldiag[modeldiag$sppocean==spp[1],], plot(decade, tssmax, type='l', col=col, ylim=c(0,1)))
for(i in 2:length(spp)) with(modeldiag[modeldiag$sppocean==spp[i],], lines(decade, tssmax, col=col))
lines(aggregate(modeldiag$tssmax, by=list(modeldiag$decade), FUN=mean, na.rm=TRUE), lwd=2)
	# accmax
with(modeldiag[modeldiag$sppocean==spp[1],], plot(decade, accmax, type='l', col=col, ylim=c(0,1)))
for(i in 2:length(spp)) with(modeldiag[modeldiag$sppocean==spp[i],], lines(decade, accmax, col=col))
lines(aggregate(modeldiag$accmax, by=list(modeldiag$decade), FUN=mean, na.rm=TRUE), lwd=2)
	# kappamax
with(modeldiag[modeldiag$sppocean==spp[1],], plot(decade, kappamax, type='l', col=col, ylim=c(0,1)))
for(i in 2:length(spp)) with(modeldiag[modeldiag$sppocean==spp[i],], lines(decade, kappamax, col=col))
lines(aggregate(modeldiag$kappamax, by=list(modeldiag$decade), FUN=mean, na.rm=TRUE), lwd=2)
	# rpb
with(modeldiag[modeldiag$sppocean==spp[1],], plot(decade, rpb, type='l', col=col, ylim=c(0,1)))
for(i in 2:length(spp)) with(modeldiag[modeldiag$sppocean==spp[i],], lines(decade, rpb, col=col))
lines(aggregate(modeldiag$rpb, by=list(modeldiag$decade), FUN=mean, na.rm=TRUE), lwd=2)

dev.off()

# On gridded data
col=rgb(0, 0, 0, 0.2)
quartz(width=12, height=9)
# pdf(width=12, height=9, file='figures/nonstationarity_projection2_bygrid.pdf')
par(mfrow=c(3,4))
	# auc
with(modeldiag[modeldiag$sppocean==spp[1],], plot(decade, auc.gr, type='l', col=col, ylim=c(0,1)))
for(i in 2:length(spp)) with(modeldiag[modeldiag$sppocean==spp[i],], lines(decade, auc.gr, col=col))
lines(aggregate(modeldiag$auc.gr, by=list(modeldiag$decade), FUN=mean, na.rm=TRUE), lwd=2)
	# tss
with(modeldiag[modeldiag$sppocean==spp[1],], plot(decade, tss.gr, type='l', col=col, ylim=c(0,1)))
for(i in 2:length(spp)) with(modeldiag[modeldiag$sppocean==spp[i],], lines(decade, tss.gr, col=col))
lines(aggregate(modeldiag$tss.gr, by=list(modeldiag$decade), FUN=mean, na.rm=TRUE), lwd=2)
	# acc
with(modeldiag[modeldiag$sppocean==spp[1],], plot(decade, acc.gr, type='l', col=col, ylim=c(0,1)))
for(i in 2:length(spp)) with(modeldiag[modeldiag$sppocean==spp[i],], lines(decade, acc.gr, col=col))
lines(aggregate(modeldiag$acc.gr, by=list(modeldiag$decade), FUN=mean, na.rm=TRUE), lwd=2)
	# sens
with(modeldiag[modeldiag$sppocean==spp[1],], plot(decade, sens.gr, type='l', col=col, ylim=c(0,1)))
for(i in 2:length(spp)) with(modeldiag[modeldiag$sppocean==spp[i],], lines(decade, sens.gr, col=col))
lines(aggregate(modeldiag$sens.gr, by=list(modeldiag$decade), FUN=mean, na.rm=TRUE), lwd=2)
	# spec
with(modeldiag[modeldiag$sppocean==spp[1],], plot(decade, spec.gr, type='l', col=col, ylim=c(0,1)))
for(i in 2:length(spp)) with(modeldiag[modeldiag$sppocean==spp[i],], lines(decade, spec.gr, col=col))
lines(aggregate(modeldiag$spec.gr, by=list(modeldiag$decade), FUN=mean, na.rm=TRUE), lwd=2)
	# kappa
with(modeldiag[modeldiag$sppocean==spp[1],], plot(decade, kappa.gr, type='l', col=col, ylim=c(0,1)))
for(i in 2:length(spp)) with(modeldiag[modeldiag$sppocean==spp[i],], lines(decade, kappa.gr, col=col))
lines(aggregate(modeldiag$kappa.gr, by=list(modeldiag$decade), FUN=mean, na.rm=TRUE), lwd=2)
	# tssmax
with(modeldiag[modeldiag$sppocean==spp[1],], plot(decade, tssmax.gr, type='l', col=col, ylim=c(0,1)))
for(i in 2:length(spp)) with(modeldiag[modeldiag$sppocean==spp[i],], lines(decade, tssmax.gr, col=col))
lines(aggregate(modeldiag$tssmax.gr, by=list(modeldiag$decade), FUN=mean, na.rm=TRUE), lwd=2)
	# accmax
with(modeldiag[modeldiag$sppocean==spp[1],], plot(decade, accmax.gr, type='l', col=col, ylim=c(0,1)))
for(i in 2:length(spp)) with(modeldiag[modeldiag$sppocean==spp[i],], lines(decade, accmax.gr, col=col))
lines(aggregate(modeldiag$accmax.gr, by=list(modeldiag$decade), FUN=mean, na.rm=TRUE), lwd=2)
	# kappamax
with(modeldiag[modeldiag$sppocean==spp[1],], plot(decade, kappamax.gr, type='l', col=col, ylim=c(0,1)))
for(i in 2:length(spp)) with(modeldiag[modeldiag$sppocean==spp[i],], lines(decade, kappamax.gr, col=col))
lines(aggregate(modeldiag$kappamax.gr, by=list(modeldiag$decade), FUN=mean, na.rm=TRUE), lwd=2)
	# rpb
with(modeldiag[modeldiag$sppocean==spp[1],], plot(decade, rpb.gr, type='l', col=col, ylim=c(0,1)))
for(i in 2:length(spp)) with(modeldiag[modeldiag$sppocean==spp[i],], lines(decade, rpb.gr, col=col))
lines(aggregate(modeldiag$rpb.gr, by=list(modeldiag$decade), FUN=mean, na.rm=TRUE), lwd=2)

dev.off()