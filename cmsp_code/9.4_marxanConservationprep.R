# Set up a Marxan with Zones run

## Set working directories
if(Sys.info()["nodename"] == "pinsky-macbookair"){
	setwd('~/Documents/Rutgers/Range projections/proj_ranges/')
	marxfolder <- '../MarZone_runs/'
	}
if(Sys.info()["nodename"] == "amphiprion.deenr.rutgers.edu"){
	setwd('~/Documents/range_projections/')
	.libPaths(new='~/R/x86_64-redhat-linux-gnu-library/3.1/') # so that it can find my old packages
	marxfolder <- 'MarZone_runs/'
	}
# could add code for Lauren's working directory here


############
## Flags
############

## choose which runs to use
## runtype refers to how the Species Distribution Models (SDMs) were fit
## projtype refers to how the SDM projections were done
#runtype <- 'test'; projtype=''
#runtype <- ''; projtype=''`
#runtype <- 'testK6noSeas'; projtype='_xreg'
runtype <- 'fitallreg'; projtype='_xreg'


####################
## helper functions
####################
require(RColorBrewer)

# normalize to 0-1
norm01 <- function(x){
	mn <- min(x)
	mx <- max(x)
	return((x-mn)/(mx-mn))
}

# parallel normalize to 0-1
# use y to guide the scale
pnorm01 <- function(x, y){
	mn <- min(y)
	mx <- max(y)
	return((x-mn)/(mx-mn))
}


############################################
## Set up a Marxan run just on 2006-2020
## use just one rcp
############################################
# choose the rcp
rcp <- 85


load(paste('data/biomassavemap_', runtype, projtype, '_rcp', rcp, '.RData', sep='')) # loads biomassavemap data.frame
	dim(biomassavemap)
load(paste('data/presmap_', runtype, projtype, '_rcp', rcp, '.RData', sep='')) # loads presmap data.frame with presence/absence information

goal <- 0.2 # proportion to capture in conservation
runname <- 'conservationtest'

# Just one region for now: NEUS Spring spring
	# pu.dat
	# planning features are each 1/4 deg square
	pus <- presmap[presmap$region=='NEFSC_NEUSSpring',c('lat', 'lon')]
	pus <- pus[!duplicated(pus),]
		dim(pus)
	pus <- pus[order(pus$lat, pus$lon),]
	pus$id <- 1:nrow(pus)
#	pus$dummycost <- runif(nrow(pus), 0, 1)
	pus$dummycost <- rep(0.1, nrow(pus))
	
#	pus <- pus[1:50,] # trim as a test
		
	pu.dat<-pus[,c('id', 'dummycost')]

	save(pus, file=paste(marxfolder, 'input', runname, '/pus.Rdata', sep=''))	
	write.csv(pu.dat, file=paste(marxfolder, 'input', runname, '/pu.dat', sep=''), row.names=FALSE, quote=FALSE)
	
	# spec.dat (takes the place of feat.dat?)
	# goal for every species present
	# not documented in 1.0.1 manual. copying format from unix example
	spps <- presmap[presmap$region=='NEFSC_NEUSSpring',c('sppocean')]
	spps <- spps[!duplicated(spps)]
		length(spps) # 123

	sppstokeep <- presmap[presmap$region=='NEFSC_NEUSSpring' & presmap$period=='2006-2020' & presmap$pres==TRUE,c('lat', 'lon', 'sppocean', 'pres')]
		dim(sppstokeep)
	sppstokeep <- merge(sppstokeep, pus[,c('lat', 'lon', 'id')]) # add pu id (and trim to focal pus)
		names(sppstokeep)[names(sppstokeep)=='id'] <- 'pu'
		dim(sppstokeep)
	
		ngrid <- aggregate(list(ngrid=sppstokeep$pu), by=list(sppocean=sppstokeep$sppocean), FUN=function(x) length(unique(x)))
		sppstokeep <- merge(sppstokeep, ngrid)
		summary(sppstokeep$ngrid) # 1 to 50

		nspps <- aggregate(list(nspp=sppstokeep$sppocean), by=list(pu=sppstokeep$pu), FUN=function(x) length(unique(x)))
		sppstokeep <- merge(sppstokeep, nspps)
		summary(sppstokeep$nspp) # 56 to 79

		sppstokeep <- sppstokeep[sppstokeep$ngrid>10,]

		length(unique(sppstokeep$sppocean)) # 114

	spps <- spps[spps %in% sppstokeep$sppocean]
		length(spps) # 121

#	spps <- data.frame(id=1:length(spps), prop=rep(goal, length(spps)), spf=rep(10000, length(spps)), name=gsub(' |_', '', spps)) #  fill spaces in species names.
	spps <- data.frame(id=1:length(spps), name=gsub(' |_', '', spps), sppocean=spps) #  fill spaces in species names.

#	spps <- spps[1:20,] # trim as a test

	spps.dat <- spps[,c('id', 'name')]
	
	save(spps, file=paste(marxfolder, 'input', runname, '/spps.Rdata', sep=''))		
	write.csv(spps.dat, file=paste(marxfolder, 'input', runname, '/spec.dat', sep=''), row.names=FALSE, quote=FALSE)

	# puvsp.dat (takes the place of puvfeat.dat?)
	# species by planning units
	# not documented in 1.0.1 manual. copying format from unix example
	puvsp <- presmap[presmap$region=='NEFSC_NEUSSpring' & presmap$period=='2006-2020',c('lat', 'lon', 'sppocean', 'pres')]
		dim(puvsp)
	puvsp <- merge(puvsp, pus[,c('lat', 'lon', 'id')]) # add pu id (and trim to pus)
		dim(puvsp)
		names(puvsp)[names(puvsp)=='id'] <- 'pu'


	puvsp$name <- gsub(' |_', '', puvsp$sppocean)
	puvsp <- merge(puvsp, spps[,c('id', 'name')])
		dim(puvsp)
		names(puvsp)[names(puvsp)=='id'] <- 'species'
	puvsp$amount <- as.numeric(puvsp$pres)
	
	length(unique(puvsp$pu))
	length(unique(puvsp$species))

	puvsp.dat <- puvsp[,c('species', 'pu', 'amount')]
		puvsp.dat <- puvsp.dat[order(puvsp.dat$pu, puvsp.dat$species),]
		puvsp.dat <- puvsp.dat[puvsp.dat$amount>0,] # trim only to presences

		table(puvsp.dat$species) # make sure all species show up in some planning units
		table(puvsp.dat$pu) # make sure all planning units have some species

		table(puvsp.dat$pu, puvsp.dat$species)

	write.csv(puvsp.dat, file=paste(marxfolder, 'input', runname, '/puvsp.dat', sep=''), row.names=FALSE, quote=FALSE)

	# zones
	zones <- data.frame(zoneid=1:2, zonename=c('available', 'conservation'))
	
	write.csv(zones, file=paste(marxfolder, 'input', runname, '/zones.dat', sep=''), row.names=FALSE, quote=FALSE)

	#costs
	costs <- data.frame(costid=1, costname='dummycost')
	
	write.csv(costs, file=paste(marxfolder, 'input', runname, '/costs.dat', sep=''), row.names=FALSE, quote=FALSE)
	
	#zone cost
	zonecost <- expand.grid(list(zoneid=zones$zoneid, costid=costs$costid))
	zonecost$multiplier <- 0
	zonecost$multiplier[zonecost$zoneid==2 & zonecost$costid==1] <- 1
	
	zonecost <- zonecost[zonecost$multiplier>0,]
	
	write.csv(zonecost, file=paste(marxfolder, 'input', runname, '/zonecost.dat', sep=''), row.names=FALSE, quote=FALSE)
	
	#boundary length
	# optional
	
	#zone boundary cost
	# optional
	
	#planning unit zone
	#optional
	
	#planning unit lock
	#optional
	
	#zone target
	# set zone-specific targets
	zonetarget <- expand.grid(list(zoneid=zones$zoneid, speciesid=spps$id))
	zonetarget$target <- 0
#	zonetarget$target[zonetarget$zoneid==1 & zonetarget$speciesid==1] <- 0.2
	zonetarget$target[zonetarget$zoneid==2] <- goal # XX proportion in conservation zones
	zonetarget$targettype <- 1 # 3: proportion of total occurrences. 1: proportion of total amount
	
	zonetarget.dat <- zonetarget[zonetarget$target>0,]
	zonetarget.dat <- zonetarget.dat[order(zonetarget.dat$zoneid, zonetarget.dat$speciesid),]
	
	write.csv(zonetarget.dat, file=paste(marxfolder, 'input', runname, '/zonetarget.dat', sep=''), row.names=FALSE, quote=FALSE)


	#input parameters
	input <- data.frame(BLM=0, PROP=0.5, RANDSEED=-1, NUMREPS=100, AVAILABLEZONE=1, NUMITNS='1000000', STARTTEMP=-1, NUMTEMP='10000', COSTTHRESH=0, THRESHPEN1=14, THRESHPEN2=10, INPUTDIR=paste('input', runname, sep=''), PUNAME='pu.dat', SPECNAME='spec.dat', PUVSPRNAME='puvsp.dat', ZONESNAME='zones.dat', COSTSNAME='costs.dat', ZONECOSTNAME='zonecost.dat', ZONETARGETNAME='zonetarget.dat', SCENNAME=runname, SAVERUN=3, SAVEBEST=3, SAVESUMMARY=3, SAVESCEN=3, SAVETARGMET=3, SAVESUMSOLN=3, SAVEPENALTY=3, SAVELOG=3, OUTPUTDIR=paste('output', runname, sep=''),  RUNMODE=1, MISSLEVEL=1, ITIMPTYPE=0, HEURTYPE=-1, CLUMPTYPE=0, VERBOSITY=3, SAVESOLUTIONSMATRIX=3, SAVEANNEALINGTRACE=0, ANNEALINGTRACEROWS=1000)


	write.table(t(input), file=paste(marxfolder, 'input.dat', sep=''), row.names=TRUE, quote=FALSE, sep=' ', col.names=FALSE)
	write.table(t(input), file=paste(marxfolder, 'input', runname, '/input.dat', sep=''), row.names=TRUE, quote=FALSE, sep=' ', col.names=FALSE)
	
	
# Go run MarZone!

##################################################################
## Set up a Marxan run on 2006-2020 and ensemble mean 2081-2100
##################################################################
load('data/biomassavemap_testK6noSeas.RData') # loads biomassavemap data.frame
	dim(biomassavemap)
load('data/presmap.RData') # loads presmap data.frame with presence/absence information

goal <- 0.2 # proportion to capture in conservation
runname <- 'conservationtest2per'

# Just one region for now: NEUS Spring spring
	# pu.dat
	# planning features are each 1/4 deg square
	pus2 <- presmap[presmap$region=='NEFSC_NEUSSpring',c('lat', 'lon')]
	pus2 <- pus2[!duplicated(pus2),]
		dim(pus2)
	pus2 <- pus2[order(pus2$lat, pus2$lon),]
	pus2$id <- 1:nrow(pus2)
#	pus2$dummycost <- runif(nrow(pus2), 0, 1)
	pus2$dummycost <- rep(0.1, nrow(pus2))
	
#	pus2 <- pus2[1:50,] # trim as a test
		
	pu2.dat<-pus2[,c('id', 'dummycost')]
	
	save(pus2, file=paste(marxfolder, 'input', runname, '/pus.Rdata', sep=''))	
	write.csv(pu2.dat, file=paste(marxfolder, 'input', runname, '/pu.dat', sep=''), row.names=FALSE, quote=FALSE)
	
	# spec.dat (takes the place of feat.dat?)
	# goal for every species present
	# not documented in 1.0.1 manual. copying format from unix example
	inds <- presmap$region=='NEFSC_NEUSSpring' & presmap$period %in% c('2006-2020', '2081-2100')
	spps2 <- paste(presmap$sppocean[inds], presmap$period[inds])
	spps2 <- spps2[!duplicated(spps2)]
		length(spps2) # 246

	sppstokeep2 <- presmap[presmap$region=='NEFSC_NEUSSpring' & presmap$period %in% c('2006-2020', '2081-2100') & presmap$pres==TRUE,c('lat', 'lon', 'period', 'sppocean', 'pres')]
		dim(sppstokeep2)
	sppstokeep2 <- merge(sppstokeep2, pus2[,c('lat', 'lon', 'id')]) # add pu id (and trim to focal pus2)
		names(sppstokeep2)[names(sppstokeep2)=='id'] <- 'pu'
		dim(sppstokeep2)
	
		ngrid2 <- aggregate(list(ngrid=sppstokeep2$pu), by=list(sppocean=sppstokeep2$sppocean, period=sppstokeep2$period), FUN=function(x) length(unique(x)))
		sppstokeep2 <- merge(sppstokeep2, ngrid2)
		summary(sppstokeep2$ngrid) # 3 to 552

		nspps2 <- aggregate(list(nspp=sppstokeep2$sppocean), by=list(pu=sppstokeep2$pu, period=sppstokeep2$period), FUN=function(x) length(unique(x)))
		sppstokeep2 <- merge(sppstokeep2, nspps2)
		summary(sppstokeep2$nspp) # 4 to 85

		sppstokeep2 <- sppstokeep2[sppstokeep2$ngrid>10,] # trim to species at least minimally common

		length(unique(paste(sppstokeep2$sppocean, sppstokeep2$period))) # 242

	spps2 <- spps2[spps2 %in% paste(sppstokeep2$sppocean, sppstokeep2$period)]
		length(spps2) # 242

#	spps2 <- data.frame(id=1:length(spps2), prop=rep(goal, length(spps2)), spf=rep(10000, length(spps2)), name=gsub(' |_', '', spps2)) #  fill spaces in species names.
	spps2 <- data.frame(id=1:length(spps2), name=gsub(' |_|-', '', spps2), sppocean=spps2) #  remove spaces and dashes in species names.

#	spps2 <- spps2[1:20,] # trim as a test

	spps2.dat <- spps2[,c('id', 'name')]
	
	save(spps2, file=paste(marxfolder, 'input', runname, '/spps.Rdata', sep=''))		
	write.csv(spps2.dat, file=paste(marxfolder, 'input', runname, '/spec.dat', sep=''), row.names=FALSE, quote=FALSE)

	# puvsp.dat (takes the place of puvfeat.dat?)
	# species by planning units
	# not documented in 1.0.1 manual. copying format from unix example
	puvsp2 <- presmap[presmap$region=='NEFSC_NEUSSpring' & presmap$period %in% c('2006-2020', '2081-2100'),c('lat', 'lon', 'period', 'sppocean', 'pres')]
		dim(puvsp2)
	puvsp2 <- merge(puvsp2, pus2[,c('lat', 'lon', 'id')]) # add pu id (and trim to focal pus2)
		dim(puvsp2)
		names(puvsp2)[names(puvsp2)=='id'] <- 'pu'


	puvsp2$name <- gsub(' |_|-', '', paste(puvsp2$sppocean, puvsp2$period))
		dim(puvsp2)
	puvsp2 <- merge(puvsp2, spps2[,c('id', 'name')]) # add species id and trim to focal species
		dim(puvsp2)
		names(puvsp2)[names(puvsp2)=='id'] <- 'species'
	puvsp2$amount <- as.numeric(puvsp2$pres)
	
	length(unique(puvsp2$pu)) # 552
	length(unique(puvsp2$species)) # 242

	puvsp2.dat <- puvsp2[,c('species', 'pu', 'amount')]
		puvsp2.dat <- puvsp2.dat[order(puvsp2.dat$pu, puvsp2.dat$species),]
		puvsp2.dat <- puvsp2.dat[puvsp2.dat$amount>0,] # trim only to presences

		table(puvsp2.dat$species) # make sure all species show up in some planning units
			range(table(puvsp2.dat$species)) # make sure all species show up in some planning units
		table(puvsp2.dat$pu) # make sure all planning units have some species
			range(table(puvsp2.dat$pu)) # make sure all planning units have some species

#		table(puvsp2.dat$pu, puvsp2.dat$species) # giant matrix of who is where

	write.csv(puvsp2.dat, file=paste(marxfolder, 'input', runname, '/puvsp.dat', sep=''), row.names=FALSE, quote=FALSE)

	# zones
	zones2 <- data.frame(zoneid=1:2, zonename=c('available', 'conservation'))
	
	write.csv(zones2, file=paste(marxfolder, 'input', runname, '/zones.dat', sep=''), row.names=FALSE, quote=FALSE)

	#costs
	costs2 <- data.frame(costid=1, costname='dummycost')
	
	write.csv(costs2, file=paste(marxfolder, 'input', runname, '/costs.dat', sep=''), row.names=FALSE, quote=FALSE)
	
	#zone cost
	zonecost2 <- expand.grid(list(zoneid=zones$zoneid, costid=costs$costid))
	zonecost2$multiplier <- 0
	zonecost2$multiplier[zonecost2$zoneid==2 & zonecost2$costid==1] <- 1
	
	zonecost2 <- zonecost2[zonecost2$multiplier>0,]
	
	write.csv(zonecost2, file=paste(marxfolder, 'input', runname, '/zonecost.dat', sep=''), row.names=FALSE, quote=FALSE)
	
	#boundary length
	# optional
	
	#zone boundary cost
	# optional
	
	#planning unit zone
	#optional
	
	#planning unit lock
	#optional
	
	#zone target
	# set zone-specific targets
	zonetarget2 <- expand.grid(list(zoneid=zones$zoneid, speciesid=spps2$id))
	zonetarget2$target <- 0
#	zonetarget2$target[zonetarget2$zoneid==1 & zonetarget2$speciesid==1] <- 0.2
	zonetarget2$target[zonetarget2$zoneid==2] <- goal # XX proportion in conservation zones
	zonetarget2$targettype <- 1 # 3: proportion of total occurrences. 1: proportion of total amount
	
	zonetarget2.dat <- zonetarget2[zonetarget2$target>0,]
	zonetarget2.dat <- zonetarget2.dat[order(zonetarget2.dat$zoneid, zonetarget2.dat$speciesid),]
	
	write.csv(zonetarget2.dat, file=paste(marxfolder, 'input', runname, '/zonetarget.dat', sep=''), row.names=FALSE, quote=FALSE)


	#input parameters
	input2 <- data.frame(BLM=0, PROP=0.5, RANDSEED=-1, NUMREPS=100, AVAILABLEZONE=1, NUMITNS='1000000', STARTTEMP=-1, NUMTEMP='10000', COSTTHRESH=0, THRESHPEN1=14, THRESHPEN2=10, INPUTDIR=paste('input', runname, sep=''), PUNAME='pu.dat', SPECNAME='spec.dat', PUVSPRNAME='puvsp.dat', ZONESNAME='zones.dat', COSTSNAME='costs.dat', ZONECOSTNAME='zonecost.dat', ZONETARGETNAME='zonetarget.dat', SCENNAME=runname, SAVERUN=3, SAVEBEST=3, SAVESUMMARY=3, SAVESCEN=3, SAVETARGMET=3, SAVESUMSOLN=3, SAVEPENALTY=3, SAVELOG=3, OUTPUTDIR=paste('output', runname, sep=''),  RUNMODE=1, MISSLEVEL=1, ITIMPTYPE=0, HEURTYPE=-1, CLUMPTYPE=0, VERBOSITY=3, SAVESOLUTIONSMATRIX=3, SAVEANNEALINGTRACE=0, ANNEALINGTRACEROWS=1000)


	write.table(t(input2), file=paste(marxfolder, 'input.dat', sep=''), row.names=TRUE, quote=FALSE, sep=' ', col.names=FALSE)
	write.table(t(input2), file=paste(marxfolder, 'input', runname, '/input.dat', sep=''), row.names=TRUE, quote=FALSE, sep=' ', col.names=FALSE)
	
	
# Go run MarZone!	


################################
# Read in results and analyze against ensemble mean
################################
runname <- 'conservationtest'
runname2 <- 'conservationtest2per'

consplan <- read.csv(paste(marxfolder, 'output', runname, '/', runname, '_best.csv', sep=''))
consplan2 <- read.csv(paste(marxfolder, 'output', runname2, '/', runname2, '_best.csv', sep=''))
load('data/rich.RData')
load('data/presmap.RData') # loads presmap data.frame with presence/absence information

# add zone to pus
pusplan <- merge(pus, consplan, by.x='id', by.y='planning_unit')
	dim(pus)
	dim(pusplan)

pusplan2 <- merge(pus2, consplan2, by.x='id', by.y='planning_unit')
	dim(pus2)
	dim(pusplan2)

# plot map of selected grids, on top of richness (#1)
	colfun <- colorRamp(rev(brewer.pal(11, 'Spectral')))
	cexs = 0.5 # to adjust
	periods <- sort(unique(rich$period))
	# quartz(width=10, height=3)
	pdf(width=10, height=3, file=paste('figures/MarZone_NEUSSpring_on_richness_', runname, '.pdf', sep=''))
	par(mfrow=c(1,length(periods)), mai=c(0.5,0.5,0.3, 0.1), las=1, mgp=c(2,1,0))
	j <- rich$region == 'NEFSC_NEUSSpring'
	for(i in 1:length(periods)){
		j2 <- rich$period == periods[i] & j
		plot(rich$lon[j2], rich$lat[j2], col=rgb(colfun(pnorm01(rich$rich[j2], rich$rich[j])), maxColorValue=255), pch=16, cex=cexs, xlab='Longitude', ylab='Latitude', main=paste('NEUSSpring', periods[i]), cex.main=0.9)

		i <- pusplan$zone == 2
		points(pusplan$lon[i], pusplan$lat[i], pch=1, col='black', cex=cexs)
	}
	legend('bottomright', legend=round(seq(min(rich$rich[j]), max(rich$rich[j]), length.out=10),2), col=rgb(colfun(norm01(seq(min(rich$rich[j]), max(rich$rich[j]), length.out=10))), maxColorValue=255), pch=16, cex=0.7, title='Species', bty='n')


	dev.off()

# plot map of selected grids, on top of richness (#1 and #2)
	colfun <- colorRamp(rev(brewer.pal(11, 'Spectral')))
	cexs = 0.5 # to adjust
	periods <- sort(unique(rich$period))
	# quartz(width=10, height=3)
	pdf(width=10, height=3, file=paste('figures/MarZone_NEUSSpring_on_richness_', runname, '&', runname2, '.pdf', sep=''))
	par(mfrow=c(1,length(periods)), mai=c(0.5,0.5,0.3, 0.1), las=1, mgp=c(2,1,0))
	j <- rich$region == 'NEFSC_NEUSSpring'
	for(i in 1:length(periods)){
		j2 <- rich$period == periods[i] & j
		plot(rich$lon[j2], rich$lat[j2], col=rgb(colfun(pnorm01(rich$rich[j2], rich$rich[j])), maxColorValue=255), pch=16, cex=cexs, xlab='Longitude', ylab='Latitude', main=paste('NEUSSpring', periods[i]), cex.main=0.9)

		i <- pusplan$zone == 2
		i2 <- pusplan2$zone == 2
		points(pusplan$lon[i], pusplan$lat[i], pch=1, col='black', cex=cexs)
		points(pusplan2$lon[i2], pusplan2$lat[i2], pch=16, col='black', cex=0.3*cexs)
	}
	legend('bottomright', legend=round(seq(min(rich$rich[j]), max(rich$rich[j]), length.out=10),2), col=rgb(colfun(norm01(seq(min(rich$rich[j]), max(rich$rich[j]), length.out=10))), maxColorValue=255), pch=16, cex=0.7, title='Species', bty='n')


	dev.off()


# evaluate # targets met in each time period
	pinds <- presmap$region == 'NEFSC_NEUSSpring' # trim to this region
	totals <- aggregate(list(total = presmap$pres[pinds]), by=list(sppocean=presmap$sppocean[pinds], period=presmap$period[pinds]), FUN=sum) # how many presences in each period
	totals <- totals[totals$sppocean %in% spps$sppocean,]
		length(unique(totals$sppocean))

	temp <- merge(presmap[pinds, ], pusplan[pusplan$zone==2,]) # only keep the conserved zones
	temp2 <- merge(presmap[pinds, ], pusplan2[pusplan2$zone==2,]) # only keep the conserved zones
	consabund <- aggregate(list(conserved = temp$pres), by=list(sppocean=temp$sppocean, period=temp$period), FUN=sum)
		dim(consabund)
	consabund2 <- aggregate(list(conserved = temp2$pres), by=list(sppocean=temp2$sppocean, period=temp2$period), FUN=sum)
		dim(consabund2)

	consabund2.1 <- merge(consabund, totals)
		dim(consabund2.1)
	consabund2.2 <- merge(consabund2, totals)
		dim(consabund2.2)

	consabund2.1$prop <- consabund2.1$conserved/consabund2.1$total
	consabund2.2$prop <- consabund2.2$conserved/consabund2.2$total

	goalsmet <- aggregate(list(nmet=consabund2.1$prop>=goal), by=list(period=consabund2.1$period), FUN=sum)
	goalsmet$mid <- sapply(strsplit(as.character(goalsmet$period), split='-'), FUN=function(x) mean(as.numeric(x)))
	goalsmet2 <- aggregate(list(nmet=consabund2.2$prop>=goal), by=list(period=consabund2.2$period), FUN=sum)
	goalsmet2$mid <- sapply(strsplit(as.character(goalsmet2$period), split='-'), FUN=function(x) mean(as.numeric(x)))

# plot goals met (solution #1)
	# quartz(width=4, height=4)
	cols = brewer.pal(4, 'Paired')
	pdf(width=4, height=4, file=paste('figures/MarZone_NEUSSpring_goalsmet_', runname, '.pdf', sep=''))

	plot(goalsmet$mid, goalsmet$nmet, xlab='Year', ylab='# Goals met', ylim=c(0, 125), type='o', pch=16, las=1, col=cols[2])
	
	dev.off()
	
# plot goals met (solution #1 and #2)
	# quartz(width=4, height=4)
	cols = brewer.pal(4, 'Paired')
	pdf(width=4, height=4, file=paste('figures/MarZone_NEUSSpring_goalsmet_', runname, '&', runname2, '.pdf', sep=''))

	plot(goalsmet$mid, goalsmet$nmet, xlab='Year', ylab='# Goals met', ylim=c(0, 125), type='o', pch=16, las=1, col=cols[2])
	points(goalsmet2$mid, goalsmet2$nmet, type='o', pch=16, col=cols[4])
	
	dev.off()
	
	

######################################################################
# Read in results and analyze against each climate model projection  #
# across each rcp                                                    #
######################################################################
runname <- 'conservationtest'
runname2 <- 'conservationtest2per'
goal <- 0.2
consplan <- read.csv(paste(marxfolder, 'output', runname, '/', runname, '_best.csv', sep=''))
consplan2 <- read.csv(paste(marxfolder, 'output', runname2, '/', runname2, '_best.csv', sep=''))
load(paste(marxfolder, 'input', runname, '/pus.Rdata', sep='')) # pus
load(paste(marxfolder, 'input', runname2, '/pus.Rdata', sep='')) # pus2
load(paste(marxfolder, 'input', runname, '/spps.Rdata', sep='')) # spps
load(paste(marxfolder, 'input', runname2, '/spps.Rdata', sep='')) # spps2

load('data/presmapbymod.RData') # loads presmap data.frame with presence/absence information from each model (slow to load)

# add zone to pus
pusplan <- merge(pus, consplan, by.x='id', by.y='planning_unit')
	dim(pus)
	dim(pusplan)

pusplan2 <- merge(pus2, consplan2, by.x='id', by.y='planning_unit')
	dim(pus2)
	dim(pusplan2)



# evaluate # targets met in each time period in each model
	pinds <- presmapbymod$region == 'NEFSC_NEUSSpring' # trim to this region
	totals <- aggregate(list(total = presmapbymod$pres[pinds]), by=list(sppocean=presmapbymod$sppocean[pinds], period=presmapbymod$period[pinds], model=presmapbymod$model[pinds]), FUN=sum) # how many presences in each period in each model
	totals <- totals[totals$sppocean %in% spps$sppocean,]
		length(unique(totals$sppocean))

	temp <- merge(presmapbymod[pinds, ], pusplan[pusplan$zone==2,]) # only keep the conserved zones
	temp2 <- merge(presmapbymod[pinds, ], pusplan2[pusplan2$zone==2,]) # only keep the conserved zones
	consabund <- aggregate(list(conserved = temp$pres), by=list(sppocean=temp$sppocean, period=temp$period, model=temp$model), FUN=sum)
		dim(consabund)
	consabund2 <- aggregate(list(conserved = temp2$pres), by=list(sppocean=temp2$sppocean, period=temp2$period, model=temp2$model), FUN=sum)
		dim(consabund2)

	intersect(names(consabund), names(totals)) # check before merging
	consabund2.1 <- merge(consabund, totals)
		dim(consabund2.1)
	consabund2.2 <- merge(consabund2, totals)
		dim(consabund2.2)

	consabund2.1$prop <- consabund2.1$conserved/consabund2.1$total
	consabund2.2$prop <- consabund2.2$conserved/consabund2.2$total

	goalsmet <- aggregate(list(nmet=consabund2.1$prop>=goal), by=list(period=consabund2.1$period, model=consabund2.1$model), FUN=sum, na.rm=TRUE) # remove NAs for species with 0 abundance
	goalsmet$mid <- sapply(strsplit(as.character(goalsmet$period), split='-'), FUN=function(x) mean(as.numeric(x)))
	goalsmet2 <- aggregate(list(nmet=consabund2.2$prop>=goal), by=list(period=consabund2.2$period, model=consabund2.2$model), FUN=sum, na.rm=TRUE)
	goalsmet2$mid <- sapply(strsplit(as.character(goalsmet2$period), split='-'), FUN=function(x) mean(as.numeric(x)))

	write.csv(consabund2.1, file=paste('output/consabund_', runname, '.csv', sep=''))
	write.csv(consabund2.2, file=paste('output/consabund_', runname2, '.csv', sep=''))
	write.csv(goalsmet, file=paste('output/goalsmet_', runname, '.csv', sep=''))
	write.csv(goalsmet2, file=paste('output/goalsmet_', runname2, '.csv', sep=''))

# compare goals met
	t.test(goalsmet$nmet[goalsmet$period=='2081-2100'], goalsmet2$nmet[goalsmet2$period=='2081-2100'])

# plot goals met (solution #1)
	# quartz(width=4, height=4)
	cols = brewer.pal(4, 'Paired')
	mods <- sort(unique(goalsmet$model))
	pdf(width=4, height=4, file=paste('figures/MarZone_NEUSSpring_goalsmetbymod_', runname, '.pdf', sep=''))

	inds <- goalsmet$model == 1
	plot(goalsmet$mid[inds], goalsmet$nmet[inds], xlab='Year', ylab='# Goals met', ylim=c(0, 125), type='o', pch=16, las=1, col=cols[1])
	for(i in 2:length(mods)){
		inds <- goalsmet$model == i
		points(goalsmet$mid[inds], goalsmet$nmet[inds], type='o', pch=16, col=cols[1])
	
	}
	ensmean <- aggregate(list(nmet=goalsmet$nmet), by=list(mid=goalsmet$mid), FUN=mean)
	lines(ensmean$mid, ensmean$nmet, col=cols[2], lwd=2)
	
	dev.off()
	
# plot goals met (solution #1 and #2)
	# quartz(width=4, height=4)
	cols = brewer.pal(4, 'Paired')
	mods <- sort(unique(goalsmet$model))
	pdf(width=4, height=4, file=paste('figures/MarZone_NEUSSpring_goalsmetbymod_', runname, '&', runname2, '.pdf', sep=''))

	inds <- goalsmet$model == 1
	plot(goalsmet$mid[inds], goalsmet$nmet[inds], xlab='Year', ylab='# Goals met', ylim=c(0, 125), type='o', pch=16, las=1, col=cols[1])
	for(i in 2:length(mods)){
		inds <- goalsmet$model == i
		points(goalsmet$mid[inds], goalsmet$nmet[inds], type='o', pch=16, col=cols[1])	
	}
	ensmean <- aggregate(list(nmet=goalsmet$nmet), by=list(mid=goalsmet$mid), FUN=mean)
	lines(ensmean$mid, ensmean$nmet, col=cols[2], lwd=2)

	inds <- goalsmet2$model == 1
	points(goalsmet2$mid[inds], goalsmet2$nmet[inds], type='o', pch=16, col=cols[3])
	for(i in 2:length(mods)){
		inds <- goalsmet2$model == i
		points(goalsmet2$mid[inds], goalsmet2$nmet[inds], type='o', pch=16, col=cols[3])
	
	}
	ensmean2 <- aggregate(list(nmet=goalsmet2$nmet), by=list(mid=goalsmet2$mid), FUN=mean)
	lines(ensmean2$mid, ensmean2$nmet, col=cols[4], lwd=2)

	
	dev.off()
	
