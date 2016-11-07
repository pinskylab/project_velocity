# Set up a Marxan with Zones run for CMSP

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

# choose the rcp (for runs using just one)
rcp <- 85

# CMSP goals
consgoal <- 0.1 # proportion of presences to capture in conservation
energygoal <- 0.2 # proportion of NPV
fishgoal <- 0.5 # proportion of biomass
amountcolnm <- 'pres' # either pres or wtcpue.proj
cost <- 0.01
zonecosts <- c(0.75, 1, 1, 1) # cost multipliers for available, conservation, fishery, and energy zones
fpf <- 10

# choose region and name this run
#myreg <- 'NEFSC_NEUSSpring'; runname1 <- 'cmsphistonly_neus'; runname2 <- 'cmsp2per_neus'
#myreg <- 'AFSC_EBS'; runname1 <- 'cmsphistonly_ebs'; runname2 <- 'cmsp2per_ebs'
#myreg <- 'DFO_NewfoundlandFall'; runname1 <- 'cmsphistonly_newf'; runname2 <- 'cmsp2per_newf'
#myreg <- 'SEFSC_GOMex'; runname1 <- 'cmsphistonly_gmex'; runname2 <- 'cmsp2per_gmex'
#myreg <- 'DFO_ScotianShelf'; runname1 <- 'cmsphistonly_scot'; runname2 <- 'cmsp2per_scot'
#myreg <- 'DFO_SoGulf'; runname1 <- 'cmsphistonly_sgulf'; runname2 <- 'cmsp2per_sgulf'
#myreg <- 'AFSC_GOA'; runname1 <- 'cmsphistonly_goa'; runname2 <- 'cmsp2per_goa'
myreg <- 'AFSC_Aleutians'; runname1 <- 'cmsphistonly_ai'; runname2 <- 'cmsp2per_ai'

# which time periods to use in the multi-period planning
planningperiods <- c('2006-2020', '2081-2100')

# folders
inputfolder1 <- paste(marxfolder, runname1, '_input', sep='')
inputfolder2 <- paste(marxfolder, runname2, '_input', sep='')
outputfolder1 <- paste(marxfolder, runname1, '_output', sep='')
outputfolder2 <- paste(marxfolder, runname2, '_output', sep='')



#####################
## Load data
#####################

load(paste('data/presmap_', runtype, projtype, '_rcp', rcp, '.RData', sep='')) # loads presmap data.frame with presence/absence information
	sort(unique(presmap$region))
windnpv <- read.csv('cmsp_data/wind_npv.csv', row.names=1)
wavenpv <- read.csv('cmsp_data/wave_npv.csv', row.names=1)
fisheryspps <- read.csv('cmsp_data/fishery_spps.csv', row.names=1) # which spp to include in fishery goal in each region

# Trim to just one region
presmap <- presmap[presmap$region==myreg,]
fisheryspps <- fisheryspps[fisheryspps$region==myreg,]
	dim(presmap)
	dim(fisheryspps)

############################################
## Set up a Marxan run just on 2006-2020
## use just one rcp
############################################

# Create directory for input and output if missing
if(!dir.exists(inputfolder1)){
	dir.create(inputfolder1)
}
if(!dir.exists(outputfolder1)){
	dir.create(outputfolder1)
}

# pu.dat
# planning features are each 1/4 deg square
	pus <- presmap[,c('lat', 'lon')]
	pus <- pus[!duplicated(pus),]
		dim(pus) # 552 (neus), (ebs), 1342 (newf), 229 (wc), 245 (gmex), 396 (scot), 186 (sgulf), 795 (goa), 231 (ai)
	pus <- pus[order(pus$lat, pus$lon),]
	pus$id <- 1:nrow(pus)
	pus$dummycost <- rep(cost, nrow(pus)) # set the same cost in each planning unit. can add separate costs for each zone.

	pu.dat<-pus[,c('id', 'dummycost')] # the version to write out for Marxan

	save(pus, file=paste(inputfolder1, '/pus.Rdata', sep=''))	
	write.csv(pu.dat, file=paste(inputfolder1, '/pu.dat', sep=''), row.names=FALSE, quote=FALSE)

# spec.dat (takes the place of feat.dat?)
# id and name for each species
# not documented in 1.0.1 manual. copying format from unix example
	spps <- presmap[presmap$pres,c('sppocean')]
	spps <- spps[!duplicated(spps)]
		length(spps) # 102 (NEUS), 214 (EBS), 111 (Newf), 118 (WC), 146 (gmex), 103 (scot), 88 (sgulf), 147 (goa), 139 (ai)

	sppstokeep <- presmap[presmap$period=='2006-2020' & presmap$pres,c('lat', 'lon', 'sppocean', 'pres')]
		dim(sppstokeep)
	sppstokeep <- merge(sppstokeep, pus[,c('lat', 'lon', 'id')]) # add pu id (and trim to focal pus)
		names(sppstokeep)[names(sppstokeep)=='id'] <- 'pu'
		dim(sppstokeep)

		ngrid <- aggregate(list(ngrid=sppstokeep$pu), by=list(sppocean=sppstokeep$sppocean), FUN=function(x) length(unique(x)))
		sppstokeep <- merge(sppstokeep, ngrid)
		summary(sppstokeep$ngrid) #

		nspps <- aggregate(list(nspp=sppstokeep$sppocean), by=list(pu=sppstokeep$pu), FUN=function(x) length(unique(x)))
		sppstokeep <- merge(sppstokeep, nspps)
		summary(sppstokeep$nspp) #

		sppstokeep <- sppstokeep[sppstokeep$ngrid> (nrow(pus)*0.05),] # trim to species found in at least 10% of grids

		length(unique(sppstokeep$sppocean)) # 67 (NEUS), 126 (EBS), 70 (Newf), 115 (WC), 140 (gmex), 48 (scot), 42 (sgulf), 115 (goa), 46 (ai)

	spps <- spps[spps %in% sppstokeep$sppocean]
		length(spps) # 67, 126, 70

	spps <- data.frame(id=1:length(spps), name=gsub(' |_', '', spps), sppocean=spps) #  fill spaces in species names.

	# set feature penalty factor
	spps$fpf <- fpf

	# add wind and wave energy feature
	spps <- rbind(spps, data.frame(id=max(spps$id)+1, name=c('energy'), sppocean=c(NA), fpf=fpf))


	spps.dat <- spps[,c('id', 'fpf', 'name')]

	save(spps, file=paste(inputfolder1, '/spps.Rdata', sep=''))		
	write.csv(spps.dat, file=paste(inputfolder1, '/spec.dat', sep=''), row.names=FALSE, quote=FALSE)

# puvsp.dat (takes the place of puvfeat.dat)
# which features are in each planning unit
# not documented in 1.0.1 manual. copying format from unix example
	# Format species data
	puvsp <- presmap[presmap$period=='2006-2020',c('lat', 'lon', 'sppocean', 'wtcpue.proj', 'pres')]
		dim(puvsp)
	puvsp <- merge(puvsp, pus[,c('lat', 'lon', 'id')]) # add pu id (and trim to focal pus)
		dim(puvsp)
		names(puvsp)[names(puvsp)=='id'] <- 'pu'


	puvsp$name <- gsub(' |_', '', puvsp$sppocean) # trim out spaces on species names
	puvsp <- merge(puvsp, spps[,c('id', 'name')]) # merge in species IDs and trim to focal species
		dim(puvsp)
		names(puvsp)[names(puvsp)=='id'] <- 'species'
	puvsp$amount <- as.numeric(puvsp[[amountcolnm]]) # use projected biomass as amount
	puvsp$amount[!puvsp$pres] <- 0 # set amount to zero where our cutoff says not present

	# Format wind and wave data
	puvenergy <- merge(windnpv, pus[,c('lat', 'lon', 'id')], all.y=TRUE)
	puvenergy <- merge(puvenergy, wavenpv)
		names(puvenergy)[names(puvenergy)=='id'] <- 'pu'
		head(puvenergy)
		dim(windnpv)
		dim(wavenpv)
		dim(puvenergy)
	puvenergy$wind_npv[puvenergy$wind_npv<0 | is.na(puvenergy$wind_npv)] <- 0 # set negative or NA NPV to 0
	puvenergy$wave_npv[puvenergy$wave_npv<0 | is.na(puvenergy$wave_npv)] <- 0
	puvenergy$amount <- puvenergy$wind_npv + puvenergy$wave_npv
	
	puvenergy$species <- spps$id[spps$name=='energy']
			
	length(unique(puvsp$pu))
	length(unique(puvenergy$pu))
	length(unique(puvsp$species))

	# Combine species, wind, and wave data for output
	puvsp.dat <- rbind(puvsp[,c('species', 'pu', 'amount')], puvenergy[,c('species', 'pu', 'amount')])
	puvsp.dat <- puvsp.dat[order(puvsp.dat$pu, puvsp.dat$species),]
	puvsp.dat <- puvsp.dat[puvsp.dat$amount>0,] # trim only to presences

		length(unique(puvsp$pu)) # 552(neus), 677 (EBS), 1342 (newf)
		length(unique(puvsp$species)) # 67 (neus), 126 (ebs), 70 (newf)
		length(unique(puvsp.dat$pu)) # 549 (neus), 677 (ebs), 1342 (newf)
		length(unique(puvsp.dat$species)) # 68 (neus), 127 (ebs), 71 (newf)
		sort(unique(table(puvsp.dat$species))) # make sure all species show up in some planning units
		sort(unique(table(puvsp.dat$pu))) # make sure all planning units have some species

		sort(unique(table(puvsp.dat$pu, puvsp.dat$species))) # should be all 0s and 1s

	write.csv(puvsp.dat, file=paste(inputfolder1, '/puvsp.dat', sep=''), row.names=FALSE, quote=FALSE)

# zones
# id and names for each zone
	zones <- data.frame(zoneid=1:4, zonename=c('available', 'conservation', 'fishery', 'energy'))

	write.csv(zones, file=paste(inputfolder1, '/zones.dat', sep=''), row.names=FALSE, quote=FALSE)

#costs
	costs <- data.frame(costid=1, costname='dummycost')

	write.csv(costs, file=paste(inputfolder1, '/costs.dat', sep=''), row.names=FALSE, quote=FALSE)

#zone cost
#to adjust the importance of each cost in each zone
	zonecost <- expand.grid(list(zoneid=zones$zoneid, costid=costs$costid))
	zonecost$multiplier <- 0
	zonecost$multiplier[zonecost$zoneid==1 & zonecost$costid==1] <- zonecosts[1] # set multiplier for available
	zonecost$multiplier[zonecost$zoneid==2 & zonecost$costid==1] <- zonecosts[2] # set multiplier for conservation
	zonecost$multiplier[zonecost$zoneid==3 & zonecost$costid==1] <- zonecosts[3] # set multiplier for fishery
	zonecost$multiplier[zonecost$zoneid==4 & zonecost$costid==1] <- zonecosts[4] # set multiplier for energy

	zonecost <- zonecost[zonecost$multiplier>0,] # trim out zeros

	write.csv(zonecost, file=paste(inputfolder1, '/zonecost.dat', sep=''), row.names=FALSE, quote=FALSE)

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

	# set conservation zone target
	consinds <- zonetarget$zoneid==zones$zoneid[zones$zonename=='conservation'] & zonetarget$speciesid %in% spps$id[spps$name != 'energy']
	zonetarget$target[consinds] <- consgoal # XX proportion
	if(amountcolnm=='pres'){
		zonetarget$targettype[consinds] <- 1 # 3: proportion of total occurrences. 1: proportion of total amount
	}
	if(amountcolnm=='wtcpue.proj'){
		zonetarget$targettype[consinds] <- 3 # 3: proportion of total occurrences. 1: proportion of total amount
	}
	

	# set fishing zone target
	fishinds <- zonetarget$zoneid==zones$zoneid[zones$zonename=='fishery'] & zonetarget$speciesid %in% spps$id[spps$sppocean %in% fisheryspps$projname]
	zonetarget$target[fishinds] <- fishgoal # XX proportion
	zonetarget$targettype[fishinds] <- 1 # 3: proportion of total occurrences. 1: proportion of total amount

	# set energy goal target
	energyinds <- zonetarget$zoneid==zones$zoneid[zones$zonename=='energy'] & zonetarget$speciesid %in% spps$id[spps$name == 'energy']
	zonetarget$target[energyinds] <- energygoal # XX proportion
	zonetarget$targettype[energyinds] <- 1 # 3: proportion of total occurrences. 1: proportion of total amount

	# format for output
	zonetarget.dat <- zonetarget[zonetarget$target>0,] # trim to only positive targets
	zonetarget.dat <- zonetarget.dat[order(zonetarget.dat$zoneid, zonetarget.dat$speciesid),] # order

	# write out	
	write.csv(zonetarget.dat, file=paste(inputfolder1, '/zonetarget.dat', sep=''), row.names=FALSE, quote=FALSE)


#input parameters
input <- data.frame(BLM=0, PROP=0.8, RANDSEED=-1, NUMREPS=100, AVAILABLEZONE=1, NUMITNS='1000000', STARTTEMP=-1, NUMTEMP='10000', COSTTHRESH=0, THRESHPEN1=0, THRESHPEN2=0, INPUTDIR=paste(runname1, '_input', sep=''), PUNAME='pu.dat', SPECNAME='spec.dat', PUVSPRNAME='puvsp.dat', ZONESNAME='zones.dat', COSTSNAME='costs.dat', ZONECOSTNAME='zonecost.dat', ZONETARGETNAME='zonetarget.dat', SCENNAME=runname1, SAVERUN=3, SAVEBEST=3, SAVESUMMARY=3, SAVESCEN=3, SAVETARGMET=3, SAVESUMSOLN=3, SAVEPENALTY=3, SAVELOG=3, OUTPUTDIR=paste(runname1, '_output', sep=''),  RUNMODE=1, MISSLEVEL=1, ITIMPTYPE=0, HEURTYPE=-1, CLUMPTYPE=0, VERBOSITY=3, SAVESOLUTIONSMATRIX=3, SAVEANNEALINGTRACE=0, ANNEALINGTRACEROWS=1000)

write.table(t(input), file=paste(marxfolder, 'input.dat', sep=''), row.names=TRUE, quote=FALSE, sep=' ', col.names=FALSE)
write.table(t(input), file=paste(inputfolder1, '/input.dat', sep=''), row.names=TRUE, quote=FALSE, sep=' ', col.names=FALSE)
	
	
# Go run MarZone!
# cd /Users/mpinsky/Documents/Rutgers/Range\ projections/MarZone_runs 
# ./MarZone_v201_Mac32 # will read in input.dat and write to the output folder



##################################################################
## Set up a Marxan run on 2006-2020 and ensemble mean 2081-2100
## This assumes that the historical-only code has been run and is loaded in memory
##################################################################
# Create directory for input and output if missing
if(!dir.exists(inputfolder2)){
	dir.create(inputfolder2)
}
if(!dir.exists(outputfolder2)){
	dir.create(outputfolder2)
}


# pu.dat
# planning features are each 1/4 deg square
	pus2 <- presmap[,c('lat', 'lon')]
	pus2 <- pus2[!duplicated(pus2),]
		dim(pus2)
	pus2 <- pus2[order(pus2$lat, pus2$lon),]
	pus2$id <- 1:nrow(pus2)
	pus2$dummycost <- rep(cost, nrow(pus2))  # set the same cost in each planning unit. can add separate costs for each zone.

	pu2.dat<-pus2[,c('id', 'dummycost')]

	save(pus2, file=paste(inputfolder2, '/pus.Rdata', sep=''))	
	write.csv(pu2.dat, file=paste(inputfolder2, '/pu.dat', sep=''), row.names=FALSE, quote=FALSE)

# spec.dat (takes the place of feat.dat?)
# document every species present
# not documented in 1.0.1 manual. copying format from unix example
	spps2 <- spps #  use the same species as in the historical-only run
		spps2$name <- as.character(spps2$name)

	sppinds <- !grepl('energy', spps$name) # don't include energy in each time period
	temp1 <- spps2[sppinds,]
	spps2$name[sppinds] <- paste(spps2$name[sppinds], gsub('-', '', planningperiods[1]), sep='')
	temp1$name <- paste(temp1$name, gsub('-', '', planningperiods[2]), sep='')
	temp1$id = temp1$id + max(spps2$id) # make sure the ids don't overlap
	spps2 <- rbind(spps2, temp1)

	# set feature penalty factor
	spps2$fpf <- fpf


	spps2.dat <- spps2[,c('id', 'fpf', 'name')]

	save(spps2, file=paste(inputfolder2, '/spps.Rdata', sep=''))		
	write.csv(spps2.dat, file=paste(inputfolder2, '/spec.dat', sep=''), row.names=FALSE, quote=FALSE)

# puvsp.dat (takes the place of puvfeat.dat?)
# table of species by planning units
# not documented in 1.0.1 manual. copying format from unix example
	puvsp2 <- presmap[presmap$period %in% c('2006-2020', '2081-2100'),c('lat', 'lon', 'period', 'sppocean', 'wtcpue.proj', 'pres')]
		dim(puvsp2)
	puvsp2 <- merge(puvsp2, pus2[,c('lat', 'lon', 'id')]) # add pu id (and trim to focal pus2)
		dim(puvsp2)
		names(puvsp2)[names(puvsp2)=='id'] <- 'pu'

	puvsp2$name <- gsub(' |_|-', '', paste(puvsp2$sppocean, puvsp2$period)) # create names, and remove spaces and dashes
		dim(puvsp2)
	puvsp2 <- merge(puvsp2, spps2[,c('id', 'name')]) # add species id and trim to focal species
		dim(puvsp2)
		names(puvsp2)[names(puvsp2)=='id'] <- 'species'
	puvsp2$amount <- as.numeric(puvsp2[[amountcolnm]]) # use projected biomass or pres as amount
	puvsp2$amount[!puvsp2$pres] <- 0 # set amount to zero where our cutoff says not present

	# Combine species, wind, and wave data for output
	puvsp2.dat <- rbind(puvsp2[,c('species', 'pu', 'amount')], puvenergy[,c('species', 'pu', 'amount')])
	puvsp2.dat <- puvsp2.dat[order(puvsp2.dat$pu, puvsp2.dat$species),]
	puvsp2.dat <- puvsp2.dat[puvsp2.dat$amount>0,] # trim only to presences

	length(unique(puvsp2$pu)) # 552 (neus), 677 (ebs), 
	length(unique(puvsp2$species)) # 134 (neus), 252 (ebs)
	length(unique(puvsp2.dat$pu)) # 552 (neus), 677 (ebs)
	length(unique(puvsp2.dat$species)) # 129 (neus), 235 (ebs)
	range(table(puvsp2.dat$species)) # make sure all species show up in some planning units
	range(table(puvsp2.dat$pu)) # make sure all planning units have some species

	write.csv(puvsp2.dat, file=paste(inputfolder2, '/puvsp.dat', sep=''), row.names=FALSE, quote=FALSE)


# zones
# id and names for each zone
	zones2 <- data.frame(zoneid=1:4, zonename=c('available', 'conservation', 'fishery', 'energy'))
	write.csv(zones2, file=paste(inputfolder2, '/zones.dat', sep=''), row.names=FALSE, quote=FALSE)

#costs
	costs2 <- data.frame(costid=1, costname='dummycost')	
	write.csv(costs2, file=paste(inputfolder2, '/costs.dat', sep=''), row.names=FALSE, quote=FALSE)


#zone cost
	zonecost2 <- expand.grid(list(zoneid=zones2$zoneid, costid=costs2$costid))
	zonecost2$multiplier <- 0
	zonecost2$multiplier[zonecost2$zoneid==1 & zonecost2$costid==1] <- zonecosts[1] # set multiplier for available
	zonecost2$multiplier[zonecost2$zoneid==2 & zonecost2$costid==1] <- zonecosts[2] # set multiplier for conservation
	zonecost2$multiplier[zonecost2$zoneid==3 & zonecost2$costid==1] <- zonecosts[3] # set multiplier for fishery
	zonecost2$multiplier[zonecost2$zoneid==4 & zonecost2$costid==1] <- zonecosts[4] # set multiplier for energy

	zonecost2 <- zonecost2[zonecost2$multiplier>0,]

	write.csv(zonecost2, file=paste(inputfolder2, '/zonecost.dat', sep=''), row.names=FALSE, quote=FALSE)

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
	zonetarget2 <- expand.grid(list(zoneid=zones2$zoneid, speciesid=spps2$id))
	zonetarget2$target <- 0

	# set conservation zone target
	consinds <- zonetarget2$zoneid==zones2$zoneid[zones2$zonename=='conservation'] & zonetarget2$speciesid %in% spps2$id[spps$name != 'energy']
	zonetarget2$target[consinds] <- consgoal # XX proportion
	if(amountcolnm=='pres'){
		zonetarget2$targettype[consinds] <- 1 # 3: proportion of total occurrences. 1: proportion of total amount
	}
	if(amountcolnm=='wtcpue.proj'){
		zonetarget2$targettype[consinds] <- 3 # 3: proportion of total occurrences. 1: proportion of total amount
	}

	# set fishing zone target
	fishinds <- zonetarget2$zoneid==zones2$zoneid[zones2$zonename=='fishery'] & zonetarget2$speciesid %in% spps2$id[spps2$sppocean %in% fisheryspps$projname]
	zonetarget2$target[fishinds] <- fishgoal # XX proportion
	zonetarget2$targettype[fishinds] <- 1 # 3: proportion of total occurrences. 1: proportion of total amount

	# set energy goal target
	energyinds <- zonetarget2$zoneid==zones2$zoneid[zones2$zonename=='energy'] & zonetarget2$speciesid %in% spps2$id[spps2$name == 'energy']
	zonetarget2$target[energyinds] <- energygoal # XX proportion
	zonetarget2$targettype[energyinds] <- 1 # 3: proportion of total occurrences. 1: proportion of total amount

	# format for output	
	zonetarget2.dat <- zonetarget2[zonetarget2$target>0,] # trim to only positive targets
	zonetarget2.dat <- zonetarget2.dat[order(zonetarget2.dat$zoneid, zonetarget2.dat$speciesid),]

	# write out	
	write.csv(zonetarget2.dat, file=paste(inputfolder2, '/zonetarget.dat', sep=''), row.names=FALSE, quote=FALSE)



#input parameters
input2 <- data.frame(BLM=0, PROP=0.5, RANDSEED=-1, NUMREPS=100, AVAILABLEZONE=1, NUMITNS='1000000', STARTTEMP=-1, NUMTEMP='10000', COSTTHRESH=0, THRESHPEN1=0, THRESHPEN2=0, INPUTDIR=paste(runname2, '_input', sep=''), PUNAME='pu.dat', SPECNAME='spec.dat', PUVSPRNAME='puvsp.dat', ZONESNAME='zones.dat', COSTSNAME='costs.dat', ZONECOSTNAME='zonecost.dat', ZONETARGETNAME='zonetarget.dat', SCENNAME=runname2, SAVERUN=3, SAVEBEST=3, SAVESUMMARY=3, SAVESCEN=3, SAVETARGMET=3, SAVESUMSOLN=3, SAVEPENALTY=3, SAVELOG=3, OUTPUTDIR=paste(runname2, '_output', sep=''), RUNMODE=1, MISSLEVEL=1, ITIMPTYPE=0, HEURTYPE=-1, CLUMPTYPE=0, VERBOSITY=3, SAVESOLUTIONSMATRIX=3, SAVEANNEALINGTRACE=0, ANNEALINGTRACEROWS=1000)


write.table(t(input2), file=paste(marxfolder, 'input.dat', sep=''), row.names=TRUE, quote=FALSE, sep=' ', col.names=FALSE)
write.table(t(input2), file=paste(inputfolder2, '/input.dat', sep=''), row.names=TRUE, quote=FALSE, sep=' ', col.names=FALSE)
	
	
# Go run MarZone!	
# cd /Users/mpinsky/Documents/Rutgers/Range\ projections/MarZone_runs 
# ./MarZone_v201_Mac32 # will read in input.dat and write to the output folder


