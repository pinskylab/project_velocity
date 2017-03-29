# Read in temperature fields and models, then make range projections

## Set working directory
if(Sys.info()["nodename"] == "pinsky-macbookair"){
	setwd('~/Documents/Rutgers/Range projections/proj_ranges/')
	projfolder = '../CEmodels_proj'
	modfolder = '../CEModels'
	climgridfolder <- '../data/'
	numcorestouse <- 2
	}
if(Sys.info()["nodename"] == "amphiprion.deenr.rutgers.edu"){
	setwd('~/Documents/range_projections/')
	projfolder = 'CEmodels_proj'
	modfolder = 'CEmodels'
	climgridfolder <- 'data/'
	numcorestouse <- 12
	# .libPaths(new='~/R/x86_64-redhat-linux-gnu-library/3.1/') # so that it can find my old packages_Muted for Jim's use
}
if(Sys.info()["user"] == "jamesmorley"){
  setwd('/Users/jamesmorley/Documents/project_velocity')
  projfolder = 'output/CEmodels_proj/'
  modfolder <- 'output/CEmodels/'
  climgridfolder <- 'data/'
}

###################
## Load libraries
###################

require(mgcv)
require(Hmisc)
# require(parallel) # for multi-core calculations

###############################################
# Choose the model fit and other flags to use
###############################################
rcp <- 85
#rcp <- 26

#runtype <- 'test'
#runtype <- 'testseason'
#runtype <- 'testK6noSeas'
runtype <- 'fitallreg_2017'
stayinregion <- FALSE


#############################
# Choose species to project #
#############################
load(paste('output/modeldiag_', runtype, '.Rdata', sep='')) # model diagnostics

#projspp <- modeldiag$sppocean[modeldiag$auc.tt >= 0.75 & !is.na(modeldiag$auc.tt)] # from Elith et al., but there may be better criteria to use

#With additional criteria: ***
projspp <- modeldiag$sppocean[modeldiag$auc.tt >= 0.75 & !is.na(modeldiag$auc.tt) & ((modeldiag$dev.pres - modeldiag$dev.pres.null > 0.05) | (modeldiag$dev.biomass - modeldiag$dev.biomass.null > 0.05))] # from Elith et al., and deviance explained by temp must be > 5% for at least one model

print(paste(sum(!is.na(modeldiag$sppocean)), 'models fit'))
print(paste(length(projspp), 'models to project')) # number of species to project to

	# look at species not selected
#	hist(modeldiag$auc.tt)
#	hist(modeldiag$auc)
#	hist(modeldiag$dev.pres - modeldiag$dev.pres.null); abline(v=0.05, col='red')
#	hist(modeldiag$dev.biomass - modeldiag$dev.biomass.null); abline(v=0.05, col='red')
#	notselected <- modeldiag[!(modeldiag$sppocean %in% projspp), c('sppocean', 'auc.tt', 'dev.pres', 'dev.pres.null', 'dev.biomass', 'dev.biomass.null')]
#		sum(notselected$auc.tt < 0.75, na.rm=TRUE)
#		with(notselected, sum((dev.pres - dev.pres.null <= 0.05) & (dev.biomass - dev.biomass.null <= 0.05), na.rm=TRUE))
#		with(notselected, sum(auc.tt < 0.75 & (dev.pres - dev.pres.null <= 0.05) & (dev.biomass - dev.biomass.null <= 0.05), na.rm=TRUE))
#
#			
#	# those that would be selected using auc instead of auc.tt
#	modeldiag$sppocean[modeldiag$auc.tt < 0.75 & modeldiag$auc >= 0.75 & !is.na(modeldiag$auc.tt) & !is.na(modeldiag$auc) & ((modeldiag$dev.pres - modeldiag$dev.pres.null > 0.05) | (modeldiag$dev.biomass - modeldiag$dev.biomass.null > 0.05))]


# find the files with these species for our chosen model fit
files <- list.files(modfolder)
files <- files[grepl(paste('_', runtype, '_', sep=''), files) & grepl(paste(gsub('/|\\(|\\)', '', projspp), collapse='|'), gsub('/|\\(|\\)', '', files))] # have to strip out parentheses and slashes from file and taxon names so that grep doesn't interpret them
length(files) # should match length of projspp

## Remove spp from planned projections IF the projection file already exists (OPTIONAL). 
#If this step is skipped, the existing files will be overwritten.
if(stayinregion) donefiles <- list.files(projfolder, pattern=paste(runtype, '_rcp', rcp, sep='')) # models made earlier
if(!stayinregion) donefiles <- list.files(projfolder, pattern=paste(runtype, '_xreg_rcp', rcp, sep='')) # models made earlier

# trim out prefix and suffix
donespp <- gsub(paste('summproj_', runtype, '_', sep=''), '', gsub('.Rdata', '', donefiles)) 
if(!stayinregion) donespp <- gsub('xreg_', '', donespp)
donespp <- gsub(paste('rcp', rcp, '_', sep=''), '', donespp)

# remove spp that we've already projected
if(length(donespp)>0){
	files <- files[!grepl(paste(gsub('/|\\(|\\)', '', donespp), collapse='|'), gsub('/|\\(|\\)', '', files))] # remove any models that we made earlier
	projspp <- projspp[!grepl(paste(gsub('/|\\(|\\)', '', donespp), collapse='|'), gsub('/|\\(|\\)', '', projspp))]
}
length(files)
length(projspp)

 
#################################
# Prep environmental data
#################################
 
load('data/projectionGrid_Feb24_2017.RData')# load projection grid to get lat/lon values
clim.grid <- proj.grid # rename as a different 'proj.grid' imported below with bathymetry
rm(proj.grid)
nrow(unique(clim.grid)) # 13,637 unique lat/lon cells for projections_same as nrow(clim.grid)

# Below any of the climate projection files can be uploaded by adjusting rcp, i, j, k, and l
rcp <- 85 
#rcp <- 26
pred.folder <- c('sst_rcp85/tos_Omon_','sst_rcp26/tos_Omon_','sbt_rcp85/temp_btm_1950_2100_','sbt_rcp26/temp_btm_1950_2100_')
modelrun <- c('bcc-csm1-1-m','bcc-csm1-1','CanESM2','CCSM4','CESM1-CAM5','CNRM-CM5','GFDL-CM3','GFDL-ESM2M','GFDL-ESM2G','GISS-E2-R','GISS-E2-H','IPSL-CM5A-LR','IPSL-CM5A-MR','MIROC-ESM','MPI-ESM-LR','NorESM1-ME')
pred.season <- c('jfm','amj','jas','ond')
pred.metric <- c('max', 'min', 'mean')

processlines <- function(x){
  if(grepl('.....', x, fixed=TRUE)){
    return(rep(as.numeric(NA), 94))
  } else {
    x <- sub('^ +', '', x) # remove leading whitespace
    return(as.numeric(unlist(strsplit(x, split=' +'))))
  }
}

# For reading in all models for a particular type of data, I’m using:
i=3; k=3; l=3 # sbt_rcp85, ond, mean
temps <- array(as.numeric(NA), dim=c(13637,94,length(modelrun)), dimnames=list(grid=1:13637, year=1:94, model=1:length(modelrun)))
for(j in 1:length(modelrun)){
  print(j)
  # THE SURFACE VS. BOTTOM TEMP FILES ARE NAMED SLIGHTLY DIFFERENT, SO THIS CONDITIONAL STATEMENT NEEDED FOR MAKING filename
  if(pred.folder[i] == 'sst_rcp85/tos_Omon_' | pred.folder[i] == 'sst_rcp26/tos_Omon_'){
    filename = paste('data/', pred.folder[i], modelrun[j], '_rcp', rcp, '_r1i1p1_1950_2100.nc_regrid.nc_2006_2100_', pred.season[k],'_', pred.metric[l], '.txt', sep="")
  } else{
    filename = paste('data/', pred.folder[i], modelrun[j], '_rcp', rcp, '_regrid.nc_2006_2100_', pred.season[k],'_', pred.metric[l], '.txt', sep="")
  }
  # filename = paste('data/', pred.folder[i], modelrun[j], '_rcp', rcp,'_regrid.nc_2006_2100_', pred.season[k],'_', pred.metric[l], '.txt', sep="")
  filein <- readLines(filename)
  temps[,,j] <- t(sapply(filein, processlines))
}
  
# Graph up some temperature projection values
library(lattice); library(RColorBrewer)
# Take a look at some temperature values_I recycle this (with the above) code to look at different combinations of model, year, etc.
summary(temps[,,12])
abc <- data.frame(cbind(temps[,,12], clim.grid))

nrow(abc[!is.na(abc$X74),]) # how many total rows of data

plot(latgrid~longrid, data=abc[!is.na(abc$X74),])
points(latgrid~longrid, col='red', data=abc[abc$X74 > 160, ])
points(latgrid~longrid, col='green', data=abc[abc$X74 < -19, ])
# another plot option
forPlot <- abc$X74# choose a given year to work with the plots_this is only necessary for establishing the cutpoints
cutpts <- c(min(forPlot, na.rm=T)-1, -10, -5, 40, 80, max(forPlot, na.rm=T)+1) # add some cutpoints to show outliers_gray is relatively normal, with two levels of both negative and positive outliers
levelplot(X74 ~ longrid * latgrid, data = abc[!is.na(abc$X74),], at = cutpts, cuts = 6, pretty = T, col.regions = (rev(brewer.pal(9, "RdBu")))) 

# ===========================================
# Malin's code for looking at empty cells
# ===========================================

# how much data do we have? (need array from previous code block)
modeldat <- apply(temps, MARGIN=c(1,3), FUN = function(x) all(!is.na(x))) # whether or not a model has data at each grid cell
colnames(modeldat) <- modelrun
t(t(colSums(modeldat))) # how many grid cells covered by each model

# how many models have data at each grid cell?
nmodeldat <- rowSums(modeldat)

hist(nmodeldat, breaks=seq(-0.5, 16.5, by=1))
sum(nmodeldat==16)
sum(nmodeldat>=13) # 9479
sum(nmodeldat>=10) # 11736

# how many models have data at each grid cell? W/OUT GFDL or Can
keep <- !grepl('GFDL', modelrun)
keep <- !grepl('GFDL|CanESM', modelrun)
sum(keep)
nmodeldat <- rowSums(modeldat[, keep])

hist(nmodeldat, breaks=seq(-0.5, 16.5, by=1))
sum(nmodeldat>=13)
sum(nmodeldat>=12) 
sum(nmodeldat>=10)

dim(temps)


# =================================================================================
# =================================================================================

load('data/ProjectionBathGrid_Feb27_2017.RData')# load projection grid 

if(!file.exists(paste(climgridfolder, 'climGrid_rcp', rcp, '.proj2_wrugos.RData', sep=''))){
	print(paste('climGrid with rugosity does not exist for RCP', rcp, '. Making it.', sep=''))
	
	load(paste(climgridfolder, 'climGrid_rcp', rcp, '.proj2.RData', sep='')) # projected temperature for each year ("clim")

	# drop unneeded columns
	clim <- clim[,!grepl('depthgrid', names(clim))] #  refer to GCM depth grids
	clim <- clim[,!grepl('bottemp.clim|surftemp.clim|delta|latgrid|longrid', names(clim))] #  the temp climatologies, deltas, and GCM lat/lon grids (1 degree)

	# add regionfact
	clim$region<- as.factor(clim$region)
	names(clim)[names(clim)=='region'] <- 'regionfact'

	# add logrugosity
	rugos <- read.csv('data/projectiongrid_latlons.1.16th_withRugosity_2015-05-06.csv')
		names(rugos)[names(rugos) == 'lon'] <- 'lon16th'
		names(rugos)[names(rugos) == 'lat'] <- 'lat16th'
		names(rugos)[names(rugos) == 'depth'] <- 'depth16th'

		gridsize=0.25 # size of grid of the climate data, in degrees
		rugos$lat <- floor(rugos$lat16th/gridsize)*gridsize + gridsize/2 # round to nearest grid center
		rugos$lon <- floor(rugos$lon16th/gridsize)*gridsize + gridsize/2

	clim <- merge(clim, rugos) # slow
	dim(clim) # 9623120 rows
	
	save(clim, file=paste(climgridfolder, 'climGrid_rcp', rcp, '.proj2_wrugos.RData', sep=''))
	
} else {
	print('climGrid with rugosity exists. Loading it')
	load(paste(climgridfolder, 'climGrid_rcp', rcp, '.proj2_wrugos.RData', sep=''))
}

############################################
## Project GAMS onto annual climate data  ##
############################################
	
options(warn=1) # print warnings as they occur

# thisprojspp <- projspp[1]

# VERY slow
doprojection <- function(thisprojspp, files, clim, projfolder, modfolder, runtype, stayinregion=TRUE){ 
	# the stayinregion flag determines whether or not we project species outside of the regions in which they were observed historically (and to which models were fit)

	# clear variables
	mods <- avemeanbiomass <- NULL
	
	# load model fits (mod and avemeanbiomass)
	fileindex <- which(grepl(gsub('/|\\(|\\)', '', thisprojspp), gsub('/|\\(|\\)', '', files)))
	print(paste(fileindex, thisprojspp, Sys.time()))

	load(paste(modfolder, '/', files[fileindex], sep='')) # loads mods and avemeanbiomass

	fitregions <- gsub('regionfact', '', grep('regionfact', names(coef(mods$mygam1)), value=TRUE)) # regions included in the model fit (if more than one region fit)
	if(length(fitregions)==0) fitregions <- names(avemeanbiomass)[which.max(avemeanbiomass)] # else, if only one region fit (and no regionfact term included), then use avemeanbiomass to pull out the name of that one region

	if(stayinregion){ # if we don't allow species to move into new regions, then we can use the average observed biomass for each region
		# add mean biomass by region
		clim$biomassmean <- 0
		clim$biomassmean[clim$regionfact %in% names(avemeanbiomass)] <- avemeanbiomass[as.character(clim$regionfact[clim$regionfact %in% names(avemeanbiomass)])] # use region to pull the correct mean biomass values
	}
	if(!stayinregion){ #else, add a standard mean biomass across all regions (use region with highest mean biomass that was in the model fit)
		clim$biomassmean <- max(avemeanbiomass[fitregions])
	}

	# smearing estimator for re-transformation bias (see Duan 1983, http://www.herc.research.va.gov/resources/faq_e02.asp)
	smear <- mean(exp(mods[['mygam2']]$residuals))

	# rows of clim to project to
	regstoproj <- names(avemeanbiomass)
	if(any(regstoproj == 'NEFSC_NEUS')){
		regstoproj <- c(regstoproj, 'NEFSC_NEUSSpring', 'NEFSC_NEUSFall') # add spring and fall surveys (how they are marked in the climatology) to the list. The models are fit to all NEUS data and don't distinguish seasons.
	}
	if(stayinregion){ # only project to regions in which this species had a climate envelope fit
		inds <- clim$regionfact %in% regstoproj
	}
	if(!stayinregion){ # project to all regions in the same ocean (Atlantic or Pacific)
		if(any(regstoproj %in% c("DFO_NewfoundlandFall", "DFO_NewfoundlandSpring", "DFO_ScotianShelf", "DFO_SoGulf", "NEFSC_NEUSFall", "NEFSC_NEUSSpring", "SEFSC_GOMex"))){ # Atlantic
			regstoproj <- c("DFO_NewfoundlandFall", "DFO_NewfoundlandSpring", "DFO_ScotianShelf", "DFO_SoGulf", "NEFSC_NEUSFall", "NEFSC_NEUSSpring", "SEFSC_GOMex")
		}
		if(any(regstoproj %in% c("AFSC_Aleutians", "AFSC_EBS", "AFSC_GOA", "AFSC_WCTri", "NWFSC_WCAnn"))){ # Pacific
			regstoproj <- c("AFSC_Aleutians", "AFSC_EBS", "AFSC_GOA", "AFSC_WCTri", "NWFSC_WCAnn")
		}
		inds <- clim$regionfact %in% regstoproj
	}
	

	# Dataframe for this species' projections
	thisproj <- clim[inds,c('regionfact', 'lat', 'lon', 'depth16th', 'year')] # dataframe to hold projections for this taxon
	for(i in 1:13) thisproj[[paste('wtcpue.proj_', i, sep='')]] <- NA 	# Add projected biomass density columns for each climate model

	# Add region factor to use in model during projection
	if(stayinregion){
		clim$regiontoproj <- as.character(clim$regionfact)
		# Adjust region names for NEUS
		if(any(regstoproj %in% 'NEFSC_NEUS')){
			clim$regiontoproj[clim$regiontoproj %in% c('NEFSC_NEUSSpring', 'NEFSC_NEUSFall')] <- 'NEFSC_NEUS'
		}
	}
	if(!stayinregion){ 
		clim$regiontoproj <- fitregions[which.max(avemeanbiomass[fitregions])] # pick the region with the highest biomass. this ensures that it had at least some presences, and so would have been in the model fit
	}

	# Calculate predictions for 2020-2100 for each model
	for(i in 1:13){ # this alone takes a long time
		print(paste(fileindex, 'model', i))
		nd <- data.frame(regionfact = clim$regiontoproj[inds], surftemp = clim[[paste('surftemp.proj_', i, sep='')]][inds], bottemp = clim[[paste('bottemp.proj_', i, sep='')]][inds], logrugosity = log(clim$rugosity[inds]+0.01), biomassmean = clim$biomassmean[inds], row.names=1:sum(inds))
		preds1 <- predict.gam(mods$mygam1, newdata = nd, type='response')
		preds2 <- exp(predict(mods$mygam2, newdata = nd, type='response'))
		preds <- preds1*preds2*smear
		preds[preds<0] <- 0
		
		pnm = paste('wtcpue.proj_', i, sep='')
		thisproj[[pnm]] = preds # can do this because clim[inds,] and thisproj are in the same order

		# set biomass projections on land to zero
		thisproj[[pnm]][thisproj$depth16th > 0] <- 0
	}
	
	# summarize by 1/4 degree grid
	summproj <- aggregate(thisproj[,grepl('wtcpue.proj', names(thisproj))], by=list(region=thisproj$regionfact, lat=thisproj$lat, lon=thisproj$lon, year=thisproj$year), FUN=mean) # also slow
	
#	print(summary(summproj))
#	print(dim(summproj))	

	thisprojspp <- gsub('/', '', thisprojspp) # would mess up saving the file if the species name had /
	outfile <- paste(projfolder, '/summproj_', runtype, '_rcp', rcp, '_', thisprojspp, '.Rdata', sep='')
	if(!stayinregion) outfile <- paste(projfolder, '/summproj_', runtype, '_xreg_rcp', rcp, '_', thisprojspp, '.Rdata', sep='') # append xreg to projections if species could be projected out of their observed regions
	save(summproj, file=outfile) # write out the projections (15MB file)
}

result <- mclapply(X= projspp, FUN=doprojection, files=files, clim=clim, projfolder=projfolder, modfolder=modfolder, runtype=runtype, stayinregion=stayinregion, mc.cores=numcorestouse) # spawn out to multiple cores. errors will be stored in results

print(result)
