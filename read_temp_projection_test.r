## Set working directory
if(Sys.info()["nodename"] == "pinsky-macbookair"){
	setwd('~/Documents/Rutgers/Range projections/project_velocity/')
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


# function to process lines of data from temperature projection file
processlines <- function(x){
	if(grepl('.....', x, fixed=TRUE)){
		return(rep(as.numeric(NA), 94))
	} else {
		x <- sub('^ +', '', x) # remove leading whitespace
		return(as.numeric(unlist(strsplit(x, split=' +'))))
	}		
}


# Malin, THIS IS CURRENTLY SET UP TO GET THE CLIMATE PROJECTION DATA WORKING_WILL MODIFY ONCE THAT IS WORKING
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


## for min and max files, chop off the column furthest to the right (this would be for 2101)
## then will have 2007-2100

## testing file reading on sbt data
i=3; k=2; j=1; l=3 # sbt_rcp85, ond, mean
filename = paste('data/', pred.folder[i], modelrun[j], '_rcp', rcp, '_regrid.nc_2006_2100_', pred.season[k],'_', pred.metric[l], '.txt', sep="")
filein <- readLines(filename)
temps <- t(sapply(filein, processlines))
row.names(temps) <- 1:nrow(temps)

	# map the available data
	modeldat <- apply(temps, MARGIN=1, FUN = function(x) all(!is.na(x))) # whether or not a model has data at each grid cell
	modeldat <- cbind(modeldat, clim.grid)
	
	plot(modeldat$longrid, modeldat$latgrid, col=c('red', 'black')[as.numeric(modeldat$modeldat)+1], cex=0.1, main=modelrun[j])
	

## test to read in a file from each model (sbt data)
i=3; k=4; l=3 # sbt_rcp85, ond, mean
temps <- array(as.numeric(NA), dim=c(13637,94,length(modelrun)), dimnames=list(grid=1:13637, year=1:94, model=modelrun))
for(j in 1:length(modelrun)){
	print(j)
	filename = paste('data/', pred.folder[i], modelrun[j], '_rcp', rcp, '_regrid.nc_2006_2100_', pred.season[k],'_', pred.metric[l], '.txt', sep="")
	filein <- readLines(filename)
	temps[,,j] <- t(sapply(filein, processlines))
}

	# save out
	save(temps, file='data/sbt_rcp85_mean_allmods.rdata')
	
	# plot available data for all
	png(file='figures/temp_projection_data_maps.png', width=14, height=14, units='in', res=300)
	par(mfrow=c(4,4))
	for(j in 1:dim(temps)[3]){	
		modeldat <- apply(temps[,,j], MARGIN=1, FUN = function(x) all(!is.na(x))) # whether or not a model has data at each grid cell
		modeldat <- cbind(modeldat, clim.grid)
	
		plot(modeldat$longrid, modeldat$latgrid, col=c('red', 'black')[as.numeric(modeldat$modeldat)+1], cex=0.01, main=modelrun[j])
	}
	dev.off()
	
# what range of temperatures in each model?
apply(temps, MARGIN=3, FUN=function(x) summary(as.numeric(x)))

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
