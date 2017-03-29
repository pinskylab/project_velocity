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
i=3; k=2; j=2; l=3 # sbt_rcp85, ond, mean
filename = paste('data/', pred.folder[i], modelrun[j], '_rcp', rcp, '_regrid.nc_2006_2100_', pred.season[k],'_', pred.metric[l], '.txt', sep="")
filein <- readLines(filename)
temps <- t(sapply(filein, processlines))



## test to read in a file from each model (sbt data)
i=3; k=4; l=3 # sbt_rcp85, ond, mean
temps <- array(as.numeric(NA), dim=c(13637,94,length(modelrun)), dimnames=list(grid=1:13637, year=1:94, model=1:length(modelrun)))
for(j in 1:length(modelrun)){
	print(j)
	filename = paste('data/', pred.folder[i], modelrun[j], '_rcp', rcp, '_regrid.nc_2006_2100_', pred.season[k],'_', pred.metric[l], '.txt', sep="")
	filein <- readLines(filename)
	temps[,,j] <- t(sapply(filein, processlines))
}

# how much data do we have? (need array from previous code block)
modeldat <- apply(temps, MARGIN=c(1,3), FUN = function(x) all(!is.na(x))) # whether or not a model has data at each grid cell
colnames(modeldat) <- modelrun
nmodeldat <- rowSums(modeldat) # how many models have data at each grid cell

hist(nmodeldat, breaks=seq(-0.5, 16.5, by=1))
	sum(nmodeldat==16)
	sum(nmodeldat>=13) # 9479
	sum(nmodeldat>=10) # 11736

dim(temps)
