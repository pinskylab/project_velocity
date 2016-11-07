# evaluate changes in richness

require(Hmisc)

## Set working directories
if(Sys.info()["nodename"] == "pinsky-macbookair"){
	setwd('~/Documents/Rutgers/Range projections/proj_ranges/')
	projfolder <- '../CEmodels_proj' # holds model projections (outside Git)
	modfolder <- '../CEModels' # holds the models (outside Git)
	}
if(Sys.info()["nodename"] == "amphiprion.deenr.rutgers.edu"){
	setwd('~/Documents/range_projections/')
	projfolder <- 'CEmodels_proj'
	modfolder <- 'CEmodels'
	}
# could add code for Lauren's working directory here

####################
## helper functions
####################
lu <- function(x) return(length(unique(x)))

# weighted mean function to use with summarize()
wmean <- function(x){ # values in col 1, weights in col 2
	inds <- !is.na(x[,1]) & !is.na(x[,2])
	return(weighted.mean(x[inds,1], x[inds,2])) 
}

# expects x to have columns 'periodmids' and 'rich' (in that order)
calcrichtrendmids <- function(x){
	mod <- lm(x[,2] ~ x[,1])
	return(mod$coefficients[2])
}

#######################################################################
## Are the same grid cells consistently low or high richness through time?
########################################################################

load('data/rich.RData') # loads rich data.frame: richness by grid cell by time period (ensemble mean)

# select lowest and highest richness cell in each region, start and end
regs <- as.character(sort(unique(rich$region)))
lowcellsstart <- vector('list', length(regs)) # one for each region
lowcellsend <- vector('list', length(regs)) # one for each region
highcellsstart <- vector('list', length(regs)) # one for each region
highcellsend <- vector('list', length(regs)) # one for each region
names(lowcells)<-regs
for(i in 1:length(regs)){
	sorted <- rich[rich$region==regs[i] & rich$period =='2006-2020',]
	sorted <- sorted[order(sorted$rich, decreasing=FALSE),]
	lowcellsstart[[i]] <- paste(sorted$lat, sorted$lon, sep=',')[1:round(0.1*(nrow(sorted)))]
	highcellsstart[[i]] <- paste(sorted$lat, sorted$lon, sep=',')[round(0.9*(nrow(sorted))):nrow(sorted)]

	sorted <- rich[rich$region==regs[i] & rich$period =='2081-2100',]
	sorted <- sorted[order(sorted$rich, decreasing=FALSE),]
	lowcellsend[[i]] <- paste(sorted$lat, sorted$lon, sep=',')[1:round(0.1*(nrow(sorted)))]
	highcellsend[[i]] <- paste(sorted$lat, sorted$lon, sep=',')[round(0.9*(nrow(sorted))):nrow(sorted)]

}

# overlap
for(i in 1:length(regs)){
	print(regs[i])
	print('low')
	print(length(lowcellsstart[[i]]))
	print(length(lowcellsend[[i]]))
	print(length(intersect(lowcellsstart[[i]], lowcellsend[[i]])))
	print('high')
	print(length(highcellsstart[[i]]))
	print(length(highcellsend[[i]]))
	print(length(intersect(highcellsstart[[i]], highcellsend[[i]])))
}