## Set working directories
if(Sys.info()["nodename"] == "pinsky-macbookair"){
	setwd('~/Documents/Rutgers/Range projections/proj_ranges/')
	}
if(Sys.info()["nodename"] == "amphiprion.deenr.rutgers.edu" & Sys.info()["login"=='malinp']){
	setwd('~/Documents/range_projections/')
	.libPaths(new='~/R/x86_64-redhat-linux-gnu-library/3.1/') # so that it can find my old packages (chron and ncdf4)
	cmip5dir <- ''
}
if(Sys.info()["nodename"] == "amphiprion.deenr.rutgers.edu" & Sys.info()["login"]=='jamesm'){
  setwd('~/Documents/range_projections/')
  .libPaths(new='~/R/x86_64-redhat-linux-gnu-library/3.1/') # so that it can find my old packages (chron and ncdf4)
  cmip5dir <-'/local/home/malinp/Documents/range_projections/'
}
# could add code for Lauren's working directory here

######################
## Useful functions
######################
#to install packages below Jim had to:
#options(repos='http://cran.rstudio.com/')
#install.packages("ncdf4")
#install.packages("chron")

require(ncdf4) # for reading ncdf files
require(chron) # for converting julian days

nyrs = function(x) as.numeric(as.character(years(x))) # return year from a chron object
mo = function(x) 1+as.numeric(strftime(x, format='%m')) # return month. strftime returns one too low... not sure why
jmo2mo = function(x, startmonth=1){ # converts from julian month (starting at 1 for startmonth in some year) to month 1:12
	temp = (x %% 12) + startmonth-1
	temp[temp>12] = temp[temp>12]-12
	temp[temp==0] = 12
	return(temp)
}
su = function(x) sort(unique(x)) # sort unique values from a vector
meanbyyear = function(x, yr) aggregate(list(x=x), by=list(yr=yr), FUN=mean)$x # take mean of x by levels of yr and return x
linearsmooth = function(x, yr=1:length(x)){ # smooth a time-series with a linear regresion
	mod = lm(x ~ yr)
	return(predict(mod, new.data = data.frame(yr = seq(min(yr), max(yr)))))
}
linearb = function(x, yr=1:length(x)){ # return slope of x vs. yrs. require at least 3 datapoints
	if(sum(!is.na(x))>2){
		mod = lm(as.numeric(x) ~ yr)
		return(coef(mod)[2])
	} else {
		return(NA)
	}
}
agmean = function(x, by) aggregate(list(x), by=by, FUN=mean)[,2] # aggregate and take the mean of a vector. NOTE: much faster to reshape so that I can use straight up mean
season = function(x){ # convert month of year to season
	out = rep(NA, length(x))
	out[x <= 3 & x>0] = 1
	out[x <= 6 & x>3] = 2
	out[x <= 9 & x>6] = 3
	out[x <= 12 & x>9] = 4
	return(out)
}
roundto = function(x,y){r = which.min(abs(x-y)); return(y[r])} # rounds x to nearest number in y


################################
## Get data about each region
################################
# load("data/trawl_allregionsforprojections_wSST_2015-06-02.RData") # loads dat to get base years for each region
load("trawl_allregionsforprojections_wSEUS_wSST_2016-05-05.RData") # loads dat to get base years for each region, using the file with SEUS

# Fix dat to match CMIP5 format
        dat = dat[!is.na(dat$lat) & !is.na(dat$lon),] # remove NA lat/lon
        dat$lon[dat$lon < 0] = dat$lon[dat$lon < 0] + 360 # fix lons to only positive to match climate data

# Get range of years, lats, and lons for each region. Use Jan 1960 as the base date (month 1)
        dat$dates = chron(dates. = paste('01', dat$month, dat$year, sep='/'), format='d/m/y') # don't have day information in dat
        baseyrs = aggregate(list(rng = dat$year), by=list(region = dat$region), FUN=range)
        basedates = aggregate(list(min = dat$dates), by=list(region = dat$region), FUN=min)
                basedates = merge(basedates, aggregate(list(max = dat$dates), by=list(region = dat$region), FUN=max))
                basedates$minmo = (nyrs(basedates$min)-1960)*12+mo(basedates$min)-1 # months since Jan 1960
                basedates$maxmo = (nyrs(basedates$max)-1960)*12+mo(basedates$max)-1
        basemonths = aggregate(list(months = dat$month), by=list(region = dat$region), FUN=su) # list the months
        baselats = aggregate(list(rng = dat$lat), by=list(region = dat$region), FUN=range)
                baselats$rng = cbind(floor(baselats$rng[,1])+0.5, ceiling(baselats$rng[,2])-0.5) # round to nearest degree center
        baselons = aggregate(list(rng = dat$lon), by=list(region = dat$region), FUN=range)
                baselons$rng = cbind(floor(baselons$rng[,1])+0.5, ceiling(baselons$rng[,2])-0.5) # round to nearest degree center
        basedepths = aggregate(list(rng = dat$depth), by=list(region = dat$region), FUN=range, na.rm=T)
                basedepths$rng[,1] = 0 # to make sure I get all the way to the surface
        regs = sort(unique(dat$region))

        rm(dat)

        # write out
        save(basedates, file='output/basedates.RData')
        save(basemonths, file='output/basemonths.RData')
        save(baselats, file='output/baselats.RData')
        save(baselons, file='output/baselons.RData')
        save(basedepths, file='output/basedepths.RData')

#################################
## Read netCDF files from IPCC ##
#################################
load('output/basedates.RData')
load('output/basemonths.RData')
load('output/baselats.RData')
load('output/baselons.RData')
load('output/basedepths.RData')
regs = sort(unique(basedates$region))

agency = c('CNRM-CERFACS', 'IPSL', 'IPSL', 'MOHC', 'MPI-M', 'MPI-M', 'MRI', 'NCAR', 'NCC', 'NCC', 'NOAA-GFDL', 'NOAA-GFDL', 'NOAA-GFDL')
model = c('CNRM-CM5', 'IPSL-CM5A-MR', 'IPSL-CM5B-LR', 'HadGem2-CC', 'MPI-ESM-LR', 'MPI-ESM-MR', 'MRI-CGCM3', 'CCSM4', 'NorESM1-M', 'NorESM1-ME', 'GFDL-CM3', 'GFDL-ESM2G', 'GFDL-ESM2M')

# Read in, average, and write out each model's historical period for each region (average)
for(j in 1:length(agency)){
        file = paste(cmip5dir, 'Froelicher_data/', agency[j], '/', model[j], '/historical/thetao_all_regrid_malin.nc', sep='')
        print(paste(j, agency[j], model[j]))
        n = nc_open(file) # Open the netCDF file
        lats = ncvar_get(n, 'LAT')
        lons = ncvar_get(n, 'LON')
        depths = ncvar_get(n, 'LEV')
        times = 1:552 # Jan 1960 to Dec 2005
        if(model[j]=='HadGem2-CC') times = 1:551 # the exception, ends Nov 2005

        # Variable to hold the climate data in the base time period for each region
        basetemps = vector('list', length(regs))
                names(basetemps) = regs

        # Loop through the regions
        for(i in 1:length(regs)){
                regtimes = basedates$minmo[i]:basedates$maxmo[i]
                xtimes = intersect(times, regtimes) # the intersecting times
                start = c(which.min(abs(lons - baselons$rng[i,1])), which.min(abs(lats - baselats$rng[i,1])), which.min(abs(depths - basedepths$rng[i,1])), xtimes[1]) # indices of the start (lon, lat, depth, time)
                end = c(which.min(abs(lons - baselons$rng[i,2])), which.min(abs(lats - baselats$rng[i,2])), which.min(abs(depths - basedepths$rng[i,2])), xtimes[length(xtimes)]) # indices of the end
                temp = ncvar_get(n, 'THETAO_REGRID', start = start, count=end - start+1) # a slice through time and space (lon, lat, depth, time)
                theselons = ncvar_get(n, 'LON', start=start[1], count= end[1]-start[1]+1)
                theselats = ncvar_get(n, 'LAT', start=start[2], count= end[2]-start[2]+1)
                thesedepths = ncvar_get(n, 'LEV', start=start[3], count= end[3]-start[3]+1)
                        print(paste(round(min(thesedepths)), 'm', 1960+floor((xtimes[1]-1)/12), '/', jmo2mo(xtimes[1]), 1960+floor((xtimes[length(xtimes)]-1)/12), '/', jmo2mo(xtimes[length(xtimes)]))) # make sure I got the depth and date conversions right. 
        
                # trim to the appropriate months
                motokeep <- unlist(basemonths$months[basemonths$region == regs[i]])
                temp <- temp[,,,jmo2mo(xtimes) %in% motokeep]
                xtimes <- xtimes[jmo2mo(xtimes) %in% motokeep]
        
                # Average across years within each lat/lon/depth
                dim(temp)
                temp2 = apply(temp, MARGIN=c(1,2,3), FUN=mean) 
                print(dim(temp2))

                basetemps[[i]] = temp2
                dimnames(basetemps[[i]]) = list(lon=theselons, lat=theselats, depth=thesedepths)
        }

        # Write out
        outfile = paste('data/tempshist_', agency[j], '_', model[j], '.RData', sep='')
        save(basetemps, file=outfile)
        
        # close
        nc_close(n)
}

# Read in and write out each model's 2006-2100 period RCP8.5 (by year and season)

for(j in 1:length(agency)){ # very slow because of apply() down below.
        print(paste(j, agency[j], model[j], Sys.time()))
        file = paste(cmip5dir, 'Froelicher_data/', agency[j], '/', model[j], '/rcp85/thetao_all_regrid_malin.nc', sep='')
        n = nc_open(file) # Open the netCDF file
        lats = ncvar_get(n, 'LAT')
        lons = ncvar_get(n, 'LON')
        depths = ncvar_get(n, 'LEV')
        startmonth = 1; startyear = 2006
        if(model[j]=='HadGem2-CC'){     startmonth = 12; startyear = 2005} # HadGEM2-CC model is different for RCP8.5 (see readme.txt), starts Dec 2005
        times = 1:length(ncvar_get(n, 'TIME')) # time, measured in months since the start month/year for this model

        # Variable to hold the climate data in the base time period for each region
        temps2100 = vector('list', length(regs))
                names(temps2100) = regs

        # times to extract for each model: Jan 2006 to Dec 2100 = 1140 months. measured in months since the model's start month/year
        keeptimes = seq(from=(2006-startyear)*12+1-startmonth+1, length.out=1140) # time, measured in months since model's start month/year
        xtimes = intersect(times, keeptimes) # the intersecting times
        if(keeptimes[1] != xtimes[1]) stop(paste('j=',j, 'did not start in January'))

        # Loop through the regions
        for(i in 1:length(regs)){
                start = c(which.min(abs(lons - baselons$rng[i,1])), which.min(abs(lats - baselats$rng[i,1])), which.min(abs(depths - basedepths$rng[i,1])), xtimes[1]) # indices of the start (lon, lat, depth, time)
                end = c(which.min(abs(lons - baselons$rng[i,2])), which.min(abs(lats - baselats$rng[i,2])), which.min(abs(depths - basedepths$rng[i,2])), xtimes[length(xtimes)]) # indices of the end

                temp = ncvar_get(n, 'THETAO_REGRID', start = start, count=end - start+1) # a slice of time: dims are lon, lat, depth, time
                theselons = ncvar_get(n, 'LON', start=start[1], count= end[1]-start[1]+1)
                theselats = ncvar_get(n, 'LAT', start=start[2], count= end[2]-start[2]+1)
                thesedepths = ncvar_get(n, 'LEV', start=start[3], count= end[3]-start[3]+1)

                # trim to the appropriate months
                motokeep <- unlist(basemonths$months[basemonths$region == regs[i]])
                temp <- temp[,,,jmo2mo(xtimes) %in% motokeep]
                thesextimes <- xtimes[jmo2mo(xtimes) %in% motokeep]

                # reshape so that month and year are in separate dimensions
                dim(temp)
                temp2 = array(temp, dim=c(dim(temp)[1:3], length(motokeep), ceiling(length(thesextimes)/length(motokeep)))) # lon, lat, depth, month, year
                if(length(temp2) > length(temp)) temp2[(length(temp)+1):length(temp2)] <- NA # if needed, fill in extra values in temp2 with NAs (R otherwise recycles the values in temp)
                dim(temp2)

                # Average across months within each lat/lon/depth/year
                dim(temp2)
                temp3 = apply(temp2, MARGIN=c(1,2,3,5), FUN=mean)
                print(dim(temp3))

                temps2100[[i]] = temp3
                dimnames(temps2100[[i]]) = list(lon=theselons, lat=theselats, depth=thesedepths,yr = 2006:(2005+ceiling(length(thesextimes)/length(motokeep))))
        }

        # Write out
        outfile = paste('data/tempsRCP85_2006-2100_', agency[j], '_', model[j], '.RData', sep='')
        save(temps2100, file=outfile)

        # close
        nc_close(n)
}

################################ Test with SEUS, only with one climate model #######################

  print(paste(1, agency[1], model[1], Sys.time()))
  file = paste('Froelicher_data/', agency[1], '/', model[1], '/rcp85/thetao_all_regrid_malin.nc', sep='')
  n = nc_open(file) # Open the netCDF file
  lats = ncvar_get(n, 'LAT')
  lons = ncvar_get(n, 'LON')
  depths = ncvar_get(n, 'LEV')
  startmonth = 1; startyear = 2006
  if(model[1]=='HadGem2-CC'){     startmonth = 12; startyear = 2005} # HadGEM2-CC model is different for RCP8.5 (see readme.txt), starts Dec 2005
  times = 1:length(ncvar_get(n, 'TIME')) # time, measured in months since the start month/year for this model
   
  # Variable to hold the climate data in the base time period for each region
  temps2100 = vector('list', length(regs))
  names(temps2100) = regs
  
  # times to extract for each model: Jan 2006 to Dec 2100 = 1140 months. measured in months since the model's start month/year
  keeptimes = seq(from=(2006-startyear)*12+1-startmonth+1, length.out=1140) # time, measured in months since model's start month/year
  xtimes = intersect(times, keeptimes) # the intersecting times
  if(keeptimes[1] != xtimes[1]) stop(paste('1=',1, 'did not start in January'))
  
  # Loop through the regions
  for(i in 1:length(regs)){
    start = c(which.min(abs(lons - baselons$rng[i,1])), which.min(abs(lats - baselats$rng[i,1])), which.min(abs(depths - basedepths$rng[i,1])), xtimes[1]) # indices of the start (lon, lat, depth, time)
    end = c(which.min(abs(lons - baselons$rng[i,2])), which.min(abs(lats - baselats$rng[i,2])), which.min(abs(depths - basedepths$rng[i,2])), xtimes[length(xtimes)]) # indices of the end
    
    temp = ncvar_get(n, 'THETAO_REGRID', start = start, count=end - start+1) # a slice of time: dims are lon, lat, depth, time
    theselons = ncvar_get(n, 'LON', start=start[1], count= end[1]-start[1]+1)
    theselats = ncvar_get(n, 'LAT', start=start[2], count= end[2]-start[2]+1)
    thesedepths = ncvar_get(n, 'LEV', start=start[3], count= end[3]-start[3]+1)
    
    # trim to the appropriate months
    motokeep <- unlist(basemonths$months[basemonths$region == regs[i]])
    temp <- temp[,,,jmo2mo(xtimes) %in% motokeep]
    thesextimes <- xtimes[jmo2mo(xtimes) %in% motokeep]
    
    # reshape so that month and year are in separate dimensions
    dim(temp)
    temp2 = array(temp, dim=c(dim(temp)[1:3], length(motokeep), ceiling(length(thesextimes)/length(motokeep)))) # lon, lat, depth, month, year
    if(length(temp2) > length(temp)) temp2[(length(temp)+1):length(temp2)] <- NA # if needed, fill in extra values in temp2 with NAs (R otherwise recycles the values in temp)
    dim(temp2)
    
    # Average across months within each lat/lon/depth/year
    dim(temp2)
    temp3 = apply(temp2, MARGIN=c(1,2,3,5), FUN=mean)
    print(dim(temp3))
    
    temps2100[[i]] = temp3
    dimnames(temps2100[[i]]) = list(lon=theselons, lat=theselats, depth=thesedepths,yr = 2006:(2005+ceiling(length(thesextimes)/length(motokeep))))
  }
  
  # Write out
  outfile = paste('tempsRCP85_2006-2100_wSEUS', agency[1], '_', model[1], '.RData', sep='')
  save(temps2100, file=outfile)
  
  # close
  nc_close(n)
  
#######################################################
  
# Read in and write out each model's 2006-2100 period RCP4.5 (by year and season)
for(j in 1:length(agency)){ # very slow because of apply() down below.
        print(paste(j, agency[j], model[j], Sys.time()))
        file = paste(cmip5dir, 'Froelicher_data/', agency[j], '/', model[j], '/rcp45/thetao_all_regrid_malin.nc', sep='')
        n = nc_open(file) # Open the netCDF file
        lats = ncvar_get(n, 'LAT')
        lons = ncvar_get(n, 'LON')
        depths = ncvar_get(n, 'LEV')
        startmonth = 1; startyear = 2006
        if(model[j]=='HadGem2-CC'){     startmonth = 12; startyear = 2005} # HadGEM2-CC model is different for RCP4.5 (see readme.txt), starts Dec 2005
        times = 1:length(ncvar_get(n, 'TIME')) # time, measured in months since the start month/year for this model

        # Variable to hold the climate data in the base time period for each region
        temps2100 = vector('list', length(regs))
                names(temps2100) = regs

        # times to extract for each model: Jan 2006 to Dec 2100 = 1140 months. measured in months since the model's start month/year
        keeptimes = seq(from=(2006-startyear)*12+1-startmonth+1, length.out=1140) # time, measured in months since model's start month/year
        xtimes = intersect(times, keeptimes) # the intersecting times
        if(keeptimes[1] != xtimes[1]) stop(paste('j=',j, 'did not start in January'))

        # Loop through the regions
        for(i in 1:length(regs)){
                start = c(which.min(abs(lons - baselons$rng[i,1])), which.min(abs(lats - baselats$rng[i,1])), which.min(abs(depths - basedepths$rng[i,1])), xtimes[1]) # indices of the start (lon, lat, depth, time)
                end = c(which.min(abs(lons - baselons$rng[i,2])), which.min(abs(lats - baselats$rng[i,2])), which.min(abs(depths - basedepths$rng[i,2])), xtimes[length(xtimes)]) # indices of the end

                temp = ncvar_get(n, 'THETAO_REGRID', start = start, count=end - start+1) # a slice of time: dims are lon, lat, depth, time
                theselons = ncvar_get(n, 'LON', start=start[1], count= end[1]-start[1]+1)
                theselats = ncvar_get(n, 'LAT', start=start[2], count= end[2]-start[2]+1)
                thesedepths = ncvar_get(n, 'LEV', start=start[3], count= end[3]-start[3]+1)

                # trim to the appropriate months
                motokeep <- unlist(basemonths$months[basemonths$region == regs[i]])
                temp <- temp[,,,jmo2mo(xtimes) %in% motokeep]
                thesextimes <- xtimes[jmo2mo(xtimes) %in% motokeep]

                # reshape so that month and year are in separate dimensions
                dim(temp)
                temp2 = array(temp, dim=c(dim(temp)[1:3], length(motokeep), ceiling(length(thesextimes)/length(motokeep)))) # lon, lat, depth, month, year
                if(length(temp2) > length(temp)) temp2[(length(temp)+1):length(temp2)] <- NA # if needed, fill in extra values in temp2 with NAs (R otherwise recycles the values in temp)
                dim(temp2)

                # Average across months within each lat/lon/depth/year
                dim(temp2)
                temp3 = apply(temp2, MARGIN=c(1,2,3,5), FUN=mean)
                print(dim(temp3))

                temps2100[[i]] = temp3
                dimnames(temps2100[[i]]) = list(lon=theselons, lat=theselats, depth=thesedepths,yr = 2006:(2005+ceiling(length(thesextimes)/length(motokeep))))
        }

        # Write out
        outfile = paste('data/tempsRCP45_2006-2100_', agency[j], '_', model[j], '.RData', sep='')
        save(temps2100, file=outfile)

        # close
        nc_close(n)
}

# Read in, average, and write out each model's control run drift (as degrees per year)

for(j in 1:length(agency)){
        file = paste(cmip5dir, 'Froelicher_data/', agency[j], '/', model[j], '/piControl/thetao_all_regrid_malin.nc', sep='')
        print(paste(agency[j], model[j]))
        n = nc_open(file) # Open the netCDF file
        lats = ncvar_get(n, 'LAT')
        lons = ncvar_get(n, 'LON')
        depths = ncvar_get(n, 'LEV')
        startmonth = 1; startyear = 2006; seasonlist = c(1,1,1,2,2,2,3,3,3,4,4,4) # seasons Jan-Mar, Apr-Jun, Jul-Sep, Oct-Dec
        if(model[j]=='HadGem2-CC'){      # HadGEM2-CC model is different for RCP8.5 (see readme.txt), starts Dec 2005
                startmonth = 12; startyear = 2005
                seasonlist = c(4,1,1,1,2,2,2,3,3,3,4,4)
        }
        times = 1:length(ncvar_get(n, 'TIME')) # time, measured in months since the start month/year for this model

        # Variable to hold the smoothed climate data for the control runs (linear regression)
        tempscontrol = vector('list', length(regs))
                names(tempscontrol) = regs

        # Loop through the regions
        for(i in 1:length(regs)){
                start = c(which.min(abs(lons - baselons$rng[i,1])), which.min(abs(lats - baselats$rng[i,1])), which.min(abs(depths - basedepths$rng[i,1])), 1) # indices of the start (lon, lat, depth, time)
                end = c(which.min(abs(lons - baselons$rng[i,2])), which.min(abs(lats - baselats$rng[i,2])), which.min(abs(depths - basedepths$rng[i,2])), start[4]-2) # indices of the end

                temp = ncvar_get(n, 'THETAO_REGRID', start = start, count=end - start+1) # a slice of time
                theselons = ncvar_get(n, 'LON', start=start[1], count= end[1]-start[1]+1)
                theselats = ncvar_get(n, 'LAT', start=start[2], count= end[2]-start[2]+1)
                thesedepths = ncvar_get(n, 'LEV', start=start[3], count= end[3]-start[3]+1)

                # trim to the appropriate months
                motokeep <- unlist(basemonths$months[basemonths$region == regs[i]])
                temp <- temp[,,,jmo2mo(xtimes) %in% motokeep]
                thesextimes <- xtimes[jmo2mo(xtimes) %in% motokeep]

                # reshape so that month and year are in separate dimensions
                dim(temp)
                temp2 = array(temp, dim=c(dim(temp)[1:3], length(motokeep), ceiling(length(thesextimes)/length(motokeep)))) # lon, lat, depth, month, year
                if(length(temp2) > length(temp)) temp2[(length(temp)+1):length(temp2)] <- NA # if needed, fill in extra values in temp2 with NAs (R otherwise recycles the values in temp)
                dim(temp2)

                # Take slope of temperature across years within each lat/lon/depth
                dim(temp2)
                temp3 = apply(temp2, MARGIN=c(1,2,3), FUN=linearb, yr=rep(1:dim(temp2)[5], each=length(motokeep)))
                print(dim(temp3))
        
                tempscontrol[[i]] = temp3
                dimnames(tempscontrol[[i]]) = list(lon=theselons, lat=theselats, depth=thesedepths)
        }

        # Write out
        outfile = paste('data/tempscontrol_', agency[j], '_', model[j], '.RData', sep='')
        save(tempscontrol, file=outfile)

        # close
        nc_close(n)
}


###################################
## Read in data from RData files ##
###################################
agency = c('CNRM-CERFACS', 'IPSL', 'IPSL', 'MOHC', 'MPI-M', 'MPI-M', 'MRI', 'NCAR', 'NCC', 'NCC', 'NOAA-GFDL', 'NOAA-GFDL', 'NOAA-GFDL')
model = c('CNRM-CM5', 'IPSL-CM5A-MR', 'IPSL-CM5B-LR', 'HadGem2-CC', 'MPI-ESM-LR', 'MPI-ESM-MR', 'MRI-CGCM3', 'CCSM4', 'NorESM1-M', 'NorESM1-ME', 'GFDL-CM3', 'GFDL-ESM2G', 'GFDL-ESM2M')

load('output/basedates.RData')

# Read in data (historical, rcp85 2006-2100, rcp45 2006-2100, and control)
histl = vector('list', length(agency)) # "l" at end so that it matches tempscontroll and temps2100l85
	for(i in 1:length(agency)){
		print(i)
		infile = paste('data/tempshist_', agency[i], '_', model[i], '.RData', sep='')
		load(infile)
		histl[[i]] = basetemps
	}
	names(histl) = paste(agency, model, sep='_')
		
tempscontroll = vector('list', length(agency)) # "l" at end so that it doesn't get overwritten by load()
	for(i in 1:length(agency)){
		print(i)
		infile = paste('data/tempscontrol_', agency[i], '_', model[i], '.RData', sep='')
		load(infile)
		tempscontroll[[i]] = tempscontrol
	}
	names(tempscontroll) = paste(agency, model, sep='_')

temps2100l85 = vector('list', length(agency)) # "l" at end so that it doesn't get overwritten by load()
	for(i in 1:length(agency)){
		print(i)
		infile = paste('data/tempsRCP85_2006-2100_', agency[i], '_', model[i], '.RData', sep='')
		load(infile)
		dim(temps2100) # check dimensions
		dimnames(temps2100)
#		temps2100 = lapply(temps2100, FUN=aperm, perm=c(1,2,3,5,4)) # transpose so that year is last dimension. may not be needed
		temps2100l85[[i]] = temps2100
	}
	names(temps2100l85) = paste(agency, model, sep='_')

temps2100l45 = vector('list', length(agency)) # "l" at end so that it doesn't get overwritten by load()
	for(i in 1:length(agency)){
		print(i)
		infile = paste('data/tempsRCP45_2006-2100_', agency[i], '_', model[i], '.RData', sep='')
		load(infile)
		dim(temps2100) # check dimensions
		dimnames(temps2100)
#		temps2100 = lapply(temps2100, FUN=aperm, perm=c(1,2,3,5,4)) # transpose so that year is last dimension. may not be needed
		temps2100l45[[i]] = temps2100
	}
	names(temps2100l45) = paste(agency, model, sep='_')

regs = names(temps2100l85[[1]])	

#####################################
## Calculate climate change deltas ##
## projections - historical        ##
#####################################

# Calculate deltas: offsets from historical periods
	# RCP8.5 deltas in each year in each region in each model
	delta.raw2100rcp85 = vector('list', length(agency)) # each item is a model
	for(i in 1:length(agency)){ # for each model
		print(i)
		delta.raw2100rcp85[[i]] = vector('list', length(histl[[i]])) # each item is a region (in a model)
		for(j in 1:length(histl[[i]])){ # for each region
			a = histl[[i]][[j]] # dims are lon, lat, depth
			dim(a) = c(dim(a), 1) # add a fake year dimension
			# make sure lon/lat/depth dimensions of histl and temps2100 match
			a = a[,,,rep(1,dim(temps2100l85[[i]][[j]])[4])] # expand histl to match dimensions of temps2100l85
			delta.raw2100rcp85[[i]][[j]] = temps2100l85[[i]][[j]] - a
		}
		names(delta.raw2100rcp85[[i]]) = names(temps2100l85[[i]])
	}
	names(delta.raw2100rcp85) = names(temps2100l85)
	rm(a) # remove our copy of histl

	# write out
	save(delta.raw2100rcp85, file='data/delta.raw2100rcp85.Rdata')

	# RCP4.5 deltas in each year in each region in each model
	delta.raw2100rcp45 = vector('list', length(agency)) # each item is a model
	for(i in 1:length(agency)){ # for each model
		print(i)
		delta.raw2100rcp45[[i]] = vector('list', length(histl[[i]])) # each item is a region (in a model)
		for(j in 1:length(histl[[i]])){ # for each region
			a = histl[[i]][[j]] # dims are lon, lat, depth
			dim(a) = c(dim(a), 1) # add a fake year dimension
			# make sure lon/lat/depth dimensions of histl and temps2100 match
			a = a[,,,rep(1,dim(temps2100l45[[i]][[j]])[4])] # expand histl to match dimensions of temps2100l45
			delta.raw2100rcp45[[i]][[j]] = temps2100l45[[i]][[j]] - a
		}
		names(delta.raw2100rcp45[[i]]) = names(temps2100l45[[i]])
	}
	names(delta.raw2100rcp45) = names(temps2100l45)
	rm(a) # remove our copy of histl

	# write out
	save(delta.raw2100rcp45, file='data/delta.raw2100rcp45.Rdata')



###########################
# Plot and compare deltas #
###########################

load('data/delta.raw2100rcp85.Rdata')
load('data/delta.raw2100rcp45.Rdata')

		# Plot raw deltas for each model in each region (averaged across all depths and lon/lats)
		require(RColorBrewer)
		cols = c('#000000', brewer.pal(12, 'Set3')) 

		#NEED TO REDO BASETEMPS ETC. FROM NEAR TOP WITH THE UPDATED 'DAT' FILE AGAIN********************
		
			# RCP8.5
		# quartz(width=11, height=8.5)
		pdf(width=11, height=8.5, file=paste('figures/delta.raw.rcp85.pdf', sep=''))
		par(mfrow=c(4,3), mai=c(0.5, 0.5, 0.3, 0.1), mgp = c(2.3,0.8,0))
		ylims = c(-1,7.5)
		xlims=c(2006,2100)
		for(i in 1:length(delta.raw2100rcp85[[1]])){ # for each region
			print(regs[i])
			plot(1,1, col='white', xlim=xlims, ylim=ylims, xlab='Year', ylab='Delta (°C)', main=names(delta.raw2100rcp85[[1]])[i])
			for(j in 1:length(delta.raw2100rcp85)){ # for each model
				y = apply(delta.raw2100rcp85[[j]][[i]], MARGIN=4, FUN=mean, na.rm=TRUE) # mean by year within region across all lon/lat/depth
				x = (2006:2100)[1:length(y)]
				lines(x, y, col=cols[j]) # for 2006-2100
			}
		}
		legend('topleft', col=cols, legend=model, cex=0.8, lty=1, ncol=2, bty='n')
	
		dev.off()

			# RCP4.5
		# quartz(width=11, height=8.5)
		pdf(width=11, height=8.5, file=paste('figures/delta.raw.rcp45.pdf', sep=''))
		par(mfrow=c(4,3), mai=c(0.5, 0.5, 0.3, 0.1), mgp = c(2.3,0.8,0))
		ylims = c(-1,7.5)
		xlims=c(2006,2100)
		for(i in 1:length(delta.raw2100rcp45[[1]])){ # for each region
			print(regs[i])
			plot(1,1, col='white', xlim=xlims, ylim=ylims, xlab='Year', ylab='Delta (°C)', main=names(delta.raw2100rcp45[[1]])[i])
			for(j in 1:length(delta.raw2100rcp45)){ # for each model
				y = apply(delta.raw2100rcp45[[j]][[i]], MARGIN=4, FUN=mean, na.rm=TRUE) # mean by year within region across all lon/lat/depth
				x = (2006:2100)[1:length(y)]
				lines(x, y, col=cols[j]) # for 2006-2100
			}
		}
		legend('topleft', col=cols, legend=model, cex=0.8, lty=1, ncol=2, bty='n')
	
		dev.off()

		# Plot control deltas for each model in each region
		require(RColorBrewer)
		cols = c('#000000', brewer.pal(12, 'Set3')) 
		ymat = sapply(tempscontroll, FUN=function(x){ sapply(x, FUN=mean, na.rm=TRUE)}) # take mean by region and model

		#quartz(width=6, height=4)
		pdf(width=6, height=4, file=paste('figures/delta.cont.pdf', sep=''))
		par(mai=c(1.3, 0.7, 0.3, 0.1), mgp = c(2.3,0.8,0))
		plot(0,0, xlim=c(0,13), ylim=c(-0.017, 0.017), col='white', xlab='', ylab='Control Drift (°C/year)', xaxt='n')
		axis(1, at=1:length(regs), labels=regs, cex.axis=0.5, las=2)
		for(i in 1:length(tempscontroll)){ # for each model
			points(1:length(regs), ymat[,i], col=cols[i])
		}
		abline(h=0,col='grey')
		legend('topleft', col=cols, legend=model, cex=0.4, pch=1, ncol=2, bty='n')
	
		dev.off()

		# Compare magnitude of raw and control deltas (both in degC per year)
			# RCP85
		compare85 = data.frame(agency = rep(agency, each=length(regs)), model = rep(model, each=length(regs)), region = rep(regs, times=length(model)), rawdelt = numeric(length(regs)*length(model)), controldelt = numeric(length(regs)*length(model)))
		for(i in 1:length(delta.raw2100rcp85[[1]])){ # for each region
			for(j in 1:length(delta.raw2100rcp85)){ # for each model
				y = apply(delta.raw2100rcp85[[j]][[i]], MARGIN=4, FUN=mean, na.rm=TRUE) # mean within region across years # is year in dim 1??
				x = (2006:2100)[1:length(y)]
				ind = compare85$region==names(delta.raw2100rcp85[[1]])[i] & paste(as.character(compare85$agency), as.character(compare85$model), sep='_') == names(delta.raw2100rcp85)[j]
				compare85$rawdelt[ind] = coef(lm(y ~ x))[2]
				compare85$controldelt[ind] = mean(tempscontroll[[j]][[i]], na.rm=TRUE)
			}
		}
		compare85$prop = compare85$rawdelt/(compare85$rawdelt + compare85$controldelt)
		compare85$ratio = compare85$controldelt/compare85$rawdelt

		write.csv(compare85, file=paste('Tables/compare_delta_control.rcp85.csv', sep='')) # write out

			# RCP45
		compare45 = data.frame(agency = rep(agency, each=length(regs)), model = rep(model, each=length(regs)), region = rep(regs, times=length(model)), rawdelt = numeric(length(regs)*length(model)), controldelt = numeric(length(regs)*length(model)))
		for(i in 1:length(delta.raw2100rcp85[[1]])){ # for each region
			for(j in 1:length(delta.raw2100rcp85)){ # for each model
				y = apply(delta.raw2100rcp85[[j]][[i]], MARGIN=4, FUN=mean, na.rm=TRUE) # mean within region across years # is year in dim 1??
				x = (2006:2100)[1:length(y)]
				ind = compare45$region==names(delta.raw2100rcp85[[1]])[i] & paste(as.character(compare45$agency), as.character(compare45$model), sep='_') == names(delta.raw2100rcp85)[j]
				compare45$rawdelt[ind] = coef(lm(y ~ x))[2]
				compare45$controldelt[ind] = mean(tempscontroll[[j]][[i]], na.rm=TRUE)
			}
		}
		compare45$prop = compare45$rawdelt/(compare45$rawdelt + compare45$controldelt)
		compare45$ratio = compare45$controldelt/compare45$rawdelt

		write.csv(compare45, file=paste('Tables/compare_delta_control.rcp45.csv', sep='')) # write out


###################################################	
# Calculate total deltas (raw delta - drift*year) #
###################################################	
load('data/delta.raw2100rcp85.Rdata')
load('data/delta.raw2100rcp45.Rdata')

	# RCP85
	require(plyr) # for aaply() function
	delta2100rcp85 = vector('list', length(agency))
	options(warn=1)
	for(i in 1:length(agency)){ # for each model
		print(paste(i, model[i]))
		delta2100rcp85[[i]] = vector('list', length(histl[[i]]))
		for(j in 1:length(histl[[i]])){ # for each region
			print(regs[j])
			#print(as.character(basedates$region[j])) # should be the same as regs[j]
			baseyear = mean(c(nyrs(basedates$min[j]), nyrs(basedates$max[j])))
			a = tempscontroll[[i]][[j]]
			dim(a) = c(dim(a), 1) # add a year dimension
			a = a[,,,rep(1,dim(temps2100l85[[i]][[j]])[4])] # expand a to match dimensions of temps2100l85	(80 years)
			a = aaply(a, .margins=c(1,2,3), .fun= function(x){return(x*((2006:2100)-baseyear))}) # multiply degC/yr by years to get degC of drift
			delta2100rcp85[[i]][[j]] = delta.raw2100rcp85[[i]][[j]] - a
		}
		names(delta2100rcp85[[i]]) = names(delta.raw2100rcp85[[i]])
	}
	names(delta2100rcp85) = names(delta.raw2100rcp85)

	print(warnings())

	save(delta2100rcp85, file=paste('data/delta2100rcp85.Rdata', sep='')) # write out

	# RCP45
	require(plyr) # for aaply() function
	delta2100rcp45 = vector('list', length(agency))
	options(warn=1)
	for(i in 1:length(agency)){ # for each model
		print(paste(i, model[i]))
		delta2100rcp45[[i]] = vector('list', length(histl[[i]]))
		for(j in 1:length(histl[[i]])){ # for each region
			print(regs[j])
			#print(as.character(basedates$region[j])) # should be the same as regs[j]
			baseyear = mean(c(nyrs(basedates$min[j]), nyrs(basedates$max[j])))
			a = tempscontroll[[i]][[j]]
			dim(a) = c(dim(a), 1) # add a year dimension
			a = a[,,,rep(1,dim(temps2100l45[[i]][[j]])[4])] # expand a to match dimensions of temps2100l85	(80 years)
			endyr <- 2006 + dim(a)[4] - 1 # needed because HadGem-CC only goes to 2099
			a = aaply(a, .margins=c(1,2,3), .fun= function(x){return(x*((2006:endyr)-baseyear))}) # multiply degC/yr by years to get degC of drift
			delta2100rcp45[[i]][[j]] = delta.raw2100rcp45[[i]][[j]] - a
		}
		names(delta2100rcp45[[i]]) = names(delta.raw2100rcp45[[i]])
	}
	names(delta2100rcp45) = names(delta.raw2100rcp45)

	print(warnings())

	save(delta2100rcp45, file=paste('data/delta2100rcp45.Rdata', sep='')) # write out



#############################################
# Plot deltas for each model in each region
#############################################

load('data/delta2100rcp85.RData')
load('data/delta2100rcp45.RData')

	# RCP85
	require(RColorBrewer)
	cols = c('#000000', brewer.pal(12, 'Set3')) 

	#quartz(width=11, height=8.5)
	pdf(width=11, height=8.5, file=paste('figures/deltas.rcp85.pdf', sep=''))
	par(mfrow=c(4,3), mai=c(0.5, 0.5, 0.3, 0.1), mgp = c(2.3,0.8,0))
	ylims = c(-1,7.5)
	xlims=c(2006,2100)
	for(i in 1:length(delta2100rcp85[[1]])){ # for each region
		plot(1,1, col='white', xlim=xlims, ylim=ylims, xlab='Year', ylab='Delta (°C)', main=names(delta2100rcp85[[1]])[i])
		for(j in 1:length(delta2100rcp85)){ # for each model
			y = apply(delta2100rcp85[[j]][[i]], MARGIN=4, FUN=mean, na.rm=TRUE)
			x = as.numeric(names(y))
			lines(x, y, col=cols[j]) # for 2006-2100
		}
	}
	legend('topleft', col=cols, legend=model, cex=0.8, lty=1, ncol=2, bty='n')

	dev.off()

	# RCP45
	require(RColorBrewer)
	cols = c('#000000', brewer.pal(12, 'Set3')) 

	#quartz(width=11, height=8.5)
	pdf(width=11, height=8.5, file=paste('figures/deltas.rcp45.pdf', sep=''))
	par(mfrow=c(4,3), mai=c(0.5, 0.5, 0.3, 0.1), mgp = c(2.3,0.8,0))
	ylims = c(-1,7.5)
	xlims=c(2006,2100)
	for(i in 1:length(delta2100rcp45[[1]])){ # for each region
		plot(1,1, col='white', xlim=xlims, ylim=ylims, xlab='Year', ylab='Delta (°C)', main=names(delta2100rcp45[[1]])[i])
		for(j in 1:length(delta2100rcp45)){ # for each model
			y = apply(delta2100rcp45[[j]][[i]], MARGIN=4, FUN=mean, na.rm=TRUE)
			x = as.numeric(names(y))
			lines(x, y, col=cols[j]) # for 2006-2100
		}
	}
	legend('topleft', col=cols, legend=model, cex=0.8, lty=1, ncol=2, bty='n')

	dev.off()


###################################	
# convert deltas to long format   #
# (for each region of each model) #
# takes lots of RAM               #
###################################
load('data/delta2100rcp85.RData')
load('data/delta2100rcp45.RData')
require(reshape2) # for melt()

# RCP85
delta2100rcp85long = vector('list', length(delta2100rcp85)) # one for each model
names(delta2100rcp85long) =names(delta2100rcp85)
an = numeric(0) # a placeholder
ac = character(0)
for(i in 1:length(delta2100rcp85)){ # for each climate model
	print(i)
	aa = data.frame(region = ac, year = an, lat = an, lon=an, depth=an, delta = an) # holds all regions for this model 2100
	for(j in 1:length(delta2100rcp85[[i]])){ # for each region
		cat(paste(j, ' ', sep=''))
		b = melt(delta2100rcp85[[i]][[j]]) # very fast
		b$region = names(delta2100rcp85[[i]])[j]
		names(b)[names(b) == 'yr'] = 'year'
		names(b)[names(b) == 'value'] = 'delta'
		aa = rbind(aa,b)	
	}


	delta2100rcp85long[[i]] = aa
	cat('\n')
}	

save(delta2100rcp85long, file=paste('data/delta2100rcp85long.RData', sep='')) # write out

# RCP45
delta2100rcp45long = vector('list', length(delta2100rcp45)) # one for each model
names(delta2100rcp45long) =names(delta2100rcp45)
an = numeric(0) # a placeholder
ac = character(0)
for(i in 1:length(delta2100rcp45)){ # for each climate model
	print(i)
	aa = data.frame(region = ac, year = an, lat = an, lon=an, depth=an, delta = an) # holds all regions for this model 2100
	for(j in 1:length(delta2100rcp45[[i]])){ # for each region
		cat(paste(j, ' ', sep=''))
		b = melt(delta2100rcp45[[i]][[j]]) # very fast
		b$region = names(delta2100rcp45[[i]])[j]
		names(b)[names(b) == 'yr'] = 'year'
		names(b)[names(b) == 'value'] = 'delta'
		aa = rbind(aa,b)	
	}


	delta2100rcp45long[[i]] = aa
	cat('\n')
}	

save(delta2100rcp45long, file=paste('data/delta2100rcp45long.RData', sep='')) # write out


