# Combine the trawl survey data together and reformat in preparation for creating climatologies and climate envelope models

library(data.table) # much of this code could be sped up by converting to data.tables
library(Hmisc)
library(stringr)
library(zoo)

# useful function: acts like sum(na.rm=T) but returns NA if all are NA
	sumna = function(x){
		if(!all(is.na(x))) return(sum(x, na.rm=T))
		if(all(is.na(x))) return(NA)
	}


###############################
## Read in and reformat data ##
###############################
setwd('/Users/mpinsky/Documents/Rutgers/Range projections')

# Read in data
	# Aleutians
	ai1 = read.csv('../NorthAmerican_survey_data/AFSC_Aleutians/2013-10-17/ai1983_2000.csv')
	ai2 = read.csv('../NorthAmerican_survey_data/AFSC_Aleutians/2013-10-17/ai2002_2012.csv')	
	ai = rbind(ai1, ai2)
	rm(ai1, ai2)

	# Eastern Bering Sea
	ebs1 = read.csv('../NorthAmerican_survey_data/AFSC_EBS/2013-10-17/ebs1982_1984.csv')
	ebs2 = read.csv('../NorthAmerican_survey_data/AFSC_EBS/2013-10-17/ebs1985_1989.csv')
	ebs3 = read.csv('../NorthAmerican_survey_data/AFSC_EBS/2013-10-17/ebs1990_1994.csv')
	ebs4 = read.csv('../NorthAmerican_survey_data/AFSC_EBS/2013-10-17/ebs1995_1999.csv')
	ebs5 = read.csv('../NorthAmerican_survey_data/AFSC_EBS/2013-10-17/ebs2000_2004.csv')
	ebs6 = read.csv('../NorthAmerican_survey_data/AFSC_EBS/2013-10-17/ebs2005_2008.csv')
	ebs7 = read.csv('../NorthAmerican_survey_data/AFSC_EBS/2013-10-17/ebs2009_2012.csv')
	ebs8 = read.csv('../NorthAmerican_survey_data/AFSC_EBS/2013-10-17/ebs2013.csv')
	ebs = rbind(ebs1, ebs2, ebs3, ebs4, ebs5, ebs6, ebs7, ebs8)
	rm(ebs1, ebs2, ebs3, ebs4, ebs5, ebs6, ebs7)

	# Gulf of Alaska
	goa1 = read.csv('../NorthAmerican_survey_data/AFSC_GOA/2013-10-17/goa1984_1987.csv')
	goa2 = read.csv('../NorthAmerican_survey_data/AFSC_GOA/2013-10-17/goa1990_1999.csv')
	goa3 = read.csv('../NorthAmerican_survey_data/AFSC_GOA/2013-10-17/goa2001_2005.csv')
	goa4 = read.csv('../NorthAmerican_survey_data/AFSC_GOA/2013-10-17/goa2007_2013.csv')
	goa = rbind(goa1, goa2, goa3, goa4)
	rm(goa1, goa2, goa3, goa4)

	# Northeast US
	load("../NorthAmerican_survey_data/NEFSC/2015-01-30/Survdat.Rdata") # survdat. comes as a data.table
	load("../NorthAmerican_survey_data/NEFSC/2015-01-30/SVSPP.RData") # species names. a data.table
	setkey(survdat, CRUISE6, STATION, STRATUM, SVSPP, CATCHSEX)
	neus <- unique(survdat) # drops length data
	neus[, c('LENGTH', 'NUMLEN') := NULL] # remove length columns
	neus = neus[,sum(BIOMASS),by=list(YEAR, EST_TOWDATE, SEASON, LAT, LON, DEPTH, SURFTEMP, BOTTEMP, CRUISE6, STATION, STRATUM, SVSPP)] # sum different sexes of same spp together
	setnames(neus, 'V1', 'wtcpue')
	spp[,c('ITISSPP', 'COMNAME', 'AUTHOR') := NULL] # remove some columns from spp data.table
	neus = merge(neus, spp, by='SVSPP') # add species names
	neus = as.data.frame(neus) # this makes the calculations less efficient... but avoids having to rewrite the code for data.tables
	rm(survdat, spp)

	# Southeast US, comes from two surveys, SEAMAP (bottom trawl survey) and MARMAP (offshore trap survey)
  # SEAMAP trawl survey
	survcatch = read.csv('/Users/abigailporay/Documents/Work/Rutgers/Data/SEAMAP data/seamap_catch_2015.csv', stringsAsFactors=FALSE)
	survhaul = read.csv('/Users/abigailporay/Documents/Work/Rutgers/Data/SEAMAP data/seamap_event_2015.csv', stringsAsFactors=FALSE) # only needed for depth
	survhaul = unique(data.frame(EVENTNAME = survhaul$EVENTNAME, DEPTHSTART = survhaul$DEPTHSTART))
	seusstrata = read.csv('/Users/abigailporay/Documents/Work/Rutgers/Data/SEAMAP data/seamap.Strata.Area.csv') # Not sure if necessary
	seus = merge(x=survcatch, y=survhaul, all.x=T, by="EVENTNAME") # Add depth data from survhaul 
	seus = cbind(seus, STRATA = as.integer(str_sub(string = seus$STATIONCODE, start = 1, end = 2))) #Create STRATA column
	seus = seus[seus$DEPTHZONE != "OUTER",] # Drop OUTER depth zone because it was only sampled for 10 years, and in specific seasons-areas
	seus = merge(x=seus, y=seusstrata, all.x=TRUE, by='STRATA') #add STRATAHECTARE to main file 

	# West Coast Trienniel (1977-2004)
	wctricatch = read.csv('../NorthAmerican_survey_data/AFSC_WestCoast/2011-12-08/CATCHWCTRIALLCOAST.csv')
	wctrihaul = read.csv('../NorthAmerican_survey_data/AFSC_WestCoast/2011-12-08/HAULWCTRIALLCOAST.csv')
	wctrispecies = read.csv('../NorthAmerican_survey_data/AFSC_WestCoast/2011-12-08/RACEBASE_SPECIES.csv')
	wctri = merge(wctricatch[,c('CRUISEJOIN', 'HAULJOIN', 'VESSEL', 'CRUISE', 'HAUL', 'SPECIES_CODE', 'WEIGHT')], wctrihaul[,c('CRUISEJOIN', 'HAULJOIN', 'VESSEL', 'CRUISE', 'HAUL', 'HAUL_TYPE', 'PERFORMANCE', 'START_TIME', 'DURATION', 'DISTANCE_FISHED', 'NET_WIDTH', 'STRATUM', 'START_LATITUDE', 'END_LATITUDE', 'START_LONGITUDE', 'END_LONGITUDE', 'STATIONID', 'BOTTOM_DEPTH', 'SURFACE_TEMPERATURE', 'GEAR_TEMPERATURE')], all.x=TRUE) # Add haul info to catch data
	wctri = merge(wctri, wctrispecies[,c('SPECIES_CODE', 'SPECIES_NAME', 'COMMON_NAME')]) #  add species names
	wctri = wctri[wctri$HAUL_TYPE==3 & wctri$PERFORMANCE==0,] # trim to standard hauls and good performance
	rm(wctricatch, wctrihaul, wctrispecies)
	
	# West Coast annual (2003-2012)
	wcannfish = read.csv('../NorthAmerican_survey_data/NWFSC/2014-02-11/wcann2003_2012fish.csv')
	wcannhaul = read.csv('../NorthAmerican_survey_data/NWFSC/2014-02-11/wcann2003_2012haul.csv')
	wcanninvert = read.csv('../NorthAmerican_survey_data/NWFSC/2014-02-11/wcann2003_2012invert.csv')
		wcanninvert$Individual.Average.Weight..kg. = NA # add a column to allow use of rbind with wcannfish
	wcanncatch = rbind(wcannfish, wcanninvert)
	wcann = merge(wcannhaul, wcanncatch)
	rm(wcannfish, wcannhaul, wcanninvert, wcanncatch)

	# Gulf of Mexico: requires a lot of prep before a single data.frame can be created
	gmexstation = read.csv('../NorthAmerican_survey_data/SEAMAP-GMex/2014-06-25/STAREC.csv')
	gmextow = read.csv('../NorthAmerican_survey_data/SEAMAP-GMex/2014-06-25/INVREC.csv')
	gmexspp = read.csv('../NorthAmerican_survey_data/SEAMAP-GMex/2014-06-25/NEWBIOCODESBIG.csv')
	gmexcruise = read.csv('../NorthAmerican_survey_data/SEAMAP-GMex/2014-06-25/CRUISES.csv')
	test = read.csv('../NorthAmerican_survey_data/SEAMAP-GMex/2014-06-25/BGSREC.csv', nrows=2) # gmexbio is a large file: only read in some columns
	biocols = c('CRUISEID', 'STATIONID', 'VESSEL', 'CRUISE_NO', 'P_STA_NO', 'GENUS_BGS', 'SPEC_BGS', 'BGSCODE', 'BIO_BGS', 'SELECT_BGS')
	colstoread = rep('NULL', ncol(test)) # NULL means don't read that column (see ?read.csv)
	colstoread[names(test) %in% biocols] = NA # NA means read that column
	gmexbio = read.csv('../NorthAmerican_survey_data/SEAMAP-GMex/2014-06-25/BGSREC.csv', colClasses=colstoread) # sped up by reading in only some columns
		# trim out young of year records (only useful for count data) and those with UNKNOWN species
		gmexbio = gmexbio[gmexbio$BGSCODE != 'T' & gmexbio$GENUS_BGS != 'UNKNOWN',]
		gmexbio = gmexbio[!duplicated(gmexbio),] # remove the few rows that are still duplicates
	newspp = data.frame(Key1 = c(503,5770), TAXONOMIC = c('ANTHIAS TENUIS AND WOODSI', 'MOLLUSCA AND UNID.OTHER #01'), CODE=c(170026003, 300000000), TAXONSIZECODE=NA, isactive=-1, common_name=c('threadnose and swallowtail bass', 'molluscs or unknown'), tsn = NA) # make two combined records where multiple species records share the same species code
	gmexspp = gmexspp[!(gmexspp$CODE %in% gmexspp$CODE[which(duplicated(gmexspp$CODE))]),] # remove the duplicates that were just combined
	gmexspp = rbind(gmexspp[,names(newspp)], newspp) # add the combined records on to the end. trim out extra columns from gmexspp

	gmex = merge(gmexbio, gmextow[gmextow$GEAR_TYPE=='ST', c('STATIONID', 'CRUISE_NO', 'P_STA_NO', 'INVRECID', 'GEAR_SIZE', 'GEAR_TYPE', 'MESH_SIZE', 'MIN_FISH', 'OP')], all.x=TRUE) # merge tow information with catch data, but only for shrimp trawl tows (ST)
	gmex = merge(gmex, gmexstation[,c('STATIONID', 'CRUISEID', 'CRUISE_NO', 'P_STA_NO', 'TIME_ZN', 'TIME_MIL', 'S_LATD', 'S_LATM', 'S_LOND', 'S_LONM', 'E_LATD', 'E_LATM', 'E_LOND', 'E_LONM', 'DEPTH_SSTA', 'TEMP_BOT', 'TEMP_SSURF', 'MO_DAY_YR', 'VESSEL_SPD', 'COMSTAT')], all.x=TRUE) # add station location and related data
	gmex = merge(gmex, gmexspp[,c('CODE', 'TAXONOMIC')], by.x='BIO_BGS', by.y='CODE', all.x=TRUE) # add scientific name
	gmex = merge(gmex, gmexcruise[,c('CRUISEID', 'VESSEL', 'TITLE')], all.x=TRUE) # add cruise title
	gmex = gmex[gmex$TITLE %in% c('Summer SEAMAP Groundfish Survey', 'Summer SEAMAP Groundfish Suvey') & as.numeric(as.character(gmex$GEAR_SIZE))==40 & gmex$MESH_SIZE == 1.63 & !is.na(gmex$MESH_SIZE) & gmex$OP %in% c(''),] # # Trim to high quality SEAMAP summer trawls, based off the subset used by Jeff Rester's GS_TRAWL_05232011.sas
	rm(gmexstation, gmextow, gmexspp, gmexcruise, test, gmexbio, newspp, colstoread)
 
	# Newfoundland (takes 3 min)
	path = '../NorthAmerican_survey_data/DFO_Newfoundland/2012-03-29/Spring and Fall/'
	files = list.files(path=path, pattern = "199[23456789]|200[0123456789]|201[012]") # only 1992 and later, since data format changes (more columns added)
	n = numeric(0); ch = numeric(0)
	newf = data.frame(recordtype = n, vessel = n, trip = n, set = n, yearl=n, monthl = n, dayl = n, settype = n, stratum = n, nafo = ch, unitarea = ch, light = n, winddir = n, windforce = n, sea = n, bottom = n, timel = n, duration = n, distance = n, operation = n, depth = n, depthmin = n, depthmax = n, depthbottom = n, surftemp = n, bottemp = n, latstart = n, lonstart = n, posmethod = n, gear = n, sppcode = n, num = n, wgt = n, latend = n, lonend = n, bottempmeth = n, geardevice = n)
	options(warn = 1)
	for(i in 1:length(files)){ # for each file
		if(i == 1) print(paste(length(files), 'files to read')) # 45 files
		cat(paste(i,' ', sep=''))
		indata = read.fwf(file=paste(path, files[i], sep=''), widths=c(
		1, # record type
		2, # vessel
		3, # trip
		3, # set
		2, # year
		2, # mo
		2, # day
		2, # set type
		3, # stratum
		2, # nafo
		3, # unit
		3, # light
		1, # winddir
		1, # wind force
		1, # sea
		1, # bottom type
		4, # time
		3, # duration
		3, # distance 
		1, # operation
		4, # depth mean
		4, # depth min
		4, # depth max
		4, # depth bottom
		3, # temp surf
		3, # temp bot
		5, # lat start
		5, # lon start
		1, # pos meth
		4, # gear
		4, # sppcode
		6, # number
		7, # wgt 
		5, # lat end
		5, # lon end
		2, # bot temp device
		2), # gear mon device
		header= FALSE, stringsAsFactors = FALSE)
		names(indata) = c('recordtype', 'vessel', 'trip', 'set', 'yearl', 'monthl', 'dayl', 'settype', 'stratum', 'nafo', 'unitarea', 'light', 'winddir', 'windforce', 'sea', 'bottom', 'timel', 'duration', 'distance', 'operation', 'depth', 'depthmin', 'depthmax', 'depthbottom', 'surftemp', 'bottemp', 'latstart', 'lonstart', 'posmethod', 'gear', 'sppcode', 'num', 'wgt', 'latend', 'lonend', 'bottempmeth', 'geardevice')
		newf = rbind(newf, indata)
	}
	spp = read.csv('../NorthAmerican_survey_data/DFO_Newfoundland/2012-03-29/Tables/GFSPCY.CODE_2012-07-05.csv', row.names=1) # add species names
	newf = merge(newf, spp, all.x=TRUE)
	newf = newf[newf$gear == 61 & !is.na(newf$gear),] # Trim to Campelen 1800 lined shrimp trawl gear
	newf = newf[newf$settype == 1,] # trim to set-type == survey
	
	
	# Scotian Shelf	
	scot = read.csv('../NorthAmerican_survey_data/DFO_ScotianShelf/2012-03-26/gscat_adj_pinsky.csv', header=TRUE)
#		names(scot)
#		dim(scot)
#		summary(scot)
	scotsetdata = read.csv('../NorthAmerican_survey_data/DFO_ScotianShelf/2012-03-26/gsinf_pinsky.csv')
#		names(scotsetdata)
#		names(scotsetdata)[names(scotsetdata)=='REMARKS'] = 'REMARKS.SET'
#		dim(scotsetdata)
#		summary(scotsetdata)
	scotspecies = read.csv('../NorthAmerican_survey_data/DFO_ScotianShelf/2012-03-26/species list.csv')
#		head(scotspecies)

	scot = merge(scot, scotsetdata, all.x=T, all.y=T, by=c('MISSION', 'SETNO')) 	# Add set data to catch data and make sure all sets are included
		dim(scot) # 185,324
	
		# Turn NAs to 0 where sets not in catch data were added
		i = which(is.na(scot$SPEC) & is.na(scot$ADJ_TOTWGT))
		scot$SPEC[i] = sort(unique(scot$SPEC))[1] # add a place-holder spp so that the row is not NA
		scot$ADJ_TOTWGT[i] = 0
		scot$ADJ_TOTNO[i] = 0

	scot = merge(scot, scotspecies, by.x='SPEC', by.y='CODE', all.x=TRUE) # Add species names
		dim(scot) # 185324
		
	# S. Gulf of St. Lawrence
	sgslcatchdata = read.csv('../NorthAmerican_survey_data/DFO_SouthernGulf/2011-12-08/southern Gulf survey data.csv')
#		names(sgslcatchdata)
#		dim(sgslcatchdata)
	sgslsetdata = read.csv('../NorthAmerican_survey_data/DFO_SouthernGulf/2011-12-08/sGSL_RV Survey sets_1971_2009.csv', na.strings='.')
#		names(sgslsetdata)
#		dim(sgslsetdata)
#		intersect(names(sgslcatchdata), names(sgslsetdata))

	names(sgslcatchdata)[names(sgslcatchdata) == 'expt'] = 'exptCatch' # so they don't conflict with names in setdata
	names(sgslcatchdata)[names(sgslcatchdata) == 'time'] = 'timeCatch'
	names(sgslcatchdata)[names(sgslcatchdata) == 'depth'] = 'depthCatch'
	intersect(names(sgslcatchdata), names(sgslsetdata)) # names that catchdata and setdata could match on
	sgsl = merge(sgslcatchdata, sgslsetdata[,c('vessel', 'cruise', 'year', 'set', 'expt', 'strat', 'month', 'day', 'time', 'depth', 't_surface', 't_bottom', 'salin_bottom')], all.x=T, all.y=T) # add temp and salinity
		dim(sgsl) # 180652 (69 more rows than with all.y=F)
	
	i = which(is.na(sgsl$latin_name) & is.na(sgsl$biomass)) # Fill in new tows that were just added
		length(i) # 69: good, matches above
		sgsl$latin_name[i] = sgsl$latin_name[1]
		sgsl$name[i] = sgsl$name[1]
		sgsl$biomass[i] = 0
		sgsl$catch[i] = 0
#		sgsl[i,]
#		sort(unique(sgsl$strat[i]))
	
 
# Create a unique haulid
	ai$haulid = paste(formatC(ai$VESSEL, width=3, flag=0), formatC(ai$CRUISE, width=3, flag=0), formatC(ai$HAUL, width=3, flag=0), sep='-')
	ebs$haulid = paste(formatC(ebs$VESSEL, width=3, flag=0), formatC(ebs$CRUISE, width=3, flag=0), formatC(ebs$HAUL, width=3, flag=0), sep='-')
	goa$haulid = paste(formatC(goa$VESSEL, width=3, flag=0), formatC(goa$CRUISE, width=3, flag=0), formatC(goa$HAUL, width=3, flag=0), sep='-')
	neus$haulid = paste(formatC(neus$CRUISE6, width=6, flag=0), formatC(neus$STATION, width=3, flag=0), formatC(neus$STRATUM, width=4, flag=0), sep='-') 
	seus$haulid = seus$EVENTNAME
	seus.shelf$haulid = seus.shelf$EVENTNAME
	wctri$haulid = paste(formatC(wctri$VESSEL, width=3, flag=0), formatC(wctri$CRUISE, width=3, flag=0), formatC(wctri$HAUL, width=3, flag=0), sep='-')
	wcann$haulid = wcann$Trawl.Id
	gmex$haulid = paste(formatC(gmex$VESSEL, width=3, flag=0), formatC(gmex$CRUISE_NO, width=3, flag=0), formatC(gmex$P_STA_NO, width=5, flag=0, format='d'), sep='-')
	newf$haulid = paste(formatC(newf$vessel, width=2, flag=0), formatC(newf$trip, width=3, flag=0), formatC(newf$set, width=3, flag=0, format='d'), sep='-')
	scot$haulid = paste(as.character(scot$MISSION), formatC(scot$SETNO, width=3, flag=0), sep='-')
	sgsl$haulid = paste(as.character(sgsl$vessel), formatC(sgsl$cruise, width=3, flag=0), formatC(sgsl$set, width=3, flag=0), sep='-')

# Extract year where needed
	wctri$year = as.numeric(substr(wctri$CRUISE, 1,4))
	wcann$year = as.numeric(gsub('Cycle ', '', wcann$Survey.Cycle))
	gmex$year = as.numeric(unlist(strsplit(as.character(gmex$MO_DAY_YR), split='-'))[seq(1,by=3,length=nrow(gmex))])
	newf$year = newf$yearl + 1900 # l stands for local date/time
	newf$year[newf$year<1950] = newf$year[newf$year<1950] + 100 # correct for 2000s
	scot$year = as.numeric(substr(as.character(scot$MISSION), 4,7))
  seus$date = as.Date(seus$DATE, "%m/%d/%y")
	seus <- cbind(seus, year = year(seus$date), month = month(seus$date)) # also made month column for seus in this step
	seus$season <- as.yearqtr(seus$date) # make a 'season' column to distinguish the spring, summer, and fall surveys
	seus$season <- factor(format(seus$season, "%q"), levels = 1:4, labels = c("winter", "spring", "summer", "fall")) # takes ~30 sec
	seus$season[seus$month == 9] <- "fall" 	#Sept EVENTS (all late-sept.) were grouped with summer, should be fall

# Extract month where needed
	ai$month = as.numeric(unlist(strsplit(as.character(ai$DATETIME), split='/'))[seq(1,length=nrow(ai), by=3)])
	ebs$month = as.numeric(unlist(strsplit(as.character(ebs$DATETIME), split='/'))[seq(1,length=nrow(ebs), by=3)])
	goa$month = as.numeric(unlist(strsplit(as.character(goa$DATETIME), split='/'))[seq(1,length=nrow(goa), by=3)])
	neus$month = as.numeric(unlist(strsplit(as.character(neus$EST_TOWDATE), split='-'))[seq(2,length=nrow(neus), by=3)])
	wctri$month = as.numeric(unlist(strsplit(as.character(wctri$START_TIME), split='/'))[seq(1,length=nrow(wctri), by=3)])
	wcann$month = as.numeric(unlist(strsplit(as.character(wcann$Trawl.Date), split='/'))[seq(1,length=nrow(wcann), by=3)])
	gmex$month = as.numeric(unlist(strsplit(as.character(gmex$MO_DAY_YR), split='-'))[seq(2,by=3,length=nrow(gmex))])
	newf$month = newf$monthl # l stands for local date/time
	scot$month = NA
	scot$month[grepl('/', scot$SDATE)] = as.numeric(unlist(strsplit(as.character(scot$SDATE[grepl('/', scot$SDATE)]), split='/'))[seq(2,by=3,length=sum(grepl('/', scot$SDATE)))]) # when date formatted with /
	scot$month[grepl('-', scot$SDATE)] = as.numeric(unlist(strsplit(as.character(scot$SDATE[grepl('-', scot$SDATE)]), split='-'))[seq(2,by=3,length=sum(grepl('-', scot$SDATE)))]) # when date formatted with -
	# already have month for sgsl

# Calculate decimal lat and lon, depth in m, where needed
	gmex$S_LATD[gmex$S_LATD == 0] = NA
	gmex$S_LOND[gmex$S_LOND == 0] = NA
	gmex$E_LATD[gmex$E_LATD == 0] = NA
	gmex$E_LOND[gmex$E_LOND == 0] = NA
	gmex$lat = rowMeans(cbind(gmex$S_LATD + gmex$S_LATM/60, gmex$E_LATD + gmex$E_LATM/60), na.rm=T) # mean of start and end positions, but allow one to be NA (missing)
	gmex$lon = -rowMeans(cbind(gmex$S_LOND + gmex$S_LONM/60, gmex$E_LOND + gmex$E_LONM/60), na.rm=T) # need negative sign since western hemisphere
	gmex$depth = gmex$DEPTH_SSTA*1.8288 # convert fathoms to meters

	newf$lat = NA
	i = newf$latstart>0 & newf$latend > 0
	newf$lat[i] = (as.numeric(substr(newf$latstart[i], 1, 2)) + as.numeric(substr(newf$latstart[i], 3, 5))/600 + as.numeric(substr(newf$latend[i], 1, 2)) + as.numeric(substr(newf$latend[i], 3, 5))/600)/2 # lat as mean of start and end lat. format of latstart and latend is DDMM.M
	i = newf$latstart>0 & newf$latend == 0 & !is.na(newf$latstart)
	newf$lat[i] = as.numeric(substr(newf$latstart[i], 1, 2)) + as.numeric(substr(newf$latstart[i], 3, 5))/600 # only use start if end is not available	
	newf$lon = NA
	i = newf$lonstart>0 & newf$lonend > 0
	newf$lon[i] = -(as.numeric(substr(newf$lonstart[i], 1, 2)) + as.numeric(substr(newf$lonstart[i], 3, 5))/600 + as.numeric(substr(newf$lonend[i], 1, 2)) + as.numeric(substr(newf$lonend[i], 3, 5))/600)/2 # lon as mean of start and end lon. format of lonstart and lonend is DDMM.M
	i = newf$lonstart>0 & newf$lonend == 0 & !is.na(newf$latstart)
	newf$lon[i] = -(as.numeric(substr(newf$lonstart[i], 1, 2)) + as.numeric(substr(newf$lonstart[i], 3, 5))/600) # only use start if end is not available
	if(class(newf$depth)=='character') newf$depth = as.numeric(newf$depth)

	scot$lat = as.numeric(substr(scot$SLAT,1,2)) + as.numeric(substr(scot$SLAT,3,4))/60
	scot$lon = -as.numeric(substr(scot$SLON,1,2)) - as.numeric(substr(scot$SLON,3,4))/60
	scot$depth = rowMeans(cbind(scot$DMIN, scot$DMAX))
 

# Calculate temperature where needed
	wcann$surftemp = NA # field is not collected, apparently (or was not provided)

	# Newfoundland one decimal place. 900 means negative
	i = newf$surftemp >= 900 & !is.na(newf$surftemp)
	newf$surftemp[i] = -(newf$surftemp[i] - 900)/10
	i = newf$surftemp < 900 & newf$surftemp > 0 & !is.na(newf$surftemp)
	newf$surftemp[i] = newf$surftemp[i]/10
	i = newf$bottemp >= 900 & !is.na(newf$bottemp)
	newf$bottemp[i] = -(newf$bottemp[i] - 900)/10
	i = newf$bottemp < 900 & newf$bottemp > 0 & !is.na(newf$bottemp)
	newf$bottemp[i] = newf$bottemp[i]/10

	# Fix -9999 to NA for SST and BT
	ai$BOT_TEMP[ai$BOT_TEMP==-9999] = NA
	ai$SURF_TEMP[ai$SURF_TEMP==-9999] = NA
	ebs$BOT_TEMP[ebs$BOT_TEMP==-9999] = NA
	ebs$SURF_TEMP[ebs$SURF_TEMP==-9999] = NA
	goa$BOT_TEMP[goa$BOT_TEMP==-9999] = NA
	goa$SURF_TEMP[goa$SURF_TEMP==-9999] = NA

	# The SST entries on Scotian Shelf in 2010 and 2011 appear suspect. There are very few (as opposed to >1000 in previous years) and are only 0 or 1. There are no entries in 2009.
	scot$SURFACE_TEMPERATURE[scot$year %in% c(2009, 2010, 2011)] = NA

	# Turn 0 values in GoMex to NA. These are outliers (way too cold) and must be mistakes.
	i = which(gmex$TEMP_SSURF == 0)
	gmex$TEMP_SSURF[i] = NA
	i = which(gmex$TEMP_BOT == 0)
	gmex$TEMP_BOT[i] = NA
	
	# 0 values in ai July and goa July are much lower than other values, seem suspect
	ai$SURF_TEMP[ai$month == 7 & ai$SURF_TEMP==0] = NA
	goa$SURF_TEMP[goa$month == 7 & goa$SURF_TEMP==0] = NA
	

	# SST histograms by month for each survey. look for outliers by month (especially 0s)
#	quartz(width=15, height=12)
#	par(mfcol=c(10,12), mai=c(0.3, 0.3, 0.1, 0.1)) # use mfcol to fill by columns
#	for(i in 1:12){ #for each month
#		if(sum(ai$month == i) >0) hist(ai$SURF_TEMP[ai$month == i], col='grey', main='ai') else plot(0,0, col='white', bty='n', xaxt='n', yaxt='n', main='')
#		if(sum(ebs$month == i) >0) hist(ebs$SURF_TEMP[ebs$month == i], col='grey', main='ebs') else plot(0,0, col='white', bty='n', xaxt='n', yaxt='n', main='')
#		if(sum(goa$month == i) >0) hist(goa$SURF_TEMP[goa$month == i], col='grey', main='goa') else plot(0,0, col='white', bty='n', xaxt='n', yaxt='n', main='')
#		if(sum(neus$month == i) >0) hist(neus$SURFTEMP[neus$month == i], col='grey', main='neus') else plot(0,0, col='white', bty='n', xaxt='n', yaxt='n', main='')
#		if(sum(wctri$month == i) >0) hist(wctri$SURFACE_TEMPERATURE[wctri$month == i], col='grey', main='wctri') else plot(0,0, col='white', bty='n', xaxt='n', yaxt='n', main='')
#		if(sum(wcann$month == i & !is.na(wcann$surftemp)) >0) hist(wcann$surftemp[wcann$month == i], col='grey', main='wcann') else plot(0,0, col='white', bty='n', xaxt='n', yaxt='n', main='')
#		if(sum(gmex$month == i) >0) hist(gmex$TEMP_SSURF[gmex$month == i], col='grey', main='gmex') else plot(0,0, col='white', bty='n', xaxt='n', yaxt='n', main='')
#		if(sum(newf$month == i & !is.na(newf$surftemp)) >0) hist(newf$surftemp[newf$month == i], col='grey', main='newf') else plot(0,0, col='white', bty='n', xaxt='n', yaxt='n', main='')
#		if(sum(scot$month == i, na.rm=TRUE) >0) hist(scot$SURFACE_TEMPERATURE[scot$month == i], col='grey', main='scot') else plot(0,0, col='white', bty='n', xaxt='n', yaxt='n', main='')
#		if(sum(sgsl$month == i) >0) hist(sgsl$t_surface[sgsl$month == i], col='grey', main='sgsl') else plot(0,0, col='white', bty='n', xaxt='n', yaxt='n', main='')
#	}

	# BT histograms by month for each survey. look for outliers by month (especially 0s)
#	quartz(width=15, height=12)
#	par(mfcol=c(10,12), mai=c(0.3, 0.3, 0.1, 0.1)) # use mfcol to fill by columns
#	for(i in 1:12){ #for each month
#		if(sum(ai$month == i) >0) hist(ai$BOT_TEMP[ai$month == i], col='grey', main='ai') else plot(0,0, col='white', bty='n', xaxt='n', yaxt='n', main='')
#		if(sum(ebs$month == i) >0) hist(ebs$BOT_TEMP[ebs$month == i], col='grey', main='ebs') else plot(0,0, col='white', bty='n', xaxt='n', yaxt='n', main='')
#		if(sum(goa$month == i) >0) hist(goa$BOT_TEMP[goa$month == i], col='grey', main='goa') else plot(0,0, col='white', bty='n', xaxt='n', yaxt='n', main='')
#		if(sum(neus$month == i) >0) hist(neus$BOTTEMP[neus$month == i], col='grey', main='neus') else plot(0,0, col='white', bty='n', xaxt='n', yaxt='n', main='')
#		if(sum(wctri$month == i) >0) hist(wctri$GEAR_TEMPERATURE[wctri$month == i], col='grey', main='wctri') else plot(0,0, col='white', bty='n', xaxt='n', yaxt='n', main='')
#		if(sum(wcann$month == i & !is.na(wcann$Temperature.At.the.Gear..degs.C.)) >0) hist(wcann$Temperature.At.the.Gear..degs.C.[wcann$month == i], col='grey', main='wcann') else plot(0,0, col='white', bty='n', xaxt='n', yaxt='n', main='')
#		if(sum(gmex$month == i) >0) hist(gmex$TEMP_BOT[gmex$month == i], col='grey', main='gmex') else plot(0,0, col='white', bty='n', xaxt='n', yaxt='n', main='')
#		if(sum(newf$month == i & !is.na(newf$bottemp)) >0) hist(newf$bottemp[newf$month == i], col='grey', main='newf') else plot(0,0, col='white', bty='n', xaxt='n', yaxt='n', main='')
#		if(sum(scot$month == i, na.rm=TRUE) >0) hist(scot$BOTTOM_TEMPERATURE[scot$month == i], col='grey', main='scot') else plot(0,0, col='white', bty='n', xaxt='n', yaxt='n', main='')
#		if(sum(sgsl$month == i) >0) hist(sgsl$t_bottom[sgsl$month == i], col='grey', main='sgsl') else plot(0,0, col='white', bty='n', xaxt='n', yaxt='n', main='')
#	}

	
# Trim out or fix speed and duration records, and other bad tows
	gmex = gmex[gmex$MIN_FISH<=60 & gmex$MIN_FISH > 0 & !is.na(gmex$MIN_FISH),] # trim out tows of 0, >60, or unknown minutes
	gmex$VESSEL_SPD[gmex$VESSEL_SPD==30] = 3 # fix typo according to Jeff Rester: 30 = 3	
	gmex = gmex[gmex$VESSEL_SPD < 5 & gmex$VESSEL_SPD > 0  & !is.na(gmex$VESSEL_SPD),] # trim out vessel speeds 0, unknown, or >5
	newf = newf[newf$operation %in% c(1,2) & newf$recordtype == 6,] # 6 is biological data, 5 is set information
	newf = newf[newf$duration<=60,]
	scot = scot[scot$TYPE==1,] # 1 is normal tows
	sgsl = sgsl[sgsl$expt %in% c(1,5),] # high quality tows: surveys and comparative tows
  seus = seus[!(seus$year == 1989 & seus$season == "spring"),] #They sampled at night this year-season

	# EFFORT fixes for the the seamap survey
  seus$EFFORT[seus$COLLECTIONNUMBER == 19900119] <- 2.78478
  seus$EFFORT[seus$COLLECTIONNUMBER == 19910105] <- 1.71273
	seus$EFFORT[seus$COLLECTIONNUMBER == 19950335] <- 0.9775
	seus$EFFORT[seus$COLLECTIONNUMBER == 19990065] <- 0.53648
	seus$EFFORT[seus$COLLECTIONNUMBER == 20070177] <- 0.99936
	seus$EFFORT[seus$COLLECTIONNUMBER == 20110393] <- 1.65726
	seus$EFFORT[seus$EVENTNAME == 1992219] <- 1.796247
  seus$EFFORT[seus$EVENTNAME == 1991423] <- 2.29
  seus$EFFORT[seus$EVENTNAME == 2001431] <- 3.18

# Add "strata" (define by lat, lon and depth bands) where needed
	stratlatgrid = floor(wctri$START_LATITUDE)+0.5 # degree bins
	stratdepthgrid = floor(wctri$BOTTOM_DEPTH/100)*100 + 50 # 100 m bins
	wctri$stratum = paste(stratlatgrid, stratdepthgrid, sep='-') # no need to use lon grids on west coast (so narrow)

	stratlatgrid = floor(wcann$Best.Latitude..dd.)+0.5 # degree bins
	stratdepthgrid = floor(wcann$Best.Depth..m./100)*100 + 50 # 100 m bins
	wcann$stratum = paste(stratlatgrid, stratdepthgrid, sep='-') # no need to use lon grids on west coast (so narrow)

	stratlatgrid = floor(gmex$lat)+0.5 # degree bins
	stratlongrid = floor(gmex$lon)+0.5 # degree bins
	stratdepthgrid = floor(gmex$depth/100)*100 + 50 # 100 m bins
	gmex$stratum = paste(stratlatgrid, stratlongrid, stratdepthgrid, sep='-')
	
	stratlatgrid = floor(seus.shelf$LATITUDESTART)+0.5 # degree bins
	stratlongrid = floor(seus.shelf$LONGITUDESTART)+0.5 # degree bins
	stratdepthgrid = floor(seus.shelf$DEPTHSTART/100)*100 + 50 # 100 m bins
	seus.shelf$stratum = paste(stratlatgrid, stratlongrid, stratdepthgrid, sep='-')
	rm(stratlatgrid, stratdepthgrid)

# Fix some lat/lon data entry errors in seus
	seus$LONGITUDESTART[seus$EVENTNAME == 1998467] <- -79.01
	seus$LONGITUDESTART[seus$EVENTNAME == 2010233] <- -81.006
	seus$LATITUDESTART[seus$EVENTNAME == 2014325] <- 34.616
  seus$LATITUDEEND[seus$EVENTNAME == 2007295] <- 30.992
  seus$LONGITUDEEND[seus$EVENTNAME == 2003403] <- -78.998
  seus$LONGITUDEEND[seus$EVENTNAME == 2011583] <- -78.997

# Correct more errors in seus with non-weighed species_mostly sea turtles
	seus$SPECIESTOTALWEIGHT <- replace(seus$SPECIESTOTALWEIGHT, seus$COLLECTIONNUMBER == 20010106 & seus$SPECIESCODE == 8713050104, 31.9) 
	seus$SPECIESTOTALWEIGHT[is.na(seus$SPECIESTOTALWEIGHT)] <- 0
  seus <- within(seus, SPECIESTOTALWEIGHT[SPECIESTOTALWEIGHT == 0 & SPECIESCODE == 5802010101] <- 1.9)
  seus$SPECIESTOTALWEIGHT <- replace(seus$SPECIESTOTALWEIGHT, seus$COLLECTIONNUMBER == 19940236 & seus$SPECIESCODE == 9002050101, 203.8) 
	seus$SPECIESTOTALWEIGHT <- replace(seus$SPECIESTOTALWEIGHT, seus$SPECIESTOTALWEIGHT == 0 & seus$SPECIESCODE == 9002040101, 46.99) 
  seus$SPECIESTOTALWEIGHT <- replace(seus$SPECIESTOTALWEIGHT, seus$COLLECTIONNUMBER == 20130188 & seus$SPECIESCODE == 9002040401, 12.77) 

# Create cpue column for seus; first must combine paired tows
	# Any columns with NA need to not get included in aggregate, and then readded, which is just the 2 temp columns.
	seus.temps <- unique(data.frame(haulid = seus$haulid, bottemp = seus$TEMPBOTTOM, surftemp = seus$TEMPSURFACE))
	seus <- aggregate(list(BIOMASS = seus$SPECIESTOTALWEIGHT), by=list(haulid = seus$haulid, stratum = seus$STRATA, stratumarea = seus$STRATAHECTARE, year = seus$year, month = seus$month, season = seus$season, 
	  lat = seus$LATITUDESTART, lon = seus$LONGITUDESTART, depth = seus$DEPTHSTART, EFFORT = seus$EFFORT, spp = seus$SPECIESSCIENTIFICNAME), FUN=sum)
	seus <- merge(x=seus, y=seus.temps, all.x=T, by="haulid")
	seus$wtcpue <- seus$BIOMASS/(seus$EFFORT*2)#yields biomass (kg) per hectare for each 'spp' and 'haulid'; EFFORT is multiplied by 2 b/c it is always identical for each of the paired tows
	rm(seus.temps)
	 
# Fix column names
	names(ai)[names(ai)=='YEAR'] = 'year'
	names(ai)[names(ai)=='LATITUDE'] = 'lat'
	names(ai)[names(ai)=='LONGITUDE'] = 'lon' 
	names(ai)[names(ai)=='STRATUM'] = 'stratum' 
	names(ai)[names(ai)=='BOT_DEPTH'] = 'depth'
	names(ai)[names(ai)=='BOT_TEMP'] = 'bottemp'
	names(ai)[names(ai)=='SURF_TEMP'] = 'surftemp'
	names(ai)[names(ai)=='SCIENTIFIC'] = 'spp'
	names(ai)[names(ai)=='WTCPUE'] = 'wtcpue'

	names(ebs)[names(ebs)=='YEAR'] = 'year'
	names(ebs)[names(ebs)=='LATITUDE'] = 'lat'
	names(ebs)[names(ebs)=='LONGITUDE'] = 'lon' # use the adjusted longitude
	names(ebs)[names(ebs)=='STRATUM'] = 'stratum' 
	names(ebs)[names(ebs)=='BOT_DEPTH'] = 'depth'
	names(ebs)[names(ebs)=='BOT_TEMP'] = 'bottemp'
	names(ebs)[names(ebs)=='SURF_TEMP'] = 'surftemp'
	names(ebs)[names(ebs)=='SCIENTIFIC'] = 'spp'
	names(ebs)[names(ebs)=='WTCPUE'] = 'wtcpue'

	names(goa)[names(goa)=='YEAR'] = 'year'
	names(goa)[names(goa)=='LATITUDE'] = 'lat'
	names(goa)[names(goa)=='LONGITUDE'] = 'lon'
	names(goa)[names(goa)=='STRATUM'] = 'stratum' 
	names(goa)[names(goa)=='BOT_DEPTH'] = 'depth'
	names(goa)[names(goa)=='BOT_TEMP'] = 'bottemp'
	names(goa)[names(goa)=='SURF_TEMP'] = 'surftemp'
	names(goa)[names(goa)=='SCIENTIFIC'] = 'spp'
	names(goa)[names(goa)=='WTCPUE'] = 'wtcpue'

	names(neus)[names(neus)=='YEAR'] = 'year'
	names(neus)[names(neus)=='SCINAME'] = 'spp'
	names(neus)[names(neus)=='LAT'] = 'lat'
	names(neus)[names(neus)=='LON'] = 'lon'
	names(neus)[names(neus)=='STRATUM'] = 'stratum'
	names(neus)[names(neus)=='DEPTH'] = 'depth'
	names(neus)[names(neus)=='BOTTEMP'] = 'bottemp'
	names(neus)[names(neus)=='SURFTEMP'] = 'surftemp'
	
# seus names were fixed in previous step that aggregated twin hauls

	names(wctri)[names(wctri)=='VESSEL'] = 'svvessel'
	names(wctri)[names(wctri) == 'START_LATITUDE'] = 'lat'
	names(wctri)[names(wctri) == 'START_LONGITUDE'] = 'lon'
	names(wctri)[names(wctri) == 'BOTTOM_DEPTH'] = 'depth'
	names(wctri)[names(wctri)=='GEAR_TEMPERATURE'] = 'bottemp'
	names(wctri)[names(wctri)=='SURFACE_TEMPERATURE'] = 'surftemp'
	names(wctri)[names(wctri) == 'SPECIES_NAME'] = 'spp'
	names(wctri)[names(wctri)=='WEIGHT'] = 'wtcpue'

	names(wcann)[names(wcann)=='Best.Latitude..dd.'] = 'lat'
	names(wcann)[names(wcann)=='Best.Longitude..dd.'] = 'lon'
	names(wcann)[names(wcann)=='Best.Depth..m.'] = 'depth'
	names(wcann)[names(wcann)=='Temperature.At.the.Gear..degs.C.'] = 'bottemp'
	names(wcann)[names(wcann)=='Species'] = 'spp'

	names(gmex)[names(gmex)=='TAXONOMIC'] = 'spp'
	names(gmex)[names(gmex)=='TEMP_BOT'] = 'bottemp'
	names(gmex)[names(gmex)=='TEMP_SSURF'] = 'surftemp'

	names(scot)[names(scot)=='MISSION'] = 'cruise'
	names(scot)[names(scot)=='STRAT'] = 'stratum'
	names(scot)[names(scot)=='SETNO'] = 'tow'
	names(scot)[names(scot)=='SURFACE_TEMPERATURE'] = 'surftemp'
	names(scot)[names(scot)=='BOTTOM_TEMPERATURE'] = 'bottemp'
	names(scot)[names(scot)=='BOTTOM_SALINITY'] = 'botsal'
	names(scot)[names(scot)=='ADJ_TOTWGT'] = 'wtcpue'
	names(scot)[names(scot)=='ADJ_TOTNO'] = 'numcpue'

	names(sgsl)[names(sgsl)=='vessel'] = 'svvessel'
	names(sgsl)[names(sgsl)=='strat'] = 'stratum'
	names(sgsl)[names(sgsl)=='set'] = 'tow'
	names(sgsl)[names(sgsl)=='latitude'] = 'lat'
	names(sgsl)[names(sgsl)=='longitude'] = 'lon' # use the adjusted longitude
	names(sgsl)[names(sgsl)=='t_surface'] = 'surftemp'
	names(sgsl)[names(sgsl)=='t_bottom'] = 'bottemp'
	names(sgsl)[names(sgsl)=='latin_name'] = 'spp'
	names(sgsl)[names(sgsl)=='biomass'] = 'wtcpue'
	names(sgsl)[names(sgsl)=='catch'] = 'numcpue'
	names(sgsl)[names(sgsl)=='salin_bottom'] = 'botsal'
	names(sgsl)[names(sgsl)=='name'] = 'common'

# Turn -9999 to NA where needed
	ai$wtcpue[ai$wtcpue==-9999] = NA
	ebs$wtcpue[ebs$wtcpue==-9999] = NA
	goa$wtcpue[goa$wtcpue==-9999] = NA

# Adjust for towed area where needed
	wctri$wtcpue = wctri$wtcpue*10000/(wctri$DISTANCE_FISHED*1000*wctri$NET_WIDTH) # weight per hectare (10,000 m2)	
	wcann$wtcpue = wcann$Haul.Weight..kg./wcann$Area.Swept.by.the.Net..hectares. # kg per hectare (10,000 m2)	
	gmex$wtcpue = gmex$SELECT_BGS /(gmex$VESSEL_SPD * 1.85200 * 1000 * gmex$MIN_FISH / 60 * gmex$GEAR_SIZE * 0.3048) # kg per m2. calc area trawled in m2: knots * 1.8 km/hr/knot * 1000 m/km * minutes * 1 hr/60 min * width of gear in feet * 0.3 m/ft # biomass per standard tow
	newf$towarea = newf$distance/10 * 1852 * 55.25 * 0.3048 # area trawled: distance in nm with one decimal * 1852 m/nm * 55.25 ft wide * 0.3048 m/ft
	newf$towarea[newf$towarea == 0] = NA
	newf$wtcpue = newf$wgt/100 /newf$towarea * 100 # in kg/100m2


# Remove a tow when paired tows exist (same lat/lon/year but different haulid, only Gulf of Mexico)
#	dups = which(duplicated(gmex[,c('year', 'lat', 'lon')]) & !duplicated(gmex$haulid)) # identify duplicate tows at same year/lat/lon
#	dupped = gmex[paste(gmex$year, gmex$lat, gmex$lon) %in% paste(gmex$year[dups], gmex$lat[dups], gmex$lon[dups]),] # all tows at these year/lat/lon
#		# sum(!duplicated(dupped$haulid)) # 26 (13 pairs of haulids)
#	gmex = gmex[!(gmex$haulid %in% unique(dupped$haulid[grep('PORT', dupped$COMSTAT)])),] # remove the port haul (this is arbitrary, but seems to be right based on the notes associated with these hauls)

# Removes rows without scientific names or with unreliable IDs
	ai = ai[ai$spp != '',]
	ebs = ebs[ebs$spp != '',]
	goa = goa[goa$spp != '',]
	neus = neus[!(neus$spp == '' | is.na(neus$spp)),]
	neus = neus[!(neus$spp %in% c('UNIDENTIFIED FISH', 'ILLEX ILLECEBROSUS EGG MOPS', 'LOLIGO PEALEII EGG MOPS')),] # remove unidentified spp
	wctri = wctri[wctri$spp != '',]
	wcann = wcann[wcann$spp != '',]
	gmex = gmex[!(gmex$spp == '' | is.na(gmex$spp)),]
	gmex = gmex[!(gmex$spp %in% c('UNID CRUSTA', 'UNID OTHER', 'UNID.FISH', 'CRUSTACEA(INFRAORDER) BRACHYURA', 'MOLLUSCA AND UNID.OTHER #01', 'ALGAE', 'MISCELLANEOUS INVERTEBR', 'OTHER INVERTEBRATES')),] # remove unidentified spp
	newf = newf[!(newf$spp == '' | is.na(newf$spp)),]
	newf = newf[!(newf$spp %in% c('EGGS, FISH(SPAWN)', 'EGGS, INVERTEBRATE', 'EGGS, SKATE CASES', 'EGGS, UNIDENTIFIED', 'OFFAL, OTHER', 'PLANT MATERIAL', 'SHELLS', 'STONE', 'UNIDENTIFIED FISH', 'UNIDENTIFIED MATERIAL')),]
	scot = scot[!(scot$spp %in% c('FISH EGGS-UNIDENTIFIED', 'SOFT CORAL UNIDENTIFIED', 'UNID REMAINS,DIGESTED', 'UNID FISH AND INVERTEBRATES', 'UNID. FISH (LARVAE,JUVENILE AND ADULTS)', 'UNIDENTIFIED', 'UNIDENTIFIED A', 'UNIDENTIFIED B', 'UNIDENTIFIED C', 'UNIDENTIFIED D', 'UNIDENTIFIED E', 'RAJA EGGS', 'SHARK (NS)', 'BUCCINIDAE EGGS', 'EGGS UNID.', 'FINFISHES (NS)', 'FOREIGN ARTICLES,GARBAGE', 'HEMITRIPTERUS AMERICANUS, EGGS', 'MARINE INVERTEBRATA (NS)', 'MYOXOCEPHALUS EGGS', 'POLYCHAETA C.,SMALL', 'PURSE LITTLE SKATE', 'PURSE WINTER SKATE', 'COTTIDAE F. UNID.')),]
	seus = seus[!(seus$spp %in% c('MISCELLANEOUS INVERTEBRATES', 'XANTHIDAE', 'MICROPANOPE NUTTINGI', 'MICROPANOPE SCULPTIPES', 'GLYPTOXANTHUS EROSUS', 'PSEUDOMEDAEUS AGASSIZII', 'ALGAE', 'DYSPANOPEUS SAYI')),]
	# nothing to remove in Southern Gulf of St. Lawrence or seus.shelf
 
# Adjust spp names for those cases where they've changed or where matching failed (GMex)
	# first convert factors to strings so that we can modify them
	i <- sapply(ai, is.factor); ai[i] <- lapply(ai[i], as.character)
	i <- sapply(ebs, is.factor); ebs[i] <- lapply(ebs[i], as.character)
	i <- sapply(goa, is.factor); goa[i] <- lapply(goa[i], as.character)
	i <- sapply(neus, is.factor); neus[i] <- lapply(neus[i], as.character)
	i <- sapply(wctri, is.factor); wctri[i] <- lapply(wctri[i], as.character)
	i <- sapply(wcann, is.factor); wcann[i] <- lapply(wcann[i], as.character)
	i <- sapply(gmex, is.factor); gmex[i] <- lapply(gmex[i], as.character)
	i <- sapply(newf, is.factor); newf[i] <- lapply(newf[i], as.character)


	ai$spp[ai$spp %in% c('Atheresthesevermanni', 'Atheresthesstomias')] = 'Atheresthessp.'
	ai$spp[ai$spp %in% c('Lepidopsettapolyxystra', 'Lepidopsettabilineata')] = 'Lepidopsettasp.'
	ai$spp[ai$spp %in% c('Myoxocephalusjaok', 'Myoxocephalusniger', 'Myoxocephaluspolyacanthocephalus', 'Myoxocephalusquadricornis', 'Myoxocephalusverrucosus')] = 'Myoxocephalussp.'
	ai$spp[ai$spp %in% c('Bathyrajaabyssicola', 'Bathyrajaaleutica', 'Bathyrajainterrupta', 'Bathyrajalindbergi', 'Bathyrajamaculata', 'Bathyrajamariposa', 'Bathyrajaminispinosa', 'Bathyrajaparmifera', 'Bathyrajasmirnovi', 'Bathyrajasp.cf.parmifera(Orretal.)', 'Bathyrajaspinosissima', 'Bathyrajataranetzi', 'Bathyrajatrachura', 'Bathyrajaviolacea')] = 'Bathyrajasp.'

	ebs$spp[ebs$spp %in% c('Atheresthes evermanni', 'Atheresthes stomias')] = 'Atheresthes sp.'
	ebs$spp[ebs$spp %in% c('Lepidopsetta polyxystra', 'Lepidopsetta bilineata')] = 'Lepidopsetta sp.'
	ebs$spp[ebs$spp %in% c('Hippoglossoides elassodon', 'Hippoglossoides robustus')] = 'Hippoglossoides sp.'
	ebs$spp[ebs$spp %in% c('Myoxocephalus jaok', 'Myoxocephalus niger', 'Myoxocephalus polyacanthocephalus', 'Myoxocephalus quadricornis', 'Myoxocephalus verrucosus', 'Myoxocephalus scorpioides')] = 'Myoxocephalus sp.'
	ebs$spp[ebs$spp %in% c('Bathyraja abyssicola', 'Bathyraja aleutica', 'Bathyraja interrupta', 'Bathyraja lindbergi', 'Bathyraja maculata', 'Bathyraja mariposa', 'Bathyraja minispinosa', 'Bathyraja parmifera', 'Bathyraja smirnovi', 'Bathyraja sp.', 'Bathyraja sp.cf.parmifera(Orretal.)', 'Bathyraja spinosissima', 'Bathyraja taranetzi', 'Bathyraja trachura', 'Bathyraja violacea')] = 'Bathyraja sp.'

	goa$spp[goa$spp %in% c('Lepidopsettapolyxystra', 'Lepidopsettabilineata')] = 'Lepidopsettasp.'
	goa$spp[goa$spp %in% c('Myoxocephalusjaok', 'Myoxocephalusniger', 'Myoxocephaluspolyacanthocephalus', 'Myoxocephalusquadricornis', 'Myoxocephalusverrucosus')] = 'Myoxocephalussp.'
	goa$spp[goa$spp %in% c('Bathyrajaabyssicola', 'Bathyrajaaleutica', 'Bathyrajainterrupta', 'Bathyrajalindbergi', 'Bathyrajamaculata', 'Bathyrajamariposa', 'Bathyrajaminispinosa', 'Bathyrajaparmifera', 'Bathyrajasmirnovi', 'Bathyrajasp.cf.parmifera(Orretal.)', 'Bathyrajaspinosissima', 'Bathyrajataranetzi', 'Bathyrajatrachura', 'Bathyrajaviolacea')] = 'Bathyrajasp.'

	wctri$spp[wctri$spp %in% c('Lepidopsetta polyxystra', 'Lepidopsetta bilineata')] = 'Lepidopsetta sp.'
	wctri$spp[wctri$spp %in% c('Bathyraja interrupta', 'Bathyraja trachura', 'Bathyraja parmifera', 'Bathyraja spinosissima')] = 'Bathyrajasp.'

	wcann$spp[wcann$spp %in% c('Lepidopsetta polyxystra', 'Lepidopsetta bilineata')] = 'Lepidopsetta sp.' # so that species match wctri
	wcann$spp[wcann$spp %in% c('Bathyraja abyssicola', 'Bathyraja aleutica', 'Bathyraja kincaidii (formerly B. interrupta)', 'Bathyraja sp. ', 'Bathyraja trachura', 'Bathyraja parmifera', 'Bathyraja spinosissima')] = 'Bathyrajasp.'

	seus$spp[seus$spp %in% c('ANCHOA HEPSETUS', 'ANCHOA MITCHILLI', 'ANCHOA CUBANA', 'ANCHOA LYOLEPIS', 'ENGRAULIS EURYSTOLE')] = 'ANCHOA'
	seus$spp[seus$spp %in% c('LIBINIA DUBIA', 'LIBINIA EMARGINATA')] = 'LIBINIA'
	

	i = gmex$GENUS_BGS == 'PELAGIA' & gmex$SPEC_BGS == 'NOCTUL'; gmex$spp[i] = 'PELAGIA NOCTILUCA'; gmex$BIO_BGS[i] = 618030201
	i = gmex$GENUS_BGS == 'MURICAN' & gmex$SPEC_BGS == 'FULVEN'; gmex$spp[i] = 'MURICANTHUS FULVESCENS'; gmex$BIO_BGS[i] = 308011501
	i = gmex$spp %in% c('APLYSIA BRASILIANA', 'APLYSIA WILLCOXI'); gmex$spp[i] = 'APLYSIA'
	i = gmex$spp %in% c('AURELIA AURITA'); gmex$spp[i] = 'AURELIA'
	i = gmex$spp %in% c('BOTHUS LUNATUS', 'BOTHUS OCELLATUS', 'BOTHUS ROBINSI'); gmex$spp[i] = 'BOTHUS'
	i = gmex$spp %in% c('CLYPEASTER PROSTRATUS', 'CLYPEASTER RAVENELII'); gmex$spp[i] = 'CLYPEASTER'
	i = gmex$spp %in% c('CONUS AUSTINI', 'CONUS STIMPSONI'); gmex$spp[i] = 'CONUS'
	i = gmex$spp %in% c('CYNOSCION ARENARIUS', 'CYNOSCION NEBULOSUS', 'CYNOSCION NOTHUS'); gmex$spp[i] = 'CYNOSCION'
	i = gmex$spp %in% c('ECHINASTER SENTUS', 'ECHINASTER SERPENTARIUS'); gmex$spp[i] = 'ECHINASTER'
	i = gmex$spp %in% c('ECHINASTER SENTUS', 'ECHINASTER SERPENTARIUS'); gmex$spp[i] = 'ECHINASTER'
	i = gmex$spp %in% c('OPISTOGNATHUS AURIFRONS', 'OPISTOGNATHUS LONCHURUS'); gmex$spp[i] = 'OPISTOGNATHUS'
	i = gmex$spp %in% c('OPSANUS BETA', 'OPSANUS PARDUS', 'OPSANUS TAU'); gmex$spp[i] = 'OPSANUS'
	i = gmex$spp %in% c('ROSSIA BULLISI'); gmex$spp[i] = 'ROSSIA'
	i = gmex$spp %in% c('SOLENOCERA ATLANTIDIS', 'SOLENOCERA NECOPINA', 'SOLENOCERA VIOSCAI'); gmex$spp[i] = 'SOLENOCERA'
	i = gmex$spp %in% c('TRACHYPENEUS CONSTRICTUS', 'TRACHYPENEUS SIMILIS'); gmex$spp[i] = 'TRACHYPENEUS'

	i = newf$spp %in% c('ARTEDIELLUS ATLANTICUS', 'ARTEDIELLUS UNCINATUS'); newf$spp[i] = 'ARTEDIELLUS  SP.'
	i = newf$spp %in% c('BUCCINUM  SP.', 'BUCCINUM TOTTENI', 'BUCCINUM UNDATUM'); newf$spp[i] = 'BUCCINIDAE' # didn't search for all included genera (many!)
	i = newf$spp %in% c('CHIONOECETES OPILIO FEMALE', 'CHIONOECETES OPILIO MALE'); newf$spp[i] = 'CHIONOECETES OPILIO'
	i = newf$spp %in% c('EUALUS GAIMARDII BELCHERI', 'EUALUS GAIMARDII GAIMARDII'); newf$spp[i] = 'EUALUS GAIMARDII'
	i = newf$spp %in% c('EUMICROTREMUS SPINOSUS VARIABILIS'); newf$spp[i] = 'EUMICROTREMUS SPINOSUS'
	i = newf$spp %in% c('GAIDROPSARUS ARGENTATUS', 'GAIDROPSARUS ENSIS'); newf$spp[i] = 'GAIDROPSARUS  SP.'
	i = newf$spp %in% c('GONATUS FABRICII'); newf$spp[i] = 'GONATUS  SP.'
	i = newf$spp %in% c('GORGONOCEPHALUS ARCTICUS', 'GORGONOCEPHALUS SP.'); newf$spp[i] = 'GORGONOCEPHALIDAE'
	i = newf$spp %in% c('HYAS ARANEUS', 'HYAS COARCTATUS'); newf$spp[i] = 'HYAS  SP.'
	i = newf$spp %in% c('LIPARIS ATLANTICUS', 'LIPARIS FABRICII', 'LIPARIS GIBBUS', 'LIPARIS LIPARIS', 'LIPARIS TUNICATUS'); newf$spp[i] = 'LIPARIDAE'
	i = newf$spp %in% c('LITHODES  SP.', 'LITHODES  SP.', 'NEOLITHODES  SP.', 'NEOLITHODES GRIMALDII'); newf$spp[i] = 'LITHODIDAE'
	i = newf$spp %in% c('LYCENCHELYS PAXILLUS', 'LYCENCHELYS SARSI', 'LYCENCHELYS VERRILLI'); newf$spp[i] = 'LYCENCHELYS  SP.'
	i = newf$spp %in% c('NOTACANTHUS NASUS'); newf$spp[i] = 'NOTACANTHIDAE'
	i = newf$spp %in% c('PANDALUS BOREALIS(FEE 1ST W/HR)', 'PANDALUS BOREALIS(FEE 1ST W/O HR)', 'PANDALUS BOREALIS(FEE 1ST)', 'PANDALUS BOREALIS(FEE OVIG.)', 'PANDALUS BOREALIS(FEE(1+) W/HR)', 'PANDALUS BOREALIS(FEE(1+) W/O HR)', 'PANDALUS BOREALIS(MALE)', 'PANDALUS BOREALIS(TRANS. W/HR)', 'PANDALUS BOREALIS(TRANS. W/O HR)'); newf$spp[i] = 'PANDALUS BOREALIS'
	i = newf$spp %in% c('PANDALUS MONTAGUI(FEE 1ST W/HR)', 'PANDALUS MONTAGUI(FEE 1ST W/O HR)', 'PANDALUS MONTAGUI(FEE OVIG.)', 'PANDALUS MONTAGUI(FEE(1+) W/HR)', 'PANDALUS MONTAGUI(FEE(1+) W/O HR)', 'PANDALUS MONTAGUI(MALE)', 'PANDALUS MONTAGUI(TRANS. W/HR)', 'PANDALUS MONTAGUI(TRANS.W/O HR))'); newf$spp[i] = 'PANDALUS MONTAGUI'
	i = newf$spp %in% c('PARALEPIS  SP.', 'PARALEPIS BREVIS (ATLANTICA)', 'PARALEPIS COREGONOIDES BOREALIS', 'ANOTOPTERIDAE', 'ANOTOPTERUS PHARAO', 'NOTOLEPIS RISSOI KROYERI'); newf$spp[i] = 'PARALEPIDIDAE' # most observations seem to be at family level # see http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=162532: Anotopteridae should be Paralepididae
	i = newf$spp %in% c('PASIPHAEA MULTIDENTATA', 'PASIPHAEA TARDA'); newf$spp[i] = 'PASIPHAEA  SP.'
	i = newf$spp %in% c('SCOPELOSAURUS  SP.'); newf$spp[i] = 'SCOPELOSAURIDAE' # according to ITIS, this should be 	Scopelosauridae
	i = newf$spp %in% c('SIMENCHELYS PARASITICUS'); newf$spp[i] = 'SIMENCHELYIDAE'
	i = newf$spp %in% c('TRIGLOPS MURRAYI', 'TRIGLOPS NYBELINI', 'TRIGLOPS PINGELI'); newf$spp[i] = 'TRIGLOPS  SP.'

	# nothing for Scotian Shelf
	# nothing for Southern Gulf of St. Lawrence

# Combine entries from spp that have now been consolidated
	ai2 = aggregate(list(wtcpue = ai$wtcpue), by = list(haulid = ai$haulid, year = ai$year, lat = ai$lat, lon = ai$lon, spp = ai$spp), FUN=sumna) # sum across entries with same haulid
		ai2 = merge(ai2, ai[!duplicated(ai$haulid) ,c('haulid', 'month', 'depth', 'stratum', 'bottemp', 'surftemp')], all.x=TRUE) # add depth, month, temperature

	ebs2 = aggregate(list(wtcpue = ebs$wtcpue), by = list(haulid = ebs$haulid, year = ebs$year, lat = ebs$lat, lon = ebs$lon, spp = ebs$spp), FUN=sumna) # sum across entries with same haulid
		ebs2 = merge(ebs2, ebs[!duplicated(ebs$haulid) ,c('haulid', 'month', 'depth', 'stratum', 'bottemp', 'surftemp')], all.x=TRUE) # add depth, month, temperature

	goa2 = aggregate(list(wtcpue = goa$wtcpue), by = list(haulid = goa$haulid, year = goa$year, lat = goa$lat, lon = goa$lon, spp = goa$spp), FUN=sumna) # use addNA() for depth so that NA values are not dropped by aggregate()
		goa2 = merge(goa2, goa[!duplicated(goa$haulid) ,c('haulid', 'month', 'depth', 'stratum', 'bottemp', 'surftemp')], all.x=TRUE) # add depth, month, temperature

	neus2 = aggregate(list(wtcpue = neus$wtcpue), by = list(haulid = neus$haulid, season = neus$SEASON, year = neus$year, lat = neus$lat, lon = neus$lon, spp = neus$spp), FUN=sumna)
		neus2 = merge(neus2, neus[!duplicated(neus$haulid) ,c('haulid', 'month', 'depth', 'stratum', 'bottemp', 'surftemp')], all.x=TRUE) # add depth, month, temperature

	seus2 = aggregate(list(wtcpue = seus$wtcpue), by = list(haulid = seus$haulid, season = seus$season, year = seus$year, lat = seus$lat, lon = seus$lon, spp = seus$spp), FUN=sumna)
		seus2 = merge(seus2, seus[!duplicated(seus$haulid) ,c('haulid', 'month', 'depth', 'stratum', 'bottemp', 'surftemp')], all.x=TRUE) # add depth, month, temperature

	wctri2 = aggregate(list(wtcpue = wctri$wtcpue), by = list(haulid = wctri$haulid, year = wctri$year, lat = wctri$lat, lon = wctri$lon, spp = wctri$spp), FUN=sumna)
		wctri2 = merge(wctri2, wctri[!duplicated(wctri$haulid) ,c('haulid', 'month', 'depth', 'stratum', 'bottemp', 'surftemp')], all.x=TRUE) # add depth, month, temperature

	wcann2 = aggregate(list(wtcpue = wcann$wtcpue), by = list(haulid = wcann$haulid, year = wcann$year, lat = wcann$lat, lon = wcann$lon, spp = wcann$spp), FUN=sumna)
		wcann2 = merge(wcann2, wcann[!duplicated(wcann$haulid) ,c('haulid', 'month', 'depth', 'stratum', 'bottemp', 'surftemp')], all.x=TRUE) # add depth, month, temperature

	gmex2 = aggregate(list(wtcpue = gmex$wtcpue), by=list(haulid = gmex$haulid, year = gmex$year, lat = gmex$lat, lon = gmex$lon, spp = gmex$spp), FUN=sumna)
		gmex2 = merge(gmex2, gmex[!duplicated(gmex$haulid) ,c('haulid', 'month', 'depth', 'stratum', 'bottemp', 'surftemp')], all.x=TRUE) # add depth, month, temperature

	newf2 = aggregate(list(wtcpue = newf$wtcpue), by=list(haulid = newf$haulid, vessel = newf$vessel, trip = newf$trip, year = newf$year, lat = newf$lat, lon = newf$lon, depth = newf$depth, spp = newf$spp), FUN=sumna)
		newf2 = merge(newf2, newf[!duplicated(newf$haulid) ,c('haulid', 'month', 'depth', 'stratum', 'bottemp', 'surftemp')], all.x=TRUE) # add depth, month, temperature
	
	# nothing for Scotian Shelf
	# nothing for Southern Gulf of St. Lawrence

# Calculate a corrected longitude for Aleutians (all in western hemisphere coordinates)
	ai2$lon[ai2$lon>0] = ai2$lon[ai2$lon>0] - 360	

# Split spring and fall surveys (NEUS, SEUS, and Newfoundland)
	# NEUS
	neusfal2 = neus2[neus2$season == "FALL",]
	neusspr2 = neus2[neus2$season == "SPRING",]
	# SEUS 
	seusspr2 = seus2[seus2$season == "spring",]
	seussum2 = seus2[seus2$season == "summer",]
	seusfal2 = seus2[seus2$season == "fall",]
	
	# Newfoundland
	surveys = read.csv('../../Princeton/Trawl Data/DFO Newfoundland/Tables/surveys_table.csv') # split spring from fall
	surveys2 = read.csv('../../Princeton/Trawl Data/DFO Newfoundland/Tables/surveys_table2009-2011.csv')
	fallseries = c(as.character(surveys$CRUISE[surveys$Series %in% c('2GH - Stratified Random Bottom Trawl Survey - Campelen 1800', 'Fall - Stratified Random Bottom Trawl Survey - Campelen 1800')]), as.character(surveys2$cruise[surveys2$season=='fall']))
	springseries = c(as.character(surveys$CRUISE[surveys$Series %in% c('Annual 3P - Stratified Random Bottom Trawl Survey - Campelen 1800', 'Spring 3LNO - Stratified Random Bottom Trawl Survey - Campelen 1800')]), as.character(surveys2$cruise[surveys2$season == 'spring']))
	cruiseid = paste(newf2$vessel, formatC(newf2$trip, width=3, flag=0), sep='') # CRUISE is concatenated vessel and trip
	newffal2 = newf2[cruiseid %in% fallseries,] # also take any fall tows in later years: need to ask Bill Brodie which cruiseids are actually correct!
	newfspr2 = newf2[cruiseid %in% springseries,] # also take any spring tows in later years: need to ask Bill Brodie which cruiseids are actually correct!


# Add a survey year that turns over July 1 (for fall surveys that cross Dec. 31)
ai2$yearsurv = ai2$year
ebs2$yearsurv = ebs2$year
goa2$yearsurv = goa2$year
neusfal2$yearsurv = neusfal2$year
neusspr2$yearsurv = neusspr2$year
wctri2$yearsurv = wctri2$year
wcann2$yearsurv = wcann2$year
gmex2$yearsurv = gmex2$year
newffal2$yearsurv = newffal2$year
	newffal2$yearsurv[newffal2$month<4] = newffal2$yearsurv[newffal2$month<4] - 1
newfspr2$yearsurv = newfspr2$year
scot$yearsurv = scot$year
sgsl$yearsurv = sgsl$year
seusspr2$yearsurv = seusspr2$year
seussum2$yearsurv = seussum2$year
seusfal2$yearsurv = seusfal2$year

# Add a region column
ai2$region = "AFSC_Aleutians"
ebs2$region = "AFSC_EBS"
goa2$region = "AFSC_GOA"
neusfal2$region = "NEFSC_NEUSFall"
neusspr2$region = "NEFSC_NEUSSpring"
wctri2$region = "AFSC_WCTri"
wcann2$region = "NWFSC_WCAnn"
gmex2$region = "SEFSC_GOMex"
newffal2$region = "DFO_NewfoundlandFall"
newfspr2$region = "DFO_NewfoundlandSpring"
scot$region = "DFO_ScotianShelf"
sgsl$region = "DFO_SoGulf"
seusspr2$region = "SCDNR_SEUSSpring"
seussum2$region = "SCDNR_SEUSSummer"
seusfal2$region = "SCDNR_SEUSFall"

# Rearrange and trim columns
nm = c('region', 'haulid', 'year', 'yearsurv', 'month', 'stratum', 'lat', 'lon', 'depth', 'surftemp', 'bottemp', 'spp', 'wtcpue')
ai2 = ai2[,nm]
ebs2 = ebs2[,nm]
goa2 = goa2[,nm]
neusfal2 = neusfal2[,nm]	
neusspr2 = neusspr2[,nm]	
wctri2 = wctri2[,nm]
wcann2 = wcann2[,nm]
gmex2 = gmex2[,nm]
newffal2 = newffal2[,nm]
newfspr2 = newfspr2[,nm]
scot2 = scot[,nm]
sgsl2 = sgsl[,nm]
seusspr2 = seusspr2[,nm]
seussum2 = seussum2[,nm]
seusfal2 = seusfal2[,nm]

# combine together
dat = rbind(ai2, ebs2, goa2, neusfal2, neusspr2, wctri2, wcann2, gmex2, newffal2, newfspr2, scot2, sgsl2) #, seusspr2, seussum2, seusfal2)
 
# examine
# summary(dat)

# Write out
save(dat, file=paste('Output/trawl_allregionsforprojections_', Sys.Date(), '.RData', sep='')) # 24.3MB as binary

setwd("/Users/abigailporay/Desktop")
dat.seus <- rbind(seusspr2, seussum2, seusfal2) # I don't rbind this with 'dat' yet b/c it may disrupt script 1.5
save(dat.seus, seusspr2, seussum2, seusfal2, file="SEUS_output_Sep23_2016.Rdata")
 