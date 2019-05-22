# Combine all the trawl survey data with the ‘trawlData’ package. 

# To install trawlData package_perhaps just for the first time:
# FIRST SET DIRECTORY TO trawlData
#library(devtools)
#devtools::install

library(trawlData)
library(ggplot2)
library(stringr)
library(zoo)

# useful function: acts like sum(na.rm=T) but returns NA if all are NA
sumna = function(x){
  if(!all(is.na(x))) return(sum(x, na.rm=T))
  if(all(is.na(x))) return(NA)
}

setwd("/Users/jim/Documents/proj_ranges")
load("Jim/trawl_allregionsforprojections_wSEUS_wSST_2016-09-26.RData") # Load 'dat' file from original data compile; to compare with trawlData method

# READ IN DATA for U.S. surveys from the trawlData package. Canadian surveys are not using trawlData, but rather using the file produced from the script 1_data_allcombine.r
ai <- copy(clean.ai)# Aleutians
ebs <- copy(clean.ebs)# Eastern Bering Sea
goa <- copy(clean.goa)# Gulf of Alaska

# Northeast US_requires some more steps and some species names get fixed here as an aggregation for 'sex' is necessary
neus <- copy(clean.neus)# Northeast US
setkey(neus, cruise, stratum, station, haulid, sex, SID)
neus <- unique(neus) # drops length data based on neus key
neus[, c('length', 'NUMLEN') := NULL] # remove length columns

# Some code to identify potential naming errors_where a raw data species name is changed to the wrong thing
abc <- unique(data.frame(ref=neus$ref, spp=neus$spp, SID=neus$SID, stringsAsFactors = F)) #list of species code combos for some QAQC
abc[is.na(abc$spp),] # all NAs for 'spp' are non species and not important
cde <- data.frame(table(abc$spp)); cde <- cde[cde$Freq>1,] # spp which have more than one ref or SID
abc[abc$spp=="Placopecten magellanicus",]
#Fixes of species names
neus[ref == "ISURUS OXYRINCHUS", spp := "Isurus oxyrinchus"]
neus[ref == "CARCHARHINUS SIGNATUS", spp := "Carcharhinus signatus"]
neus[ref == "CAULOLATILUS CYANOPS", spp := "Caulolatilus cyanops"]
neus[ref == "OCYURUS CHRYSURUS", spp := "Ocyurus chrysurus"]
neus[ref == "CITHARICHTHYS GYMNORHINUS", spp := "Citharichthys gymnorhinus"]
neus[ref == "SYNODUS INTERMEDIUS", spp := "Synodus intermedius"]
neus[ref == "HOLACANTHUS CILIARIS", spp := "Holacanthus ciliaris"]
neus[ref == "MYCTEROPERCA BONACI", spp := "Mycteroperca bonaci"]
rm(abc, cde)

neus <- neus[,sum(weight),by=list(year, haulid, datetime, season, lat, lon, depth, stemp, btemp, cruise, station, stratum, spp, ref)] # sum different sexes of same spp together
setnames(neus, 'V1', 'wtcpue') # Here wtcpue is just the weight in the catch, so assumes the trawl is an effort of '1' every time.

# West Coast Trienniel (1977-2004)
wctri <- copy(clean.wctri)
wctri <- wctri[HAUL_TYPE==3 & PERFORMANCE==0] # trim to standard hauls and good performance
# West Coast annual (2003-2012)
wcann <- copy(clean.wcann)

# Gulf of Mexico
gmex <- copy(clean.gmex)
# Can combine this code later on, I was going through step by step to learn the survey better
gmex <- gmex[geartype=="ST"]# only use Shrimp Trawls
gmex <- gmex[gearsize==40]# only use the 40' trawl, drop some other random nets, but mostly drops a bunch of tows from Texas, which uses a 20' net
gmex <- gmex[meshsize==1.63]# Most common used mesh size by far
gmex <- gmex[towduration > 9 & towduration < 61] # gmex tow durations are all over the place, this drops the handful of major outliers
gmex <- gmex[!is.na(towduration)] # remove hauls with no towduration info
# Reduce to summer and fall surveys and then conform some name as they changed the name and mispelled it etc. over time
gmex <- gmex[!(survey.name=="Comparative Tow" | survey.name=="Fall SEAMAP Plankton Survey" | survey.name=="SEAMAP Comparitive Tow" | survey.name=="SEAMAP Comparitve Tow" | survey.name=="Winter SEAMAP Groundfish Survey" | survey.name=="Winter SEAMAP Trawl Survey" | survey.name=="Spring SEAMAP Trawl Survey" | survey.name=="Spring SEAMAP Groundfish Survey" | survey.name=="Spring Groundfish Survey" | survey.name=="Spring Groundfish/Plankton Survey")]
gmex[survey.name == "Fall SEAMAP Groundfish Suvey", survey.name := "Fall SEAMAP Groundfish Survey"]
gmex[survey.name == "SEAMAP Fall Groundfish Survey", survey.name := "Fall SEAMAP Groundfish Survey"]
gmex[survey.name == "Fall Groundfish Survey", survey.name := "Fall SEAMAP Groundfish Survey"]
gmex[survey.name == "SEAMAP Summer Groundfish Survey", survey.name := "Summer SEAMAP Groundfish Survey"]
gmex <- gmex[OP == ""] # Trim to high quality hauls
gmex$haulid = paste(formatC(gmex$vessel, width=3, flag=0), formatC(gmex$CRUISE_NO, width=3, flag=0), formatC(gmex$P_STA_NO, width=5, flag=0, format='d'), sep='-')# Make a haulid (the one in trawlData is different in the 'cruise' column used)

# Southeast US_here I upload my own file b/c trawlData is not up to date for seus
survcatch = read.csv('/Users/jim/Documents/Work/Rutgers/Data/seamap_catch_2015.csv', stringsAsFactors=FALSE)
survhaul = read.csv('/Users/jim/Documents/Work/Rutgers/Data/seamap_event_2015.csv', stringsAsFactors=FALSE) # only needed for depth
survhaul = unique(data.frame(EVENTNAME = survhaul$EVENTNAME, DEPTHSTART = survhaul$DEPTHSTART))
seus = merge(x=survcatch, y=survhaul, all.x=T, by="EVENTNAME") # Add depth data from survhaul 
seus = cbind(seus, STRATA = as.integer(str_sub(string = seus$STATIONCODE, start = 1, end = 2))) #Create STRATA column
seus = seus[seus$DEPTHZONE != "OUTER",] # Drop OUTER depth zone because it was only sampled for 10 years, and in specific seasons-areas
rm(survcatch, survhaul)
seus$haulid = seus$EVENTNAME
#extract year and month and season
seus$date = as.Date(seus$DATE, "%m/%d/%y")
seus <- cbind(seus, year = year(seus$date), month = month(seus$date)) # also made month column for seus in this step
seus$season <- as.yearqtr(seus$date) # make a 'season' column to distinguish the spring, summer, and fall surveys
seus$season <- factor(format(seus$season, "%q"), levels = 1:4, labels = c("winter", "spring", "summer", "fall")) # takes ~30 sec
seus$season[seus$month == 9] <- "fall" 	#Sept EVENTS (all late-sept.) were grouped with summer, should be fall

# Northeast inshore (NEAMAP)
neamap.catch <- data.table(read.csv("Jim/NEAMAP_Catch.csv", stringsAsFactors=F))
neamap.spp <- data.table(read.csv("Jim/NEAMAP_sppList.csv", stringsAsFactors=F))
neamap <- merge(neamap.catch, neamap.spp, all.x=T, by="VIMSCODE", sort=F)
neamap.hauls <- data.table(read.csv("Jim/NEAMAP_Station.csv", stringsAsFactors=F))
neamap <- merge(neamap, neamap.hauls, all.x=T, by="STATION", sort=F)
rm(neamap.catch, neamap.hauls, neamap.spp)
#Extract date and season
neamap$year <- as.integer(str_sub(string = neamap$DATE, start=1, end=4))
neamap$month <- as.integer(str_sub(string = neamap$DATE, start=6, end=7))
neamap$season <- ifelse(neamap$month < 6, "spring", "fall")
neamap$haulid <- neamap$STATION
# Bring in water quality data_this file has water column CTD casts, so need to extract surface and bottom temperature measures; first into separate files
neamap.wq <- data.table(read.csv("Jim/neamap_waterQual.csv", stringsAsFactors=F))
neamap.st <- neamap.wq[PARAMETER ==  "WT" & LAYER == "S"]
neamap.st <- data.table(STATION = neamap.st$STATION, surftemp = neamap.st$VALUE)
neamap.bt <- neamap.wq[PARAMETER ==  "WT" & LAYER == "B"]
neamap.bt <- data.table(STATION = neamap.bt$STATION, bottemp = neamap.bt$VALUE)
neamap <- merge(neamap, neamap.bt, all.x=T, by="STATION", sort=F)
neamap <- merge(neamap, neamap.st, all.x=T, by="STATION", sort=F)
rm(neamap.wq, neamap.bt, neamap.st)

# Extract month where needed
ai$month <- as.integer(str_sub(string=ai$datetime, start=6, end=7))
ebs$month <- as.integer(str_sub(string=ebs$datetime, start=6, end=7))
goa$month <- as.integer(str_sub(string=goa$datetime, start=6, end=7))
neus$month <- as.integer(str_sub(string=neus$datetime, start=6, end=7))
wctri$month <- as.integer(str_sub(string=wctri$datetime, start=6, end=7))
wcann$month <- as.integer(str_sub(string=wcann$datetime, start=6, end=7))
gmex$month <- as.integer(str_sub(string=gmex$datetime, start=6, end=7))

# Calculate temperature where needed and change values that are clearly off into NA
wcann$surftemp = NA # field is not collected, apparently (or was not provided)
# check outlier temps and remove wonky values_wcann, wctri, ebs, seus, neamap, Newfoundland,  Scotian Shelf, and SoGulf are good-to-go
goa[btemp > 20] <- NA; goa[stemp > 25] <- NA; goa[stemp==0] <- NA
ai[stemp == 0] <- NA
gmex[btemp > 35] <- NA; gmex[stemp < 12] <- NA
neus[stemp==0] <- NA#a bunch of '0' for stemp; many have 'NA' for btemp so they don't really matter, almost all of the rest seem suspicious given the date &/or bottom temp
neus$btemp[neus$haulid=="196721- 309-1160"] <- NA

#     save(ai, dat, ebs, gmex, goa, neamap, neus, seus, wcann, wctri, file="Jim/projections_Nov4_2016.rdata")

# for gmex remove towspeeds that are NA or less than 1
gmex = gmex[gmex$VESSEL_SPD < 5 & gmex$VESSEL_SPD > 0  & !is.na(gmex$VESSEL_SPD),] # trim out vessel speeds 0, unknown, or >5
seus = seus[!(seus$year == 1989 & seus$season == "spring"),] #They sampled at night this year-season


#====================================================
#====================================================

# Add columns to calculate cpue#### 
neamap$distance <- neamap$TowDistTrack # This is the actual tow distance, ie not the straight line between points, but there are some NAs
neamap$distance[is.na(neamap$distance)] <- neamap$TOWDIST[is.na(neamap$distance)] # Where NAs are present for TowDistTrack, use the straightline distance provided (it is a near 1:1 relationhip anyway)
neamap$sampleAreaM2 <- neamap$distance * neamap$NetWidth
neamap$hectares <- neamap$sampleAreaM2/10000
neamap$wtcpue <- neamap$TotWght/neamap$hectares 

# Combine categories of blue crab into one species
neamap$spp[neamap$spp %in% c('Callinectes sapidus, adult fem', 'Callinectes sapidus, juv fem', 'Callinectes sapidus, male')] = "Callinectes sapidus"
neamap <- neamap[!VIMSCODE==9001] # This species code has NA for names, there is only 5 incidence, so no big deal if it actually is a species I lack the name for

