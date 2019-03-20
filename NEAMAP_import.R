library(data.table)
library(ggplot2)
library(stringr)

setwd("/Users/jim/Documents/proj_ranges")
neamap.catch <- data.table(read.csv("Jim/NEAMAP_Catch.csv", stringsAsFactors=F))
neamap.spp <- data.table(read.csv("Jim/NEAMAP_sppList.csv", stringsAsFactors=F))
neamap <- merge(neamap.catch, neamap.spp, all.x=T, by="VIMSCODE", sort=F)
neamap.hauls <- data.table(read.csv("Jim/NEAMAP_Station.csv", stringsAsFactors=F))
# length(unique(neamap.hauls$STATION)) # = 2720; so each haul has one unique row
# length(unique(neamap.catch$STATION)) # = 2720; all hauls represented in catch data
neamap <- merge(neamap, neamap.hauls, all.x=T, by="STATION", order=F)

# Make a species list for QAQC purposes ####
# spp <- unique(data.frame(neamap$VIMSCODE, neamap$COMMON, neamap$spp))
# setorder(spp)
# abc <- data.frame(table(spp$neamap.VIMSCODE)); abc[abc$Freq!=1,]
# abc <- data.frame(table(spp$neamap.COMMON)); abc[abc$Freq!=1,]
# abc <- data.frame(table(spp$neamap.spp)); abc[abc$Freq!=1,]
# spp[spp$neamap.spp=="Cnidaria",] # The only instance where redundancy in species naming occurs, but is non-species level and will be dropped later
# rm(abc, spp)

# Add columns to calculate cpue#### 
neamap$distance <- neamap$TowDistTrack # This is the actual tow distance, ie not the straight line between points, but there are some NAs
neamap$distance[is.na(neamap$distance)] <- neamap$TOWDIST[is.na(neamap$distance)] # Where NAs are present for TowDistTrack, use the straightline distance provided (it is a near 1:1 relationhip anyway)
neamap$sampleAreaM2 <- neamap$distance * neamap$NetWidth
neamap$hectares <- neamap$sampleAreaM2/10000
neamap$wtcpue <- neamap$TotWght/neamap$hectares 

# Combine categories of blue crab into one species
neamap$spp[neamap$spp %in% c('Callinectes sapidus, adult fem', 'Callinectes sapidus, juv fem', 'Callinectes sapidus, male')] = "Callinectes sapidus"
neamap <- neamap[!VIMSCODE==9001] # This species code has NA for names, there is only 5 incidence, so no big deal if it actually is a species I lack the name for
 
#Extract date and season
neamap$year <- as.integer(str_sub(string = neamap$DATE, start=1, end=4))
neamap$month <- as.integer(str_sub(string = neamap$DATE, start=6, end=7))
neamap$season <- ifelse(neamap$month < 6, "spring", "fall")

  

