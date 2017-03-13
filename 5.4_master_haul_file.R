library(data.table) 
library(ggplot2)
library(lattice)
library(maps)
setwd('/Users/jamesmorley/Documents/project_velocity')
  
# Import master hauls file and sediment data
load('data/haulInfo_course_Feb8_2017.RData') # load all haul info data_with 1/20 resolution for bathymetry
sediment <- readRDS('data/benthic_hab_hauls.rds') # updated hauls data with sediment variables and Halpern (hard/soft) designation

# Merge in sediment variables to master hauls file
hauls <- merge(haulsTest, sediment[,c('haulid', 'GRAINSIZE', 'GRAVEL', 'SAND', 'MUD', 'habitat')], by='haulid', all.x=T, sort=F)
rm(sediment, haulsTest)
  
# salvage some hauls back with sediment data from the two old sediment files_I only do this for GRAINSIZE, the variable of interest
abc <- hauls[is.na(hauls$GRAINSIZE),]
sedOld1 <- readRDS('data/sediment.all.out.rds') # For sediment variables
sedOld1 <- sedOld1[sedOld1$haulid %in% abc$haulid,]
sedOld1 <- sedOld1[!is.na(sedOld1$GRAINSIZE),]
hauls <- merge(hauls, sedOld1[,c('haulid', 'GRAINSIZE')], by='haulid', all.x=T, sort=F)
hauls$GRAINSIZE <- ifelse(is.na(hauls$GRAINSIZE.x), hauls$GRAINSIZE.y, hauls$GRAINSIZE.x)
hauls$GRAINSIZE.x <- hauls$GRAINSIZE.y <- NULL
rm(abc, sedOld1)
 
# Create column based on 'habitat' that is just hard vs. soft_as currently the other aspect of this is indicating depth (i.e. 'subt', 'shelf', 'slope')
hauls$habitatFact <- as.character(hauls$habitat)
hauls$habitatFact[hauls$habitatFact %in% c("rocky_reef", "shelf_hard", "slope_hard", "deep_hard")] <- "hard"
hauls$habitatFact[hauls$habitatFact %in% c("subt_soft", "slope_soft", "shelf_soft", "deep_soft")] <- "soft" # 1.4% of hauls are designated as 'hard'_this is after removing the SCDNR_SEUSReef survey, which is mostly hard
hauls$habitatFact <- factor(hauls$habitatFact, levels=c('soft', 'hard'), ordered=FALSE) 
 
# Bring in SODA climate variables
soda <- read.csv('data/haulInfo_Dec26_2016_soda_tlf.csv', stringsAsFactors = FALSE) # Bring in the SODA climate variables
hauls <- merge(hauls, soda[,c('haulid', 'SST.actual', 'SST.seasonal.mean', 'SST.min', 'SST.max', 'SBT.actual', 'SBT.seasonal', 'SBT.min', 'SBT.max')], all=TRUE, by='haulid', sort=FALSE)
hauls <- hauls[!hauls$region=='SCDNR_SEUSReef',] # remove the seus trap survey
rm(soda)
 
# Below some plots/regressions to examine how effectively SODA matches observed temps, and also looks at covariation between predictors
plot(surftemp_orig~SST.actual, cex=.1, data=hauls)
summary(lm(surftemp_orig~SST.actual, data=hauls)) # slope=0.865(.0009); p<.0000...1 (14 zeros); r2=0.89;F-statistic: 8.379e+05 on 1 and 102048 DF 
plot(bottemp~SBT.actual, cex=.1, data=hauls)
summary(lm(bottemp~SBT.actual, data=hauls)) # slope=0.976(.0012); p<.0000...1 (14 zeros); r2=0.85; F-statistic: 6.787e+05 on 1 and 120859 DF 
plot(surftemp_orig~SST.seasonal.mean, cex=.1, data=hauls)
summary(lm(surftemp_orig~SST.seasonal.mean, data=hauls)) # slope=0.906(.0009); p<.0000...1 (15 zeros); r2=0.902;F-statistic: 9.379e+05 on 1 and 102048 DF
plot(bottemp~SBT.seasonal, cex=.1, data=hauls)
summary(lm(bottemp~SBT.seasonal, data=hauls)) # slope=1.003(.0012); p<.0000...1 (14 zeros); r2=0.86; F-statistic: 7.41e+05 on 1 and 120859 DF
plot(SST.min~SBT.min, cex=.1, data=hauls)
summary(lm(SST.min~SBT.min, data=hauls)) # slope=1.027(.0011); p<.0000...1 (15 zeros); r2=0.874;F-statistic: 9.422e+05 on 1 and 136126 DF 
plot(SST.max~SBT.max, cex=.1, data=hauls) # MAX BOTTOM TEMPS ARE GENERALLY CONFORMRED TO NOT EXCEED SURFACE TEMPS, SEEMS REASONABLE
summary(lm(SST.max~SBT.max, data=hauls)) # slope=0.698180(0.001237); p<.0000...1 (14 zeros); r2=0.70; F-statistic: 3.184e+05 on 1 and 136126 DF
plot(SST.seasonal.mean~SBT.seasonal, cex=.1, data=hauls) # SEASONAL BOTTOM TEMPS ARE GENERALLY CONFORMRED TO NOT EXCEED SURFACE TEMPS, SEEMS REASONABLE
summary(lm(SST.seasonal.mean~SBT.seasonal, data=hauls)) # slope=0.927217(0.001598); p<.0000...1 (14 zeros); r2=0.71; F-statistic: 3.369e+05 on 1 and 136126 DF
# From all this, it seems that we should drop surface minimum temps as they are highly correlated with bottom minimums (r2 = 0.87)
# Clearly there are some other covariation concerns, but the minimum temps is by far the worst, and makes the most sense considering winter destratification occurs in most/all surveyed areas

# Below are some changes to prepare for models
hauls$survey <- hauls$region #keeps all the seasonal surveys separate in "survey")_for calcuation of mean biomass
hauls$surveyfact <- as.factor(hauls$survey)

# Combine different seasonal surveys within regions
hauls$region[hauls$region %in% c("NEFSC_NEUSFall","NEFSC_NEUSSpring")] <- "NEFSC_NEUS" #NEFSC as one region in "region"
hauls$region[hauls$region %in% c("DFO_NewfoundlandFall","DFO_NewfoundlandSpring")] <- "DFO_Newfoundland"
hauls$region[hauls$region %in% c("DFO_ScotianShelfFall","DFO_ScotianShelfSpring", "DFO_ScotianShelfSummer")] <- "DFO_ScotianShelf"
hauls$region[hauls$region %in% c("SCDNR_SEUSFall","SCDNR_SEUSSpring", "SCDNR_SEUSSummer")] <- "SCDNR_SEUS"
hauls$region[hauls$region %in% c("SEFSC_GOMexFall","SEFSC_GOMexSummer")] <- "SEFSC_GOMex"
hauls$region[hauls$region %in% c("VIMS_NEAMAPFall","VIMS_NEAMAPSpring")] <- "VIMS_NEAMAP"
hauls$regionfact <- as.factor(hauls$region) # the version to use for models

hauls$rugosity <- log(hauls$tri + 1)

# drop 1932 rows with missing predictor values_gams drop rows missing predictor values
hauls <- hauls[complete.cases(hauls$rugosity, hauls$SBT.seasonal, hauls$habitatFact, hauls$GRAINSIZE),] # the multiple sediment and SODA variables have the same NA rows (within each of those two categories) 

hauls$ocean[hauls$region %in% c("AFSC_Aleutians", "AFSC_EBS", "AFSC_GOA", "AFSC_WCTri", "NWFSC_WCAnn")] <- "Pac"
hauls$ocean[hauls$region %in% c("DFO_Newfoundland", "DFO_ScotianShelf","DFO_SoGulf","NEFSC_NEUS", "SCDNR_SEUS", "SEFSC_GOMex", "VIMS_NEAMAP")] <- "Atl"

save(hauls, file='data/master_hauls_March7_2017.RData')

# SCRIPT ENDS HERE

# ==========================================================================================================================================
# Below produces 'haulsTrim' which is used for estimating annual mean cpue, which we no longer use.
# If we bring this predictor back in, probably should go over this again (esp. for DFO surveys), b/c we added new sediment data to fill in some holes since the first time
# ==========================================================================================================================================

# The following year-survey combos grossly undersampled and thus not reliable for annual biomass predictor
hauls <- hauls[!(hauls$survey== 'DFO_NewfoundlandFall' & hauls$year == 2008),] # only 7 hauls
hauls <- hauls[!(hauls$survey== 'SEFSC_GOMexFall' & hauls$year == 1986),] # Only 24 hauls, all in very small area

# Check out spatial effort history for each survey_Need to delineate consistent sampling footprint of each survey to estimate annual cpue
# This code I recycled with each region (as needed) to determine how to delineate each survey region
# ggplot(haulsTrim[haulsTrim$survey == 'AFSC_WCTri',], aes(x=lon, y=lat)) + geom_point() + facet_wrap(~year)
# with(hauls[hauls$survey== 'AFSC_WCTri' & hauls$year == 1983,], range(lon))
# hauls[(hauls$survey== 'NEFSC_NEUSFall' & hauls$year == 2014 & hauls$lat < 34.5),]

# Need to drop the 'inshore' strata for NEUS and the 'offshore' strata for SEUS
inshore <- c('3030','3040','3060','3070','3090','3100','3120','3130','3150','3160','3180','3190','3210','3220','3240','3250','3270','3280','3300','3310','3330','3340','3360',
             '3370','3390','3400','3420','3430','3550','3580','7510') 
offshore <- c('24', '26', '28', '30', '32', '34', '36', '38', '40', '42', '44', '46', '48', '50', '52', '54', '56', '60', '62', '64', '66', '68')
 
# Create separate file to calculate annual average cpue for species 'i'_trims out hauls from inconsistent sampling region
haulsTrim <- hauls[!(hauls$survey == 'DFO_ScotianShelfFall' & hauls$lat > 46.1),]
haulsTrim <- haulsTrim[!(haulsTrim$survey == 'DFO_ScotianShelfFall' & haulsTrim$lat > 43.1 & haulsTrim$lon < -66),]
haulsTrim <- haulsTrim[!(haulsTrim$survey == 'DFO_ScotianShelfFall' & haulsTrim$lat > 44 & haulsTrim$lon < -64.85),]
haulsTrim <- haulsTrim[!(haulsTrim$survey == 'DFO_ScotianShelfSpring' & haulsTrim$lat > 43.1 & haulsTrim$lon < -66),]
haulsTrim <- haulsTrim[!(haulsTrim$survey == 'DFO_ScotianShelfSpring' & haulsTrim$lat > 44.6 & haulsTrim$lon < -64.3),]
haulsTrim <- haulsTrim[!(haulsTrim$survey == 'DFO_NewfoundlandFall' & haulsTrim$lat > 55.33),]
haulsTrim <- haulsTrim[!(haulsTrim$survey == 'NEFSC_NEUSFall' & haulsTrim$lat < 34.5),]
haulsTrim <- haulsTrim[!(haulsTrim$survey == 'NEFSC_NEUSFall' & haulsTrim$lon > -65.5),]
haulsTrim <- haulsTrim[!(haulsTrim$survey == 'NEFSC_NEUSSpring' & haulsTrim$lat < 34.5),]
haulsTrim <- haulsTrim[!(haulsTrim$survey == 'NEFSC_NEUSSpring' & haulsTrim$lon > -65.5),]
haulsTrim <- haulsTrim[!haulsTrim$stratum %in% inshore,] # These strata values only found in NEUS
haulsTrim <- haulsTrim[!(haulsTrim$region == 'SEFSC_GOMex' & haulsTrim$lon > -87.5),]
haulsTrim <- haulsTrim[!(haulsTrim$stratum %in% offshore & haulsTrim$region == 'SCDNR_SEUS'),] # removes the offshore strata for SEUS in fall and spring
haulsTrim <- haulsTrim[!(haulsTrim$survey == 'AFSC_WCTri' & haulsTrim$lat > 48.5),]
haulsTrim <- haulsTrim[!(haulsTrim$survey == 'AFSC_WCTri' & haulsTrim$lon > -121.75),]
rm(inshore, offshore)
