library(data.table) 
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
# This was not ultimately used in the habitat models_hard bottom is too rare, which caused some issues
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
hauls$survey <- hauls$region #keeps all the seasonal surveys separate in "survey"_for calcuation of mean biomass
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
hauls <- hauls[complete.cases(hauls$rugosity, hauls$SBT.seasonal, hauls$GRAINSIZE),] # the multiple sediment and SODA variables have the same NA rows (within each of those two categories) 

hauls$ocean[hauls$region %in% c("AFSC_Aleutians", "AFSC_EBS", "AFSC_GOA", "AFSC_WCTri", "NWFSC_WCAnn")] <- "Pac"
hauls$ocean[hauls$region %in% c("DFO_Newfoundland", "DFO_ScotianShelf","DFO_SoGulf","NEFSC_NEUS", "SCDNR_SEUS", "SEFSC_GOMex", "VIMS_NEAMAP")] <- "Atl"

save(hauls, file='data/master_hauls_March7_2017.RData')

# Run some diagnostics among predictor variables to look more closely at colinearity_above I only compared SST with SBT
haulsE <- hauls[hauls$longrid > -100,]
haulsW <- hauls[hauls$longrid < -100,]
predsE <- data.frame(cbind(rugosity=haulsE$rugosity, GRAINSIZE=haulsE$GRAINSIZE, SST.seasonal.mean=haulsE$SST.seasonal.mean, SST.max=haulsE$SST.max, SBT.seasonal=haulsE$SBT.seasonal, SBT.min=haulsE$SBT.min, SBT.max=haulsE$SBT.max))
predsW <- data.frame(cbind(rugosity=haulsW$rugosity, GRAINSIZE=haulsW$GRAINSIZE, SST.seasonal.mean=haulsW$SST.seasonal.mean, SST.max=haulsW$SST.max, SBT.seasonal=haulsW$SBT.seasonal, SBT.min=haulsW$SBT.min, SBT.max=haulsW$SBT.max))
preds <- data.frame(cbind(rugosity=hauls$rugosity, GRAINSIZE=hauls$GRAINSIZE, SST.seasonal.mean=hauls$SST.seasonal.mean, SST.max=hauls$SST.max, SBT.seasonal=hauls$SBT.seasonal, SBT.min=hauls$SBT.min, SBT.max=hauls$SBT.max))

cor(preds) # SBT.seasonal highly correlated with all temperature metrics, all > 0.805, the other two SBT metrics are both > 0.919 correlation with SBT.seasonal
# all other potential temperature combinations have a correlation of > 0.75
# Rugosity and grainsize are not strongly correlated with any of the temp variables or eachother.
# The 'cor' function is not the same as an R2 value, it is higher
cor(predsW) # correlations are a bit weaker on the west coast for a number of combinations
