# This script takes the output from the projections for parameter uncertainty, with all the parameter resampling
  # Mainly it creates 'uncert_final' for each species, which is used in one of the primary figures.
  # There are also some basic figures on here, but not used in the MS. Some of the code here is repeated on the main script for the figure.
  # AT THE BOTTOM OF THIS SCRIPT IS CODE TO CREATE THE giffs FOR VISUALIZING PARAMETER UNCERTAINTY

setwd('/Users/jamesmorley/Documents/project_velocity')
library(Hmisc)
library(data.table)
        
RCP <- c(26,45,85)
season <- 'jas'
modelrun <- c('bcc-csm1-1-m','bcc-csm1-1','CanESM2','CCSM4','CESM1-CAM5','CNRM-CM5','GFDL-CM3','GFDL-ESM2M','GFDL-ESM2G','GISS-E2-R','GISS-E2-H','HadGEM2-ES','IPSL-CM5A-LR','IPSL-CM5A-MR','MIROC-ESM','MIROC5','MPI-ESM-LR','NorESM1-ME')
nicheMod <- c('GLM', 'GAM')
nicheMethod <- c('PA', 'biom')
years <- c("2007-2020", "2021-2040", "2041-2060", "2061-2080", "2081-2100")
 
species <- 'sebastes alutus_Pac'
projfolder <- paste('output/CEmodels_Proj_Uncertainty_2018', species, sep='_')
#load('data/EEZ_grid_east.RData') # load the EEZ for U.S.
load('data/EEZ_grid_west.RData') # load the EEZ for U.S.

filename <- paste('output/CEmodels_proj_Nov2017/', species, '_rcp26_jas_prediction_AGG.RData', sep='')
load(filename)
#pred.agg$EEZ <- EEZ_east
pred.agg$EEZ <- EEZ_west

pred.agg <- pred.agg[,c(2,3,20)] # reference file for the U.S. EEZ
pred.agg <- unique(pred.agg)
rm(EEZ_east, EEZ_west)

# Create a grid area correction file ======================
library(raster)

filename <- paste(projfolder, '/', species, '_rcp26_bcc-csm1-1-m_GAM_iter40.RData', sep='') # to get lat value boundaries
load(filename) 

grid_rast <- raster(ymn=min(biom_final$latitude) - 4, ymx=max(biom_final$latitude) + 6, xmn=min(biom_final$longitude), xmx=max(biom_final$longitude), resolution=1/20)
grid_rast[] <- area(grid_rast)[] # in km2
grid_area <- data.frame(rasterToPoints(grid_rast), stringsAsFactors = FALSE) # all the lat/lon expanded and all three values in the rasterStack included
grid_area <- unique(data.frame(lat = grid_area$y, area = grid_area$layer))

area_mod <- lm(area~lat + I(lat^2), grid_area) # use this model (r2 value of 1.0) to get grid-cell areas
save(area_mod, file='grid_area_correction_lm.RData')
rm(grid_rast, grid_area, PA_final, biom_final)

# Begin loop =========================================================
  
uncert <- data.frame(nicheMod=character(), rcp=numeric(), modelrun=character(), years=character(), iter=integer(), centroid_lat=numeric(), centroid_lon=numeric(),  habitat_south=numeric(), habitat_north=numeric(), centroidPA_lat=numeric(), centroidPA_lon=numeric(), habitat_southPA=numeric(), habitat_northPA=numeric(), stringsAsFactors = F)
uncert_final <- NULL
   
for(i in 1:3){ # RCP
  rcp <- RCP[i]
  print(paste('Beginning rcp', rcp, '   ', Sys.time()))
  for(j in 1:18){ # GCMs
    gcm <- modelrun[j]
    print(paste('model run:', j, '  ', gcm))
    for(k in 1:2){ #nicheMod
      nicModel <- nicheMod[k]
      # Load the projection
      filename <- paste(projfolder, '/', species, '_rcp', rcp, '_', gcm, '_', nicModel, '_iter40.RData', sep='') 
      load(filename) 
      
      # Carve out areas of non interest (e.g., GMex, Canada, etc.) & U.S. EEZ
      PA_final <- merge(PA_final, pred.agg, by=c('latitude', 'longitude'), all.x=T) # merge east coast EEZ file
      biom_final <- merge(biom_final, pred.agg, by=c('latitude', 'longitude'), all.x=T) # merge east coast EEZ file
      
      # I recycled these plot lines for a number of things throughout the script
      #plot(latitude~longitude, cex=.1, PA_final[PA_final$year_bin=='2021-2040',])
      #points(latitude~longitude, cex=.1, col='green', PA_final[PA_final$region=='South',])
      
      #PA_final <- PA_final[PA_final$EEZ == T,]
      #biom_final <- biom_final[biom_final$EEZ == T,]
      # DROP GULF OF MEXICO
      if(grepl('_Atl', species)){
        PA_final <- PA_final[!PA_final$longitude < -82,]
        PA_final <- PA_final[!(PA_final$longitude < -80.75 & PA_final$latitude < 27),]
        biom_final <- biom_final[!biom_final$longitude < -82,]
        biom_final <- biom_final[!(biom_final$longitude < -80.75 & biom_final$latitude < 27),]
      }
      # LOAD THE GRID CELL AREA CORRECTIONS HERE
      biom_final$area <- coef(area_mod)[1] + coef(area_mod)[2]*biom_final$latitude + coef(area_mod)[3]*(biom_final$latitude^2)
      PA_final$area <- coef(area_mod)[1] + coef(area_mod)[2]*PA_final$latitude + coef(area_mod)[3]*(PA_final$latitude^2)
      PA_final$area_rel <- PA_final$area/max(PA_final$area)
      
      # Identify south and north areas_this depends on species
      if(species == 'paralichthys dentatus_Atl' | species == 'centropristis striata_Atl'){ # For a Cape Hatteras split (BSB)
        PA_final$region <- ifelse(PA_final$latitude > 35.25, 'North', 'South')
        PA_final$region[PA_final$EEZ == FALSE] <- NA # drop non U.S. waters
        biom_final$region <- ifelse(biom_final$latitude > 35.25, 'North', 'South')
        biom_final$region[biom_final$EEZ == FALSE] <- NA # drop non U.S. waters
      }
      if(species == 'homarus americanus_Atl'){ # For a Gulf of Maine/Georges Bank versus Southern NE/MAB split (Lobster)
        PA_final$region <- ifelse(PA_final$longitude > -70, 'North', 'South')
        PA_final$region[PA_final$region == 'South' & PA_final$latitude > 41.75] <- 'North'
        PA_final$region[PA_final$latitude < 35.25] <- NA # drop south of Cape Hatteras
        PA_final$region[PA_final$EEZ == FALSE] <- NA # drop non U.S. waters
      
        biom_final$region <- ifelse(biom_final$longitude > -70, 'North', 'South')
        biom_final$region[biom_final$region == 'South' & biom_final$latitude > 41.75] <- 'North'
        biom_final$region[biom_final$latitude < 35.25] <- NA # drop south of Cape Hatteras
        biom_final$region[biom_final$EEZ == FALSE] <- NA # drop non U.S. waters
        }
      if(species == 'anoplopoma fimbria_Pac' | species == 'hippoglossus stenolepis_Pac' | species == 'sebastes alutus_Pac'){ # For an Alaska-U.S. split, won't include Canada for habitat analysis
        PA_final$region <- ifelse(PA_final$latitude > 50, 'North', 'South')
        PA_final$region[PA_final$EEZ == FALSE] <- NA # drop non U.S. waters
        biom_final$region <- ifelse(biom_final$latitude > 50, 'North', 'South')
        biom_final$region[biom_final$EEZ == FALSE] <- NA # drop non U.S. waters
      }
      if(species == 'doryteuthis opalescens_Pac'){ # For U.S. west coast vs. West Canada split
        PA_final$region <- ifelse(PA_final$longitude > -130 & PA_final$EEZ == TRUE, 'South', NA)
        PA_final$region[PA_final$latitude > 45 & PA_final$latitude < 58 & PA_final$EEZ == FALSE] <- 'North'
        biom_final$region <- ifelse(biom_final$longitude > -130 & biom_final$EEZ == TRUE, 'South', NA)
        biom_final$region[biom_final$latitude > 45 & biom_final$latitude < 58 & biom_final$EEZ == FALSE] <- 'North'
      }
      
      for(l in 1:5){ # Time period
        time_bin <- years[l]
        PA_sub <- PA_final[PA_final$year_bin == time_bin,]
        biom_sub <- biom_final[biom_final$year_bin == time_bin,]
        for(m in 1:40){ # iteration
          # Calculate habitat area and centroid for each region and each iteration, using the grid area scalar correction
          habitat_south <- sum(biom_sub[,m+2][!is.na(biom_sub$region) & biom_sub$region=='South'] * biom_sub$area[!is.na(biom_sub$region) & biom_sub$region=='South'])
          habitat_north <- sum(biom_sub[,m+2][!is.na(biom_sub$region) & biom_sub$region=='North'] * biom_sub$area[!is.na(biom_sub$region) & biom_sub$region=='North'])
          habitat_southPA <- sum(PA_sub[,m+2][!is.na(biom_sub$region) & PA_sub$region=='South'] * PA_sub$area_rel[!is.na(biom_sub$region) & PA_sub$region=='South'])
          habitat_northPA <- sum(PA_sub[,m+2][!is.na(biom_sub$region) & PA_sub$region=='North'] * PA_sub$area_rel[!is.na(biom_sub$region) & PA_sub$region=='North'])
          
          centroid_lat <- wtd.mean(biom_sub$latitude, weights = (biom_sub[,m+2] * biom_sub$area))
          centroid_lon <- wtd.mean(biom_sub$longitude, weights = (biom_sub[,m+2] * biom_sub$area))
          centroidPA_lat <- wtd.mean(PA_sub$latitude, weights = (PA_sub[,m+2] * PA_sub$area_rel))
          centroidPA_lon <- wtd.mean(PA_sub$longitude, weights = (PA_sub[,m+2] * PA_sub$area_rel))
          
          uncert[m,] <- data.frame(nicheMod=nicModel, rcp=rcp, modelrun=gcm, years=time_bin, iter=m, centroid_lat=centroid_lat, centroid_lon=centroid_lon, habitat_south=habitat_south, habitat_north=habitat_north, centroidPA_lat=centroidPA_lat, centroidPA_lon=centroidPA_lon, habitat_southPA=habitat_southPA, habitat_northPA=habitat_northPA, stringsAsFactors = F)
        }
        uncert_final <- rbind(uncert_final, uncert)
      }
    }
  }
}
  
filename <- paste('output/Uncert_summary_', species, '_Nov2018.RData', sep='')
save(uncert_final, file=filename)

 
# SUMMARY GRAPHS AND ANALYSES =======================================================
# BEGINNING WITH SQUID, ADD A %CHANGE FOR THE BIOMASS FIGURES_IN ADDITION TO THE ABSOLUTE VALUES 
    # THEN PERHAPS JUST ADD TO MAIN FIGURES SCRIPT_JUST ADAPT FROM uncert_final
 
# Make an x-axis
uncert_final$xAxis <- 1 
uncert_final$xAxis[uncert_final$years=='2021-2040'] <- 2
uncert_final$xAxis[uncert_final$years=='2041-2060'] <- 3
uncert_final$xAxis[uncert_final$years=='2061-2080'] <- 4
uncert_final$xAxis[uncert_final$years=='2081-2100'] <- 5

# Raw data plots
plot(habitat_north~xAxis, ylim=c(min(uncert_final$habitat_north), max(uncert_final$habitat_north)), uncert_final)
plot(habitat_south~xAxis,  ylim=c(min(uncert_final$habitat_south), max(uncert_final$habitat_south)), uncert_final)
plot(habitat_northPA~xAxis, ylim=c(min(uncert_final$habitat_northPA), max(uncert_final$habitat_northPA)), uncert_final)
plot(habitat_southPA~xAxis,  ylim=c(min(uncert_final$habitat_southPA), max(uncert_final$habitat_southPA)), uncert_final)
plot(centroid_lat~xAxis, ylim=c(min(uncert_final$centroid_lat), max(uncert_final$centroid_lat)), uncert_final)
plot(centroidPA_lat~xAxis,  ylim=c(min(uncert_final$centroidPA_lat), max(uncert_final$centroidPA_lat)), uncert_final)


# Partition the variance with a general lm ================================================================================
# first do centroid model
  
library(reshape2)
 
uncert_cent <- uncert_final[,c(1:6,10,14)]
uncert_cent <- melt(uncert_cent, id.vars=c('nicheMod', 'rcp', 'modelrun', 'years', 'iter', 'xAxis'), measure.vars=c('centroid_lat', 'centroidPA_lat'), variable.name='method', value.name='centroid')
# Make factor columns of everything...just to be sure everything works right
uncert_cent$nicheModel <- as.factor(paste(uncert_cent$nicheMod, uncert_cent$method, sep='_'))
uncert_cent$rcp_fact <- as.factor(uncert_cent$rcp)
uncert_cent$modelrun_fact <- as.factor(uncert_cent$modelrun)
uncert_cent$years_fact <- as.factor(uncert_cent$years) # for regression tree only

# REGRESSION TREE to see if any factors are contributing most of variance and should be excluded as 'wonky' =============================================================================================
library(rpart) # A package that does regression trees
tree_cent <- rpart(centroid~nicheModel + rcp_fact + modelrun_fact + years_fact, method='anova', data=uncert_cent)
complex.param <- tree_cent$cptable[which.min(tree_cent$cptable[,"xerror"]),"CP"] # calculate a penalty value for overfitting
tree_cent2 <- prune(tree_cent, cp=complex.param) # refit the model with the penalty value

plot(tree_cent2)
# text(tree_cent2, cex=.7, use.n=TRUE, all=TRUE)
plot(tree_cent2, uniform=T) # use this and add an inset with the actual branch lengths
text(tree_cent2, cex=.7, use.n=TRUE, all=TRUE)

#printcp(tree_cent2)
#plotcp(tree_cent2)
rsq.rpart(tree_cent2)
print(tree_cent2) # explains how everything was split_what went where
summary(tree_cent2) #includes relative importance of the variables_among a lot of other info

 
# Fit linear models in order to partition the SS ====================================================================

# I THINK I CAN DO A DOMINANCE ANALYSIS_BUT STATE THE ASSUMPTION THAT ALL THE VARIATION IS ACCOUNTED FOR IN THIS 'UNIVERSE' OF PROJECTIONS
  # STATE THAT A SIMPLE PARTITIONING OF THE SS WAS IMPOSSIBLE DUE SIGNIFICANT INTERACTION TERMS
  # EACH BAR IN GRAPH WILL SHOW THE TOTAL SS TO SHOW HOW OVERALL VARIATION CHANGES THROUGH TIME
  # CALCULATE PERCENTAGE INFLUENCE OF PARAMETER UNCERTAINTY BASED ON RESIDUAL SS AND TOTAL SS.
  # THE REMAINING SS IN THE BARS WILL BE PARTITIONED BASED ON RELATIVE IMPORTANCE OF PREDICTORS, RESCALED TO THE REMAINING SS (AFTER RESIDS FOR PARAM UNCER)

# HERE I USE THE FUNCTION FROM SCRIPT: Dominance analysis.R

mod_cent2020 <- dom_centroid('2007-2020', uncert_cent)
mod_cent2040 <- dom_centroid('2021-2040', uncert_cent)
mod_cent2060 <- dom_centroid('2041-2060', uncert_cent)
mod_cent2080 <- dom_centroid('2061-2080', uncert_cent)
mod_cent2100 <- dom_centroid('2081-2100', uncert_cent)

SStable <- data.table(rbind(mod_cent2020, mod_cent2040, mod_cent2060, mod_cent2080, mod_cent2100))
SStable_trans <- transpose(SStable)
colnames(SStable_trans) <- SStable$years
SStable_trans <- SStable_trans[2:5,]
rownames(SStable_trans) <- names(SStable)[2:5]
SStable_trans <- data.matrix(SStable_trans)
barplot(SStable_trans, main="", ylab='Variation (Sum of Squares)', col=c("dodgerblue2","orange","green3","red"), legend = c('Niche model','GCM','RCP','Parameter'), args.legend=list(x='topleft')) 
box()


# VIOLIN PLOTS ========================================================================
# uncert_cent_Mod <- uncert_cent[!uncert_cent$nicheModel=='GLM_centroid_lat',]
uncert_cent_Mod <- uncert_cent # just so I don't need multiple formats below

cent_mean <- aggregate(list(mean_cent = uncert_cent_Mod$centroid), by=list(years = uncert_cent_Mod$years), FUN=mean) 
cent_mean$xAxis <- c(1:5)

ggplot(uncert_cent_Mod, aes(x=years, y=centroid)) + geom_violin() + 
  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_line(data=cent_mean, aes(x=xAxis, y=mean_cent)) + 
  geom_violin(data=uncert_cent_Mod, aes(x=years, y=centroid), fill='gray80') + geom_point(data=cent_mean, aes(x=xAxis, y=mean_cent), shape=21, size=3, fill='white')

# Save files used to make the graphs
filename <- paste('output/CEmodels_Uncertainty_2018/', species, '_figures_centroid.RData', sep='')
save(uncert_final, SStable, uncert_cent_Mod, cent_mean, uncert_cent, file=filename)

 

# =======================================================================================
# Now do habitat North ==================================================================
# =======================================================================================

# FOR BIOMASS, THE PA AND BIOMASS MODELS AREN'T ON THE SAME PLAYING FIELD
 
uncert_cent <- uncert_final[,c(1:5,9,13,14)]
uncert_cent <- melt(uncert_cent, id.vars=c('nicheMod', 'rcp', 'modelrun', 'years', 'iter', 'xAxis'), measure.vars=c('habitat_north', 'habitat_northPA'), variable.name='method', value.name='habitat_north')
uncert_cent$nicheModel <- as.factor(paste(uncert_cent$nicheMod, uncert_cent$method, sep='_'))
uncert_cent$rcp_fact <- as.factor(uncert_cent$rcp)
uncert_cent$modelrun_fact <- as.factor(uncert_cent$modelrun)
uncert_cent$years_fact <- as.factor(uncert_cent$years) # for regression tree only

# REGRESSION TREE to see if any factors are contributing most of variance and should be excluded as 'wonky' =============================================================================================
tree_cent <- rpart(habitat_north~nicheModel + rcp_fact + modelrun_fact + years_fact, method='anova', data=uncert_cent)
complex.param <- tree_cent$cptable[which.min(tree_cent$cptable[,"xerror"]),"CP"] # calculate a penalty value for overfitting
tree_cent2 <- prune(tree_cent, cp=complex.param) # refit the model with the penalty value

plot(tree_cent2)
# text(tree_cent2, cex=.7, use.n=TRUE, all=TRUE)
plot(tree_cent2, uniform=T) # use this and add an inset with the actual branch lengths
text(tree_cent2, cex=.7, use.n=TRUE, all=TRUE)

#printcp(tree_cent2)
#plotcp(tree_cent2)
rsq.rpart(tree_cent2)
print(tree_cent2) # explains how everything was split_what went where
summary(tree_cent2) #includes relative importance of the variables_among a lot of other info


# Linear mixed effects models ==================================================
library(lme4)
library(MuMIn)

colnames(uncert_cent)[8] <- 'habitat'
# fit linear models_WITH THE RANDOM EFFECTS 
mod_cent2020 <- dom_habitat('2007-2020', uncert_cent)
mod_cent2040 <- dom_habitat('2021-2040', uncert_cent)
mod_cent2060 <- dom_habitat('2041-2060', uncert_cent)
mod_cent2080 <- dom_habitat('2061-2080', uncert_cent)
mod_cent2100 <- dom_habitat('2081-2100', uncert_cent)

SStable <- data.table(rbind(mod_cent2020, mod_cent2040, mod_cent2060, mod_cent2080, mod_cent2100))
SStable_trans <- transpose(SStable)
colnames(SStable_trans) <- SStable$years
SStable_trans <- SStable_trans[2:5,]
rownames(SStable_trans) <- names(SStable)[2:5]
SStable_trans <- data.matrix(SStable_trans)
barplot(SStable_trans, main="", ylab='Variation (Sum of Squares)', col=c("dodgerblue2","orange","green3","red"), legend = c('Niche model','GCM','RCP','Parameter'), args.legend=list(x='topleft')) 
box()


# VIOLIN PLOTS ========================================================================
uncert_cent_Mod <- uncert_cent[uncert_cent$method=='habitat_north',]
uncert_cent_ModPA <- uncert_cent[uncert_cent$method=='habitat_northPA',]

cent_mean <- aggregate(list(habitat = uncert_cent_Mod$habitat), by=list(years = uncert_cent_Mod$years), FUN=mean) 
cent_mean$xAxis <- c(2,4,6,8,10)
cent_meanPA <- aggregate(list(habitat = uncert_cent_ModPA$habitat), by=list(years = uncert_cent_ModPA$years), FUN=mean) 
cent_meanPA$xAxis <- c(1,3,5,7,9)

#ggplot(uncert_cent_Mod, aes(x=years, y=sqrt(habitat_north))) + geom_violin() + geom_line(data=cent_mean, aes(x=xAxis, y=habitat_north)) + 
 # geom_violin(data=uncert_cent_Mod, aes(x=years, y=habitat_north)) + geom_point(data=cent_mean, aes(x=xAxis, y=habitat_north), shape=1, size=3) 
#ggplot(uncert_cent_ModPA, aes(x=years, y=habitat_north)) + geom_violin() + geom_line(data=cent_meanPA, aes(x=xAxis, y=habitat_north)) + 
 # geom_violin(data=uncert_cent_ModPA, aes(x=years, y=habitat_north)) + geom_point(data=cent_meanPA, aes(x=xAxis, y=habitat_north), shape=1, size=3) 

# combined violin plot for both niche model types 
# NEED TO CHANGE THE X-AXIS, AND see if I can specify the tick marks on the right to include some smaller values
# lobster: 1.72
# fluke: 1.62
# BSB: 1.72
# halibut: 1.465
val <- 1.7
ggplot(uncert_cent_ModPA, aes(x=as.factor(xAxis-0.05), y=habitat)) + geom_violin() +
  geom_violin(data = uncert_cent_Mod, aes(x=as.factor(xAxis+.05), y=habitat^(1/val))) + # adjust the exponent denominator here and elsewhere in block to adjust scaling of the second plot and axis
  geom_line(data=cent_meanPA, aes(x=xAxis, y=habitat)) +
  geom_line(data=cent_mean, aes(x=xAxis, y=habitat^(1/val))) +
  geom_violin(data=uncert_cent_ModPA, aes(x=as.factor(xAxis-0.05), y=habitat)) +
  geom_violin(data = uncert_cent_Mod, aes(x=as.factor(xAxis+.05), y=habitat^(1/val)), fill='grey75') +
  geom_point(data=cent_meanPA, aes(x=xAxis, y=habitat), shape=1, size=2) +
  geom_point(data=cent_mean, aes(x=xAxis, y=habitat^(1/val)),size=2, fill="black", colour="black") +
  scale_y_continuous(expand = c(0, 0),  limits=c(0,2500), sec.axis = ~.^val) # make sure this matches denominators above, to scale secon axis correctly

# Save files used to make the graphs
filename <- paste('output/CEmodels_Uncertainty_2018/', species, '_figures_habNorth.RData', sep='')
save(uncert_final, SStable, uncert_cent_Mod, uncert_cent_ModPA, cent_mean, cent_meanPA, uncert_cent, file=filename)


# =====================================================================================================================
# Try the habitat analysis, but using the %change between periods as the variable of interest==========================

uncert_cent <- uncert_cent[,c(2,3,4,5,8,9)]
niche <- c('GAM_habitat_north','GAM_habitat_northPA','GLM_habitat_north','GLM_habitat_northPA')
hab_change <- data.frame(nicheMod=character(), rcp=character(), modelrun=character(), iter=integer(), time1=numeric(), time2=numeric(), time3=numeric(), time4=numeric(), stringsAsFactors = F)
hab_chang_fin <- NULL

for(i in 1:3){ # RCP
  rcp <- RCP[i]
  for(j in 1:18){ # GCMs
    gcm <- modelrun[j]
    for(k in 1:4){ #nicheMod
      nicModel <- niche[k]
      uncert_cent_sub <- uncert_cent[uncert_cent$rcp==rcp & uncert_cent$modelrun==gcm & uncert_cent$nicheModel==nicModel,]
      uncert_cent_cast <- dcast(uncert_cent_sub, rcp+modelrun+iter+nicheModel~years, value.var='habitat')  
      time1 <- ((uncert_cent_cast$'2021-2040' - uncert_cent_cast$'2007-2020')/uncert_cent_cast$'2007-2020')*100
      time2 <- ((uncert_cent_cast$'2041-2060' - uncert_cent_cast$'2007-2020')/uncert_cent_cast$'2007-2020')*100
      time3 <- ((uncert_cent_cast$'2061-2080' - uncert_cent_cast$'2007-2020')/uncert_cent_cast$'2007-2020')*100
      time4 <- ((uncert_cent_cast$'2081-2100' - uncert_cent_cast$'2007-2020')/uncert_cent_cast$'2007-2020')*100
      
      hab_change <- data.frame(nicheMod=nicModel, rcp=rcp, modelrun=gcm, iter=1:40, time1=time1, time2=time2, time3=time3, time4=time4, stringsAsFactors = F)
      hab_chang_fin <- rbind(hab_chang_fin, hab_change)
    }
  }
}

hab_chang_fin <- melt(hab_chang_fin, id.vars=c('nicheMod', 'rcp', 'modelrun', 'iter'), measure.vars=c('time1','time2','time3','time4'), variable.name='years', value.name='habitat')

hab_chang_fin$nicheModel <- as.factor(hab_chang_fin$nicheMod)
hab_chang_fin$rcp_fact <- as.factor(hab_chang_fin$rcp)
hab_chang_fin$modelrun_fact <- as.factor(hab_chang_fin$modelrun)
hab_chang_fin$years_fact <- as.factor(hab_chang_fin$years) # for regression tree only

# ======================================================================================================================
# REGRESSION TREE to see if any factors are contributing most of variance and should be excluded as 'wonky' =============================================================================================
tree_cent <- rpart(habitat~nicheModel + rcp_fact + modelrun_fact + years_fact, method='anova', data=hab_chang_fin)
complex.param <- tree_cent$cptable[which.min(tree_cent$cptable[,"xerror"]),"CP"] # calculate a penalty value for overfitting
tree_cent2 <- prune(tree_cent, cp=complex.param) # refit the model with the penalty value

plot(tree_cent2)
# text(tree_cent2, cex=.7, use.n=TRUE, all=TRUE)
plot(tree_cent2, uniform=T) # use this and add an inset with the actual branch lengths
text(tree_cent2, cex=.7, use.n=TRUE, all=TRUE)

#printcp(tree_cent2)
#plotcp(tree_cent2)
rsq.rpart(tree_cent2)
print(tree_cent2) # explains how everything was split_what went where
summary(tree_cent2) #includes relative importance of the variables_among a lot of other info


# Linear mixed effects models ==================================================
# fit linear models
mod_cent2020 <- dom_percChan('time1', hab_chang_fin)
mod_cent2040 <- dom_percChan('time2', hab_chang_fin)
mod_cent2060 <- dom_percChan('time3', hab_chang_fin)
mod_cent2080 <- dom_percChan('time4', hab_chang_fin)

SStable <- data.table(rbind(mod_cent2020, mod_cent2040, mod_cent2060, mod_cent2080))
SStable_trans <- transpose(SStable)
colnames(SStable_trans) <- SStable$years
SStable_trans <- SStable_trans[2:5,]
rownames(SStable_trans) <- names(SStable)[2:5]
SStable_trans <- data.matrix(SStable_trans)
barplot(SStable_trans, main="", ylab='Variation (Sum of Squares)', col=c("dodgerblue2","orange","green3","red"), legend = c('Niche model','GCM','RCP','Parameter'), args.legend=list(x='topleft')) 
box()


# VIOLIN PLOTS ========================================================================
cent_mean <- aggregate(list(habitat = hab_chang_fin$habitat), by=list(years = hab_chang_fin$years), FUN=mean) 
cent_mean$xAxis <- c(1,2,3,4)

# presently I have this logged so probably need to take off for next species
ggplot(hab_chang_fin, aes(x=years_fact, y=log(habitat))) + geom_violin() + 
  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_line(data=cent_mean, aes(x=xAxis, y=log(habitat))) + 
  geom_violin(data=hab_chang_fin, aes(x=years_fact, y=log(habitat)), fill='gray80') + geom_point(data=cent_mean, aes(x=xAxis, y=log(habitat)), shape=21, size=3, fill='white')


# Save files used to make the graphs
filename <- paste('output/CEmodels_Uncertainty_2018/', species, '_figures_habNorth_percChan.RData', sep='')
save(uncert_final, SStable, hab_chang_fin, cent_mean, uncert_cent, file=filename)


# =======================================================================================
# Now do habitat South ==================================================================
# =======================================================================================
 
uncert_cent <- uncert_final[,c(1:5,8,12,14)]
uncert_cent <- melt(uncert_cent, id.vars=c('nicheMod', 'rcp', 'modelrun', 'years', 'iter', 'xAxis'), measure.vars=c('habitat_south', 'habitat_southPA'), variable.name='method', value.name='habitat_south')
uncert_cent$nicheModel <- as.factor(paste(uncert_cent$nicheMod, uncert_cent$method, sep='_'))
uncert_cent$rcp_fact <- as.factor(uncert_cent$rcp)
uncert_cent$modelrun_fact <- as.factor(uncert_cent$modelrun)
uncert_cent$years_fact <- as.factor(uncert_cent$years) # for regression tree only

# REGRESSION TREE to see if any factors are contributing most of variance and should be excluded as 'wonky' =============================================================================================
tree_cent <- rpart(habitat_south~nicheModel + rcp_fact + modelrun_fact + years_fact, method='anova', data=uncert_cent)
complex.param <- tree_cent$cptable[which.min(tree_cent$cptable[,"xerror"]),"CP"] # calculate a penalty value for overfitting
tree_cent2 <- prune(tree_cent, cp=complex.param) # refit the model with the penalty value

plot(tree_cent2)
# text(tree_cent2, cex=.7, use.n=TRUE, all=TRUE)
plot(tree_cent2, uniform=T) # use this and add an inset with the actual branch lengths
text(tree_cent2, cex=.7, use.n=TRUE, all=TRUE)

#printcp(tree_cent2)
#plotcp(tree_cent2)
rsq.rpart(tree_cent2)
print(tree_cent2) # explains how everything was split_what went where
summary(tree_cent2) #includes relative importance of the variables_among a lot of other info


# fit linear models_WITH THE RANDOM EFFECTS ==========================================================
colnames(uncert_cent)[8] <- 'habitat'
# fit linear models_WITH THE RANDOM EFFECTS 
mod_cent2020 <- dom_habitat('2007-2020', uncert_cent)
mod_cent2040 <- dom_habitat('2021-2040', uncert_cent)
mod_cent2060 <- dom_habitat('2041-2060', uncert_cent)
mod_cent2080 <- dom_habitat('2061-2080', uncert_cent)
mod_cent2100 <- dom_habitat('2081-2100', uncert_cent)

SStable <- data.table(rbind(mod_cent2020, mod_cent2040, mod_cent2060, mod_cent2080, mod_cent2100))
SStable_trans <- transpose(SStable)
colnames(SStable_trans) <- SStable$years
SStable_trans <- SStable_trans[2:5,]
rownames(SStable_trans) <- names(SStable)[2:5]
SStable_trans <- data.matrix(SStable_trans)
barplot(SStable_trans, main="", ylab='Variation (Sum of Squares)', col=c("dodgerblue2","orange","green3","red"), legend = c('Niche model','GCM','RCP','Parameter'), args.legend=list(x='topleft')) 
box()


# VIOLIN PLOTS ========================================================================================
uncert_cent_Mod <- uncert_cent[uncert_cent$method=='habitat_south',]
uncert_cent_ModPA <- uncert_cent[uncert_cent$method=='habitat_southPA',]

cent_mean <- aggregate(list(habitat = uncert_cent_Mod$habitat), by=list(years = uncert_cent_Mod$years), FUN=mean) 
cent_mean$xAxis <- c(2,4,6,8,10)
cent_meanPA <- aggregate(list(habitat = uncert_cent_ModPA$habitat), by=list(years = uncert_cent_ModPA$years), FUN=mean) 
cent_meanPA$xAxis <- c(1,3,5,7,9)

#ggplot(uncert_cent_Mod, aes(x=years, y=habitat_south)) + geom_violin() + geom_line(data=cent_mean, aes(x=xAxis, y=habitat_south)) + 
 #geom_violin(data=uncert_cent_Mod, aes(x=years, y=habitat_south)) + geom_point(data=cent_mean, aes(x=xAxis, y=habitat_south), shape=1, size=3) 
#ggplot(uncert_cent_ModPA, aes(x=years, y=habitat_south)) + geom_violin() + geom_line(data=cent_meanPA, aes(x=xAxis, y=habitat_south)) + 
 #geom_violin(data=uncert_cent_ModPA, aes(x=years, y=habitat_south)) + geom_point(data=cent_meanPA, aes(x=xAxis, y=habitat_south), shape=1, size=3) 
# lobster: 1.62
# fluke: 1.15
# combined violin plot for both niche model types 
# NEED TO CHANGE THE X-AXIS, AND see if I can specify the tick marks on the right to include some smaller values
val <- 1.75
ggplot(uncert_cent_ModPA, aes(x=as.factor(xAxis-0.05), y=habitat)) + geom_violin() +
  geom_violin(data = uncert_cent_Mod, aes(x=as.factor(xAxis+.05), y=habitat^(1/val))) + # adjust the exponent denominator here and elsewhere in block to adjust scaling of the second plot and axis
  geom_line(data=cent_meanPA, aes(x=xAxis, y=habitat)) +
  geom_line(data=cent_mean, aes(x=xAxis, y=habitat^(1/val))) +
  geom_violin(data=uncert_cent_ModPA, aes(x=as.factor(xAxis-0.05), y=habitat)) +
  geom_violin(data = uncert_cent_Mod, aes(x=as.factor(xAxis+.05), y=habitat^(1/val)), fill='grey75') +
  geom_point(data=cent_meanPA, aes(x=xAxis, y=habitat), shape=1, size=2) +
  geom_point(data=cent_mean, aes(x=xAxis, y=habitat^(1/val)),size=2, fill="black", colour="black") +
  scale_y_continuous(expand = c(0, 0),  limits=c(0,2400), sec.axis = ~.^val) # make sure this matches denominators above, to scale secon axis correctly


# Save files used to make the graphs
filename <- paste('output/CEmodels_Uncertainty_2018/', species, '_figures_habSouth.RData', sep='')
save(uncert_final, SStable, uncert_cent_Mod, uncert_cent_ModPA, cent_mean, cent_meanPA, uncert_cent, file=filename)



# =====================================================================================================================
# Try the habitat analysis, but using the %change between periods as the variable of interest==========================

uncert_cent <- uncert_cent[,c(2,3,4,5,8,9)]
niche <- c('GAM_habitat_south','GAM_habitat_southPA','GLM_habitat_south','GLM_habitat_southPA')
hab_change <- data.frame(nicheMod=character(), rcp=character(), modelrun=character(), iter=integer(), time1=numeric(), time2=numeric(), time3=numeric(), time4=numeric(), stringsAsFactors = F)
hab_chang_fin <- NULL

for(i in 1:3){ # RCP
  rcp <- RCP[i]
  for(j in 1:18){ # GCMs
    gcm <- modelrun[j]
    for(k in 1:4){ #nicheMod
      nicModel <- niche[k]
      uncert_cent_sub <- uncert_cent[uncert_cent$rcp==rcp & uncert_cent$modelrun==gcm & uncert_cent$nicheModel==nicModel,]
      uncert_cent_cast <- dcast(uncert_cent_sub, rcp+modelrun+iter+nicheModel~years, value.var='habitat')  
      time1 <- ((uncert_cent_cast$'2021-2040' - uncert_cent_cast$'2007-2020')/uncert_cent_cast$'2007-2020')*100
      time2 <- ((uncert_cent_cast$'2041-2060' - uncert_cent_cast$'2007-2020')/uncert_cent_cast$'2007-2020')*100
      time3 <- ((uncert_cent_cast$'2061-2080' - uncert_cent_cast$'2007-2020')/uncert_cent_cast$'2007-2020')*100
      time4 <- ((uncert_cent_cast$'2081-2100' - uncert_cent_cast$'2007-2020')/uncert_cent_cast$'2007-2020')*100
      
      hab_change <- data.frame(nicheMod=nicModel, rcp=rcp, modelrun=gcm, iter=1:40, time1=time1, time2=time2, time3=time3, time4=time4, stringsAsFactors = F)
      hab_chang_fin <- rbind(hab_chang_fin, hab_change)
    }
  }
}

hab_chang_fin <- melt(hab_chang_fin, id.vars=c('nicheMod', 'rcp', 'modelrun', 'iter'), measure.vars=c('time1','time2','time3','time4'), variable.name='years', value.name='habitat')

hab_chang_fin$nicheModel <- as.factor(hab_chang_fin$nicheMod)
hab_chang_fin$rcp_fact <- as.factor(hab_chang_fin$rcp)
hab_chang_fin$modelrun_fact <- as.factor(hab_chang_fin$modelrun)
hab_chang_fin$years_fact <- as.factor(hab_chang_fin$years) # for regression tree only

# ======================================================================================================================
# REGRESSION TREE to see if any factors are contributing most of variance and should be excluded as 'wonky' =============================================================================================
tree_cent <- rpart(habitat~nicheModel + rcp_fact + modelrun_fact + years_fact, method='anova', data=hab_chang_fin)
complex.param <- tree_cent$cptable[which.min(tree_cent$cptable[,"xerror"]),"CP"] # calculate a penalty value for overfitting
tree_cent2 <- prune(tree_cent, cp=complex.param) # refit the model with the penalty value

plot(tree_cent2)
# text(tree_cent2, cex=.7, use.n=TRUE, all=TRUE)
plot(tree_cent2, uniform=T) # use this and add an inset with the actual branch lengths
text(tree_cent2, cex=.7, use.n=TRUE, all=TRUE)

#printcp(tree_cent2)
#plotcp(tree_cent2)
rsq.rpart(tree_cent2)
print(tree_cent2) # explains how everything was split_what went where
summary(tree_cent2) #includes relative importance of the variables_among a lot of other info


# Linear mixed effects models ==================================================
# fit linear models
mod_cent2020 <- dom_percChan('time1', hab_chang_fin)
mod_cent2040 <- dom_percChan('time2', hab_chang_fin)
mod_cent2060 <- dom_percChan('time3', hab_chang_fin)
mod_cent2080 <- dom_percChan('time4', hab_chang_fin)

SStable <- data.table(rbind(mod_cent2020, mod_cent2040, mod_cent2060, mod_cent2080))
SStable_trans <- transpose(SStable)
colnames(SStable_trans) <- SStable$years
SStable_trans <- SStable_trans[2:5,]
rownames(SStable_trans) <- names(SStable)[2:5]
SStable_trans <- data.matrix(SStable_trans)
barplot(SStable_trans, main="", ylab='Variation (Sum of Squares)', col=c("dodgerblue2","orange","green3","red"), legend = c('Niche model','GCM','RCP','Parameter'), args.legend=list(x='topleft')) 
box()
 

# VIOLIN PLOTS ========================================================================
cent_mean <- aggregate(list(habitat = hab_chang_fin$habitat), by=list(years = hab_chang_fin$years), FUN=mean) 
cent_mean$xAxis <- c(1,2,3,4)

ggplot(hab_chang_fin, aes(x=years_fact, y=habitat)) + geom_violin() + 
  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_line(data=cent_mean, aes(x=xAxis, y=habitat)) + 
  geom_violin(data=hab_chang_fin, aes(x=years_fact, y=habitat), fill='gray80') + geom_point(data=cent_mean, aes(x=xAxis, y=habitat), shape=21, size=3, fill='white')


# Save files used to make the graphs
filename <- paste('output/CEmodels_Uncertainty_2018/', species, '_figures_habSouth_percChan.RData', sep='')
save(uncert_final, SStable, hab_chang_fin, cent_mean, uncert_cent, file=filename)


 
# ====================================================================================
# INVESTIGATE PARAMETER UNCERTAINTY BY LOOKING AT PROJECTIONS ===================
# ====================================================================================
# Based on this, may go back up and then redo certain elements if some of the habitat models look bad

library(latticeExtra)
library(maps)

gcm <- modelrun[7] # use GFDL-CM3 as the standard
rcp <- 85

nicModel <- nicheMod[2] # Change this to '2' and then rerun the entire code section below

filename <- paste(projfolder, '/', species, '_rcp', rcp, '_', gcm, '_', nicModel, '_iter40.RData', sep='') 
load(filename) 
 
# Carve out Gmex and northern areas (course carving, not for manuscript but for a PPT)
biom_final <- biom_final[!biom_final$longitude < -82,]
biom_final <- biom_final[!(biom_final$longitude < -80.75 & biom_final$latitude < 27),]
biom_final <- biom_final[!biom_final$latitude > 47,]
biom_final <- biom_final[!biom_final$latitude < 31,]
PA_final <- PA_final[!PA_final$longitude < -82,]
PA_final <- PA_final[!(PA_final$longitude < -80.75 & PA_final$latitude < 27),]
PA_final <- PA_final[!PA_final$latitude > 47,]
PA_final <- PA_final[!PA_final$latitude < 31,]

xlimit = c(min(biom_final$longitude),max(biom_final$longitude))
ylimit = c(min(biom_final$latitude), max(biom_final$latitude))
#scale85 = seq(0, log(max(biom_final[,1:40][biom_final$year_bin=='2081-2100',])+1.1), length.out=20)
scale85 = seq(0, log(max(biom_final[,1:40][biom_final$year_bin=='2081-2100',])+1.01), length.out=20)
cols = colorRampPalette(colors = c('gray90', 'blue', 'dark blue', 'black'))

Projmap <- map('world', xlim=xlimit, ylim=ylimit, plot=F) 
Projmap <- data.frame(lon=Projmap$x, lat=Projmap$y)

 
# plot up all 40 iterations to use in making a gif
filename <- paste('figures/Uncertainty_mods/Param_uncert/', species, '_paramUncert_rcp',rcp, '_', gcm, '_', nicModel, 'LogBiomass_Nov2018.pdf', sep='')
pdf(height=4, width=3.5, bg='white', file=filename)

levelplot(log(biom_final[,1][biom_final$year_bin=='2081-2100']+1) ~ longitude*latitude, data=biom_final[biom_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(log(biom_final[,2][biom_final$year_bin=='2081-2100']+1) ~ longitude*latitude, data=biom_final[biom_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(log(biom_final[,3][biom_final$year_bin=='2081-2100']+1) ~ longitude*latitude, data=biom_final[biom_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(log(biom_final[,4][biom_final$year_bin=='2081-2100']+1) ~ longitude*latitude, data=biom_final[biom_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(log(biom_final[,5][biom_final$year_bin=='2081-2100']+1) ~ longitude*latitude, data=biom_final[biom_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(log(biom_final[,6][biom_final$year_bin=='2081-2100']+1) ~ longitude*latitude, data=biom_final[biom_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(log(biom_final[,7][biom_final$year_bin=='2081-2100']+1) ~ longitude*latitude, data=biom_final[biom_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(log(biom_final[,8][biom_final$year_bin=='2081-2100']+1) ~ longitude*latitude, data=biom_final[biom_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(log(biom_final[,9][biom_final$year_bin=='2081-2100']+1) ~ longitude*latitude, data=biom_final[biom_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(log(biom_final[,10][biom_final$year_bin=='2081-2100']+1) ~ longitude*latitude, data=biom_final[biom_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(log(biom_final[,11][biom_final$year_bin=='2081-2100']+1) ~ longitude*latitude, data=biom_final[biom_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(log(biom_final[,12][biom_final$year_bin=='2081-2100']+1) ~ longitude*latitude, data=biom_final[biom_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(log(biom_final[,13][biom_final$year_bin=='2081-2100']+1) ~ longitude*latitude, data=biom_final[biom_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(log(biom_final[,14][biom_final$year_bin=='2081-2100']+1) ~ longitude*latitude, data=biom_final[biom_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(log(biom_final[,15][biom_final$year_bin=='2081-2100']+1) ~ longitude*latitude, data=biom_final[biom_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(log(biom_final[,16][biom_final$year_bin=='2081-2100']+1) ~ longitude*latitude, data=biom_final[biom_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(log(biom_final[,17][biom_final$year_bin=='2081-2100']+1) ~ longitude*latitude, data=biom_final[biom_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(log(biom_final[,18][biom_final$year_bin=='2081-2100']+1) ~ longitude*latitude, data=biom_final[biom_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(log(biom_final[,19][biom_final$year_bin=='2081-2100']+1) ~ longitude*latitude, data=biom_final[biom_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(log(biom_final[,20][biom_final$year_bin=='2081-2100']+1) ~ longitude*latitude, data=biom_final[biom_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(log(biom_final[,21][biom_final$year_bin=='2081-2100']+1) ~ longitude*latitude, data=biom_final[biom_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(log(biom_final[,22][biom_final$year_bin=='2081-2100']+1) ~ longitude*latitude, data=biom_final[biom_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(log(biom_final[,23][biom_final$year_bin=='2081-2100']+1) ~ longitude*latitude, data=biom_final[biom_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(log(biom_final[,24][biom_final$year_bin=='2081-2100']+1) ~ longitude*latitude, data=biom_final[biom_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(log(biom_final[,25][biom_final$year_bin=='2081-2100']+1) ~ longitude*latitude, data=biom_final[biom_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(log(biom_final[,26][biom_final$year_bin=='2081-2100']+1) ~ longitude*latitude, data=biom_final[biom_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(log(biom_final[,27][biom_final$year_bin=='2081-2100']+1) ~ longitude*latitude, data=biom_final[biom_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(log(biom_final[,28][biom_final$year_bin=='2081-2100']+1) ~ longitude*latitude, data=biom_final[biom_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(log(biom_final[,29][biom_final$year_bin=='2081-2100']+1) ~ longitude*latitude, data=biom_final[biom_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(log(biom_final[,30][biom_final$year_bin=='2081-2100']+1) ~ longitude*latitude, data=biom_final[biom_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(log(biom_final[,31][biom_final$year_bin=='2081-2100']+1) ~ longitude*latitude, data=biom_final[biom_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(log(biom_final[,32][biom_final$year_bin=='2081-2100']+1) ~ longitude*latitude, data=biom_final[biom_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(log(biom_final[,33][biom_final$year_bin=='2081-2100']+1) ~ longitude*latitude, data=biom_final[biom_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(log(biom_final[,34][biom_final$year_bin=='2081-2100']+1) ~ longitude*latitude, data=biom_final[biom_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(log(biom_final[,35][biom_final$year_bin=='2081-2100']+1) ~ longitude*latitude, data=biom_final[biom_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(log(biom_final[,36][biom_final$year_bin=='2081-2100']+1) ~ longitude*latitude, data=biom_final[biom_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(log(biom_final[,37][biom_final$year_bin=='2081-2100']+1) ~ longitude*latitude, data=biom_final[biom_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(log(biom_final[,38][biom_final$year_bin=='2081-2100']+1) ~ longitude*latitude, data=biom_final[biom_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(log(biom_final[,39][biom_final$year_bin=='2081-2100']+1) ~ longitude*latitude, data=biom_final[biom_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(log(biom_final[,40][biom_final$year_bin=='2081-2100']+1) ~ longitude*latitude, data=biom_final[biom_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')

dev.off()


# Now with PA models ===============================================
scale85 = seq(0, max(PA_final[,1:40][PA_final$year_bin=='2081-2100',])+.01, length.out=20)

filename <- paste('figures/Uncertainty_mods/Param_uncert/', species, '_paramUncert_rcp',rcp, '_', gcm, '_', nicModel, 'PresAbs_Nov2018.pdf', sep='')
pdf(height=4, width=3.5, bg='white', file=filename)

levelplot(PA_final[,1][PA_final$year_bin=='2081-2100'] ~ longitude*latitude, data=PA_final[PA_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(PA_final[,2][PA_final$year_bin=='2081-2100'] ~ longitude*latitude, data=PA_final[PA_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(PA_final[,3][PA_final$year_bin=='2081-2100'] ~ longitude*latitude, data=PA_final[PA_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(PA_final[,4][PA_final$year_bin=='2081-2100'] ~ longitude*latitude, data=PA_final[PA_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(PA_final[,5][PA_final$year_bin=='2081-2100'] ~ longitude*latitude, data=PA_final[PA_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(PA_final[,6][PA_final$year_bin=='2081-2100'] ~ longitude*latitude, data=PA_final[PA_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(PA_final[,7][PA_final$year_bin=='2081-2100'] ~ longitude*latitude, data=PA_final[PA_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(PA_final[,8][PA_final$year_bin=='2081-2100'] ~ longitude*latitude, data=PA_final[PA_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(PA_final[,9][PA_final$year_bin=='2081-2100'] ~ longitude*latitude, data=PA_final[PA_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(PA_final[,10][PA_final$year_bin=='2081-2100'] ~ longitude*latitude, data=PA_final[PA_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(PA_final[,11][PA_final$year_bin=='2081-2100'] ~ longitude*latitude, data=PA_final[PA_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(PA_final[,12][PA_final$year_bin=='2081-2100'] ~ longitude*latitude, data=PA_final[PA_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(PA_final[,13][PA_final$year_bin=='2081-2100'] ~ longitude*latitude, data=PA_final[PA_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(PA_final[,14][PA_final$year_bin=='2081-2100'] ~ longitude*latitude, data=PA_final[PA_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(PA_final[,15][PA_final$year_bin=='2081-2100'] ~ longitude*latitude, data=PA_final[PA_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(PA_final[,16][PA_final$year_bin=='2081-2100'] ~ longitude*latitude, data=PA_final[PA_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(PA_final[,17][PA_final$year_bin=='2081-2100'] ~ longitude*latitude, data=PA_final[PA_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(PA_final[,18][PA_final$year_bin=='2081-2100'] ~ longitude*latitude, data=PA_final[PA_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(PA_final[,19][PA_final$year_bin=='2081-2100'] ~ longitude*latitude, data=PA_final[PA_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(PA_final[,20][PA_final$year_bin=='2081-2100'] ~ longitude*latitude, data=PA_final[PA_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(PA_final[,21][PA_final$year_bin=='2081-2100'] ~ longitude*latitude, data=PA_final[PA_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(PA_final[,22][PA_final$year_bin=='2081-2100'] ~ longitude*latitude, data=PA_final[PA_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(PA_final[,23][PA_final$year_bin=='2081-2100'] ~ longitude*latitude, data=PA_final[PA_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(PA_final[,24][PA_final$year_bin=='2081-2100'] ~ longitude*latitude, data=PA_final[PA_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(PA_final[,25][PA_final$year_bin=='2081-2100'] ~ longitude*latitude, data=PA_final[PA_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(PA_final[,26][PA_final$year_bin=='2081-2100'] ~ longitude*latitude, data=PA_final[PA_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(PA_final[,27][PA_final$year_bin=='2081-2100'] ~ longitude*latitude, data=PA_final[PA_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(PA_final[,28][PA_final$year_bin=='2081-2100'] ~ longitude*latitude, data=PA_final[PA_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(PA_final[,29][PA_final$year_bin=='2081-2100'] ~ longitude*latitude, data=PA_final[PA_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(PA_final[,30][PA_final$year_bin=='2081-2100'] ~ longitude*latitude, data=PA_final[PA_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(PA_final[,31][PA_final$year_bin=='2081-2100'] ~ longitude*latitude, data=PA_final[PA_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(PA_final[,32][PA_final$year_bin=='2081-2100'] ~ longitude*latitude, data=PA_final[PA_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(PA_final[,33][PA_final$year_bin=='2081-2100'] ~ longitude*latitude, data=PA_final[PA_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(PA_final[,34][PA_final$year_bin=='2081-2100'] ~ longitude*latitude, data=PA_final[PA_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(PA_final[,35][PA_final$year_bin=='2081-2100'] ~ longitude*latitude, data=PA_final[PA_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(PA_final[,36][PA_final$year_bin=='2081-2100'] ~ longitude*latitude, data=PA_final[PA_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(PA_final[,37][PA_final$year_bin=='2081-2100'] ~ longitude*latitude, data=PA_final[PA_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(PA_final[,38][PA_final$year_bin=='2081-2100'] ~ longitude*latitude, data=PA_final[PA_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(PA_final[,39][PA_final$year_bin=='2081-2100'] ~ longitude*latitude, data=PA_final[PA_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
levelplot(PA_final[,40][PA_final$year_bin=='2081-2100'] ~ longitude*latitude, data=PA_final[PA_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
          at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')

dev.off()

 

# =====================================================================================================
# =====================================================================================================

# Below is code that produces tandem plots_still needs work
#abc <- levelplot(log(biom_final[,40][biom_final$year_bin=='2081-2100']+1) ~ longitude*latitude, data=biom_final[biom_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
 #         at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
#cde <- levelplot(PA_final[,40][PA_final$year_bin=='2081-2100'] ~ longitude*latitude, data=PA_final[PA_final$year_bin=='2081-2100',], xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
 #         at = scale85, col.regions=cols, colorkey = FALSE) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')
#update(c(abc, cde, x.same = TRUE), layout = c(2, 1), ylab = NULL) 


# I'm not using this one any more =================================================================================================

# BELOW IS AN ALTERNATE TO THE BAR CHART
SStableb_melt <- data.frame(SStableb)
SStableb_melt$SS_source <- rownames(SStableb)
SStableb_melt <- data.frame(cbind(SS_source = rownames(SStableb), SStableb), stringsAsFactors = F)
rownames(SStableb_melt) <- NULL
colnames(SStableb_melt) <- c('2007-2020','2021-2040','2041-2060','2061-2080','2081-2100','SS_source')

library(ggplot2)

SStableb_melt <- melt(SStableb_melt, id.vars=c('SS_source'), value.name='SS')
SStableb_melt$xAxis <- as.numeric(SStableb_melt$variable)
SStableb_melt$categ <- rep(c('b','c','d','a','e'), 5)
ggplot(SStableb_melt, aes(x=xAxis, y=SS, fill=categ)) + geom_area()


# DO ANOTHER STACKED BAR CHART THAT DOES NOT SUM TO 1 (i.e. absolute variance values), that are centered on the mean centroid value
# Similar to Cheung et al. 2016

uncert_cent_Mod <- uncert_cent[!(uncert_cent$nicheModel=='GLM_centroid_lat'),]
hist(uncert_cent_Mod$centroid[uncert_cent_Mod$years=='2007-2020'])

# need to remake the SS table to sum to .9 and drop the residual error - Adjust if I switch to violins, then I'll probably want to sum to 1 (but still drop resid error)?
# May not need this anymore
SStablec <- data.frame(anova(mod_cent2020b)[1:4,2] / sum(anova(mod_cent2020b)[1:4,2]), anova(mod_cent2040b)[1:4,2] / sum(anova(mod_cent2040b)[1:4,2]), anova(mod_cent2060b)[1:4,2] / sum(anova(mod_cent2060b)[1:4,2]), 
                       anova(mod_cent2080b)[1:4,2] / sum(anova(mod_cent2080b)[1:4,2]), anova(mod_cent2100b)[1:4,2] / sum(anova(mod_cent2100b)[1:4,2]))
SStablec <- data.matrix(SStablec) 
colnames(SStablec) <- c('2007-2020', '2021-2040', '2041-2060', '2061-2080', '2081-2100')
SStablec <- SStablec * 0.9
rownames(SStablec) <- c('Niche_mod', 'model_run', 'rcp', 'iter')

# Need a loop to form the basis for this graph
cent_mean <- data.frame(rcp_low=numeric(), rcp_high=numeric(), modelrun_low=numeric(), modelrun_high=numeric(), niche_low=numeric(), niche_high=numeric(), iter_low=numeric(), iter_high=numeric())
for(i in 1:5){ # 5 time blocks
  time_bin <- years[i]
  rcp_quants <- c(0.5 - SStablec['rcp',time_bin]/2, 0.5 + SStablec['rcp',time_bin]/2)
  rcp_vals <- quantile(uncert_cent_Mod$centroid[uncert_cent_Mod$years==time_bin], probs=rcp_quants)
  gcm_quants <- c(rcp_quants[1] - SStablec['model_run',time_bin]/2, rcp_quants[2] + SStablec['model_run',time_bin]/2)
  gcm_vals <- quantile(uncert_cent_Mod$centroid[uncert_cent_Mod$years==time_bin], probs=gcm_quants)
  niche_quants <- c(gcm_quants[1] - SStablec['Niche_mod',time_bin]/2, gcm_quants[2] + SStablec['Niche_mod',time_bin]/2)
  niche_vals <- quantile(uncert_cent_Mod$centroid[uncert_cent_Mod$years==time_bin], probs=niche_quants)
  iter_quants <- c(niche_quants[1] - SStablec['iter',time_bin]/2, niche_quants[2] + SStablec['iter',time_bin]/2)
  iter_vals <- quantile(uncert_cent_Mod$centroid[uncert_cent_Mod$years==time_bin], probs=iter_quants)
  cent_mean[i,] <- data.frame(rcp_low=as.numeric(rcp_vals[1]), rcp_high=as.numeric(rcp_vals[2]), modelrun_low=as.numeric(gcm_vals[1]), modelrun_high=as.numeric(gcm_vals[2]), niche_low=as.numeric(niche_vals[1]), niche_high=as.numeric(niche_vals[2]), iter_low=as.numeric(iter_vals[1]), iter_high=as.numeric(iter_vals[2]))
}

cent_mean <- cbind(cent_mean, cent_MEAN)

plot(mean_cent~xAxis, type='l', ylim=c(32.25,39), ylab='Centroid range - 90% of projections', xlab='Time period', cent_mean)
arrows(x0=cent_mean$xAxis, y0=cent_mean$iter_low, x1=cent_mean$xAxis, y1=cent_mean$iter_high, length=0, col='red', lwd=25)
arrows(x0=cent_mean$xAxis, y0=cent_mean$niche_low, x1=cent_mean$xAxis, y1=cent_mean$niche_high, length=0, col='dodgerblue2', lwd=25)
arrows(x0=cent_mean$xAxis, y0=cent_mean$modelrun_low, x1=cent_mean$xAxis, y1=cent_mean$modelrun_high, length=0, col='orange', lwd=25)
arrows(x0=cent_mean$xAxis, y0=cent_mean$rcp_low, x1=cent_mean$xAxis, y1=cent_mean$rcp_high, length=0, col='green3', lwd=25)
points(mean_cent~xAxis, type='l', cent_mean)



# =====================================================================================
# uneeded stuff for estimating lat-cell areas
# now get the latitude values to match up as sig-figs altered in raster
#lats <- unique(biom_final$latitude)
#grid_area$latitude <- NA

#for(i in 1:nrow(grid_area)){
#  lat <- grid_area$lat[i]
#  lat_adj <- data.frame(lats=lats, diff=lats-lat)
#  lat_adj$min <- abs(lat_adj$diff)
#  grid_area$latitude[i] <- lat_adj$lats[lat_adj$min==min(lat_adj$min)][1]
#}
#grid_area$dups <- duplicated(grid_area$latitude)
#grid_area <- grid_area[!grid_area$dups == T,]
#grid_area$lat <- NULL 
#grid_area$dups <- NULL
#rownames(grid_area) <- NULL

#grid_area$latitude <- round(grid_area$latitude, digits=3)# to prevent rounding error in the loop
#grid_area <- rbind(grid_area, c(14.413, 62.125))# this latitude row got knocked off in raster
#rm(grid_rast,lats, lat, lat_adj)
